#' @title corr
#' @param object, mbSet class
#' @param parallel, logical
#' @param ncore, integer number of core,default is 4
#' @docType methods
#' @importFrom psych corr.test
#' @return the correlation data table
#' @examples
#' library(data.table)
#' path <- system.file('data',package = 'mbOmic')
#' load(file.path(path,'metabolites_and_genera.rda'))
#' object <-
#'        mbSet(
#'              m = metabolites,
#'              b = genera
#'              )
#' res <- corr(object, method = 'spearman')
#' @export
setGeneric('corr', function(object, parallel = FALSE, ncore=4, ...){
        standardGeneric('corr')
})


#' @description genetic methods to perform the correlation test
#' @title corr
#' @importFrom magrittr %>%
#' @importFrom parallel makeCluster stopCluster parLapply clusterEvalQ clusterExport
#' @docType methods
setMethod(f = 'corr', signature(object = 'mbSet'),
          definition = function(object, parallel = FALSE, ncore=4, ...){
                  m <- m(object)
                  mn <- m$rn
                  samples <- mbSamples(object)
                  m <- m[,..samples]
                  m <- t(m)
                  colnames(m) <- mn

                  b <- b(object)
                  bn <- b$rn
                  b <- b[,..samples]
                  b <- t(b)
                  colnames(b) <- bn

                  if (!parallel) {
                          options(warn=-1)
                          res <- psych::corr.test(b, m, ...)
                          options(warn=0)
                          res <- cor2df(res)
                  } else {
                          cat("Using parallel with ", ncore, "cores!\n")
                          cl <- makeCluster(ncore)
                          on.exit(stopCluster(cl))
                          clusterEvalQ(cl,
                                       {
                                               library(magrittr,quietly = TRUE)
                                               library(data.table,quietly = TRUE)
                                       })
                          clusterExport(cl, list('m','b'),
                                        envir = environment())

                          res <-
                                  parLapply(cl, colnames(b),
                                            fun = function(x){
                                                    options(warn=-1)
                                                    res <- psych::corr.test(setNames(data.frame(b[,x]),
                                                                              x),
                                                                     m, ...)
                                                    options(warn=0)
                                                    cor2df(res)
                                            }) %>%
                                  do.call('rbind', .)
                          }
                  res[,
                      `:=`("b" = sapply(
                              strsplit(
                                      as.character(res$pair), split = ' : '
                              ), `[`, 1),
                           "m" = sapply(
                                   strsplit(
                                           as.character(res$pair), split = ' : '
                                   ), `[`, 2))]

                  res <- res[, .(b, m, rho, p, padj)]
                  res

          })

#' @title cor2df
#' @return data table
cor2df <-
        function(res){
                padj <- as.data.table(res$p.adj,
                                      keep.rownames = TRUE) %>%
                        melt(id = 'rn')
                padj <- padj[,
                             list(pair=interaction(padj$rn,padj$variable,sep = ' : '),
                                  padj = padj$value)
                ]
                p <- as.data.table(res$p,
                                   keep.rownames = TRUE) %>%
                        melt(id = 'rn')
                p <- p[,
                       list(pair=interaction(p$rn,p$variable,sep = ' : '),
                            p = p$value)
                ]
                r <- as.data.table(res$r,
                                   keep.rownames = TRUE) %>%
                        melt(id = 'rn')
                r <- r[,
                       list(pair=interaction(r$rn,r$variable,sep = ' : '),
                            rho = r$value)
                ]
                setkey(r,'pair')
                setkey(p, 'pair')
                setkey(padj, 'pair')
                padj[p][r]

        }
