#' corr
#'
#' Calculating spearman correlation between each OTU and
#' each metabolites.
#'
#'
#' @docType methods
#' @importFrom psych corr.test
#' @return the correlation data table
#' @examples
#' library(data.table)
#' path <- system.file('extdata','metabolites_and_genera.rda',package = 'mbOmic')
#' load(path)
#' object <-
#'        mbSet(
#'              m = metabolites,
#'              b = genera
#'              )
#' res <- corr(object, method = 'spearman')
#' @export
#' @param object mbSet class
#' @param parallel logical
#' @param ncore integer number of core,default is 4
#' @param method character default is spearman
setGeneric('corr', function(object, method = 'spearman', parallel = FALSE, ncore=4){
        standardGeneric('corr')
})


#' corr
#'
#' @description genetic methods to perform the correlation test.
#'
#' @importFrom parallel makeCluster stopCluster parLapply clusterEvalQ clusterExport
#' @importFrom stats setNames
#' @importFrom data.table tstrsplit
#' @docType methods
#' @param object mbSet class
#' @param parallel logical
#' @param ncore integer number of core,default is 4
#' @param method character default is spearman
#' @return correlation matrix
setMethod(f = 'corr', signature(object = 'mbSet'),
          definition = function(object, method = 'spearman',
                                parallel = FALSE, ncore=4){
                  m <- m(object)
                  mn <- m$rn
                  samples <- mbSamples(object)
                  m <- t(m[, samples, with = FALSE])
                  colnames(m) <- mn
                  b <- b(object)
                  bn <- b$rn
                  b <- t(b[, samples, with = FALSE])
                  colnames(b) <- bn
                  quiteCorr <- function(a, b){
                          options(warn=-1)
                          res <- psych::corr.test(a, b, method = method)
                          options(warn=0)
                          res <- cor2df(res)
                  }
                  if (!parallel) {
                          res <- quiteCorr(b, m)
                  } else {
                          cat("Using parallel with ", ncore, "cores!\n")
                          cl <- makeCluster(ncore)
                          on.exit(stopCluster(cl))
                          clusterEvalQ(cl,
                                       {
                                               library(magrittr, quietly = TRUE)
                                               library(data.table, quietly = TRUE)
                                       })
                          clusterExport(cl, list('m','b'), envir = environment())

                          res <-
                                  parLapply(cl, colnames(b),
                                            fun = function(x){
                                                    quiteCorr(setNames(data.frame(b[,x]),x),
                                                              m)
                                            })
                          res <- do.call('rbind', res)
                  }
                  res[,c("b", "m") := tstrsplit(pair, split = " : ")]
                  res[, c('b', 'm', 'rho', 'p', 'padj'), with = FALSE]

          })

#' cor2df
#'
#' Coverting the correlation result to data table class.
#'
#' @return data table
#' @param res output of corr.test
cor2df <-
        function(res){
                padj <- as.data.table(res$p.adj,
                                      keep.rownames = TRUE)  |>
                        melt(id = 'rn')
                padj <- padj[,
                             list(pair=interaction(padj$rn,padj$variable,sep = ' : '),
                                  padj = padj$value)
                ]
                p <- as.data.table(res$p,
                                   keep.rownames = TRUE) |>
                        melt(id = 'rn')
                p <- p[,
                       list(pair=interaction(p$rn,p$variable,sep = ' : '),
                            p = p$value)
                ]
                r <- as.data.table(res$r,
                                   keep.rownames = TRUE) |>
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
