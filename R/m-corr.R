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
#' names(genera)[1] <- 'rn'
#' names(metabolites)[1] <- 'rn'
#' b <- bSet(b = genera)
#' m <- mSet(m = metabolites)
#' res <- corr(m, b, method = 'spearman')
#' @export
#' @param m mSet class
#' @param b bSet class
#' @param parallel logical
#' @param ncore integer number of core,default is 4
#' @param method character default is spearman
setGeneric('corr', function(m, b, method = 'spearman', parallel = FALSE, ncore=4){
        standardGeneric('corr')
})


#' corr
#'
#' @description genetic methods to perform the correlation test.
#'
#' @importFrom parallel makeCluster stopCluster parLapply clusterEvalQ clusterExport
#' @importFrom stats setNames p.adjust
#' @importFrom data.table tstrsplit
#' @importFrom doParallel registerDoParallel
#' @docType methods
#' @param m mSet class
#' @param b bSet class
#' @param parallel logical
#' @param ncore integer number of core,default is 4
#' @param method character default is spearman
#' @return correlation matrix
setMethod(f = 'corr', signature = c(m='mSet', b='bSet'),
          definition = function(m, b, method = 'spearman',
                                parallel = FALSE, ncore=4){
                  m_dt <- m@dt
                  mn <- features(m)
                  m_samples <- samples(m)
                  b_dt <- b@dt
                  bn <- features(b)
                  b_samples <- samples(b)
                  stopifnot("metabolites and OTU abundance matrix has no consistant samples" = identical(sort(m_samples), sort(b_samples)))

                  samp <- intersect(m_samples, b_samples)
                  m_dt <- t(m_dt[, samp, with = FALSE])
                  colnames(m_dt) <- mn

                  b_dt <- t(b_dt[, samp, with = FALSE])
                  colnames(b_dt) <- bn

                  quiteCorr <- function(a, b){
                          options(warn=-1)
                          res <- psych::corr.test(a, b, method = method)
                          options(warn=0)
                          res <- cor2df(res)
                  }
                  if (!parallel) {
                          res <- quiteCorr(b_dt, m_dt)
                  } else {
                          message("Using parallel with ", ncore, "cores!\n")
			  registerDoParallel(cores = ncore)
                          cl <- makeCluster(ncore)
                          on.exit(stopCluster(cl))
                          clusterExport(cl, list('m_dt','b_dt', 'cor2df'), envir = environment())

                          res <-
                                  parLapply(cl, colnames(b_dt),
                                            fun = function(x){
                                                    quiteCorr(setNames(data.frame(b_dt[,x]),x),
                                                              m_dt)
                                            })
                          res <- do.call('rbind', res)
                  }
                  res[,c("b", "m") := tstrsplit(pair, split = " : ")]
                  res$padj <- p.adjust(res$p, method = 'BH')
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
                padj <- data.table::as.data.table(res$p.adj,
                                      keep.rownames = TRUE)  |>
                        data.table::melt(id = 'rn')
                padj <- padj[,
                             list(pair=interaction(padj$rn,padj$variable,sep = ' : '),
                                  padj = padj$value)
                ]
                p <- data.table::as.data.table(res$p,
                                   keep.rownames = TRUE) |>
                        data.table::melt(id = 'rn')
                p <- p[,
                       list(pair=interaction(p$rn,p$variable,sep = ' : '),
                            p = p$value)
                ]
                r <- data.table::as.data.table(res$r,
                                   keep.rownames = TRUE) |>
                        data.table::melt(id = 'rn')
                r <- r[,
                       list(pair=interaction(r$rn,r$variable,sep = ' : '),
                            rho = r$value)
                ]
                data.table::setkey(r,'pair')
                data.table::setkey(p, 'pair')
                data.table::setkey(padj, 'pair')
                padj[p][r]

        }
