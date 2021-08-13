#' @title show
#' @export
#' @return print
setMethod(f = "show",
          signature = "mbSet",
          definition = function(object){
                  cat("Samples: ", object@samples,"\n")
                  cat("number of OTU: ", object@nb, "\n")
                  cat("number of metabolites: ", object@nm, "\n")
          })


#' @rdname cor.test
#' @title cor.test
#' @param object, mbSet class
#' @param parallel, logical
#' @param ncore, integer number of core
#' @docType methods
#' @importFrom psych corr.test
#' @return the correlation data table
#' @examples
#' res <- cor.test(object, method = 'spearman')
#' @export
setGeneric('cor.test', function(object, parallel = FALSE, ncore=4, ...){
        standardGeneric('cor.test')
})


#' @description genetic methods to perform the correlation test
#' @rdname cor.test
#' @importFrom psych corr.test
#' @importFrom magrittr %>%
#' @importFrom parallel makeCluster parLapply clusterEvalQ clusterExport
#' @docType methods
setMethod(f = 'cor.test', signature(object = 'mbSet'),
          definition = function(object, parallel = FALSE, ncore=4, ...){
                  m <- t(m.extra(object))
                  colnames(m) <- object@ms
                  b <- t(b.extra(object))
                  colnames(b) <- object@bs
                  if (!parallel){
                          res <- corr.test(b,m, ...)
                          res <- cor2df(res)
                  } else {
                          cl <- makeCluster(ncore)
                          on.exit(stopCluster(cl))
                          clusterEvalQ(cl,
                                       {
                                               library(magrittr,quietly = TRUE)
                                               library(data.table,quietly = TRUE)
                                               library(psych,quietly = TRUE)
                                        })
                          clusterExport(cl, list('m','b'),
                                        envir = environment())

                          res <-
                                  parLapply(cl, object@bs,
                                           fun = function(x){
                                                   res <- corr.test(setNames(data.frame(b[,x]),
                                                                             x),
                                                                           m, ...)
                                                   cor2df(res)
                                           }) %>%
                                  do.call('rbind', .)
                  }
                  res[, `:=`("b" = sapply(strsplit(as.character(res$pair), split = ' : '), `[`, 1),
                              "m" = sapply(strsplit(as.character(res$pair), split = ' : '), `[`, 2))]
                  res <- res[,.(b,m,rho,p, padj)]
                  res

          })

#' wgcna
#' @title wgcna
#' @docType methods
#' @rdname wgcna
#' @author Congcong Gong
#' @param object, mbSet class
#' @param minN, integer the minimum number of sample
#' @param power, integer, if the pickSoftThreshold function (WGCNA) can find appropriate power, this param is invalid
#' @examples
#' net <- wgcna(object, minN = 2, power = 9,message = FALSE)
#' @export
#' @return network

setGeneric(name = 'wgcna',
           def =
                   function(object, group, minN, power = NULL,powers = 1:30,
                            threshold.d = 0.05, threshold = 0.8,message = FALSE, ...) {
                           standardGeneric('wgcna')
                   })

#' wgcna
#' @title wgcna
#' @importFrom WGCNA goodSamplesGenes blockwiseModules cor
#' @importFrom psych corr.test
#' @return network

setMethod('wgcna', 'mbSet', function(object,group,minN,power = NULL,powers = 1:30,
                                     threshold.d = 0.05, threshold = 0.8,message =FALSE, ...){
        m <- m.extra(object)
        m <- m.filter(m,minN)
        cat("\nOne: detect the power for softThreshold\n")

        invisible(
                capture.output(
                                {
                                        pw <- pickST(m, plot = TRUE,
                                              threshold.d=threshold.d,
                                              threshold=threshold,powers=powers)
                                }
                        )
                )
        if (is.na(pw) & is.null(power)) stop("No power detected! pls set the power parameter")
        if (is.na(pw) | !is.null(power)) pw <- power
        # id <- featureNames(mExpSet(object))
        cat("\nusing the powser: ", pw, "to constructe net!\n")
        cat("\nTwo: Automatic network construction and module detection by one-step method\n")
        if(message){
                net <-  blockwiseModules(datExpr = t(m),
                                         power = pw,...)
        } else {
                invisible(capture.output({
                        net <-  blockwiseModules(datExpr = t(m),
                                                 power = pw,...)
                }))
        }

        #structure(net$colors, MEs = net$MEs, class = 'mb.module')
        names(net$colors) <- object@ms
        net$ms <- object@ms
        net$bs <- object@bs
        net
})

