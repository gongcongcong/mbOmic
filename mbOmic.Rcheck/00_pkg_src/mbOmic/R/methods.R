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

#' @rdname clean_analytes
#' @title clean_analytes
#' @param object, mbSet class
#' @param m_thres, integer, Removal of metabolites only measured in m_thres of samples
#' @param b_thres, integer, Removal of OTU only measured in b_thres of samples
#' @docType methods
#' @examples
#' mb <- clean_analytes(mbSet, 2, 2)
#' @return mbSet
#' @import data.table
#' @export
setGeneric('clean_analytes',function(object, m_thres=2, b_thres=2){
        standardGeneric('clean_analytes')
})
#' @rdname clean_analytes
#' @description Removal of analytes only measured in few samples
#' @docType methods
setMethod('clean_analytes', 'mbSet', function(object, m_thres=2, b_thres=2){
        m <- object@m
        b <- object@b
        samples <- intersect(colnames(m), colnames(b))
        index_m <- apply(m[,..samples],1, function(x)sum(x!=0)>m_thres)
        m <- m[index_m]
        index_b <- apply(b[,..samples],1, function(x)sum(x!=0)>b_thres)
        b <- b[index_b]
        mbSet(m, b)
})

#' @rdname corr
#' @title corr
#' @param object, mbSet class
#' @param parallel, logical
#' @param ncore, integer number of core
#' @docType methods
#' @importFrom psych corr.test
#' @return the correlation data table
#' @examples
#' res <- corr(object, method = 'spearman')
#' @export
setGeneric('corr', function(object, parallel = FALSE, ncore=4, ...){
        standardGeneric('corr')
})


#' @description genetic methods to perform the correlation test
#' @rdname corr
#' @importFrom psych corr.test
#' @importFrom magrittr %>%
#' @importFrom parallel makeCluster parLapply clusterEvalQ clusterExport
#' @docType methods
setMethod(f = 'corr', signature(object = 'mbSet'),
          definition = function(object, parallel = FALSE, ncore=4, ...){
                  m <- t(m.extra(object))
                  colnames(m) <- object@ms
                  b <- t(b.extra(object))
                  colnames(b) <- object@bs
                  if (!parallel) {
                          res <- corr.test(b, m, ...)
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
                   function(object, group, minN=2, power = NULL,powers = 1:30,
                            threshold.d = 0.05, threshold = 0.8,message = FALSE, ...) {
                           standardGeneric('wgcna')
                   })

#' wgcna
#' @title wgcna
#' @importFrom WGCNA goodSamplesGenes blockwiseModules cor
#' @importFrom psych corr.test
#' @return network

setMethod('wgcna', 'mbSet', function(object,group,minN=2,power = NULL,powers = 1:30,
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

