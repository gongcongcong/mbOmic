#' @title coExpress
#' @docType methods
#' @author Congcong Gong
#' @param object, mbSet class
#' @param minN, integer the minimum number of sample
#' @param power, integer, if the pickSoftThreshold function (WGCNA) can find appropriate power, this param is invalid
#' @export
#' @return network

setGeneric(name = 'coExpress',
           def =
                   function(object, power = NULL, powerVec = 1:30,
                            threshold = 0.8, message = TRUE, ...) {
                           standardGeneric('coExpress')
                   })

#' @title coExpress
#' @importFrom WGCNA goodSamplesGenes blockwiseModules cor
#' @return network
setMethod('coExpress', 'mbSet', function(object,power = NULL, powerVec = 1:30,
                                         threshold = 0.8, message = TRUE, ...){

        if (!message){
                sink(tempfile())
                on.exit(sink())
        }
        object <- clean_analytes(object)
        m <- m(object)
        mN <- m$rn
        m <- m[,rn:=NULL]
        cat("\nOne: detect the power for softThreshold\n")
        pw <- pickST(
                m,
                plot = FALSE,
                threshold.d = 0.05,
                threshold = threshold,
                powers = powerVec
        )

        if (is.na(pw)) {
                if (is.null(power)) {
                        stop("No power detected! pls set the power parameter\n")
                } else(
                        pw <- power
                )
        }
        # id <- featureNames(mExpSet(object))
        cat("\nusing the power: ", pw, "to constructe net!\n")
        cat("\nTwo: Eetwork construction and module detection was done\n")
        net <-  quiteRun(blockwiseModules(datExpr = t(m),
                                          power = pw, ...))
        res <- net$colors
        res <- res[res!='grey']
        res <- data.frame(table(res))
        names(res) <- c("Module","Size")
        cat("There are ", nrow(res), "modules were constructed: \n")
        apply(res, 1, cat, '\n')
        #structure(net$colors, MEs = net$MEs, class = 'mb.module')
        names(net$colors) <- mN
        net

})

#' @title pickST
#' @description Automatic network construction and module detection by one-step method
#' @importFrom WGCNA pickSoftThreshold
#' @return power

pickST <-
        function(m, threshold.d = 0.05, threshold = 0.8, plot = TRUE,powers =NULL){
                if (is.null(powers)) powers <-  c(c(1:10), seq(from = 12, to=20, by=2))
                suppressWarnings(
                        sft <- quiteRun(pickSoftThreshold(m,dataIsExpr = TRUE, powerVector = powers, verbose = 5))
                )
                powers <- sft$fitIndices[,1]
                sft <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
                res <- which(abs(diff(sft)) <= threshold.d)+1
                if (sum(sft[res]>=threshold)>0) {
                        res <- res[which(sft[res]>=threshold)][1]
                        cat("\n Scale Free Topology Model Fit,signed R^2: ", sft[res], '\n')
                } else {
                        cat("\nNo proper softThreshold was detected!!\n")
                        res <- NA_real_
                }
                if(plot) {
                        plot(powers, sft, xlab = "Soft Threshold (power)",
                             ylab = expression(paste("Scale Free Topology Model Fit,signed ",
                                                     R^2)
                             ),
                             main = "Scale independence",
                             ylim = c(-2,1)
                        )
                        points(res, sft[res], pch = 16, col = 'darkred')
                        text(res, sft[res]-0.1, labels = res, col = 'red')
                        text(1,threshold+0.2,labels=threshold,col = 'blue')
                        text(20,-1,labels=paste(
                                "Threshold: ", threshold,
                                "\ntol: ", threshold.d,
                                sep = ''
                        ))
                        abline(h = threshold, col = 'blue', lty =2)
                }
                res
        }


