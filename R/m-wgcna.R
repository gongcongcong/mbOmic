#' coExpress
#'
#' `coExpress` identify the co-expression metabolites by
#' performing the metabolites adjant matrix basing their
#' relative abundace. This process was mainly implemented
#' using the WGCNA package.
#'
#' @docType methods
#' @author Congcong Gong
#' @export
#' @return network
#' @param object mbSet class
#' @param power integer if the pickSoftThreshold function (WGCNA) can find appropriate power, this param is invalid
#' @param powerVec vector was passed to PickST function to get the power value
#' @param threshold numeric as the threshold to filter power value
#' @param message logical whether to show verbose info
#' @param plot logical whether plot in PickST function
#' @param ... additional arguments passed to WGCNA
#' @examples
#' library(data.table)
#' path <- system.file('extdata', 'metabolites_and_genera.rda', package = 'mbOmic')
#' load(path)
#' names(genera)[1] <- 'rn'
#' names(metabolites)[1] <- 'rn'
#' b <- bSet(b = genera)
#' m <- mSet(m = metabolites)
#' res <- corr(m, b, method = 'spearman')
#' net <- coExpress(m, minN = 2, power = 9, message = FALSE)
setGeneric(name = 'coExpress',
           def =
                   function(object, power = NULL, powerVec = seq_len(30),
                            threshold = 0.8, message = TRUE, plot = FALSE, ...) {
                           standardGeneric('coExpress')
                   })

#' coExpress
#'
#' @importFrom WGCNA goodSamplesGenes blockwiseModules cor
#' @param object mbSet class
#' @param power integer, if the pickSoftThreshold function (WGCNA) can find appropriate power, this param is invalid
#' @param powerVec vector was passed to PickST function to get the power value
#' @param threshold numeric as the threshold to filter power value
#' @param message logical, whether to show verbose info
#' @param plot logical, whether plot in PickST function
#' @param ... args passed to WGCNA
#' @return network
setMethod('coExpress', 'Set', function(object, power = NULL,
                                       powerVec, threshold = 0.8,
                                       message = TRUE, plot = FALSE,
                                       ...){

        stopifnot("object should be mSet" = inherits(object, 'mSet'))
        if (!message){
                sink(tempfile())
                on.exit(sink())
        }
        object <- clean_analytes(object)
        dt <- object@dt
        mN <- features(object)
        message("\nOne: detect the power for softThreshold\n")
        if (missing(powerVec)) {
                powerVec <-  c(seq_len(10), seq(from = 12, to=20, by=2))
                }
        if (is.null(power)) {
                pw <- quiteRun(
                        pickST(
                                dt,
                                plot = plot,
                                threshold.d = 0.05,
                                threshold = threshold,
                                powers = powerVec
                        )

                )
        } else pw <- power
        if (is.na(pw)) {
                stop("No power detected! pls set the power parameter\n")
        }
        # id <- featureNames(mExpSet(object))
        message("\nusing the power: ", pw, "to constructe net!\n")
        message("\nTwo: Network construction and module detection was done\n")
        net <-  quiteRun(blockwiseModules(datExpr = t(dt),
                                          power = pw, ...))
        res <- net$colors
        res <- res[res!='grey']
        res <- data.frame(table(res))
        names(res) <- c("Module","Size")
        message("====> There are ", nrow(res), " modules were constructed: \n")
        apply(res, 1, function(x) message("====|| ", x, "\n"))
        # structure(net$colors, MEs = net$MEs, class = 'mb.module')
        names(net$colors) <- mN
        net

})

#' pickST
#'
#' Picking up the Soft threshold of metabolites abundance matrix.
#'
#' @importFrom WGCNA pickSoftThreshold
#' @return integer
#' @importFrom graphics points text abline
#' @param m data.table metabolites abundance
#' @param threshold.d numeric rol
#' @param threshold numeric threshold
#' @param plot logical whether to plot
#' @param powers vector

pickST <-
        function(m, threshold.d = 0.05, threshold = 0.8, plot = TRUE,powers = NULL){
                quiteRun(
                        sft <- quiteRun(pickSoftThreshold(m,dataIsExpr = TRUE, powerVector = powers, verbose = 5))
                )
                powers <- sft$fitIndices[, 1]
                sft <- -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2]
                res <- which(abs(diff(sft)) <= threshold.d) + 1
                if (sum(sft[res] >= threshold) > 0) {
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
                        abline(h = threshold, col = 'blue', lty = 2)
                }
                res
        }


