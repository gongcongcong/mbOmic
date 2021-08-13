#'@title print
#'@export
#'@return logical
print.mbSet <-
        function(obj){
        cat("Samples: ", obj@samples,"\n")
        cat("number of OTU: ", obj@nb, "\n")
        cat("number of metabolites: ", obj@nm, "\n")
}

#' obtain the m express matrix from mbSet
#' create the mbSet class
#' @title mbSet
#' @import data.table
#' @param m, the dataframe of metabolites abundance profile
#' @param b, the dataframe of OTU abundance profile
#' @examples
#' library(data.table)
#' path <- system.file('data',package = 'mbOmic')
#' b <- readxl::read_xlsx(file.path(path,
#'                     '41598_2021_85433_MOESM2_ESM.xlsx'),
#'                     sheet = 'Table S2',skip = 1)
#' m <- readxl::read_xlsx(file.path(path,
#'                                  '41598_2021_85433_MOESM2_ESM.xlsx'),
#'                     sheet = 'Table S4',skip = 1)
#' setDT(b)
#' setDT(m)
#' cn <- c('Name_des',grep(x = names(m), 'BS|SS.*', value = TRUE))
#' m <- m[,..cn]
#' object <-
#'        mbSet(
#'              m = m,
#'              b = b
#'              )
#' @export
#' @return mbSet class
mbSet <-
        function(m, b, pData = NULL){
                m <- as.data.table(m)
                b <- as.data.table(b)
                nb <- nrow(b)
                nm <- nrow(m)
                samples <- intersect(names(m), names(b))
                m.id <- setdiff(names(m), samples)
                b.id <- setdiff(names(b), samples)
                rownames(m) <- as.vector(unlist(m[,..m.id]))
                rownames(b) <- as.vector(unlist(b[,..b.id]))
                ms <- rownames(m)
                bs <- rownames(b)
                if(is.null(pData)){
                        pData <- data.table()
                }
                new('mbSet', m = m , b = b,
                    pData = pData, nb = nb, nm = nm,
                    samples = samples, ms = factor(ms), bs = factor(bs))
        }



#' obtain the b slot from mbSet
#' @title b.extra
#' @return the otu abundance data table
b.extra <- function(obj) {s <- obj@samples; return(obj@b[,..s])}
#' obtain the m slot from mbSet
#' @title m.extra
#' @return the metabolites abundance data table
m.extra <- function(obj) {s <- obj@samples; return(obj@m[,..s])}
#' obtain the samples
#' @title samples
#' @return a vector containing sample names
samples.extra <- function(obj) return(obj@samples)

#' @title nb
#' @return number of otu
nb.extra <- function(obj) return(obj@nb)
#' @title nm
#' @return number of metabolites
nm.extra <- function(obj) return(obj@nm)

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

#' test the metabolite profile
#' @importFrom WGCNA goodGenes
#' @return filter bad metabolites
m.filter <-
        function(m, minNSamples=2){
                suppressMessages({
                        goodMetabolites <- goodGenes(m, verbose=2,
                                                     minNSamples = minNSamples,
                                                     minNGenes = 1, weights = NULL, minFraction = 1/2,
                                                     tol = NULL, minRelativeWeight = 0.1, indent = 0)
                })
                if (sum(!goodMetabolites)>0){
                               cat("\nRemoving metabolites(occured in less than ",minNSamples,
                                   ") : ",
                                   paste(names(m)[!sum(!goodMetabolites)>0], sep = ", \n")
                                   )
                } else {
                        cat("\nAll metabolites are ok!\n")
                }
                m[, ..goodMetabolites]
        }

#' Automatic network construction and module detection by one-step method
#' @importFrom WGCNA pickSoftThreshold
#' @return power

pickST <-
        function(m, threshold.d = 0.05, threshold = 0.8,plot = TRUE,powers =NULL){
                if (is.null(powers)) powers <-  c(c(1:10), seq(from = 12, to=20, by=2))
                invisible(capture.output(
                                sft <-pickSoftThreshold(m,dataIsExpr = TRUE, powerVector = powers, verbose = 5)
                                ))
                powers <- sft$fitIndices[,1]
                sft <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
                res <- which(abs(diff(sft)) <= threshold.d)+1
                if (sum(sft[res]>=threshold)>0) {
                        res <- res[which(sft[res]>=threshold)][1]
                        cat("\n Scale Free Topology Model Fit,signed R^2: ", sft[res], '\n')
                } else {
                        cat("\nNo proper softThreshold was detected!!")
                        res <- NA_real_
                }
                if(plot) {
                        plot(powers, sft, xlab = "Soft Threshold (power)",
                             ylab = expression(paste("Scale Free Topology Model Fit,signed ",
                                                     R^2)
                             ),
                             main = "Scale independence",
                             ylim = c(-1,1)
                        )
                        points(res, sft[res], pch = 16, col = 'darkred')
                        text(res, sft[res]-0.1, labels = res, col = 'red')
                        text(1,threshold+0.1,labels=threshold,col = 'blue')
                        text(17,0,labels=paste(
                                "Threshold: ", threshold,
                                "\ntol: ", threshold.d,
                                sep = ''
                        ))
                        abline(h = threshold, col = 'blue', lty =2)
                }
                res
        }



