#' @title dist_jsd
#'
#' Calculating Jensen-Shannon divergence basing the formulate
#' \eqn{JSD(P || Q) = \sqrt{0.5 \times KLD (P||\frac{P + Q}{2}) + 0.5 \times KLD(Q||\frac{P+Q}{2})} }.
#'
#' @title
#' @param b bacterial abundance matrix
#' @examples
#' data(iris)
#' dist <- dist_jsd(iris[, 1:4])
#' @return Jensen-Shannon divergence data frame
#' @export
dist_jsd <- function(b) {
        KLD <- function(x, y)
                sum(x * log(x / y))
        JSD <-
                function(x, y)
                        sqrt(0.5 * KLD(x, (x + y) / 2) + 0.5 * KLD(y, (x + y) / 2))
        ## function apply jsd  to vector and matrix
        dist_ <- function(x, d, f) {
                lsD <- d
                if (!inherits(d, "data.frame")) lsD <- as.data.frame(lsD)
                vapply(
                        as.list(lsD),
                        FUN = function(i) {
                                f(i, x)
                        },
                        numeric(1)
                )
        }
        dist <- apply(b, 2, dist_, f = JSD, d = b)
        dist <- as.data.frame(dist)
        class(dist) <- c("data.frame", "jsd")
        dist
}

#' cluster_jsd
#'
#' `cluster_jsd` cluster the sample using the pam based the JSD.
#'
#' @param dist distance matrix, if miss distance matrix the
#' bacterial abundance matrix `b` was used to calculate the
#' Jensen-Shannon divergence
#' @param b bacterial abundance matrix, if the distance matrix
#'  `dist` was given, it is useless
#' @param k number of cluster
#' @examples
#'
#' data(iris)
#' ## using the matrix to cluster samples.
#' res_cluster <- cluster_jsd(b=t(iris[,1:4]), k = 3)
#' ## or using the resoult of dist_jsd
#' ## dist <- dist_jsd(iris[, 1:4])
#' ## res_cluster <- cluster_jsd(dist = iris[,1:4], k = 3)
#' table(res_cluster, iris[, 5])
#'
#' @return cluster vector
#' @export
cluster_jsd <- function(dist, b, k) {
        if (missing(dist)) dist <- dist_jsd(b)
        if (!inherits(dist, "jsd")) stop("distance should been Jensen-Shannon divergence; the dist matrix need has the jsd class")
        as.vector(cluster::pam(dist, k = k, diss = TRUE)$clustering)
}

#' estimate_k
#'
#' To estimate the optimal cluster number ,
#' `estimate_k` takes advantage of two measures,
#' Calinski-Harabasz (CH) Index and silhouette coefficients.
#'
#' @param b bacterial abundance matrix
#' @param verK number vector of cluster
#' @return list
#' @examples
#' data(iris)
#' ret <- estimate_k(b = t(iris[, 1:4]), 1:10)
#' @export
estimate_k <- function(b, verK = 2:10) {
        if (min(verK) <= 1) verK <- verK[verK >1]
        verK <- sort(verK) #verctor, number of cluster
        dist <- dist_jsd(b) #matrix, jsd
        verCluster <- integer(ncol(b)) #cluster verctor using the optimal number
        optK <- 0 #optimal number of cluster
        optCHI <- 0 #CHI with optimal number ofcluster
        Silhouette <- 0 #silhouette coefficients with optimal number ofcluster
        ret <- lapply(verK, function(x) {
                verTmpCluster <- cluster_jsd(dist = dist, k = x)
                CHI <- clusterSim::index.G1(t(b), verTmpCluster, d = dist, centrotypes = "medoids")
                silhouette <- cluster::silhouette(verTmpCluster, dist)[,3] |> mean()# silhouette coefficients
                if (CHI > optCHI) {
                        optCHI <<- CHI
                        optK <<- x
                        verCluster <<- verTmpCluster
                        Silhouette <<- silhouette
                        }
                list(CHI = CHI, silhouette = silhouette)
                })
        verCHI = vapply(ret, function(x) x[['CHI']], numeric(1)) #CHI verctor
        verSilhouette = vapply(ret, function(x) x[['silhouette']], numeric(1)) #vector silhouette coefficients
        names(ret) <- verK
        names(verCluster) <- colnames(b)
        ret <- list(verCHI = verCHI,
                    verK = verK,
                    verSilhouetteCoef = verSilhouette,
                    jsd = dist,
                    optK = optK,
                    optCHI = optCHI,
                    Silhouette = Silhouette,
                    verOptCluster = verCluster)
        class(ret) <- "verCHI"
        ret
}

#' print.verCHI
#'
#' print the optialmal numer of cluster, max CHI, and Silhouette
#' and plot the Calinski-Harabasz (CH) Index and silhouette coefficients simultaneously.
#'
#' @param x estimate_k output
#' @param verbose logical
#' @param plotting logical
#' @param cluster cluster result
#' @param ... print
#' @examples
#' data(iris)
#' ret <- estimate_k(b = t(iris[,1:4]), 1:10)
#' ret
#' @return no return
#' @export
#' @importFrom graphics lines text par
#' @importFrom stats quantile

print.verCHI <- function(x, ..., verbose = TRUE, plotting = TRUE, cluster) {
        verCHI <- x
        if (!verbose){
                sink(tempfile())
                on.exit(sink())
        }
        par(mfrow = c(1, 2))
        n <- length(verCHI$verK)
        cat("optimal number of cluster: ", verCHI$optK, "\n")
        cat("Max CHI: ", verCHI$optCHI, "\n")
        cat("Silhouette: ", verCHI$Silhouette, "\n")
        if (!missing(cluster)){
                cat("\nConfusion Matrix: \n")
                print(table(cluster, verCHI$verOptCluster))
        }
        if (plotting) {
                #plot CHI
                plot(verCHI$verK, verCHI$verCHI, type = 'h',
                     xlab = "Number of cluster",
                     ylab = "Calinski-Harabasz (CH) index")
                lines(verCHI$optK, verCHI$optCHI, type = "h", col = "red")
                text(quantile(verCHI$verK, 0.65), quantile(verCHI$verCHI, 0.8),
                     label = paste0("Optimal number of cluster is ", verCHI$optK, "\n",
                                    "Max CHI is ", round(verCHI$optCHI, 2)), cex = 0.4)
                #plot silhouette coefficients
                plot(verCHI$verK, verCHI$verSilhouetteCoef, type = 'h',
                     xlab = "Number of cluster", ylab = "Silhouette coefficients")
                lines(verCHI$optK, verCHI$optCHI, type = "h", col = "red")

        }
        invisible()
}

#' enterotyping
#'
#' identifing enterotype basing the OTU abundance using
#' the cluster anlysis refer to the XXX works.
#'
#' @param b bacterial abundance matrix
#' @param cluster cluster vector
#' @param threshold abundance threshold
#' @return list
#' @export
#' @examples
#'
#' dat <- read.delim('http://enterotypes.org/ref_samples_abundance_MetaHIT.txt')
#' dat <- impute::impute.knn(as.matrix(dat), k = 100)
#' rest <- read.table(system.file('extdata', 'enterotype.txt', package = 'mbOmic'))
#' rest <- rest[colnames(dat$data), ]
#' res2 <- estimate_k(dat$data+0.001)
#' table(res2$verOptCluster, rest$ET)
#'
#' ret=enterotyping(dat$data+0.01, res2$verOptCluster, threshold = 0.02)
#' ret
#'
enterotyping <- function(b, cluster, threshold = 0.25) {
        vecBacter <- rownames(b)
        grp <- split(colnames(b), cluster)
        # load("../data/lsGenus.Rda")
        # Firmicutes <- lsGenus$Firmicutes
        # FirmicutesIdx <- vecBacter %in% Firmicutes
        FirmicutesIdx <- grepl("firmicutes", vecBacter, ignore.case = TRUE)
        BacteroidesIdx <- grepl("bacteroides", vecBacter, ignore.case = TRUE)
        PrevotellaIdx <- grepl("prevotella", vecBacter, ignore.case = TRUE)
        dfMeans <- lapply(grp, function(x){
                c(Bacteroides = mean(b[BacteroidesIdx, x], na.rm = TRUE),
                     Prevotella = mean(b[PrevotellaIdx, x], na.rm = TRUE),
                     Firmicutes = mean(b[FirmicutesIdx, x], na.rm = TRUE))

                }) |> as.data.frame()
        ret <- apply(dfMeans, 1, function(x){
                maxValueIdx <- which.max(x)
                if(length(maxValueIdx) != 1){
                        return(NA_real_)
                }
                if (x[maxValueIdx] <= threshold)  maxValueIdx <- NA_real_
                maxValueIdx
        })
        vecCluster <- unique(cluster)
        vecClusterLeft <- setdiff(vecCluster, ret)
        ret[is.na(ret)] <- vecClusterLeft[which.max(as.numeric(table(cluster))[vecClusterLeft])]
        tmp <- names(ret)
        names(tmp) <- ret
        list(enterotypes = ret,
             data = data.frame(samples = names(cluster),
                               cluster=cluster,
                               enterotype = tmp[as.character(cluster)]))
}

