#' @title dist_jsd
#'
#' Calculating Jensen-Shannon divergence basing the formulate
#' \eqn{JSD(P || Q) =
#'       \sqrt{0.5 \times KLD (P||\frac{P + Q}{2}) + 0.5 \times KLD(Q||\frac{P+Q}{2})} }.
#'
#' @title
#' @param b a `bSet` class
#' @examples
#' data(iris)
#' b <- data.table::as.data.table(iris[1:6, 1:4])
#' b <- bSet(b = b,
#'           Features = letters[1:6],
#'           Samples = LETTERS[1:4])
#' dist <- dist_jsd(b)
#' @return Jensen-Shannon divergence data frame
#' @export
dist_jsd <- function(b) {
        stopifnot("b must be bSet class" = inherits(b, 'bSet'))
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
        dist <- apply(b@dt, 2, dist_, f = JSD, d = b@dt)
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
#' b <- bSet(b = data.table::as.data.table(t(iris[1:6, 1:4])),
#'           Features = letters[1:4],
#'           Samples = LETTERS[1:6])
#' res_cluster <- cluster_jsd(b=b, k = 3)
#' ## or using the resoult of dist_jsd
#' ## dist <- dist_jsd(iris[, 1:4])
#' ## res_cluster <- cluster_jsd(dist = iris[,1:4], k = 3)
#' table(res_cluster, iris[1:6, 5])
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
#' b <- bSet(b = data.table::as.data.table(t(iris[1:6, 1:4])),
#'           Features = letters[1:4],
#'           Samples = LETTERS[1:6])
#' ret <- estimate_k(b = b, 2:3)
#' @export
estimate_k <- function(b, verK = 2:10) {
        stopifnot("b must be bSet class" = inherits(b, 'bSet'))
        if (min(verK) <= 1) verK <- verK[verK > 1]
        verK <- sort(verK) #verctor, number of cluster
        dist <- dist_jsd(b) #matrix, jsd
        verCluster <- integer(ncol(b)) #cluster verctor using the optimal number
        optK <- 0 #optimal number of cluster
        optCHI <- 0 #CHI with optimal number ofcluster
        Silhouette <- 0 #silhouette coefficients with optimal number ofcluster
        ret <- lapply(verK, function(x) {
                verTmpCluster <- cluster_jsd(dist = dist, k = x)
                CHI <- clusterSim::index.G1(t(b@dt), verTmpCluster, d = dist, centrotypes = "medoids")
                silhouette <- cluster::silhouette(verTmpCluster, dist)[,3] |> mean()# silhouette coefficients
                if (CHI > optCHI) {
                        optCHI <<- CHI
                        optK <<- x
                        verCluster <<- verTmpCluster
                        Silhouette <<- silhouette
                        }
                list(CHI = CHI, silhouette = silhouette)
                })
        verCHI <-  vapply(ret, function(x) x[['CHI']], numeric(1)) #CHI verctor
        verSilhouette  <-  vapply(ret, function(x) x[['silhouette']], numeric(1)) #vector silhouette coefficients
        names(ret) <- verK
        names(verCluster) <- samples(b)
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
#' b <- bSet(b = data.table::as.data.table(t(iris[1:6, 1:4])),
#'           Features = letters[1:4],
#'           Samples = LETTERS[1:6])
#' ret <- estimate_k(b =b, 2:3)
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
        message("optimal number of cluster: ", verCHI$optK, "\n")
        message("Max CHI: ", verCHI$optCHI, "\n")
        message("Silhouette: ", verCHI$Silhouette, "\n")
        if (!missing(cluster)){
                message("\nConfusion Matrix: \n")
                message(table(cluster, verCHI$verOptCluster))
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
#'dat <- read.delim('http://enterotypes.org/ref_samples_abundance_MetaHIT.txt')
# dat <- impute::impute.knn(as.matrix(dat), k = 100)
# dat <- as.data.frame(dat$data+0.001)
# setDT(dat, keep.rownames = TRUE)
# dat <- bSet(b =  dat)
# rest <- read.table(system.file('extdata', 'enterotype.txt', package = 'mbOmic'))
# rest <- rest[samples(dat),]
# res2 <- estimate_k(dat)
# ret <- enterotyping(dat, res2$verOptCluster)
# ret

enterotyping <- function(b, cluster, threshold = 0.02) {
        stopifnot("b must be bSet class" = inherits(b, 'bSet'))
        vecBacter <- features(b)
        grp <- split(samples(b), cluster)
        samples <- names(cluster)
        dt <- b@dt
        rownames(dt) <- features(b)
        colnames(dt) <- samples(b)
        b <- dt
        # load("../data/lsGenus.Rda")
        # Firmicutes <- lsGenus$Firmicutes
        # FirmicutesIdx <- vecBacter %in% Firmicutes
        RuminococcusIdx <- grepl("^ruminococcus$", vecBacter, ignore.case = TRUE)
        BacteroidesIdx <- grepl("^bacteroides$", vecBacter, ignore.case = TRUE)
        PrevotellaIdx <- grepl("^prevotella$", vecBacter, ignore.case = TRUE)
        dfMeans <- lapply(grp, function(x){
                c(Bacteroides = mean(as.numeric(b[BacteroidesIdx, x, with = FALSE]), na.rm = TRUE),
                  Prevotella = mean(as.numeric(b[PrevotellaIdx, x, with = FALSE]), na.rm = TRUE),
                  Ruminococcus = mean(as.numeric(b[RuminococcusIdx, x, with = FALSE]), na.rm = TRUE))
                }) |> as.data.frame()
        ret <- data.frame(
                          max = apply(dfMeans, 1, max),
                          which = apply(dfMeans, 1, which.max)
                          )
        ret$cluster <- NA
        tmpDupIdx <- ret$which == ret$which[duplicated(ret$which)]
        if (sum(tmpDupIdx) == 3) {
                ret$cluster[which.max(ret$max)] <- paste0('cluster ', ret$which[which.max(ret$max)])
        } else if (sum(tmpDupIdx) == 2) {
                ret$cluster[!tmpDupIdx] <- paste0('cluster ', ret$which[!tmpDupIdx])
                nm <- rownames(ret)[tmpDupIdx][which.max(ret$max[tmpDupIdx])]
                ret[nm, 'cluster'] <- paste0('cluster ', ret[nm, 'which'])
        } else {
                ret$cluster <- paste0('cluster ', ret$which)
        }
        ret$cluster <- ifelse(ret$max < threshold, NA_character_, ret$cluster)
        setDT(ret, keep.rownames = "Enterotype")
        cluster <- data.table(Samples = names(cluster),
                              which = cluster)
        out <- merge(cluster, ret, on = 'which')
        out$Enterotype <- ifelse(is.na(out$cluster), NA_character_, out$Enterotype)
        uncluster <- samples[!samples %in% out$Samples]
        list(enterotypes = ret,
             data = out[!is.na(out$Enterotype), c('Samples', 'Enterotype', 'cluster'), with = FALSE],
             UnIdentifiedSamples = uncluster
             )
}

