#' check function for mbSet
#' @title check_mbSet
#' @return logical
check_mbSet <-
        function(object) {
                errors <- character()
                mSN <- sort(setdiff(names(object@m),'rn'))
                bSN <- sort(setdiff(names(object@b),'rn'))
                anno.m <- object@anno.m
                anno.b <- object@anno.b
                m <- object@m
                b <- object@b
                nm <- nrow(m)
                nb <- nrow(b)
                if (nrow(anno.m) != nm | nrow(anno.b) != nb){
                        msg <- "Features number is not same between abundance data and annotation data!"
                        errors <- c(errors, msg)
                }

                if (!identical(sort(m$rn), sort(anno.m$rn)) | !identical(sort(b$rn), sort(anno.b$rn))){
                        msg <- "Features names is not identical!"
                        errors <- c(errors,msg)
                }
                if (!identical(mSN, bSN)){
                        msg <- "Samples is not match between metabolites and microbiotas"
                        errors <- c(errors, msg)
                }
                if (!identical(mSN, object@samples) | !identical(bSN, object@samples)){
                        msg <- "Samples names is wrong!"
                        errors <- c(errors, msg)
                }
                if (length(errors) == 0) TRUE else errors
        }

#' @slot b  abundance of OTU, column response to samples, row response to microbiota
#' @slot m abundance abundance, column response to samples, row response to metabolites
#' @slot anno.m metabolites features
#' @slot anno.b microbiota features
#' @slot samples samples names
#' @name mbSet
#' @title mbSet
#' @docType class
#' @export
setClass(
        Class = 'mbSet',
        representation = representation(
                m = 'data.table',
                b = 'data.table',
                anno.m = 'data.table',
                anno.b = 'data.table',
                samples = "character"
        ),
        validity = check_mbSet,
        prototype = prototype(
                m = data.table(),
                b = data.table(),
                anno.m = data.table(),
                anno.b = data.table(),
                samples = vector('character')
        )
)


#' @title create a mbSet class
#' @import data.table
#' @param m, the dataframe of metabolites abundance profile
#' @param b, the dataframe of OTU abundance profile
#' @examples
#' path <- system.file('data',package = 'mbOmic')
#' load(file.path(path,'metabolites_and_genera.rda'))
#' object <-
#'        mbSet(
#'              m = metabolites,
#'              b = genera
#'              )
#' @export
#' @return mbSet object
mbSet <-
        function(m, b, anno.m, anno.b, m.f, b.f){
                m <- as.data.table(m,keep.rownames = TRUE)
                b <- as.data.table(b,keep.rownames = TRUE)
                samples <- intersect(names(m), names(b))
                if (length(samples) == 0 | identical(samples, 'rn')) {
                        stop("There are no samples in common of m and b!")
                }
                if (!'rn' %in% names(m)) {
                        if (missing(m.f)) {
                                m.f <- setdiff(names(m), samples)
                                if (length(m.f) != 1) {
                                        stop("The m does not has rownames and the m.f not setted!")
                                }
                        }
                        names(m)[names(m)==m.f] <- 'rn'
                }
                if (!'rn' %in% names(b)) {
                        if (missing(b.f)) {
                                b.f <- setdiff(names(b), samples)
                                if (length(b.f) != 1) {
                                        stop("The b does not has rownames and the m.f not setted!")
                                }
                        }
                        names(b)[names(b)==b.f] <- 'rn'
                }
                samples <- setdiff(samples, 'rn')
                samples <- sort(samples)
                id <- c('rn', samples)
                m <- m[, ..id]
                b <- b[, ..id]
                if (missing(anno.m)) {
                        anno.m <- data.table(rn=m$rn)
                }
                if (missing(anno.b)) {
                        anno.b <- data.table(rn=b$rn)
                }
                new('mbSet', m = m , b = b,
                    anno.b = anno.b, anno.m = anno.m,
                    samples = samples)
        }

#' @title show
#' @return print
setMethod(f = "show",
          signature = "mbSet",
          definition = function(object){
                  cat("Samples(",length(mbSamples(object)),"): ", mbSamples(object),"\n")
                  cat("number of OTU: ", bNum(object), "\n")
                  cat("number of metabolites: ", mNum(object), "\n")
          })



