#' @title check_mbSet
#'
#' `check_mbSet` check the mbSet class.
#'
#' @param object mbSet
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
                        msg <- "Features number is't same between
                                abundance and annotation data"
                        errors <- c(errors, msg)
                }
                sorted_m_rn <- identical(sort(m$rn), sort(anno.m$rn))
                sorted_b_rn <- identical(sort(b$rn), sort(anno.b$rn))
                if (!sorted_m_rn | !sorted_b_rn){
                        msg <- "Features names is not identical!"
                        errors <- c(errors,msg)
                }
                if (!identical(mSN, bSN)){
                        msg <- "Samples is't match between
                                metabolites and microbiotas"
                        errors <- c(errors, msg)
                }
                identical_m_sam <- identical(mSN, object@samples)
                identical_b_sam <- identical(bSN, object@samples)
                if (!identical_m_sam | !identical_b_sam){
                        msg <- "Samples names is wrong!"
                        errors <- c(errors, msg)
                }
                if (length(errors) == 0) TRUE else errors
        }

#' mbSet
#'
#' `mbSet` is the main S4 class of mbOmic. It contains the
#'  metabolites (`m`) and OTU (`b`) abundance matrix.
#'  The `anno.m` and `anno.b` are the annotation of the `m`
#'  and `b`, respectially, which need to implement on the later.
#'
#' @slot b  OTU, column response to samples, row response to microbiota
#' @slot m metabolites, column response to samples, row response to metabolites
#' @slot anno.m metabolites features
#' @slot anno.b microbiota features
#' @slot samples samples names
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


#' mbSet
#' @description create a mbSet class.
#' @import data.table
#' @param m the dataframe of metabolites abundance profile
#' @param b the dataframe of OTU abundance profile
#' @param anno.m dataframe for annotaing metabolites
#' @param anno.b dataframe for annotating OTU
#' @param m.f character th column contains metabolite ID
#' @param b.f character the column contains OTU ID
#' @examples
#' path <- system.file('extdata',
#'                     'metabolites_and_genera.rda',
#'                      package = 'mbOmic')
#' load(path)
#' object <-
#'        mbSet(
#'              m = metabolites,
#'              b = genera
#'              )
#' @export
#' @return mbSet object
#' @importFrom methods new
#' @importFrom stats setNames
mbSet <-
        function(m, b, anno.m, anno.b, m.f, b.f){
                m <- as.data.table(m, keep.rownames = TRUE)
                b <- as.data.table(b, keep.rownames = TRUE)
                samples <- intersect(names(m), names(b))
                if (length(samples) == 0 | identical(samples, 'rn')) {
                        stop("There are no samples
                             in common of m and b!")
                }
                if (!'rn' %in% names(m)) {
                        if (missing(m.f)) {
                                m.f <- setdiff(names(m), samples)
                                if (length(m.f) != 1) {
                                        stop("The m does not has rownames
                                              and the m.f not setted!")
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
                m <- m[, id, with = FALSE]
                b <- b[, id, with = FALSE]
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

#' @description show.
#' @title show
#' @return print
#' @param object mbSet
setMethod(f = "show", signature = "mbSet", definition = function(object){
        cat("Samples(",length(mbSamples(object)),"): ", mbSamples(object),"\n")
        cat("number of OTU: ", bNum(object), "\n")
        cat("number of metabolites: ", mNum(object), "\n")
          })



