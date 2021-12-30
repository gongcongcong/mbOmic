#' mNum
#'
#' Calculate the metabolites number from the mbSet class object.
#'
#' @param object mbSet class
#' @examples
#' library(data.table)
#' path <- system.file('extdata', 'metabolites_and_genera.rda', package = 'mbOmic')
#' load(path)
#' object <-
#'        mbSet(
#'              m = metabolites,
#'              b = genera
#'              )
#' mNum(object)
#' @export
#' @return metabolites number
mNum <-
        function(object) {
                nrow(m(object))
        }
#' bNum
#'
#' Calculate the OTU number from the mbSet class object.
#'
#' @param object mbSet class
#' @export
#' @examples
#' library(data.table)
#' path <- system.file('extdata', 'metabolites_and_genera.rda', package = 'mbOmic')
#' load(path)
#' object <-
#'        mbSet(
#'              m = metabolites,
#'              b = genera
#'              )
#' bNum(object)
#' @return OTU number
bNum <-
        function(object) {
                nrow(b(object))
        }
#' @title m
#'
#' Extract the metabolites abundance matrix from the mbSet class object.
#'
#' @param object mbSet class
#' @return data.table metabolites aundance matrix
#' @export
#' @examples
#' library(data.table)
#' path <- system.file('extdata', 'metabolites_and_genera.rda', package = 'mbOmic')
#' load(path)
#' object <-
#'        mbSet(
#'              m = metabolites,
#'              b = genera
#'              )
#' m(object)
m <-
        function(object) {
                object@m
        }
#' @title b
#'
#' Extract the OTU abundance matrix from the mbSet class object.
#'
#' @param object mbSet class
#' @return data.table OTU aundance matrix
#' @export
#' @examples
#' library(data.table)
#' path <- system.file('extdata', 'metabolites_and_genera.rda', package = 'mbOmic')
#' load(path)
#' object <-
#'        mbSet(
#'              m = metabolites,
#'              b = genera
#'              )
#' b(object)
b <-
        function(object) {
                object@b
        }
#' @title mbSamples
#' extract the samples names from a mbSet class object.
#' @param object mbSet class
#' @export
#' @return samples number
#' @examples
#'
#' library(data.table)
#' path <- system.file('extdata', 'metabolites_and_genera.rda', package = 'mbOmic')
#' load(path)
#' object <-
#'        mbSet(
#'              m = metabolites,
#'              b = genera
#'              )
#' mbSamples(object)
mbSamples <-
        function(object){
                object@samples
        }


#' clean_analytes
#'
#' Remove unqualified metabolites and OTU which only
#' measured in few samples.
#'
#' @param object, mbSet class
#' @param m_thres, integer, Removal of metabolites only measured in m_thres of samples
#' @param b_thres, integer, Removal of OTU only measured in b_thres of samples
#' @import data.table
#' @export
#' @return mbSet object
#' @examples
#' library(data.table)
#' path <- system.file('extdata', 'metabolites_and_genera.rda', package = 'mbOmic')
#' load(path)
#' object <-
#'        mbSet(
#'              m = metabolites,
#'              b = genera
#'              )
#' clean_analytes(object)
clean_analytes <-
        function(object, m_thres=2, b_thres=2){
                m <- m(object)
                b <- b(object)
                samples <- mbSamples(object)
                index_m <- apply(m[, samples, with = FALSE], 1, function(x) sum(x != 0) > m_thres)
                cat(length(m$rn[!index_m]), "metabolites Removed: ", m$rn[!index_m], '\n')
                index_b <- apply(b[, samples, with = FALSE], 1, function(x) sum(x != 0) > b_thres)
                cat(length(b$rn[!index_b]), "microbiota Removed: ", b$rn[!index_b], '\n')
                b <- b[index_b]
                m <- m[index_m]
                mbSet(m, b)
}


#' @title quiteRun
#'
#' Run expression in quite mode.
#'
#' @param x expression
#' @return no return

quiteRun <- function(x) {
        sink(tempfile())
        on.exit(sink())
        invisible(force(x))
}




