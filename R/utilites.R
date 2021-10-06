#' @rdname extract
#' @title mNum
#' @export
mNum <-
        function(object) {
                nrow(m(object))
        }
#' @rdname extract
#' @title bNum
#' @export
bNum <-
        function(object) {
                nrow(b(object))
        }
#' @title m
#' @rdname extract
#' @export
m <-
        function(object) {
                object@m
        }
#' @title b
#' @rdname extract
#' @export
b <-
        function(object) {
                object@b
        }
#' @title mbSamples
#' @rdname extract
#' @export
mbSamples <-
        function(object){
                object@samples
        }

#' @title clean_analytes
#' @param object, mbSet class
#' @param m_thres, integer, Removal of metabolites only measured in m_thres of samples
#' @param b_thres, integer, Removal of OTU only measured in b_thres of samples
#' @import data.table
#' @export
clean_analytes <-
        function(object, m_thres=2, b_thres=2){
                m <- m(object)
                b <- b(object)
                samples <- mbSamples(object)
                index_m <- apply(m[, ..samples], 1, function(x) sum(x != 0) > m_thres)
                cat(length(m$rn[!index_m]), "metabolites Removed: ", m$rn[!index_m], '\n')
                index_b <- apply(b[, ..samples], 1, function(x) sum(x != 0) > b_thres)
                cat(length(b$rn[!index_b]), "microbiota Removed: ", b$rn[!index_b], '\n')
                b <- b[index_b]
                m <- m[index_m]
                mbSet(m, b)
}


#' @title quiteRun

quiteRun <- function(x) {
        sink(tempfile())
        on.exit(sink())
        invisible(force(x))
}




