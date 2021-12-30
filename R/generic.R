#' dim.mbSet
#'
#' Outputing dimension of `m` and `b`.
#'
#' @export
#' @return dimension
#' @param x mbSet class
setMethod('dim','mbSet',
          function(x) {
                  m.dim <- dim(m(x))
                  b.dim <- dim(b(x))
                  return(list(m=m.dim, b=b.dim))
          }
          )


