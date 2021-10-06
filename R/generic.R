#' @title dim
#' @docType method
#' @export
setMethod('dim','mbSet',
          function(x) {
                  m.dim <- dim(m(x))
                  b.dim <- dim(b(x))
                  return(list(m=m.dim, b=b.dim))
          }
          )


