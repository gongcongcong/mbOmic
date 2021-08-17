#' check function for mbSet
#' @title check_mbSet
#' @return logical
check_mbSet <-
        function(object) {
                errors <- character()
                mSN <- sort(names(object@m))
                bSN <- sort(names(object@b))
                if (nrow(object@pData) != 0){
                        pSN <- names(object@pData)
                        if (!all(mSN %in% pSN) | !all(bSN %in% pSN))  {
                                msg <- "The names of m or b is not all in pData's!"
                                errors <- c(errors, msg)
                        }
                }

                if (length(errors) == 0) TRUE else errors
        }
#'  An S4 class to represent a bank account.
#' @slot b  abundance of OTU
#' @slot m abundance abundance
#' @slot pData phenodata
#' @name mbSet
#' @rdname mbSet
#' @export
mbSet <-
        setClass(
                Class = 'mbSet',
                representation = representation(m = 'data.table',
                                             b = 'data.table',
                                             pData = 'data.table',
                                             samples = "character",
                                             nb = 'integer',
                                             nm = 'integer',
                                             ms = 'factor',
                                             bs = 'factor'),
                validity = check_mbSet
        )






