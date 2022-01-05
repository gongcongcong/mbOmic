#' Set
#'
#' `Set` is a virtual class as the base class to
#' extend the `mSet` and `bSet`.
#'
#' @slot Samples character a character vector contains the samples
#' @slot Features character a character vector contains the features
#' @slot dt data.table
#' @import data.table
#' @docType class
#' @export
#' @return no return
#' @examples
#' ##`Set` is virtual class.
setClass("Set",
         representation = representation(Samples = 'character',
                                         Features = "character",
                                         dt = "data.table",
                                         'VIRTUAL'),
         validity = function (object) {
                 Sam <- unique(object@Samples)
                 len.Sam <- length(Sam)
                 Fea <- unique(object@Features)
                 len.Fea <- length(Fea)
                 nr <- base::nrow(object@dt)
                 nc <- base::ncol(object@dt)
                 msg <- c('')
                 valid <- TRUE
                 if (len.Sam < 1 | len.Fea < 1)  {
                         msg <- c(msg, "Samples and Features should not be NULL")
                         valid <- FALSE
                 }
                 if (nr != len.Fea | nc != len.Sam) {
                         msg <- c(msg, "Dimension of data is not consistent with the length of samples and features")
                         valid <- FALSE
                 }
                 if (valid) {
                         TRUE
                 } else {
                         msg
                 }
         }
         )

#' mSet
#'
#' `mSet` is a S4 class extended from the virtual `Set` used as
#' object to store metabolites abundance matrix.
#'
#' @slot Samples character a character vector contains the samples
#' @slot Features character a character vector contains the features
#' @slot dt data.table metabolties abundance matrix
#' @importClassesFrom data.table data.table
#' @export
#' @examples
#' m.path <- system.file("extdata", "metabolites_and_genera.rda", package = "mbOmic")
#' load(m.path)
#' names(metabolites)[1] <- 'rn'
#' bSet(b = metabolites)
#' @return S4 class
#' @docType class
setClass("mSet", slots = c(Samples = 'character',
                           Features = "character",
                           dt = "data.table"),
         contains = "Set")


#' bSet
#'
#' `bSet` class is similar to the `mSet` class but it store
#' the OTU abundance matrix rather than the metabolite abundance.
#'
#' @slot Samples character a character vector contains the samples
#' @slot Features character a character vector contains the features
#' @slot dt data.table OTU abundance matrix
#' @importClassesFrom data.table data.table
#' @export
#' @examples
#' b.path <- system.file("extdata", "metabolites_and_genera.rda", package = "mbOmic")
#' load(b.path)
#' names(genera)[1] <- 'rn'
#' bSet(b = genera)
#' @return S4 class
#' @docType class
setClass("bSet", slots = c(Samples = 'character',
                           Features = "character",
                           dt = "data.table"),
         contains = 'Set')

