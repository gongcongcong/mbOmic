#'
#' vmh
#'
#' `vmh` is a API for vmh
#'
#' @slot type a character, "microbe", "metabolites", "gene", or "food items"
#' @slot net a data.table class contain the network of microbe, metabolites, gene, and food items
#' @importClassesFrom data.table data.table
#' @export
#' @examples
#' b.path <- system.file("extdata", "metabolites_and_genera.rda", package = "mbOmic")
#' load(b.path)
#' names(genera)[1] <- 'rn'
#' bSet(b = genera)
#' @return vmh class
#' @docType class
setClass("vmh", slots = c(type = 'character',
                          net = "data.table"),
         validity = function(object) {
                 errors = character()
                 types = c("microbe", "metabolites", "gene", "food items")
                 if (!identical(names(object@net), types)) {
                         errors = c(errors, "net shape wrong!")
                 }
                 if (type %in% types) {
                         errors = c(errors, "type wrong!")
                 }
                 if (length(errors)==0) {
                         TRUE
                 } else {
                         errors
                 }
         }
         )

