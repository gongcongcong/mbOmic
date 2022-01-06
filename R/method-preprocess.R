#' clean_analytes
#'
#' Remove unqualified metabolites and OTU which only
#' measured in few samples.
#'
#' @param object a `bSet` or `mSet` class
#' @param fea_num, integer, Removal of features only measured in `fea_num` of samples
#' @import data.table
#' @export
#' @return return a quantified `bSet` or `mSet` object
#' @examples
#' library(data.table)
#' path <- system.file('extdata', 'metabolites_and_genera.rda', package = 'mbOmic')
#' load(path)
#' names(genera)[1] <- 'rn'
#' names(metabolites)[1] <- 'rn'
#' b <- bSet(b = genera)
#' m <- mSet(m = metabolites)
#' clean_analytes(b)
clean_analytes <-
        function(object, fea_num = 2){
                stopifnot("object should be the mSet or bSet class" = inherits(object, c("mSet", "bSet")))
                samples <- samples(object)
                dt <- object@dt
                index <- apply(dt,
                               1,
                               function(x) sum(x != 0) > fea_num)
                message(sum(!index), "features removed: ", features(object)[!index], '\n')

                object[index]
        }
