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


#' @title mb_initialize
#'
#' initialized function for `mSet` and `bSet` classs.
#'
#' @param type character `mSet` or `bSet`
#' @param dt data.table abundance matrix
#' @param Samples character samples names
#' @param Features character Features names
#' @return object
#' @import data.table
#'
mb_initialize <-
        function(type, dt, Samples, Features) {
                stopifnot("dt should be data.table" = inherits(dt, 'data.table'))
                if (missing(Features)) {
                        stopifnot("Could not find the features column!" = 'rn' %in% names(dt))
                        Features <- dt$rn
                        dt <- dt[, !"rn"]
                } else {
                        stopifnot("Features should be character" = inherits(Features, 'character'))
                }
                if (missing(Samples)) {
                        Samples <- setdiff(names(dt), 'rn')
                } else {
                        stopifnot("Samples should be data.table" = inherits(Samples, 'character'))
                }
                .Object <- new(type, dt = dt, Samples = Samples, Features = Features)
                validObject(.Object)
                .Object
        }
#' mSet
#'
#' Function to return `mSet` class.
#'
#' @export
#' @return a `mSet` class
#' @param m data.table, metabolites aundance matrix. if `rn` column is not contained
#'  in this data.table, the `Features` parameter should be given by character vector.
#' @param ... `Samples` or `Features`. if the `Samples` not given, the colnames of
#' `m` data.table will be taken as the Samples names.
#' @importFrom methods new
#' @examples
#' path <- system.file('extdata', 'metabolites_and_genera.rda', package = 'mbOmic')
#' load(path)
#' names(metabolites)[1] <- 'rn'
#' mSet(m = metabolites)
mSet <- function(m, ...){
        mb_initialize("mSet", dt = m, ...)
}
#' bSet
#'
#' Function to return `bSet` class.
#'
#' @export
#' @return a `bSet` class
#' @param b data.table, metabolites aundance matrix. if `rn` column is not contained
#'  in this data.table, the `Features` parameter should be given by character vector.
#' @param ... `Samples` or `Features`. if the `Samples` not given, the colnames of
#' `b` data.table will be taken as the Samples names.
#' @importFrom methods new
#' @examples
#' path <- system.file('extdata', 'metabolites_and_genera.rda', package = 'mbOmic')
#' load(path)
#' names(genera)[1] <- 'rn'
#' bSet(b = genera)
bSet <- function(b, ...){
        mb_initialize("bSet", dt = b, ...)
}

#' nrow
#'
#' nrow method for `Set` to obtain the number of features.
#' nrow and ncol return the number of features or samples of a `Set` class.
#'
#' @param x a `Set` Object
#' @export
#' @return an `integer` of length 1
setMethod('nrow', 'Set', function(x){
        return(length(x@Features))
})
#' ncol
#'
#' ncol method for `Set` to obtain the number of samples.
#' nrow and ncol return the number of features or samples of a `Set` class.
#'
#' @param x a `Set` Object
#' @export
#' @return an `integer` of length 1
setMethod('ncol', 'Set', function(x){
        return(length(x@Samples))
})

#' features
#' @rdname features
#'
#' @description `features` get or set features names vector.
#' @usage
#' features(object)
#' features(object) <- value
#' @param object a `Set` object
#' @export
#' @importFrom methods validObject
#' @return `features` return a character vector or update the features names
#' @examples
#' library(data.table)
#' path <- system.file('extdata', 'metabolites_and_genera.rda', package = 'mbOmic')
#' load(path)
#' names(metabolites)[1] <- 'rn'
#' m <- mSet(m = metabolites)
#' features(m)
features <- function(object){
        object@Features
}

#' features
#' @rdname features
#' @description `features` get or set features names vector.
#' @usage
#' features(object)
#' features(object) <- value
#' @export
#' @return `features` return a character vector or update the features names
#' @param object a `Set` object
#' @param value character vector
#' @importFrom methods validObject
"features<-" <- function(object, value){
        object@Features <- value
        validObject(object)
        object
}

#' samples
#' @rdname samples
#'
#' @description `samples` get or set samples names vector.
#' @usage
#' samples(object)
#' samples(object) <- value
#' @param object a `Set` object
#' @export
#' @return `samples` return a character vector or update the samples names
#' @examples
#' library(data.table)
#' path <- system.file('extdata', 'metabolites_and_genera.rda', package = 'mbOmic')
#' load(path)
#' names(metabolites)[1] <- 'rn'
#' m <- mSet(m = metabolites)
#' features(m)
samples <- function(object){
        object@Samples
}
#' samples
#' @rdname samples
#' @description `samples` get or set samples names vector.
#' @usage
#' samples(object)
#' samples(object) <- value
#' @export
#' @return `samples` return a character vector or update the samples names
#' @param object a `Set` object
#' @param value character vector
"samples<-" <- function(object, value){
        object@Samples <- value
        validObject(object)
        object
}
#' show
#'
#' show method for the `Set` class.
#'
#' @param object a `Set` class
#' @export
#' @return none return but print the main info of `Set` class
#' @importFrom methods show
setMethod('show', 'Set', function(object){
        nc <- ncol(object)
        nr <- nrow(object)
        if (nr > 5) {
                cat("1. Features(", nc, "): \n\t", features(object)[seq_len(5)], '...\n')
        } else {
                cat("1. Features(", nc, "): \n\t", features(object), '\n')
        }
        if (nc > 5) {
                cat("2. Samples(", nr, "): \n\t", samples(object)[seq_len(5)], '...\n' )
                cat("3. Top 5 Samples data:\n")
                print(object@dt[, seq_len(5)])
        } else {
                cat("2. Samples(", nr, "): \n\t", samples(object), '\n')
                cat("3. data:\n")
                print(object@dt)
        }
})


#' subset_Set
#'
#' subset `Set` class.
#'
#' @param x a `mSet` or `bSet` class
#' @param i the data.table i
#' @param j the data.table j
#' @return return a subset object
#' @export
#' @importFrom methods new
"[.Set" <- function(x, i, j){
        f <- features(x)[i]
        if(missing(j)) {
                dt <- x@dt[i]
        } else {
                dt <- x@dt[i, j, with = FALSE]
        }
        s <- names(dt)
        new(class(x)[1], dt = dt, Features = f, Samples = s)
}

