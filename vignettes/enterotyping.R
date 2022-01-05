## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(mbOmic)
library(data.table)

## -----------------------------------------------------------------------------
dat <- read.delim('http://enterotypes.org/ref_samples_abundance_MetaHIT.txt')
dat <- impute::impute.knn(as.matrix(dat), k = 100)
dat <- as.data.frame(dat$data+0.001) 
setDT(dat, keep.rownames = TRUE)
dat <- bSet(b =  dat)
rest <- read.table(system.file('extdata', 'enterotype.txt', package = 'mbOmic'))
rest <- rest[samples(dat),]
res2 <- estimate_k(dat)

## -----------------------------------------------------------------------------
table(res2$verOptCluster, rest$ET)

## -----------------------------------------------------------------------------
ret=enterotyping(dat, 
                 res2$verOptCluster) 
ret

## -----------------------------------------------------------------------------
devtools::session_info()

