## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  tidy = FALSE,
  warning = FALSE,
  message = FALSE
)


## ----setup--------------------------------------------------------------------
library(data.table)
library(mbOmic)

## -----------------------------------------------------------------------------
knitr::include_graphics('./img/mbOmic-workflow.svg')

## -----------------------------------------------------------------------------
path <- system.file('data',package = 'mbOmic')
load(file.path(path,'metabolites_and_genera.rda'))

## -----------------------------------------------------------------------------
mb <-
  mbSet(
      m = metabolites,
    b = genera
  )

## -----------------------------------------------------------------------------
samples.extra(mb)

## -----------------------------------------------------------------------------
nb <- nb.extra(mb)
nm <- nm.extra(mb)
cat("The mb object contains", nb, "generas and", nm, "metabolites\n")

## -----------------------------------------------------------------------------
mb <- clean_analytes(mb,m_thres = 2,b_thres = 2)

## ---- fig.width=8-------------------------------------------------------------
net <- try({
  wgcna(mb,message = FALSE,threshold.d = 0.02, threshold = 0.8)
})
class(net)

## ----fig.width=8--------------------------------------------------------------
net <- wgcna(mb,message = FALSE,threshold.d = 0.02, threshold = 0.8, power = 9)

## ---- fig.width=8-------------------------------------------------------------
plot_wgcna(net)

## -----------------------------------------------------------------------------
res <- corr(mb, method = 'spearman')
head(res)

## ---- fig.width=8, fig.height=8-----------------------------------------------
plot_network(net, res[abs(rho)>=0.85])

