## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  tidy = FALSE,
  warning = FALSE,
  message = FALSE
)


## ----setup--------------------------------------------------------------------
library(mbOmic)

## -----------------------------------------------------------------------------
# knitr::include_graphics(system.file('extdata', 'intro.png', 'mbOmic'))

## -----------------------------------------------------------------------------
path <- system.file('extdata', 'metabolites_and_genera.rda', package = 'mbOmic')

load(path)

## -----------------------------------------------------------------------------
names(metabolites)[1] <- 'rn'
m <- mSet(m = metabolites)
m

## -----------------------------------------------------------------------------
samples(m)
#[1] "BS1" "BS2" "BS3" "BS4" "BS5" "BS6" "SS1" "SS2" "SS3" "SS4" "SS5" "SS6"

## -----------------------------------------------------------------------------
m <- clean_analytes(m, fea_num = 2)

## ---- fig.width= 7, fig.height= 6---------------------------------------------
net <- try({
  coExpress(m,message = TRUE,threshold.d = 0.02, threshold = 0.8, plot = TRUE)
})
class(net)

## -----------------------------------------------------------------------------
net <- coExpress(m,message = TRUE,threshold.d = 0.02, threshold = 0.8, power = 9)

## ---- fig.align='center', fig.width= 7, fig.height= 6-------------------------
plot_coExpress(net)

## -----------------------------------------------------------------------------
b <- genera
names(b)[1] <- 'rn'
b <- bSet(b=b)
spearm <- corr(m = m, b = b, method = 'spearman')
# head(spearm)
spearm[p<=0.001]

## ---- fig.align='center', fig.width= 7, fig.height= 6-------------------------
plot_network(net, spearm[abs(rho) >= 0.75 & p <= 0.001], threshold = 0.75, show_text = FALSE)

## ---- fig.align='center', fig.width= 7, fig.height= 6-------------------------
plot_network(net, spearm[abs(rho) >= 0.75 & padj <= 0.05], threshold = 0.75)

## -----------------------------------------------------------------------------
devtools::session_info()

