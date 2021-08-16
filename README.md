# mbOmic

<!-- badges: start -->

<!-- badges: end -->

The `mbOmic` package contains a set of analysis functions for microbiomics and metabolomics data, designed to analyze the inter-omic correlation between microbiology and metabolites.

## Installation

You can install the released version of mbOmic from [github](https://github.com/gongcongcong/mbOmic.git) with:

``` r
library(devtools)
install_github('gongcongcong/mbSet')
```

## WorkFlow

The `mbOmic` package contains a set of analysis functions for microbiomics and metabolomics data, designed to analyze the inter-omic correlation between microbiology and metabolites, referencing the workflow of Jonathan Braun et al[^McHardy].

[^McHardy]: McHardy, I. H., Goudarzi, M., Tong, M., Ruegger, P. M., Schwager, E., Weger, J. R., Graeber, T. G., Sonnenburg, J. L., Horvath, S., Huttenhower, C., McGovern, D. P., Fornace, A. J., Borneman, J., & Braun, J. (2013). Integrative analysis of the microbiome and metabolome of the human intestinal mucosal surface reveals exquisite inter-relationships. *Microbiome*, *1*(1), 17. <https://doi.org/10.1186/2049-2618-1-17>

![mbSetWF](https://github.com/gongcongcong/mbOmic/blob/master/vignettes/img/mbOmic-workflow.svg)

## Example

Load metabolites and OTU abundance data of plant.[^Huang] The OTU had been binned into genera level and were save as the metabolites_and_genera.rda file

[^Huang]: Huang, W., Sun, D., Chen, L., & An, Y. (2021). Integrative analysis of the microbiome and metabolome in understanding the causes of sugarcane bitterness. Scientific Reports, 11(1), 1-11.

``` r
path <- system.file('data',package = 'mbOmic')
load(file.path(path,'metabolites_and_genera.rda'))
```

### Construct the `mbSet` object.

`mbSet` is S4 class containing the metabolites and OTU abundance.

We can use `mbSet` function to create directrly `mbSet` class.

``` r
mb <-
  mbSet(
      m = metabolites,
    b = genera
  )
```

Extract the samples names from `mbSet` class by `samples.extra` function.

``` r
samples.extra(mb)
```

``` r
nb <- nb.extra(mb)
nm <- nm.extra(mb)
cat("The mb object contains", nb, "generas and", nm, "metabolites\n")
```

### Remove bad analytes (OTU and metatoblites)

Removal of analytes only measured in \<2 of samples can perform by `clean_analytes`.

``` r
mb <- clean_analytes(mb,m_thres = 2,b_thres = 2)
```

### Generate metabolite module

The WGCNA package is used to generate metabolite modules. The first step is to pick up the softthreshold. This processs is enclosed into `wgcna` function of `mbOmic`. The `threshold.d` and `threshold` parameters are used to detect whether is $R^2$ changing and appropriate.

``` r
net <- try({
  wgcna(mb,message = FALSE,threshold.d = 0.02, threshold = 0.8)
})
class(net)
```

If you can't get a good scale-free topology index no matter how high set the soft-thresholding power, you can set directly the power by the parameter `power`. The appropriate soft-thresholding power can be chosen based on the number of samples as in the table below (recommend by the authors of `WGCNA` package).

| **Number of samples** | **Unsigned and signed hybrid networks** | **Signed networks** |
|:---------------------:|:---------------------------------------:|:-------------------:|
|         \<20          |                    9                    |         18          |
|        20\~30         |                    8                    |         16          |
|        30\~40         |                    7                    |         14          |
|         \>40          |                    6                    |         12          |

``` r
net <- wgcna(mb,message = FALSE,threshold.d = 0.02, threshold = 0.8, power = 9)
```

Result Visualization is performed by function `plotDendroAndColors` of `WGCNA` package.

``` r
plot_wgcna(net)
```

### Calculate the Spearman metabolite-genera correlation

you can calculate the Spearman correlation between metabolites and OTUs by `corr` function. It return a data table containing `rho`, `p value`, and `adjust p value`.

``` r
res <- corr(mb, method = 'spearman')
head(res)
```

### plot the network

Finally, you can vaisulize the network by `plot_network` function, taking the `wgcna` and `corr` output. The orange nodes correspondes to OTU(genera)).

``` r
plot_network(net, res[abs(rho)>=0.85])
```
