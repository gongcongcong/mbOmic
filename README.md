# mbOmic

The `mbOmic` package contains a set of analysis functions for microbiomics and metabolomics data, designed to analyze the inter-omic correlation between microbiology and metabolites.

## Installation

You can install the released version of mbOmic from [Github](https://github.com/gongcongcong/mbOmic.git) with:

``` r
library(devtools)
install_github('gongcongcong/mbOmic')
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
ls()
#[1] "genera"      "metabolites" "path" 
```

### Construct `mbSet` object.

`mbSet` is S4 class containing the metabolites and OTU abundance.

We can use `mbSet` function to directly create `mbSet` class.

``` r
mb <-
  mbSet(
    m = metabolites,
    b = genera
  )
mb
# Samples( 12 ):  BS1 BS2 BS3 BS4 BS5 BS6 SS1 SS2 SS3 SS4 SS5 SS6 
# number of OTU:  18 
# number of metabolites:  247
```

There are some extract function to extract information of a `mbSet`, such as `mbSamples`, `b`, and `m`.

Extract the samples names from `mbSet` class by function `mbSamples`.

``` r
mbSamples(mb)
#[1] "BS1" "BS2" "BS3" "BS4" "BS5" "BS6" "SS1" "SS2" "SS3" "SS4" "SS5" "SS6"
```

``` r
bacteria <- b(mb)
metabolites <- m(mb)
cat("The mb object contains", nrow(bacteria), 
    "generas and", nrow(metabolites), "metabolites\n")
#The mb object contains 18 generas and 247 metabolites
```

### Remove bad analytes (OTU and metatoblites)

Removal of analytes only measured in \<2 of samples can perform by `clean_analytes`.

``` r
mb <- clean_analytes(mb,m_thres = 2,b_thres = 2)
```

### Generate metabolite module

`mbOmic` can generate metabolite module by `coExpress` function. The `coExpress` function is the encapsulation of one-step network construction and module detection of `WGCNA` package. The `coExpress` function firstly pick up the soft-threshold. The `threshold.d` and `threshold` parameters are used to detect whether is $R^2$ changing and appropriate.

``` r
net <- try({
  coExpress(mb,message = TRUE,threshold.d = 0.02, threshold = 0.8)
})
class(net)
```

If you can't get a good scale-free topology index no matter how high set the soft-thresholding power, you can directly set the power value by the parameter `power`. The appropriate soft-thresholding power can be chosen based on the number of samples as in the table below (recommend by `WGCNA` package).

| **Number of samples** | **Unsigned and signed hybrid networks** | **Signed networks** |
|:---------------------:|:---------------------------------------:|:-------------------:|
|         \<20          |                    9                    |         18          |
|        20\~30         |                    8                    |         16          |
|        30\~40         |                    7                    |         14          |
|         \>40          |                    6                    |         12          |

``` r
net <- coExpress(mb,message = TRUE,threshold.d = 0.02, threshold = 0.8, power = 9)

```

Result Visualization is performed by function `plot_coExpress`.

``` r
plot_coExpress(net)
```

### Calculate the Spearman metabolite-genera correlation

you can calculate the correlation between metabolites and OTUs by `corr` function. It return a data table containing `rho`, `p value`, and `adjust p value`.

``` r
res <- corr(mb, method = 'spearman')
head(res)
```

### plot the network

Finally, you can vaisulize the network by `plot_network` function, taking the `coExpress`and `corr` output. The orange nodes correspondes to OTU(genera)).

``` r
plot_network(net, res[abs(rho)>=0.85])
```
