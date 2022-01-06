# mbOmic

The `mbOmic` package contains a set of analysis functions for microbiomics and metabolomics data, designed to analyze the inter-omic correlation between microbiology and metabolites.

## Installation

You can install the released version of mbOmic from [Github](https://github.com/gongcongcong/mbOmic.git) with:

``` r
library(devtools)
install_github('gongcongcong/mbOmic')
```

## WorkFlow

The `mbOmic` package contains a set of analysis functions for microbiomics and metabolomics data, designed to analyze the inter-omic correlation between microbiology and metabolites, referencing the workflow of Jonathan Braun et al.[^1].

[^1]: McHardy, I. H., Goudarzi, M., Tong, M., Ruegger, P. M., Schwager, E., Weger, J. R., Graeber, T. G., Sonnenburg, J. L., Horvath, S., Huttenhower, C., McGovern, D. P., Fornace, A. J., Borneman, J., & Braun, J. (2013). Integrative analysis of the microbiome and metabolome of the human intestinal mucosal surface reveals exquisite inter-relationships. *Microbiome*, *1*(1), 17. <https://doi.org/10.1186/2049-2618-1-17>

![](https://raw.githubusercontent.com/gongcongcong/mbOmic/master/inst/doc/workflow.svg)

## Example

Load metabolites and OTU data of plant.[^2] The OTU had been binned into genera level and were save as the metabolites_and_genera.rda file

[^2]: Huang, W., Sun, D., Chen, L., & An, Y. (2021). Integrative analysis of the microbiome and metabolome in understanding the causes of sugarcane bitterness. Scientific Reports, 11(1), 1-11.

``` r
library(data.table)
library(mbOmic)
path <- system.file('extdata', 'metabolites_and_genera.rda',package = 'mbOmic')
load(path)
ls()
#[1] "genera"      "metabolites" "path" 
```

### Construct `mSet` and `bSet` object.

`mSet` and `bSet` are S4 classes containing the metabolites and OTU, respectively. Them could be created by `mSet` and `bSet` function, respectively. A `data.table` storing the metabolites or OTU matrix could be used to construct `mSet` or `bSet` class. If it contains `rn` column, the `rn` column will be used as featur e. Otherwise, a character vector must be given to `mSet` or `bSet` by the `Features` parameter. Meanwhile, if the parameter `Samples` is missing, the colnames of `data.table` excepting the `rn` will be used as names of samples.

``` r
names(genera)[1] <- 'rn'
names(metabolites)[1] <- 'rn'
b <- bSet(b = genera)
## or
## b <- bSet(b = genera[, !"rn"], Features = genera$rn, Samples = names(genera)[-1])
b
m <- mSet(m = metabolites)
## or
## m <- bSet(m = metabolites[, !"rn"], Features = metabolites, Samples = names(metabolites)[-1])
m
```

Some extract function to get information of a `mSet` and `bSet`, such as `samples`, `features`, `nrow`, and `ncol` could get the sample names, features, sample numbers, and features numbers.

Extract the samples names from `mbSet` class by function `mbSamples`.

``` r
samples(b) #Samples(m)
#[1] "BS1" "BS2" "BS3" "BS4" "BS5" "BS6" "SS1" "SS2" "SS3" "SS4" "SS5" "SS6"
```

What's more, names of samples or features can be set using `samples<-` or `features<-`.

``` r
Samples(b) <- your_samples_names_vector
```

### Remove bad analytes (OTU and metatoblites)

Removal of analytes (metabolites or OTU) only measured in \<2 of samples can perform by `clean_analytes`.

``` r
b <- clean_analytes(b, 2) 
m <- clean_analytes(m, 2)
```

### Generate metabolite module

`mbOmic` can generate metabolite module by `coExpress` function. The `coExpress` function is the encapsulation of one-step network construction and module detection of `WGCNA` package. The `coExpress` function firstly pick up the soft-threshold. The `threshold.d` and `threshold` parameters are used to detect whether is $R^2$ changing and appropriate.

``` r
net <- try({
  coExpress(m, message = TRUE,threshold.d = 0.02, threshold = 0.8)
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
net <- coExpress(m,message = TRUE,threshold.d = 0.02, threshold = 0.8, power = 9)
# 0 features removed:   
# 
# One: detect the power for softThreshold
# 
# using the power:  9 to constructe net!
# 
# Two: Network construction and module detection was done
# ====> There are  6 modules were constructed: 
# ====||  blue 47 
# ====||  brown 46 
# ====||  green 27 
# ====||  red 22 
# ====||  turquoise 70 
# ====||  yellow 35 
```

Result Visualization is performed by function `plot_coExpress`.

``` r
plot_coExpress(net)
```
![](https://raw.githubusercontent.com/gongcongcong/mbOmic/master/inst/doc/cluster_dendrogram.svg)

### Calculate the Spearman metabolite-genera correlation

you can calculate the correlation between metabolites and OTUs by `corr` function. It return a data table containing `rho`, `p value`, and `adjust p value`.

``` r
corr_spearman <- corr(m, b, method = 'spearman')
head(corr_spearman)
#                         b           m        rho           p       padj
# 1:            Acidibacter Delavirdine  0.6853147 0.013905969 0.03269484
# 2:           Acidothermus Delavirdine  0.7272727 0.007355029 0.02440333
# 3:       Anaeromyxobacter Delavirdine -0.7762238 0.002992864 0.01657070
# 4: Candidatus_Udaeobacter Delavirdine -0.6993007 0.011374199 0.02997610
# 5:           Chujaibacter Delavirdine  0.6223776 0.030675895 0.05060669
# 6:           Conexibacter Delavirdine  0.5804196 0.047855977 0.06912530
```

### plot the network

Finally, you can vaisulize the network by `plot_network` function, taking the `coExpress`and `corr` output. The orange nodes correspondes to OTU(genera).

``` r
plot_network(net, corr_spearman[abs(rho)>=0.85])
```
![](https://raw.githubusercontent.com/gongcongcong/mbOmic/master/inst/doc/network.svg)

### identification of enterotype

Construct `bSet` class using the OTU abundance matrix in genera level.

``` r
dat <- read.delim('http://enterotypes.org/ref_samples_abundance_MetaHIT.txt')
dat <- impute::impute.knn(as.matrix(dat), k = 100)
dat <- as.data.frame(dat$data+0.001) 
setDT(dat, keep.rownames = TRUE)
dat <- bSet(b =  dat)
dat

# 1. Features( 278 ): 
# 	 Bacteroides Prevotella Eubacterium Faecalibacterium Alistipes ...
# 2. Samples( 51 ): 
# 	 MH0277 MH0087 MH0156 MH0444 MH0333 ...
# 3. Top 5 Samples data:
# [1] 1 2 3 4 5
```


Then estimating the approate numbers of cluster can implement by `estimate_k` function.
``` r
res2 <- estimate_k(dat)
res2
# optimal number of cluster:  4 
# Max CHI:  164.6422 
# Silhouette:  0.1814455 
```
![](https://raw.githubusercontent.com/gongcongcong/mbOmic/master/inst/doc/verCHI.svg)

Enterotype of samples validates the result of `estimate_k`.

``` r
rest <- read.table(system.file('extdata', 'enterotype.txt', package = 'mbOmic'))
rest <- rest[samples(dat),]
table(res2$verOptCluster, rest$ET)
#   ET_B ET_F ET_P
# 1    0   21   19
# 2   67    5    0
# 3    0    0   40
# 4    3  123    0
```

`enterotyping` function can estimate the enterotype using the `bSet` class.

``` r
ret=enterotyping(dat, res2$verOptCluster) 
ret

# $enterotypes
#      Enterotype        max which   cluster
# 1:  Bacteroides 0.36724946     2 cluster 2
# 2:   Prevotella 0.29692944     3 cluster 3
# 3: Ruminococcus 0.02416713     4 cluster 4
# 
# $data
#      Samples   Enterotype   cluster
#   1:  MH0087  Bacteroides cluster 2
#   2:  MH0156  Bacteroides cluster 2
#   3:  MH0444  Bacteroides cluster 2
#   4:  MH0333  Bacteroides cluster 2
#   5:  MH0233  Bacteroides cluster 2
#  ---                               
# 234:  MH0012 Ruminococcus cluster 4
# 235:  MH0415 Ruminococcus cluster 4
# 236:  MH0457 Ruminococcus cluster 4
# 237:  MH0442 Ruminococcus cluster 4
# 238:  MH0448 Ruminococcus cluster 4
# 
# $UnIdentifiedSamples
#  [1] "MH0277" "MH0161" "MH0046" "MH0175" "MH0152" "MH0104" "MH0151" "MH0189"
#  [9] "MH0030" "MH0157" "MH0063" "MH0075" "MH0141" "MH0169" "MH0050" "MH0286"
# [17] "MH0096" "MH0053" "MH0217" "MH0098" "MH0009" "MH0197" "MH0065" "MH0173"
# [25] "MH0168" "MH0070" "MH0077" "MH0288" "MH0200" "MH0031" "MH0183" "MH0132"
# [33] "MH0144" "MH0124" "MH0430" "MH0276" "MH0407" "MH0428" "MH0126" "MH0447"
```
