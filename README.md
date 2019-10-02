
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RILStEp

<!-- badges: start -->

<!-- badges: end -->

RILStEp(Recombinant Inbred Lines Stepwise Epistasis analysis) package is
the epistasis analysis tool. This package enables to detect the
epistatic relationships between SNPs for RIL population by comparison of
2 models based on bayes factor.<br> first model : only QTL effects.
second model : QTL effects & epistasis effects.

## Installation

You can install RILStEp from
[GitHub](https://github.com/slt666666/RILStEp) with:

``` r
# install.packages("devtools")
devtools::install_github("slt666666/RILStEp")
```

## Example

This is a basic example script:

``` r
library(RILStEp)

### loading dataset
loaded_data <- load_data("phenotype.csv", "genotype.csv", "trait1")

### check all combinations of 2 SNPs
result1 <- rilstep(loaded_data, "result1", core_num = 8)

### using 1 SNP in each 500 SNPs
result2 <- rilstep(loaded_data, "result2", interval = 500)

### specify QTL-like SNPs by user.
result3 <- rilstep(loaded_data, "result3", qtls = c("chr08_19928351", "chr09_3909046"))

### using SNPs in specific regions
result4 <- rilstep(loaded_data, "result4", regions = c("chr03_2132221:chr10_9330401", "chr03_2132221:chr10_9330401"))
```

option `region=c()` can specify regions. rilstep check all combinations
of SNPs in specified
    regions.

    regions = c("chr03_1234")                                    # chr03_1234 x all SNPs
    regions = c("chr03_1234:chr10_5678")                         # chr03_1234 ~ chr10_5678 x all SNPs
    regions = c("chr03_1234", "chr10_5678")                      # chr03_1234 x chr10_5678
    regions = c("chr03_1234:chr10_5678", "chr04_9012:chr7_3456") # chr03_1234~chr10_5678 x chr04_9012~chr7_3456
