mvNUR Workflow 2: mvBIMBAM
================
Lacey W. Heinsberg



# Copyright

Copyright 2023, University of Pittsburgh. All Rights Reserved.  
License: CC BY-SA 3.0
([link](https://creativecommons.org/licenses/by-sa/3.0/))

# Overview

Please note that this code was adapted from previous work led by Jerry
Zhang, Lacey Heinsberg, and Daniel Weeks:

Zhang JZ, Heinsberg LW, Krishnan M, Hawley NL, Major TJ, Carlson JC,
Harré Hindmarsh J, Watson H, Qasim M, Stamp LK, Dalbeth N, Murphy R, Sun
G, Cheng H, Naseri T, Reupena MS, Kershaw EE, Deka R, McGarvey ST,
Minster RL, Merriman TR, Weeks DE. Multivariate analysis of a missense
variant in CREBRF reveals associations with measures of adiposity in
people of Polynesian ancestries. Genet Epidemiol. 2023
Feb;47(1):105-118. doi: 10.1002/gepi.22508. Epub 2022 Nov 9. PMID:
36352773; PMCID: PMC9892232.

[GitHub Repository](https://github.com/lwheinsberg/mvCREBRF)

Note this code is dependent upon `00_mvNUR_DataPrep.Rmd` which creates
the mvBIMBAM input files detailed below.

## mvBIMBAM

mvBIMBAM is a command line program for multivariate genetic association
analysis of multiple related phenotypes. In the mvBIMBAM framework, a
global null model representing no association between phenotypes and
genotypes is compared with an exhaustive combination of alternative
models in which all different combinations of phenotype-genotype
associations are considered (see `mph2` details below). For the
alternative models, the mvBIMBAM methodology splits phenotype-genotype
associations into all possible partitions of U, D, and I, each
representing ‘unassociated’, ‘directly’, and ‘indirectly’ associated.
Details about the analyses and interpretation can be found below.
Installation and troubleshooting instructions can be found in the README
file.

# Load libraries

Load the libraries needed to run the code.

``` r
library(tidyverse)
select = dplyr::select
library(ggplot2)
# The 'preprocessCore' is a Bioconductor package
library(preprocessCore)
library(pander)
library(stringr)
```

# Data

Most of our data (created in `00_mvNUR_DataPrep.Rmd`) are simply called
by the system (see below) and not formerly loaded into R. Here though,
we will store a list of phenotypes and variants of interest. See
README.md for information about phenotypes and SNPs of interest.

``` r
# Define the phenotypes of interest
traits <- c("EMO_tscore", "bdito", "FAT_tscore", "paohcif", "EPSscore", "pain")
trait_mapping <- c("Anxiety", "Depression", "Fatigue", "Cognitive function", "Sleepiness", "Pain")

# Define variants of interest 
genes <- c("rs4880", "rs5746136", "rs1041740", "rs10432782", "rs4135225", "rs7522705")
```

# mvBIMBAM analyses

We will now apply mvBIMBAM to the synthetic data set we are working
with. Rather than navigating the analyses in the terminal, we can make
our code reproducible in R and call mvBIMBAM through the R using
`call system()`. More information about the program can be found under
the original [mvBIMBAM
documentation](https://github.com/heejungshim/mvBIMBAM/blob/master/doc/mvBIMBAM_Manual.pdf)
on GitHub. In brief, the command we need contains the following
parameters:

-g: location of the genotype data (for our example the location is:
./inputs/synthetic_geno_bimbam.txt) -p: location of the phenotype data
(for our example the location is: ./inputs/synthetic_pheno_bimbam.txt)
-o: prefix for output files (for our example we will use: bimbam_out)
-f: number of phenotypes to analyze (length(traits), for our example
this number is 6) -A: indicate multiple values for sigma (see mvBIMBAM
manual for more details; for our example we will use a range of values
from 0.05 to 0.4) -mph1 vs. mph2: analysis settings (see mvBIMBMAM
manual for more details and our note below; for this example we will run
mph2)

A note about `mph1` vs. `mph2`: In brief, `-mph1` performs a simple
multivariate test of the null hypothesis vs a general multivariate
alternative. Because it does not consider all partitions of the
alternative hypothesis, it is quite fast, can handle a moderately large
number of phenotypes, and be run on a whole-genome scale. However, it
was originally intended to be used only as needed in an initial
filtering step before running the more computationally-intensive
`-mph 2` analyses. Because of the way `-mph1` functions, it does not
allow you to ask interesting questions such as which phenotypes are
associated with each SNP. To do this more detailed analysis use the
`-mph 2` option. Note that if you have only a handful of phenotypes, you
may be able to run the `-mph 2` analysis genome wide on all variants and
skip `-mph1` all together.

A quick warning from the creators of mvBIMBAM: the prior used in
mvBIMBAM is based on the idea that if a variant is associated with one
phenotype, then it is likely associated with multiple phenotypes. That
is, it is relatively permissive of associations with multiple
phenotypes, and does not attempt to be “skeptical” of additional
associations. In this sense it is more suited to hypothesis generation
than of hypothesis testing. The prior can be changed to be more
skeptical, but it requires that you dig into the last output file,
below. This adjustment is beyond the scope of this guide.

    ## [1] "bimbam -g ./inputs/synthetic_geno_bimbam.txt -p ./inputs/synthetic_pheno_bimbam.txt -o bimbam_out -mph 2 -f 6 -A 0.05 -A 0.1 -A 0.2 -A 0.4"

    ## [1] "-bimbam: file 0 has 763 individual and 6 snps"   
    ## [2] "-bimbam: read file 0 again "                     
    ## [3] "-bimbam: number of phenotypes = 6"               
    ## [4] "total = 729"                                     
    ## [5] "output/bimbam_out.mph.txt has been created."     
    ## [6] "output/bimbam_out.mph.prob.txt has been created."
    ## [7] "output/bimbam_out.mph.BFs.txt has been created." 
    ## attr(,"status")
    ## [1] 1

    bimbam -g ./inputs/synthetic_geno_bimbam.txt -p ./inputs/synthetic_pheno_bimbam.txt 
    -o bimbam_out -mph 2 -f 6 -A 0.05 -A 0.1 -A 0.2 -A 0.4

Once we run the chunk above, we can then find our output at
./output/bimbam_out.mph.BFs.txt.

# Results

## Bayes factors (BF)

The evidence against the null hypothesis is the sum of Bayes factors
(BF) (log10 scale) of all partitions weighted by a diffuse prior. Let’s
read in our BF factor results file to take a look.

``` r
# Initialize an empty matrix to store the values
s <- readLines("./output/bimbam_out.mph.BFs.txt")
m <- matrix(nrow = length(s), ncol = length(traits)+2)
snp_names <- matrix(nrow = length(s), ncol = 1)

# Process each line and clean up result matrix
for (i in 1:length(s)) {
  # Split the line by spaces
  split_values <- strsplit(s[i], "\\s+")[[1]]
  # Convert the remaining values to numeric format
  numeric_values <- as.numeric(split_values[-1])
  # Store the numeric values in the matrix
  m[i, 1:(length(numeric_values))] <- numeric_values
  snp_names[i] <- split_values[1]
}
rownames(m) <- snp_names
# Store results in m1 object and update column names
m1 <- m
colnames(m1) <- c("BF", "BF_All_Partition", trait_mapping)
```

``` r
# View results
pander(m1, digits = 4, caption = "Bayes Factors for Synthetic Example Data Set")
```

|                |  BF   | BF_All_Partition | Anxiety | Depression | Fatigue |
|:--------------:|:-----:|:----------------:|:-------:|:----------:|:-------:|
|   **rs4880**   | 11.64 |        12        | -0.2122 |   0.9347   | -0.405  |
| **rs5746136**  | 8.751 |      9.662       | 0.1084  |  -0.3045   | 0.6435  |
| **rs1041740**  | 13.69 |       14.4       |  12.52  |   0.6702   |  4.629  |
| **rs10432782** | 4.878 |      5.472       |  1.619  |    1.85    | -0.3186 |
| **rs4135225**  | 32.48 |      33.08       | 0.6325  |   1.603    |  25.11  |
| **rs7522705**  | 1.672 |      2.006       | -0.2782 |  -0.2978   | -0.2944 |

Bayes Factors for Synthetic Example Data Set (continued below)

|                | Cognitive function | Sleepiness |  Pain   |
|:--------------:|:------------------:|:----------:|:-------:|
|   **rs4880**   |       2.041        |  0.02955   |  3.437  |
| **rs5746136**  |       0.2681       |   6.421    | 0.6853  |
| **rs1041740**  |       6.156        |  0.00889   | -0.2024 |
| **rs10432782** |      0.05509       |  -0.05979  | -0.1757 |
| **rs4135225**  |       5.586        |   4.684    |  11.63  |
| **rs7522705**  |       0.5971       |  -0.3565   | 0.2123  |

Interpretation: The above table presents the overall log10 BF for the
network, and individual log10 BFs for each trait. Strong evidence of
association is defined as log10 BF \> 5; suggestive evidence is defined
as 1.5 \< log10 BF \< 5; and negligible evidence is defined as log10 BF
\< 1.5.

In this synthetic data set, we see that all SNPs (rs7522705 less so and
rs4135225 more so) appear to play a strong role in this network. rs4880
has suggestive associations with cognitive function and pain; rs5746136
has significant associations with sleepiness; rs1041740 has significant
associations with anxiety and cognitive function and a suggestive
association with fatigue, etc.

## Bayesian posterior probabilities of association

In addition to log10 BFs, probabilities for no association, direct
association, and indirect association are given as output while marginal
posterior probabilities of association (MPPA) are calculated by summing
the marginal posterior probabilities of direct and indirect association.

``` r
# Read in second result file
s <- readLines("./output/bimbam_out.mph.prob.txt")

# Initialize an empty list to store the matrix for each rs#
matrix_list <- vector("list", length(s))

# Process each line
for (i in 1:length(s)) {
  # Split the line by spaces
  split_values <- strsplit(s[i], "\\s+")[[1]]
  # Store the rs# value
  rs <- split_values[1]
  
  m2 <- matrix(na.omit(as.numeric(str_split(s, " ")[[i]])), nrow = 2)
  m2 <- rbind(m2, 1 - (m2[1,] + m2[2,]))
  
  # Set the column and row names
  colnames(m2) <- trait_mapping
  rownames(m2) <- c(paste0(rs, "_Unassociated"), paste0(rs, "_Directly"), paste0(rs, "_Indirectly"))

  # Store the matrix in the list, using rs# as the name
  matrix_list[[i]] <- list(rs = rs, m2 = m2)
}

# Combine all matrices into a single matrix
m2 <- do.call(rbind, lapply(matrix_list, function(x) x$m))
for (i in 1:ncol(m2)) {
  m2[, i] <- as.numeric(gsub("\"", "", m2[, i]))
}
```

``` r
pander(m2, digits = 4, caption = "Bayesian Probabilities for Synthetic Example Data Set")
```

|                             | Anxiety | Depression | Fatigue |
|:---------------------------:|:-------:|:----------:|:-------:|
|   **rs4880_Unassociated**   | 0.3751  |   0.0881   | 0.3875  |
|     **rs4880_Directly**     | 0.3415  |   0.9115   | 0.2468  |
|    **rs4880_Indirectly**    | 0.2834  |  0.00035   | 0.3657  |
| **rs5746136_Unassociated**  | 0.1947  |   0.3647   | 0.07658 |
|   **rs5746136_Directly**    | 0.8037  |   0.5952   | 0.9207  |
|  **rs5746136_Indirectly**   | 0.00154 |  0.04015   | 0.00272 |
| **rs1041740_Unassociated**  |    0    |  0.07304   |    0    |
|   **rs1041740_Directly**    |    1    |   0.9007   | 0.4676  |
|  **rs1041740_Indirectly**   |    0    |  0.02628   | 0.5324  |
| **rs10432782_Unassociated** | 0.00504 |  0.00493   | 0.3862  |
|   **rs10432782_Directly**   | 0.9831  |   0.9789   | 0.3531  |
|  **rs10432782_Indirectly**  | 0.01188 |  0.01618   | 0.2607  |
| **rs4135225_Unassociated**  | 0.08721 |  0.01237   |    0    |
|   **rs4135225_Directly**    | 0.9128  |   0.532    |    1    |
|  **rs4135225_Indirectly**   |    0    |   0.4556   |    0    |
| **rs7522705_Unassociated**  | 0.4168  |    0.54    | 0.4684  |
|   **rs7522705_Directly**    |  0.458  |   0.4192   | 0.5191  |
|  **rs7522705_Indirectly**   | 0.1253  |  0.04077   | 0.01252 |

Bayesian Probabilities for Synthetic Example Data Set (continued below)

|                             | Cognitive function | Sleepiness |  Pain   |
|:---------------------------:|:------------------:|:----------:|:-------:|
|   **rs4880_Unassociated**   |      0.00077       |   0.2343   | 0.00035 |
|     **rs4880_Directly**     |       0.9992       |   0.358    | 0.9997  |
|    **rs4880_Indirectly**    |         0          |   0.4077   |    0    |
| **rs5746136_Unassociated**  |       0.1642       |     0      | 0.05097 |
|   **rs5746136_Directly**    |       0.8279       |   0.9999   | 0.5256  |
|  **rs5746136_Indirectly**   |      0.00791       |   9e-05    | 0.4234  |
| **rs1041740_Unassociated**  |         0          |   0.2551   | 0.2203  |
|   **rs1041740_Directly**    |       0.9235       |   0.7228   | 0.3356  |
|  **rs1041740_Indirectly**   |       0.0765       |  0.02212   | 0.4441  |
| **rs10432782_Unassociated** |       0.3096       |   0.2886   | 0.3648  |
|   **rs10432782_Directly**   |       0.6901       |   0.4881   | 0.4547  |
|  **rs10432782_Indirectly**  |      0.00032       |   0.2233   | 0.1805  |
| **rs4135225_Unassociated**  |         0          |     0      |    0    |
|   **rs4135225_Directly**    |       0.361        |   0.3242   | 0.9953  |
|  **rs4135225_Indirectly**   |       0.639        |   0.6758   | 0.00466 |
| **rs7522705_Unassociated**  |       0.1034       |   0.5502   | 0.3012  |
|   **rs7522705_Directly**    |       0.8893       |   0.272    | 0.6018  |
|  **rs7522705_Indirectly**   |      0.00732       |   0.1778   | 0.09697 |

Interpretation: The numbers above can be interpreted as the probability
of no association, direct association, or indirect association, which
together sum to 1 (i.e., 100%). Both directly and indirectly associated
phenotypes are associated with a given genotype, but indirectly
associated phenotypes are conditionally independent of the genotype
given the presence of a directly associated phenotype in the model.

For example, in this synthetic data set, there is a suggested 92.4%
probability that there is a direct association between rs1041740 and
anxiety. In addition, while there is a 0% chance of no association
between rs4135225 and sleep, the likelihood of direct vs. indirect
effects is more split at 32.4% and 67.6%, respectively. While mvBIMBAM
doesn’t inform of us what variable indirect associations are conditional
on, we can look back at our bnlearn network for more information.

# Save data and results

To conclude, save the Bayes factors from the mvBIMBAM results.

## mvBIMBAM BFs

``` r
save(m1, m2, table1, file = "data/mvBimBam_Results.RDdata")
```

# Conclusion

And with that, we conclude our `mvBIMBAM` tutorial! We hope that this
code is documented in enough detail so that you can easily adapt it for
your own projects, but feel free to reach out with any questions! -Lacey

# Session information

``` r
sessionInfo()
```

    ## R version 4.2.1 (2022-06-23)
    ## Platform: aarch64-apple-darwin20 (64-bit)
    ## Running under: macOS Ventura 13.4
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] pander_0.6.5          preprocessCore_1.60.2 lubridate_1.9.2      
    ##  [4] forcats_1.0.0         stringr_1.5.0         dplyr_1.1.2          
    ##  [7] purrr_1.0.1           readr_2.1.4           tidyr_1.3.0          
    ## [10] tibble_3.2.1          ggplot2_3.4.2         tidyverse_2.0.0      
    ## [13] knitr_1.42           
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.10      pillar_1.9.0     compiler_4.2.1   tools_4.2.1     
    ##  [5] digest_0.6.31    timechange_0.2.0 evaluate_0.20    lifecycle_1.0.3 
    ##  [9] gtable_0.3.3     pkgconfig_2.0.3  rlang_1.1.0      cli_3.6.1       
    ## [13] rstudioapi_0.14  yaml_2.3.7       xfun_0.38        fastmap_1.1.1   
    ## [17] withr_2.5.0      generics_0.1.3   vctrs_0.6.1      hms_1.1.3       
    ## [21] grid_4.2.1       tidyselect_1.2.0 glue_1.6.2       R6_2.5.1        
    ## [25] fansi_1.0.4      rmarkdown_2.21   tzdb_0.3.0       magrittr_2.0.3  
    ## [29] scales_1.2.1     htmltools_0.5.5  colorspace_2.1-0 utf8_1.2.3      
    ## [33] stringi_1.7.12   munsell_0.5.0
