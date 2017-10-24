# CREAM (Clustering of genomic REgions Analysis Method) #
[![Travis Build Status](https://travis-ci.org/bhklab/CREAM.svg?branch=master)](https://travis-ci.org/bhklab/CREAM) [![codecov](https://codecov.io/gh/bhklab/CREAM/branch/master/graph/badge.svg)](https://codecov.io/gh/bhklab/CREAM) [![AppVeyor Build status](https://ci.appveyor.com/api/projects/status/gvxbin36u3yqx50s?svg=true)](https://ci.appveyor.com/project/kofiav/cream-l3j9o)

Overview
--------

CREAM is a new method for identification of clusters of genomic regions within chromosomes. Primarily, it is used for calling Clusters of cis-Regulatory Elements (COREs). CREAM uses genome-wide maps of genomic regions in the tissue or cell type of interest, such as those generated from chromatin-based assays including DNaseI, ATAC or ChIP-Seq.

CREAM considers proximity of the elements within chromosomes of a given sample to identify COREs in the following steps:

1. It identifies window size or the maximum allowed distance between the elements within each CORE 
 
2. It identifies number of elements which should be clustered as a CORE 
 
3. It calls COREs 
 
4. It filters the COREs with lowest order which does not pass the threshold considered in the approach.

Installation
------------

``` r
# Installing the development version from GitHub:
# install.packages("devtools")
devtools::install_github("bhklab/CREAM")
```

Usage
-----

``` r
# Identify COREs using CREAM
CREAM( in_path = "inst/extdata/A549_Chr21.bed", out_path = "A549_Chr21_COREs.bed", MinLength = 1000, peakNumMin = 2 )
```

Getting help
------------

Contact us by filing an issue in the CREAM [issues](https://github.com/bhklab/CREAM/issues) page.
