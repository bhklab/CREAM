# CREAM #

Overview
--------

CREAM is a new method for identification of clusters of functional regions (COREs) within chromosomes. It uses genome-wide maps of functional regions in the tissue or cell type of interest, such as those generated from chromatin-based assays including DNaseI, ATAC or ChIP-Seq.

CREAM considers proximity of the elements within chromosomes of a given sample to identify COREs in the following steps:

1. It identifies window size or the maximum allowed distance between the elements within each CORE 
 
2. It identifies number of elements which should be clustered as a CORE 
 
3. It calls COREs 
 
4. It filters the COREs with lowest order which does not pass the threshold considered in the appraoch.

Installation
------------

``` r
# Installing the development version from GitHub:
# install.packages("devtools")
devtools::install_github("bhklab/CREAM")
```

Usage
-----

Under construction

Getting help
------------

Contact us by filing an issue in the CREAM [issues](https://github.com/bhklab/CREAM/issues) page.