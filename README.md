# bcbioR

<!-- badges: start -->

[![R-CMD-check](https://github.com/bcbio/bcbioR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bcbio/bcbioR/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

The goal of `bcbioR` is to create guidelines for NGS data interpretation based on the experience of the Harvard Chan Bioinformatics Core and everybody who contributes to this package.

## Installation

You can install the development version of bcbioR from [GitHub](https://github.com/) with:

```         
# install.packages("devtools")
devtools::install_github("bcbio/bcbioR")
devtools::install_github("bcbio/bcbioR",ref = "devel")
```

## Quick start

### Set base project

use `setwd()` to set your current directory to the place where you want to work. The bcbioR functions will automatically write to whatever directory you have set.

```
setwd("/path/to/analysis/folder")
```

The following code will pop up a Rmd template will populate that folder with HCBC data structure guidelines

```
use_bcbio_projects(path,pipeline="nf-core/rnaseq")
use_bcbio_projects(path,pipeline="singlecell")
```

### Set RNAseq report folder

This code will populate the folder with HCBC data structure guidelines and Rmd code: **You do not need to create a reports folder prior to running this code. This will create and populate the reports folder.**

``` 
bcbio_templates(type="rnaseq", outpath="/path/to/analysis/folder/reports")
```

### Discover more…

Go to the vignette to know more `vignette("bcbioR_quick_start",package="bcbioR")`

## How to Contribute

### Open an issue

-   If you find a bug
-   If you want a new feature
-   If you want to add code to the templates

### Modify the code

-   Clone the repository
-   Make sure you are in the `devel` branch
-   Create a new branch `git checkout -b feature1`
-   Modify you code, add and commit
-   Push to GitHub the new branch
-   Create a PR from your branch to `devel`
-   Assignt the PR to me or Alex

Some best practices when developing:

-   install `devtools`
-   Use `usethis::use_import_from("stringr","str_replace_all")` to add a new function you are using in the code.

### Contributors

-   Lorena Pantano
-   Alex Bartlett
-   Emma Berdan
-   Heather Wick
-   James Billingsley
-   Zhu Zhuo
-   Elizabeth Partan
