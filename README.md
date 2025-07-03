
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ProTransDeconv
ProTransDeconv performs data transformation, coefficient of variation summarization, and identification of cell-type-specific proteins using protein quantification data derived from Mass Spectrometry (MS).

## WorkFlow

This is the workflow of the package：

<figure>
<img src="images/Workflow.png" alt="Workflow Diagram" />
</figure>

<!-- badges: start -->
<!-- badges: end -->

## Installation

You can install the ProTransDeconv like so:

``` r
# install devtools if necessary
install.packages('devtools')

# install the ProTransDeconv package
devtools::install_github('HuangLabAtUAB/ProTransDeconv')

# load the ProTransDeconv package
library(ProTransDeconv)
```

## Example

This is a basic example which shows you how to run ProTransDeconv:

``` r
library(ProTransDeconv)
results <- ProTransDeconv(
    data = Protein_Quantification,
    type = "intensity", run_bmind = FALSE,
    cell_proportion = Cell_Type_Proportion
)
```

Result Explanation:

``` r
### Output: `results`
The result is a list containing the following components:
- `transformed_list`: transformation

summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```


You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
