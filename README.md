
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ProTransDeconv
ProTransDeconv performs data transformation, coefficient of variation (CV) summarization, and identification of cell-type-specific proteins using protein quantification data derived from Mass Spectrometry (MS).

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
### Output: results
The result is a list containing the following elements:
transformed_list: Protein expression matrices processed using different transformation methods.
cv_summary: A summary of the CV under different transformation methods, including mean, median, proportion greater than 0.25, and so on.
cv_plot_data: Gene-wise CV values across various transformation methods, used for visualization.
EDec, rodeo, csSAM, bMIND:  If cell type proportions are provided, these contain the estimated cell-type-specific expression profiles and specificity scores using respective deconvolution methods under different transformations.
gene_cell_correlation: If cell type proportions are provided, this shows the correlation between cell type proportions and bulk protein expression across samples for each gene.
gold_standard_markers: If cell type proportions are provided, this includes the inferred gold-standard marker genes for each cell type based on correlation and significance thresholds.

```


You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
