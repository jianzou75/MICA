## Mutual Information Concordance Analysis (MICA)

MICA is an R package for detecting biomarkers with concordant multi-class expression pattern across multiple omics studies from an information theoretic perspective.

### Installation

The package can be downloaded and installed from GitHub.

```r
install.packages('devtools')
devtools::install_github('jianzou75/MICA')
```


### Quick start

Users can call the below function to perform a complete analysis.

```r
metabolism_result <- mica.full(metabolism_data, metabolism_label)
```
