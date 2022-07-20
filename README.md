## Multi-Study Multi-Class Concordance Analysis (MSCC)
 
MSCC is an R package for calculating the concordance among multiple classes across multiple studies. 
For example, if we have multiple studies of gene expression matrices showing the expression patterns across 
multiple classes, MSCC could then identify the genes sharing same expression pattern across those studies.

### Installation

The package can be downloaded and installed from GitHub.

```
install.packages('devtools')
devtools::install_github('jianzou75/CASCAM')
```


### Quick start
Users can call the below function to perform a complete analysis.

```
agemap_result <- mscc.permute.test(study.data.list = agemap_data_list, study.label.list = agemap_label_list)
```
