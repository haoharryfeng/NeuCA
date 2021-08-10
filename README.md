# NeuCA: Neural-network based Cell Annotation tool

-------------------
We developed NEUral-network based Cell Annotation, `NeuCA`, a R/Bioconductor tool for cell type annotation using single-cell RNA-seq data. It is a supervised cell label assignment method that uses existing scRNA-seq data with known labels to train a neural network-based classifier, and then predict cell labels in single-cell RNA-seq data of interest.

### `NeuCA` is currently available on Bioconductor: <https://bioconductor.org/packages/NeuCA/>


# 1. Introduction
The fast advancing single cell RNA sequencing (scRNA-seq) technology enables transcriptome study in heterogeneous tissues at a single cell level. The initial important step of analyzing scRNA-seq data is to accurately annotate cell labels. We present a neural-network based cell annotation method NeuCA. When closely correlated cell types exist, NeuCA uses the cell type tree information through a hierarchical structure of neural networks to improve annotation accuracy. Feature selection is performed in hierarchical structure to further improve classification accuracy. When cell type correlations are not high, a feed-forward neural network is adopted.
![workflow](vignettes/workflow.png)

NeuCA depends on the following packages in R/Bioconductor:

- **keras**, for neural-network interface in R
- **limma**, for linear model framework and testing markers
- **SingleCellExperiment**, for data organization formatting
- **e1071**, for probability and predictive functions.


# 2. Preparing NeuCA input files: `SingleCellExperiment` class

The scRNA-seq data input for NeuCA must be objects of the Bioconductor `SingleCellExperiment`. You may need to read corresponding vignettes on how to create a SingleCellExperiment from your own data. An example is provided here to show how to do that, but please note this is not a comprehensive guidance for SingleCellExperiment.

**Step 1**: Load in example scRNA-seq data.
```
library(NeuCA)
#Baron_scRNA is the training scRNA-seq dataset
#Seg_scRNA is the testing scRNA-seq dataset
data("Baron_scRNA")
data("Seg_scRNA")
```

**Step 2a**: Prepare training data as a SingleCellExperiment object
```
Baron_anno = data.frame(Baron_true_cell_label, row.names = colnames(Baron_counts))
Baron_sce = SingleCellExperiment(
    assays = list(normcounts = as.matrix(Baron_counts)),
    colData = Baron_anno
    )
# use gene names as feature symbols
rowData(Baron_sce)$feature_symbol <- rownames(Baron_sce)
# remove features with duplicated names
Baron_sce <- Baron_sce[!duplicated(rownames(Baron_sce)), ]
```

**Step 2b**: Similarly, prepare testing data as a SingleCellExperiment object. Note the true cell type labels are not necessary (and of course often not available).
```
Seg_anno = data.frame(Seg_true_cell_label, row.names = colnames(Seg_counts))
Seg_sce <- SingleCellExperiment(
    assays = list(normcounts = as.matrix(Seg_counts)),
    colData = Seg_anno
)
# use gene names as feature symbols
rowData(Seg_sce)$feature_symbol <- rownames(Seg_sce)
# remove features with duplicated names
Seg_sce <- Seg_sce[!duplicated(rownames(Seg_sce)), ]
```

# 3. NeuCA training and prediction
**Step 3**: with both training and testing data as objects in SingleCellExperiment class, now we can train the classifier in NeuCA and predict testing dataset’s cell types. This process can be achieved with one line of code:
```
predicted.label = NeuCA(train = Baron_sce, test = Seg_sce, 
                        model.size = "big", verbose = FALSE)
```
NeuCA can detect whether highly-correlated cell types exist in the training dataset, and automatically determine if a general neural-network model will be adopted or a marker-guided hierarchical neural-network will be adopted for classification.

Users have the option to determine the complexity of the neural-network used in NeuCA by specifying the desired `model.size` argument. Here, “big”, “medium” and “small” are 3 potential choices, reflecting large, medium and small number of nodes and layers in neural-network, respectively.


# 4. Predicted cell types
`predicted.label` is a vector of the same length with the number of cells in the testing dataset, containing all cell’s predicted cell type. It can be viewed directly:
```
head(predicted.label)
## [1] "alpha" "gamma" "gamma" "gamma" "gamma" "alpha"
table(predicted.label)
## predicted.label
##       alpha        beta       delta      ductal endothelial       gamma 
##         331         109          56          65           9         132
```

[**Optional**] If you have the true cell type labels for the testing dataset, you may evaluate the predictive performance by a confusion matrix:
```
table(predicted.label, Seg_true_cell_label)
##                Seg_true_cell_label
## predicted.label alpha beta delta ductal endothelial gamma
##     alpha         328    0     0      0           0     3
##     beta            0  109     0      0           0     0
##     delta           0    0    56      0           0     0
##     ductal          1    0     0     64           0     0
##     endothelial     0    0     0      0           9     0
##     gamma           0    0     0      0           0   132
```
You may also draw a Sankey diagram to visualize the prediction accuracy:
![workflow](vignettes/sankey1.png)


# Session info

```
## R version 4.0.3 (2020-10-10)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Big Sur 10.16
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] networkD3_0.4               NeuCA_0.1.0                
##  [3] SingleCellExperiment_1.12.0 SummarizedExperiment_1.20.0
##  [5] Biobase_2.50.0              GenomicRanges_1.42.0       
##  [7] GenomeInfoDb_1.26.2         IRanges_2.24.1             
##  [9] S4Vectors_0.28.1            BiocGenerics_0.36.0        
## [11] MatrixGenerics_1.2.1        matrixStats_0.58.0         
## [13] e1071_1.7-6                 limma_3.46.0               
## [15] keras_2.4.0                 BiocStyle_2.18.1           
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.6             XVector_0.30.0         compiler_4.0.3        
##  [4] BiocManager_1.30.10    zlibbioc_1.36.0        bitops_1.0-6          
##  [7] base64enc_0.1-3        class_7.3-17           tools_4.0.3           
## [10] zeallot_0.1.0          digest_0.6.27          jsonlite_1.7.2        
## [13] evaluate_0.14          lattice_0.20-41        pkgconfig_2.0.3       
## [16] rlang_0.4.10           igraph_1.2.6           Matrix_1.2-18         
## [19] DelayedArray_0.16.1    yaml_2.2.1             xfun_0.21             
## [22] GenomeInfoDbData_1.2.4 stringr_1.4.0          knitr_1.31            
## [25] htmlwidgets_1.5.3      generics_0.1.0         grid_4.0.3            
## [28] reticulate_1.18        R6_2.5.0               rmarkdown_2.7         
## [31] bookdown_0.21          magrittr_2.0.1         whisker_0.4           
## [34] tfruns_1.5.0           htmltools_0.5.1.1      tensorflow_2.4.0      
## [37] stringi_1.5.3          proxy_0.4-25           RCurl_1.98-1.2
```
