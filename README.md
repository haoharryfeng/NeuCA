# NeuCA
## Neural-network based Cell Annotation tool

We developed NEUral-network based Cell Annotation, `NeuCA`, a tool for cell type annotation using single-cell RNA-seq data. It is a supervised cell label assignment method that uses existing scRNA-seq data with known labels to train a neural network-based classifier, and then predict cell labels in single-cell RNA-seq data of interest.

# 1. Introduction
The fast advancing single cell RNA sequencing (scRNA-seq) technology enables transcriptome study in heterogeneous tissues at a single cell level. The initial important step of analyzing scRNA-seq data is to accurately annotate cell labels. We present a neural-network based cell annotation method NeuCA. When closely correlated cell types exist, NeuCA uses the cell type tree information through a hierarchical structure of neural networks to improve annotation accuracy. Feature selection is performed in hierarchical structure to further improve classification accuracy. When cell type correlations are not high, a feed-forward neural network is adopted.

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
