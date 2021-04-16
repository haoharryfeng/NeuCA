\name{NeuCA}
\alias{NeuCA}

\title{
A NEUral-network based Cell Annotation (NeuCA) tool for cell type annotation using single-cell RNA-seq data.
}
\description{
NeuCA is a supervised cell label assignment method that uses existing scRNA-seq data with known labels to train a neural network-based classifier, and then predict cell labels in data of interest.
}
\usage{
NeuCA(train, test, model.size = "big", verbose = FALSE)
}

\arguments{
  \item{train}{
A training scRNA-seq dataset where cell labels are already known. Must be an object of SingleCellExperiment class. Must contain cell labels as the first column in its colData.
}
  \item{test}{
A testing scRNA-seq dataset where cell labels are unknown. Must be an object of SingleCellExperiment class.
}
  \item{model.size}{
an ordinal variable indicating the complexity of the neural-network. Must be one of the following: "big", "medium" or "small"
}
  \item{verbose}{
A Boolean variable (TRUE/FALSE) indicating whether additional information about the training and testing process will be printed.
}
}
\details{
When closely correlated cell types exist, NeuCA uses the cell type tree information through a hierarchical structure of neural networks to improve annotation accuracy. Feature selection is performed in hierarchical structure to further improve classification accuracy. When cell type correlations are not high, a feed-forward neural network is adopted.
}
\value{
NeuCA returns a vector of predicted cell types. The output vector has the same length with the number of cells in the testing dataset.
}
\author{
Hao Feng <hxf155@case.edu>
}
\note{
The input single-cell RNA-seq data, for both training and testing, should be objects of SingleCellExperiment class. The true cell type labels, for the training dataset, should be stored as the first column in its SingleCellExperiment "colData"" object.
}

\examples{
#1. Load in example scRNA-seq data
#Baron_scRNA is the training scRNA-seq dataset
#Seg_scRNA is the testing scRNA-seq dataset
data("Baron_scRNA")
data("Seg_scRNA")

#2. Create SingleCellExperiment object as the input for NeuCA (if data are not already in SingleCellExperiment format)
Baron_anno = data.frame(Baron_true_cell_label, row.names = colnames(Baron_counts))
Baron_sce = SingleCellExperiment(
    assays = list(normcounts = as.matrix(Baron_counts)),
    colData = Baron_anno
    )
# use gene names as feature symbols
rowData(Baron_sce)$feature_symbol <- rownames(Baron_sce)
# remove features with duplicated names
Baron_sce <- Baron_sce[!duplicated(rownames(Baron_sce)), ]

#similarly for Seg data
Seg_anno = data.frame(Seg_true_cell_label, row.names = colnames(Seg_counts))
Seg_sce <- SingleCellExperiment(
    assays = list(normcounts = as.matrix(Seg_counts)),
    colData = Seg_anno
)
# use gene names as feature symbols
rowData(Seg_sce)$feature_symbol <- rownames(Seg_sce)
# remove features with duplicated names
Seg_sce <- Seg_sce[!duplicated(rownames(Seg_sce)), ]


#3. NeuCA training and cell type prediction
predicted.label = NeuCA(train = Baron_sce, test = Seg_sce, model.size = "big", verbose = FALSE)
head(predicted.label)
#Seg_sce have its ground true cell type stored, compare the predicted vs. the truth.
sum(predicted.label==colData(Seg_sce)[,1])/length(predicted.label)
}
