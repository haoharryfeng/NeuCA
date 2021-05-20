\name{Seg_counts}
\alias{Seg_counts}
\alias{Seg_true_cell_label}
\docType{data}
\title{
Single-cell RNA-seq example dataset: Seg data
}
\description{
Seg_counts is a matrix of scRNA-seq data. Each row represents one gene. Each column represents one cell. Seg_true_cell_label is a vector of the same length as the columns of the matrix, containing the true cell labels for each cell in the same order. 
}
\usage{data(Seg_scRNA)}
\format{
Seg_counts is a matrix of gene expression values. 
Seg_true_cell_label is a vector of true cell labels for each cell. 
}
\examples{
data(Seg_scRNA)
dim(Seg_counts)
Seg_counts[1:5,1:5]
length(Seg_true_cell_label)
head(Seg_true_cell_label)
}
\keyword{datasets}
