\name{Baron_counts}
\alias{Baron_counts}
\alias{Baron_true_cell_label}
\docType{data}
\title{
Single-cell RNA-seq example dataset: Baron data
}
\description{
Baron_counts is a matrix of scRNA-seq data. Each row represents one gene. Each column represents one cell. Baron_true_cell_label is a vector of the same length as the columns of the matrix, containing the true cell labels for each cell in the same order. 
}
\usage{data(Baron_scRNA)}
\format{
Baron_counts is a matrix of gene expression values. 
Baron_true_cell_label is a vector of true cell labels for each cell. 
}
\examples{
data(Baron_scRNA)
dim(Baron_counts)
Baron_counts[1:5,1:5]
length(Baron_true_cell_label)
head(Baron_true_cell_label)
}
\keyword{datasets}
