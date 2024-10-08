% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/11_auxiliary_functions.R
\name{extract_hierarchy}
\alias{extract_hierarchy}
\title{Extract Hierarchy from Hierarchical Clustering Object}
\usage{
extract_hierarchy(hclust_obj)
}
\arguments{
\item{hclust_obj}{An object of class \code{hclust} produced by hierarchical clustering.The leafs of this object have to be labeled with strings.}
}
\value{
A list where each element represents a level of the hierarchy and contains the clusters at that level.
}
\description{
This function extracts the hierarchy of clusters from an object of class \code{hclust}.
}
\examples{
dist_matrix <- as.dist(matrix(runif(100), nrow = 10, dimnames = list(LETTERS[1:10], LETTERS[1:10])))
hc <- hclust(dist_matrix, method = "ward.D2")
hierarchy <- extract_hierarchy(hc)
print(hierarchy)
}
