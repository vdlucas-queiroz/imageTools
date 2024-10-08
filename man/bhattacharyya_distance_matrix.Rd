% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/11_auxiliary_functions.R
\name{bhattacharyya_distance_matrix}
\alias{bhattacharyya_distance_matrix}
\title{Calculate Bhattacharyya Distance Matrix for All Pairs of Classes}
\usage{
bhattacharyya_distance_matrix(
  df_samples,
  classes,
  p_function = "gaussian",
  L = NULL
)
}
\arguments{
\item{df_samples}{A data frame containing samples with a column named 'Class' indicating the class of each sample.}

\item{classes}{A vector of class labels to include in the distance matrix.}

\item{p_function}{A character string specifying the probability distribution function to be used ('gaussian' or 'gamma'). Default is 'gaussian'.}

\item{L}{An optional numeric value used in the calculation of the Bhattacharyya distance for the 'gamma' distribution.}
}
\value{
A distance matrix representing the Bhattacharyya distances between all pairs of classes.
}
\description{
This function calculates the Bhattacharyya distance matrix for all pairs of classes in a data frame of samples using either a Gaussian or Gamma distribution.
}
\examples{
df_samples <- data.frame(Class = rep(c("Class1", "Class2", "Class3"), each = 100),
                         Feature1 = rnorm(300), Feature2 = rnorm(300))
classes <- c("Class1", "Class2", "Class3")
distance_matrix_gaussian <- bhattacharyya_distance_matrix(df_samples, classes, p_function = 'gaussian')

df_samples_gamma <- data.frame(Class = rep(c("Class1", "Class2", "Class3"), each = 100),
                               Feature1 = rgamma(300, shape = 2, rate = 1))
distance_matrix_gamma <- bhattacharyya_distance_matrix(df_samples_gamma, classes, p_function = 'gamma', L = 2)
}
