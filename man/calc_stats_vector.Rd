% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/13_confusion_matrix.R
\name{calc_stats_vector}
\alias{calc_stats_vector}
\title{Calculate Statistics for a Numeric Vector}
\usage{
calc_stats_vector(vec)
}
\arguments{
\item{vec}{A numeric vector.}
}
\value{
A data frame containing the mean, standard deviation, and 95% confidence interval (lower and upper) for the vector.
}
\description{
This function calculates the mean, standard deviation, and 95% confidence interval for a numeric vector.
}
\examples{
vec <- rnorm(100)
stats <- calc_stats_vector(vec)
print(stats)
}
