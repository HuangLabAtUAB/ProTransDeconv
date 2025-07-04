% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_transformations.R
\name{Min_max_normalize}
\alias{Min_max_normalize}
\title{Min-Max Normalization}
\usage{
Min_max_normalize(data, a = 0, b = 1)
}
\arguments{
\item{data}{A numeric matrix or data frame.}

\item{a}{Lower bound of the target range.}

\item{b}{Upper bound of the target range.}
}
\value{
A normalized matrix.
}
\description{
Scales data to a specified range (default 0 to 1) column-wise.
}
