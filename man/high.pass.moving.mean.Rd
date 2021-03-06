% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rutils.R
\name{high.pass.moving.mean}
\alias{high.pass.moving.mean}
\title{Returns moving average(mean) of a vector where values lower than a given threshold are filtered}
\usage{
high.pass.moving.mean(v, threshold = 0, weight = 14, na.rm = T)
}
\arguments{
\item{v}{numeric or integer: Given vector for which the filtered moving mean is computed}

\item{threshold}{numeric (scalar): value specifying the filtering threshold (default is \code{0})}

\item{weight}{integer (scalar): Inclusion window size for computing the moving average}

\item{na.rm}{logical: Should missing values be filtered? (default is TRUE)}
}
\value{
numeric (scalar): Filtered mean computed
}
\description{
Returns moving average(mean) of a vector where values lower than a given threshold are filtered
}
\examples{
high.pass.moving.mean(1:100, weight = 10, threshold = 95)
[1] 98
# (96 + 97 + 98 + 99 + 100)/5 = 98

}
