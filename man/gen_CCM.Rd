% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_CCM.R
\name{gen_CCM}
\alias{gen_CCM}
\title{Positive Semidefinite Correlation Matrix}
\usage{
gen_CCM(d)
}
\arguments{
\item{d}{the dimension of matrix.}
}
\value{
a correlation coefficients matrix.
}
\description{
Generate a positive semidefinite correlation coefficients matrix
}
\examples{
X = gen_CCM(4)
print(X)
}
