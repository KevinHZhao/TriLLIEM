% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MASlrt.R
\name{MaSlrt}
\alias{MaSlrt}
\title{Mating asymmetry likelihood ratio test}
\usage{
MaSlrt(dat)
}
\arguments{
\item{dat}{A dataframe formatted with the same columns as \link{example_dat4R}}
}
\value{
A p-value for the Likelihood Ratio Test.
}
\description{
\code{MaSlrt()} runs a Likelihood Ratio Test for a trios model with Mating Asymmetry.
}
\examples{
MaSlrt(example_dat4R)
}
