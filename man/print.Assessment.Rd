% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Assessment-class.R
\name{print.Assessment}
\alias{print.Assessment}
\title{Print Assessment Objects}
\usage{
\method{print}{Assessment}(x, ...)
}
\arguments{
\item{x}{An object of class \code{Assessment} and of either subclass \code{DataMap} or subclass \code{Results}.}

\item{...}{Further printing parameters.}
}
\value{
Invisibly returns the input object \code{x}
}
\description{
The \code{print} method for \code{Assessment} objects
}
\details{
If \code{x} is of subclass \code{DataMap}, the length of the genome is printed along with any supplied identifying
information for the genome.

If \code{x} is of subclass \code{Results}, the number of genes in each category and the accuracy scores are printed out
along with any supplied identifying information.
}
\examples{

print(readRDS(system.file("extdata",
                          "MGAS5005_PreSaved_DataMapObj.rds",
                          package = "AssessORF")))

print(readRDS(system.file("extdata",
                          "MGAS5005_PreSaved_ResultsObj_Prodigal.rds",
                          package = "AssessORF")))

}
\seealso{
\code{\link{Assessment-class}}
}
