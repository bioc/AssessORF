% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CompareAssessmentResults.R
\name{CompareAssessmentResults}
\alias{CompareAssessmentResults}
\title{Compare Assessment Results}
\usage{
CompareAssessmentResults(obj1,
                         obj2,
                         printSummary = TRUE,
                         returnDetails = FALSE)
}
\arguments{
\item{obj1, obj2}{Objects of class \code{Assessment} and subclass \code{Results} to compare against each other.
Alternatively, either \code{obj1} or \code{obj2} (or both) can be a two-element character vector that specifies one
of such objects from \code{AssessORFData}. The first element in the vector should be the strain identifier, and the
second element should be the gene source identifier. Both objects should have been generated from the same mapping
object.}

\item{printSummary}{Logical indicating whether or not to print out a summary of the comparison analysis.}

\item{returnDetails}{Logical indicating whether or not to return a list of details from the comparison analysis.
See the next section for what items are in the list.}
}
\value{
If \code{returnDetails} is true, the function returns a 11-item list. Otherwise, the function invisibly
returns object 1.
}
\description{
Compare two objects of class Assessment. subclass Results to determine how their gene sets and the
corresponding category assignments vary
}
\details{
Since the same mapping object (an object class \code{Assessment} and subclass \code{DataMap}) can be used to assess
multiple sets of genes for genome, it is meaningful to compare how those gene sets and their category assignments
from \code{AssessGenes} vary from one another. To make describing this function easier, let us assume that one set
of genes consists of the complete set of predictions made by a gene-finding program on a particular strain's genome
and that the other set of genes consists of the complete set of predictions from a second gene-finding program.

When gene-finding programs predict genes for a genome, they make a decision on which regions of the genome code for
proteins. There is (usually) only one option for the stop codon that ends a particular coding region, but there are
typically multiple options available for the start codon that will mark the beginning of a region. It is therefore
useful to find out which (general) coding regions the two programs agree on by determining which stops are found in
both sets of predicted genes. From there, the starts each program picked for those shared coding regions can be
compared to see whether they agree or not. If the same start is chosen by both programs for a particular shared
stop / coding region, then that is an example of a gene predicted by both programs. If the starts chosen by the two
programs for a particular shared stop / coding region are different, then that it is an example of both programs
agreeing that that particular region of the genome codes for protein but disagreeing on where in the genome that
region starts. It would be interesting to see what category was assigned to each start, especially if one start has
evolutionary conservation and the other does not.

This function compares the set of genes and corresponding category assignments in the object specified by
\code{obj1} (object 1) to the set of genes and corresponding category assignments in the object specified by
\code{obj2} (object 2). It then reports the results of the comparison analysis in the format specified by the
logical parameters \code{printSummary} and \code{returnDetails}.

If \code{printSummary} is true, the function prints out the following information: the number of shared coding
regions (i.e., the number of stops in both gene sets), the number of shared genes (i.e., the number of times both a
start and its corresponding stop are found in both sets), and the number of instances where a stop is found in both
gene sets but the corresponding starts in each set disagree. For the shared stop - different start set, the function
also prints the number of instances where the start from one object has conservation evidence while the
corresponding start in the other object does not.

If \code{returnDetails} is true, the function returns a 11-item list. Each item of the list is described below. The
contents of the object 1 and object 2 gene vectors correspond to the ordering of the genes inside object 1 or object
2, respectively. For the category assignment matrix for shared stop - different start set, it is possible for the
gene in object 1 to be assigned to the same category as the corresponding gene from object 2, and the table reflects
that.
\itemize{
 \item{"StrainID"}{Same as the strain identifier inside \code{obj1} and \code{obj2}}
 \item{Species}{Same as the strain identifier inside objects 1 and 2}
 \item{Obj1_GeneSource}{Same as the gene source identifier inside \code{obj1}}
 \item{Obj2_GeneSource}{Same as the gene source identifier inside \code{obj2}}
 \item{Obj1_Genes_SharedCodingRegions}{The genes from object 1 that share a stop with a gene from the object 2}
 \item{Obj2_Genes_SharedCodingRegions}{The genes from object 2 that share a stop with a gene from the object 1}
 \item{Obj1_Genes_SharedGenes}{The genes from object 1 that share a start and stop with a gene from the object 2}
 \item{Obj2_Genes_SharedGeness}{The genes from object 2 that share a start and stop with a gene from the object 1}
 \item{Obj1_Genes_SharedStopDiffStart}{The genes from object 1 that share a stop with a gene from the object 2 but
 have a different start from the corresponding object 2 gene}
 \item{Obj2_Genes_SharedStopDiffStart}{The genes from object 2 that share a stop with a gene from the object 1 but
 have a different start from the corresponding object 1 gene}
 \item{CategoryTable_SharedStopDiffStart}{A 12-by-12 matrix describing the number of times the gene from object 1
 was assigned one category and the gene from object 2 was assigned some other category for the
 shared stop - different start set}
}

Please ensure that \code{obj1} and \code{obj2} come from the same strain / mapping object. The function will do its
best to make sure the identifying information for \code{obj1} and \code{obj2} match.

\code{printSummary} and \code{returnDetails} cannot both be FALSE.
}
\examples{

## Example showing how to use the function with the AssessORFData package:

\dontrun{
compare1 <- CompareAssessmentResults(obj1 = c("MGAS5005", "Prodigal"),
                                     obj2 = c("MGAS5005", "GeneMarkS2"),
                                     printSummary = TRUE,
                                     returnDetails = TRUE)
}

resObj1 <- readRDS(system.file("extdata",
                               "MGAS5005_PreSaved_ResultsObj_Prodigal.rds",
                               package = "AssessORF"))
                               
resObj2 <- readRDS(system.file("extdata",
                               "MGAS5005_PreSaved_ResultsObj_GeneMarkS2.rds",
                               package = "AssessORF"))
                               
compare2 <- CompareAssessmentResults(obj1 = resObj1,
                                     obj2 = resObj2,
                                     printSummary = TRUE,
                                     returnDetails = TRUE)

}
\seealso{
\code{\link{Assessment-class}}, \code{\link{AssessGenes}}
}
