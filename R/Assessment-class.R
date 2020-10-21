#' @export
#' 
#' @method as.matrix Assessment
#'
#' @title Tabulate the Category Assignments for Assessment Results Objects
#' @description The \code{as.matrix} method for \code{Assessment} and subclass \code{Results} objects
#'
#' @param x An object of class \code{Assessment} and subclass \code{Results}.
#'
#' @param ... Additional arguments.
#'
#' @details
#' \code{as.matrix.Assessment} tabulates and returns the number of times each category appears in the \code{CategoryAssignments}
#' vector within the given \code{Results} object. If the number of genes for any the 14 main gene / ORF categories is zero, a
#' count (of zero) will still be included for that category.
#'
#' @return A one-row matrix with the counts for the number of genes/ORFs that fall into each category. The corresponding
#' category codes serve as the column names, and the name of the row is the strain ID.
#'
#' @seealso \code{\link{Assessment-class}}
#' 
#' @examples
#'
#' as.matrix(readRDS(system.file("extdata",
#'                               "MGAS5005_PreSaved_ResultsObj_Prodigal.rds",
#'                               package = "AssessORF")))
#'
as.matrix.Assessment <- function(x, ...) {
  if (class(x)[1] != "Assessment") {
    stop("'x' must be an object of class 'Assessment'.")
  }
  
  if (class(x)[2] == "Results") {
    speciesName <- x$Species
    strainID <- x$StrainID
    
    if ((speciesName != "") && (strainID != "")) {
      genomeID <- paste(speciesName, strainID, sep = "_")
    } else if ((speciesName != "")) {
      genomeID <- speciesName
    } else if ((strainID != "")) {
      genomeID <- strainID
    } else {
      genomeID <- ""
    }
    
    geneSource <- x$GeneSource
    
    if ((genomeID != "") && (geneSource != "")) {
      rowTitle <- paste(genomeID, geneSource)
    } else if (genomeID != "") {
      rowTitle <- paste(genomeID, "UnknownSource")
    } else if (geneSource != "") {
      rowTitle <- paste("UnknownGenome", geneSource)
    } else {
      rowTitle <- "UnknownGenome UnknownSource"
    }
    
    ## --------------------------------------------------------------------------------------------------------------- ##
    
    ## Get the number of genes in each category.
    catSumTable <- table(x$CategoryAssignments)
    
    allCatSums <- integer(14L)
    
    names(allCatSums) <- c("Y CS+ PE+", "Y CS+ PE-", "Y CS- PE+", "Y CS- PE-",
                           "Y CS< PE!", "Y CS- PE!", "Y CS! PE+", "Y CS! PE-",
                           "Y CS> PE+", "Y CS> PE-", "Y CS< PE+", "Y CS< PE-",
                           "N CS< PE+", "N CS- PE+")
    
    allCatSums[names(catSumTable)] <- catSumTable
    
    allCatSums["N CS< PE+"] <- sum(nrow(x$'N_CS<_PE+_ORFs'), na.rm = TRUE)
    allCatSums["N CS- PE+"] <- sum(nrow(x$'N_CS-_PE+_ORFs'), na.rm = TRUE)
    
    return(matrix(allCatSums, nrow = 1, ncol = length(allCatSums),
                  dimnames = list(rowTitle, names(allCatSums))))
  } else {
    stop("x must be of subclass 'Results'.")
  }
}

#' @export
#' 
#' @method print Assessment
#'
#' @title Print Assessment Objects
#' @description The \code{print} method for \code{Assessment} objects
#'
#' @param x An object of class \code{Assessment} and of either subclass \code{DataMap} or subclass \code{Results}.
#'
#' @param ... Further printing parameters.
#'
#' @details
#' If \code{x} is of subclass \code{DataMap}, the length of the genome is printed along with any supplied identifying
#' information for the genome.
#'
#' If \code{x} is of subclass \code{Results}, the number of genes in each category and the accuracy scores are printed out
#' along with any supplied identifying information.
#' 
#' @return Invisibly returns the input object \code{x}
#'
#' @seealso \code{\link{Assessment-class}}
#' 
#' @examples
#'
#' print(readRDS(system.file("extdata",
#'                           "MGAS5005_PreSaved_DataMapObj.rds",
#'                           package = "AssessORF")))
#' 
#' print(readRDS(system.file("extdata",
#'                           "MGAS5005_PreSaved_ResultsObj_Prodigal.rds",
#'                           package = "AssessORF")))
#'
print.Assessment <- function(x, ...) {
  if (class(x)[1] != "Assessment") {
    stop("'x' must be an object of class 'Assessment'.")
  }
  
  if (class(x)[2] == "DataMap") {
    speciesName <- x$Species
    strainID <- x$StrainID
    
    if ((speciesName != "") && (strainID != "")) {
      genomeID <- paste(speciesName, strainID)
    } else if ((speciesName != "")) {
      genomeID <- speciesName
    } else if ((strainID != "")) {
      genomeID <- strainID
    } else {
      genomeID <- ""
    }
    
    printOut <- "An Assessment object that maps proteomics data and evolutionary data\n"
    
    if (genomeID != "") {
      printOut <- paste(printOut, genomeID, "\n", sep = "")
    }
    
    cat(paste(printOut, "Genome Length: ", x$GenomeLength, "\n", sep = ""))
    
  } else if (class(x)[2] == "Results") {
    
    speciesName <- x$Species
    strainID <- x$StrainID
    
    if ((speciesName != "") && (strainID != "")) {
      genomeID <- paste(speciesName, strainID)
    } else if ((speciesName != "")) {
      genomeID <- speciesName
    } else if ((strainID != "")) {
      genomeID <- strainID
    } else {
      genomeID <- ""
    }
    
    printOut <- "An Assessment object that categorizes predicted genes\n"
    
    if (genomeID != "") {
      printOut <- paste(printOut, "Strain: ", genomeID, "\n", sep = "")
    }
    
    if (x$GeneSource != "") {
      printOut <- paste(printOut, "Number of Genes Provided from ", x$GeneSource, ": ", x$NumGenes, "\n\n", sep = "")
    } else {
      printOut <- paste(printOut, "Number of Genes Provided: ", x$NumGenes, "\n\n", sep = "")
    }
    
    ## --------------------------------------------------------------------------------------------------------------- ##
    
    ## Get the number of genes in each category.
    catSumTable <- table(x$CategoryAssignments)
    
    allCatSums <- integer(14L)
    
    names(allCatSums) <- c("Y CS+ PE+", "Y CS+ PE-", "Y CS- PE+", "Y CS- PE-",
                           "Y CS< PE!", "Y CS- PE!", "Y CS! PE+", "Y CS! PE-",
                           "Y CS> PE+", "Y CS> PE-", "Y CS< PE+", "Y CS< PE-",
                           "N CS< PE+", "N CS- PE+")
    
    allCatSums[names(catSumTable)] <- catSumTable
    allCatSums["N CS< PE+"] <- sum(nrow(x$'N_CS<_PE+_ORFs'), na.rm = TRUE)
    allCatSums["N CS- PE+"] <- sum(nrow(x$'N_CS-_PE+_ORFs'), na.rm = TRUE)
    
    ## --------------------------------------------------------------------------------------------------------------- ##
    
    printOut <- paste(printOut, "Score Using All Evidence: ",
                      round(ScoreAssessmentResults(x), 4), "\n", sep = "")
    
    printOut <- paste(printOut, "Score Using Only Proteomic Evidence: ",
                      round(ScoreAssessmentResults(x, "p"), 4), "\n", sep = "")
    
    printOut <- paste(printOut, "Score Using Only Conserved Start Evidence: ",
                      round(ScoreAssessmentResults(x, "c"), 4), "\n", sep = "")
    
    ## --------------------------------------------------------------------------------------------------------------- ##
    
    catIndex <- names(allCatSums) == "Y CS+ PE+"
    printOut <- paste0(printOut, names(allCatSums[catIndex]),
                       " (correct with strong evidence): ", allCatSums[catIndex], "\n")
    
    catIndex <- (names(allCatSums) == "Y CS+ PE-") | (names(allCatSums) == "Y CS- PE+")
    printOut <- paste0(printOut, paste0(names(allCatSums[catIndex]),
                                        " (correct with some evidence): ",
                                        allCatSums[catIndex], "\n" , collapse = ""))
    
    catIndex <- names(allCatSums) == "Y CS- PE-"
    printOut <- paste0(printOut, names(allCatSums[catIndex]),
                       " (no evidence): ", allCatSums[catIndex], "\n" )
    
    catIndex <- names(allCatSums) == "Y CS< PE!"
    printOut <- paste0(printOut, names(allCatSums[catIndex]),
                       " (definitely incorrect): ", allCatSums[catIndex], "\n")
    
    catIndex <- (names(allCatSums) == "Y CS- PE!") | (names(allCatSums) == "Y CS! PE+") |
      (names(allCatSums) == "Y CS! PE-")
    printOut <- paste0(printOut, paste0(names(allCatSums[catIndex]),
                                        " (likely incorrect): ",
                                        allCatSums[catIndex], "\n" , collapse = ""))
    
    catIndex <- (names(allCatSums) == "Y CS> PE+") | (names(allCatSums) == "Y CS> PE-") |
      (names(allCatSums) == "Y CS< PE+") | (names(allCatSums) == "Y CS< PE-")
    printOut <- paste0(printOut, paste0(names(allCatSums[catIndex]),
                                        " (potentially incorrect): ",
                                        allCatSums[catIndex], "\n" , collapse = ""))
    
    catIndex <- names(allCatSums) == "N CS< PE+"
    printOut <- paste0(printOut, names(allCatSums[catIndex]),
                       " (likely missing genes): ", allCatSums[catIndex], "\n")
    
    catIndex <- names(allCatSums) == "N CS- PE+"
    printOut <- paste0(printOut, names(allCatSums[catIndex]),
                       " (potentially missing genes): ", allCatSums[catIndex], "\n")
    
    ## --------------------------------------------------------------------------------------------------------------- ##
    
    printOut <- paste(printOut, "\nNumber of Genes with Supporting Evidence: ",
                      sum(allCatSums[c("Y CS+ PE+", "Y CS+ PE-", "Y CS- PE+")]), sep = "")
    
    printOut <- paste(printOut, "\nNumber of Genes with Contradictory Evidence: ",
                      sum(allCatSums[c("Y CS> PE+", "Y CS> PE-", "Y CS< PE+", "Y CS< PE-",
                                       "Y CS! PE+", "Y CS! PE-", "Y CS< PE!", "Y CS- PE!")]),
                      sep = "")
    
    printOut <- paste(printOut, "\nNumber of ORFs with Protein Evidence and No Given Start: ",
                      sum(allCatSums[c("N CS< PE+", "N CS- PE+")]), "\n", sep = "")
    
    cat(printOut)
    
  } else {
    stop("'x' is an object of unrecognized subclass '", class(x)[2], "'.")
  }
  
  invisible(x)
}

#' @export
#' @import graphics
#' @importFrom  grDevices gray rgb
#' 
#' @method plot Assessment
#'
#' @title Plot Assessment Objects
#' @description The \code{plot} method for \code{Assessment} objects
#'
#' @param x An object of class \code{Assessment} and of either subclass \code{DataMap} or subclass \code{Results}.
#'
#' @param y An optional object of class \code{Assessment} and of either subclass \code{DataMap} or subclass \code{Results}. Its
#' subclass must be different than the subclass of \code{x}
#'
#' @param minConCovRatio_GV Minimum value of the conservation to coverage ratio needed to call a start conserved. Must range
#' from 0 to 1. Lower values allow more conserved starts through. Only used with the genome viewer. Default value is recommended.
#'
#' @param interactive_GV Logical specifying whether or not the genome viewer plot should be interactive. Default is TRUE.
#'
#' @param rangeStart_GV,rangeEnd_GV Optional positive integer values that can be specified when generating a genome viewer plot in
#' order to have the plot zoom into the range of genomic positions between those values. Both values must be within the bounds of
#' the genome. Omitted (NA) by default, which results in a plot spanning the whole genome.
#'
#' @param ... Further plotting parameters.
#'
#' @details
#' If out of \code{x} and \code{y} only \code{x} is specified and \code{x} is of subclass \code{Results}, a bar chart describing
#' the number of genes in each category is plotted. For the predicted gene categories, bars are colored by the correctness of
#' that category, where dark green represents "definitely correct", light green represents "likely correct", white represents
#' "no evidence", dark red represents "definitely incorrect", light red represents "likely incorrect", and grey represents
#' "potentially incorrect". For the two categories that come from ORFs without predicted genes, dark blue represents
#' "likely missing" and light blue represents "potentially missing".
#' 
#' If out of \code{x} and \code{y} only \code{x}} is specified and \code{x} is of subclass \code{DataMap}, a genome viewer plot
#' showing how the proteomics data and evolutionary conservation data maps to the central genome is generated.
#'
#' If both \code{x} and \code{y} are specified, each of a different subclass, a genome viewer plot showing how the proteomics
#' data, evolutionary conservation data, and set of genes map to the central genome is generated.
#' 
#' @section Genome viewer plot:
#' 
#' In the genome viewer plot, predicted starts are magenta lines, predicted stops are cyan lines, genome stops are yellow lines,
#' conserved starts are gray lines, and proteomic hits are blue / red / green blocks.
#' 
#' If \code{interactive_GV} is set to FALSE, a static genome viewer plot is generated. If \code{interactive_GV} is set to TRUE,
#' a genome viewer plot that can be interacted with using the \code{locator} is generated. In order to interact with the plot,
#' the user needs to click on the graphics window one or more times and then terminate the locator. One click will scroll the
#' viewer either to the left or the right (based on which side is closer to the click). Two clicks will zoom the viewer into the
#' horizontal range between the two click points. Three clicks will zoom out 10-fold, and four clicks will zoom out completely to
#' the entire genome. To stop interaction with the locator, click zero times then terminate the locator. Depending on the
#' graphical device, terminating the locator can either done by pressing the 'Finish' / 'Stop' button, hitting the 'Esc' key, or
#' right-clicking the graphics device.
#' 
#' By default, the genome viewer will cover all positions in the genome. If instead both \code{rangeStart_GV} and
#' \code{rangeEnd_GV} are validly specified, the genome viewer will only span the range of genomic positions between those two
#' values (at least initially).
#' 
#' @return Invisibly returns the input object \code{x}
#'
#' @seealso \code{\link{Assessment-class}}, \code{\link{locator}}
#' 
#' @examples
#'
#' currMapObj <- readRDS(system.file("extdata",
#'                                   "MGAS5005_PreSaved_DataMapObj.rds",
#'                                   package = "AssessORF"))
#' 
#' currResObj <- readRDS(system.file("extdata",
#'                                   "MGAS5005_PreSaved_ResultsObj_Prodigal.rds",
#'                                   package = "AssessORF"))
#'
#' plot(currMapObj)
#' 
#' plot(currResObj)
#' 
#' plot(currMapObj, currResObj)
#' 
#' plot(currResObj, currMapObj)
#'
plot.Assessment <- function(x, y = NULL,
                            minConCovRatio_GV = 0.8, interactive_GV = TRUE,
                            rangeStart_GV = NA_integer_, rangeEnd_GV = NA_integer_, ...) {
  
  if (class(x)[1] != "Assessment") {
    stop("'x' must be an object of class 'Assessment'.")
  }
  
  if (!is.null(y)) {
    if (class(y)[1] != "Assessment") {
      stop("'y' must be an object of class 'Assessment'.")
    }
    
    if (class(x)[2] == "DataMap") {
      cMapObj <- x
      
      if (class(y)[2] != "Results") {
        stop("'y' must be an object of subclass 'Results' when x is an object of subclass 'DataMap'.")
      }
      
      cResObj <- y
    } else if (class(x)[2] == "Results") {
      cResObj <- x
      
      if (class(y)[2] != "DataMap") {
        stop("'y' must be an object of subclass 'DataMap' when x is an object of subclass 'Results'.")
      }
      
      cMapObj <- y
    } else {
      stop("'x' is an unrecognized object of class '", class(x)[2], "'.")
    }
    
    PlotAssessmentMapping(mapObj = cMapObj, resObj = cResObj, minConCovStart = minConCovRatio_GV,
                          interactivePlot = interactive_GV,
                          initialPos1 = rangeStart_GV, initialPos2 = rangeEnd_GV)
    
  } else if (class(x)[2] == "DataMap") {
    
    PlotAssessmentMapping(mapObj = x, resObj = NULL, minConCovStart = minConCovRatio_GV,
                          interactivePlot = interactive_GV,
                          initialPos1 = rangeStart_GV, initialPos2 = rangeEnd_GV)
    
  } else if (class(x)[2] == "Results") {
    geneSource <- x$GeneSource
    
    part2 <- "Gene Category Assignments"
    
    if (geneSource != "") {
      part2 <- paste(geneSource, part2)
    }
    
    ## --------------------------------------------------------------------------------------------------------------- ##
    
    speciesName <- x$Species
    strainID <- x$StrainID
    
    if ((speciesName != "") && (strainID != "")) {
      plotTitle <- bquote(italic(.(speciesName))~.(strainID)~.(part2))
    } else if ((speciesName != "")) {
      plotTitle <- bquote(italic(.(speciesName))~.(part2))
    } else if ((strainID != "")) {
      plotTitle <- bquote(~.(strainID)~.(part2))
    } else {
      plotTitle <- part2
    }
    
    ## --------------------------------------------------------------------------------------------------------------- ##
    
    ## Get the number of genes in each category
    catSumTable <- table(x$CategoryAssignments)
    
    nonCollapseIdxs <- which((names(catSumTable) != "Y CS! PE+") &
                               (names(catSumTable) != "Y CS! PE-") &
                               (names(catSumTable) != "Y CS> PE+") &
                               (names(catSumTable) != "Y CS> PE-") &
                               (names(catSumTable) != "Y CS< PE+") &
                               (names(catSumTable) != "Y CS< PE-"))
    
    allCatSums <- integer(11L)
    
    names(allCatSums) <-  c("Y CS+ PE+", "Y CS+ PE-", "Y CS- PE+", "Y CS- PE-",
                            "Y CS< PE!", "Y CS- PE!", "Y CS! PE\u00B1",
                            "Y CS> PE\u00B1", "Y CS< PE\u00B1",
                            "N CS< PE+", "N CS- PE+")
    
    allCatSums[names(catSumTable)[nonCollapseIdxs]] <- catSumTable[nonCollapseIdxs]
    
    allCatSums["Y CS! PE\u00B1"] <- sum(catSumTable[which((names(catSumTable) == "Y CS! PE+") |
                                                            (names(catSumTable) == "Y CS! PE-"))], na.rm = TRUE)
    
    allCatSums["Y CS> PE\u00B1"] <- sum(catSumTable[which((names(catSumTable) == "Y CS> PE+") |
                                                            (names(catSumTable) == "Y CS> PE-"))], na.rm = TRUE)
    
    allCatSums["Y CS< PE\u00B1"] <- sum(catSumTable[which((names(catSumTable) == "Y CS< PE+") |
                                                            (names(catSumTable) == "Y CS< PE-"))], na.rm = TRUE)
    
    allCatSums["N CS< PE+"] <- sum(nrow(x$'N_CS<_PE+_ORFs'), na.rm = TRUE)
    allCatSums["N CS- PE+"] <- sum(nrow(x$'N_CS-_PE+_ORFs'), na.rm = TRUE)
    
    ## --------------------------------------------------------------------------------------------------------------- ##
    
    tenFactor <- 10 ^ floor(log10(max(allCatSums)))
    digitMult <- max(allCatSums) / tenFactor
    lastDigit <- floor(digitMult)
    
    if ((digitMult - lastDigit) <= 0.5) {
      currYMax <- (lastDigit + 0.5) * tenFactor
    } else {
      currYMax <- (lastDigit + 1) * tenFactor
    }
    
    ## --------------------------------------------------------------------------------------------------------------- ##
    
    catColors <- c("darkgreen", "lightgreen", "lightgreen", "white",
                   "darkred", "indianred1", "indianred1",
                   "gray", "gray",
                   "darkblue", "lightblue", "darkred")
    
    par(mfrow=c(1,1), mar = c(6, 4, 6, 2))
    
    barplot(allCatSums, las = 2, ylim = c(0, currYMax),
            col = catColors, col.lab = "black",
            family = "mono")
    
    title(main = plotTitle, line = 5)
    
    legend(x = ceiling(par("usr")[1]), y = 1.2 * currYMax, bty = "n", ncol = 3, xpd = TRUE,
           fill = c("darkgreen", "lightgreen", "white", "darkred", "indianred1", "gray", "darkblue", "lightblue"),
           legend = c("Definitely Correct", "Likely Correct", "No Evidence",
                      "Definitely Incorrect", "Likely Incorrect", "Potentially Incorrect",
                      "Likely Missing", "Potentially Missing"))
    
  } else {
    stop("'x' is an object of unrecognized subclass '", class(x)[2], "'.")
  }
  
  invisible(x)
}

#' @export
#' @importFrom graphics mosaicplot
#' @importFrom stats quantile
#' 
#' @method mosaicplot Assessment
#'
#' @title Plot Genes by Category and Length
#' @description The \code{mosaicplot} method for \code{Assessment} object
#'
#' @param x An object of class \code{Assessment} and subclass \code{Results}.
#'
#' @param ... Further \code{mosaicplot} parameters.
#'
#' @details
#' \code{mosaicplot.Assessment} plots all the genes in the given \code{Results} object by category and length. This set of genes
#' includes both the supplied predicted genes as well as open reading frames with proteomics evidence but no predicted start.
#' 
#' The set of genes are separated into ten quantile bins based on the length of the gene/open reading frame. The genes are then
#' plotted by length bin and category in a mosaic format, with each column representing a length bin and each row/block
#' representing a category.
#' 
#' @return Invisibly returns the input object \code{x}
#'
#' @seealso \code{\link{Assessment-class}}
#' 
#' @examples
#'
#' mosaicplot(readRDS(system.file("extdata",
#'                                "MGAS5005_PreSaved_ResultsObj_Prodigal.rds",
#'                                package = "AssessORF")))
#'
mosaicplot.Assessment <- function(x, ...) {
  if (class(x)[1] != "Assessment") {
    stop("'x' must be an object of class 'Assessment'.")
  }
  
  if (class(x)[2] == "Results") {
    geneSource <- x$GeneSource
    
    part2 <- "Gene Category Assignments by Nucleotide Length"
    
    if (geneSource != "") {
      part2 <- paste(geneSource, part2)
    }
    
    ## --------------------------------------------------------------------------------------------------------------- ##
    
    speciesName <- x$Species
    strainID <- x$StrainID
    
    if ((speciesName != "") && (strainID != "")) {
      plotTitle <- bquote(italic(.(speciesName))~.(strainID)~.(part2))
    } else if ((speciesName != "")) {
      plotTitle <- bquote(italic(.(speciesName))~.(part2))
    } else if ((strainID != "")) {
      plotTitle <- bquote(~.(strainID)~.(part2))
    } else {
      plotTitle <- part2
    }
    
    ## --------------------------------------------------------------------------------------------------------------- ##
    
    geneLengths <- x$GeneRightPos - x$GeneLeftPos + 1
    names(geneLengths) <- x$CategoryAssignments
    
    if (!is.null(nrow(x$'N_CS<_PE+_ORFs'))){
      cat1BLens <- x$'N_CS<_PE+_ORFs'[, 'Length']
      names(cat1BLens) <- rep("N CS< PE+", nrow(x$'N_CS<_PE+_ORFs'))
      
      geneLengths <- c(geneLengths, cat1BLens)
    }
    
    if (!is.null(nrow(x$'N_CS-_PE+_ORFs'))){
      cat1ALens <- x$'N_CS-_PE+_ORFs'[, 'Length']
      names(cat1ALens) <- rep("N CS- PE+", nrow(x$'N_CS-_PE+_ORFs'))
      
      geneLengths <- c(geneLengths, cat1ALens)
    }
    
    quantBins <- quantile(geneLengths, seq(0, 1, 0.1))
    
    catByLenTable <- table(cut(geneLengths, quantBins, dig.lab = 5),
                           names(geneLengths))
    
    catCodeNames <- c("Y CS+ PE+", "Y CS+ PE-", "Y CS- PE+", "Y CS- PE-",
                      "Y CS< PE!", "Y CS- PE!", "Y CS! PE+", "Y CS! PE-",
                      "Y CS> PE+", "Y CS> PE-", "Y CS< PE+", "Y CS< PE-",
                      "N CS< PE+", "N CS- PE+")
    
    catColors <- c("darkgreen", "lightgreen", "lightgreen", "white",
                   "darkred", "indianred1", "indianred1", "indianred1",
                   "gray", "gray", "gray", "gray",
                   "darkblue", "lightblue", "darkred")
    
    catCodeNames <- catCodeNames[catCodeNames %in% colnames(catByLenTable)]
    
    usedCats <- which(colSums(catByLenTable[, catCodeNames], na.rm = TRUE) != 0)
    
    catCodeNames <- catCodeNames[usedCats]
    
    catColors <- catColors[usedCats]
    
    par(col.lab = "black", family = "mono")
    plot(catByLenTable[, catCodeNames], main = plotTitle, col = catColors, las = 1)
    
  } else {
    stop("'x' must be of subclass 'Results'.")
  }
  
  invisible(x)
}
