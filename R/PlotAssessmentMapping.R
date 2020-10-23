PlotAssessmentMapping <- function(mapObj, resObj = NULL, minConCovStart = 0.8, interactivePlot = TRUE,
                                  initialPos1 = NA_integer_, initialPos2 = NA_integer_) {
  
  ## Get the length of the genome.
  genomeLength <- mapObj$GenomeLength
  
  ## ------------------------------------------------------------------------ ##
  
  ## Error check 'minConStart' and 'interactivePlot'.
  
  if ((!is.numeric(minConCovStart)) || (anyNA(minConCovStart))) {
    stop("Must specify a valid real number for the start conservation to coverage ratio threshold.")
  }
  
  if (length(minConCovStart) != 1) {
    stop("Must specify exactly one number for the start conservation to coverage ratio threshold.")
  }
  
  if ((minConCovStart <= 0) || (minConCovStart >= 1)) {
    stop("Must specify a number between 0 and 1 (inclusive) for the start conservation to coverage ratio threshold..")
  }
  
  if ((!is.logical(interactivePlot)) || (anyNA(interactivePlot)) || (length(interactivePlot) != 1)) {
    stop("Argument specifying whether the genome viewer should be interactive",
         " must be of type logical, be either TRUE or FALSE, and consist of only 1 element.")
  }
  
  ## ------------------------------------------------------------------------ ##
  
  ## Specifying a pair of genome positions to use as the starting range for the
  ## plot is optional. Error check accordingly.
  
  areBothNA <- as.integer(is.na(initialPos1)) + as.integer(is.na(initialPos2))
  areBothNull <- as.integer(is.null(initialPos1)) + as.integer(is.null(initialPos2))
  
  if ((areBothNA == 2) || (areBothNull == 2)) {
    useInitialRange <- FALSE
  } else {
    if ((areBothNA == 1) || (areBothNull == 1)) {
      stop("Must specifiy 2 positions in order for the the Genome Viewer to zoom into an initial range.")
    }
    
    if ((length(initialPos1) != 1) || (!is.numeric(initialPos1)) || (initialPos1 %% 1 != 0) ||
        (length(initialPos2) != 1) || (!is.numeric(initialPos2)) || (initialPos2 %% 1 != 0)) {
      stop( "Initial range to plot with the Genome Viewer must be given as exactly two separate integer numbers.")
    }
    
    if ((initialPos1 <= 0L) || (initialPos1 > genomeLength) || (initialPos2 <= 0L) || (initialPos2 > genomeLength)) {
      stop("Initial range to plot with the Genome Viewer must be within the bounds of the genome.")
    }
    
    useInitialRange <- TRUE
  }
  
  ## ------------------------------------------------------------------------ ##
  
  ## Get the location of the stops in the genome.
  stops <- mapObj$StopsByFrame

  ## Get information about the proteomics evidence.
  plotProt <- mapObj$HasProteomics
  fwdProt <- mapObj$FwdProtHits
  revProt <- mapObj$RevProtHits
  
  ## Get information about the conserved starts. For each reading direction, 
  ## the ratio of start codon conservation to coverage serves as a measure
  ## of evolutionary conservation of the corresponding start in the genome.
  plotConStarts <- mapObj$HasConservation
  
  fwdCov <- mapObj$FwdCoverage
  fwdConStart <- mapObj$FwdConStarts
  
  revCov <- mapObj$RevCoverage
  revConStart <- mapObj$RevConStarts

  ## If a results object is provided, get information about the location of
  ## predicted genes & plot them. Otherwise, skip plotting predicted genes.
  addPred <- FALSE

  if (!is.null(resObj)) {
    addPred <- TRUE

    predLeft <- resObj$GeneLeftPos
    predRight <- resObj$GeneRightPos
    predStrand <- resObj$GeneStrand
  }
  
  ## ------------------------------------------------------------------------ ##

  ## Get the species name and strain ID that the genome corresponds to.
  speciesName <- mapObj$Species
  strainID <- mapObj$StrainID

  ## Generate the plot title based on that species name and strain ID.
  if ((speciesName != "") && (strainID != "")) {
    plotTitle <- bquote(italic(.(speciesName))~.(strainID)~"Genome Viewer")
  } else if ((speciesName != "")) {
    plotTitle <- bquote(italic(.(speciesName))~"Genome Viewer")
  } else if ((strainID != "")) {
    plotTitle <- bquote(~.(strainID)~"Genome Viewer")
  } else {
    plotTitle <- "Genome Viewer"
  }
  
  ## ------------------------------------------------------------------------ ##

  ## This function determines the color of proteomic hit blocks and calculates
  ## the saturation of that color based on their scores.
  chooseColor <- function(color, protScores, start, end, max=20) {

    stopifnot(length(start)==length(end),
              color > 1 && color < 5,
              color==floor(color),
              all(start <= length(protScores)),
              all(end <= length(protScores)),
              all(start > 0),
              all(end > 0),
              all(end >= start))

    maxes <- integer(length(start))

    for (i in seq_along(start)) {
      maxes[i] <- max(protScores[start[i]:end[i]])
      if (maxes[i] > max)
        maxes[i] <- max
      if (maxes[i] < 0)
        maxes[i] <- 0
    }

    if (color==2) {
      colors <- rgb(1 - maxes/max, 1 - maxes/max, 1)
    } else if (color==3) {
      colors <- rgb(1, 1 - maxes/max, 1 - maxes/max)
    } else if (color==4) {
      colors <- rgb(1 - maxes/max, 1, 1 - maxes/max)
    }

    return(colors)
  }
  
  ## ------------------------------------------------------------------------ ##

  ## This function plots all conserved starts, proteomic hits, genome stops, and
  ## predicted genes within the given range of forward strand positions.
  plotRange <- function(fwdRange) {
    ## Determine the complementary set of reverse strand positions to plot. 
    revRange <- genomeLength - fwdRange + 1L

    ## Lay out the plot grid. There will be one "plot" for each strand. The
    ## forward strand plot needs to be above the reverse strand plot but 
    ## generated second so its x-values will be used with the locator.
    layout(matrix(c(2, 1)))

    ## Loop through the six reading frames, going backwards. This means that the
    ## reverse strand reading frames will be plotted first and that the third
    ## frame of each strand will be the first one plotted for that strand.
    for (frame in c(6, 5, 4, 3, 2, 1)) {
      
      ## For the third reading frame of each strand, generate the plot for that
      ## strand and then plot the predicted genes across the 3 reading frames.
      if (frame == 3L) {
        ## Generate the plot area for the forward strand (frames 1, 2, & 3).
        ## Since the forward strand plot is the top plot, include the title.
        par(mar=c(2.1, 4.1, 2.1, 0.5))

        plot(NA,
             xlim = range(fwdRange),
             ylim = c(0, 3),
             xlab = "",
             ylab = "Reading frame (forward)",
             main = plotTitle,
             xaxs = "i",
             yaxs = "i",
             yaxt = "n")
        
        ## Plot the y-axis, which will describe the reading frames.
        axis(2, 0:2 + 0.5, c(3, 2, 1))

        if (addPred) {
          ## If a results object was provided, plot the forward strand genes.
          predFwd <- which(predStrand == "+")
          abline(v = predLeft[predFwd], col = "magenta", lwd = 1.25)
          abline(v = predRight[predFwd], col = "cyan")
        }
        
        ## Add lines separating the reading frames.
        abline(h = 0:3)
        
      } else if (frame == 6L) {
        ## Generate the plot area for the reverse strand (frames 4, 5, & 6).
        ## Since the forward strand plot is the top plot, include the title.
        par(mar=c(4.1, 4.1, 0.1, 0.5))

        plot(NA,
             xlim = rev(range(revRange)),
             ylim = c(0, 3),
             xlab = "Position",
             ylab = "Reading frame (reverse)",
             xaxs = "i",
             yaxs = "i",
             yaxt = "n")
        
        ## Plot the y-axis, which will describe the reading frames.
        axis(2, 0:2 + 0.5, c(6, 5, 4))

        if (addPred) {
          ## If a results object was provided, plot the reverse strand genes.
          predRev <- which(predStrand == "-")
          abline(v = genomeLength - predLeft[predRev] + 1, col="cyan")
          abline(v = genomeLength - predRight[predRev] + 1, col="magenta", lwd = 1.25)
        }
        
        ## Add lines separating the reading frames.
        abline(h = 0:3)
      }
      
      ## Next, for all reading frames, plot the proteomic hit blocks, genome
      ## stops, and conserved starts.
      if (frame <= 3) {
        ## This is a forward strand frame. Use the corresponding data.
        range <- fwdRange
        
        ## Determine the top and bottom y-values for current reading frame.
        yBot <- 3 - frame
        yTop <- 4 - frame

        ## If there is proteomics evidence, plot it.
        if (plotProt) {
          ## Find all the proteomic hit blocks in the current reading frame and
          ## then determine the start and end position for each one.
          protBlocks <- fwdProt[[frame]][[seq]][range] > 0
  
          protBlockBegs <- which(diff(c(0, protBlocks)) == 1)
          protBlockEnds <- which(diff(c(protBlocks, 0)) == -1)
  
          ## Plot the proteomic hit blocks, determining the color & saturation
          ## based on the current reading frame and proteomic hit strengths.
          if (length(protBlockBegs) > 0) {
            rect(range[protBlockBegs] - 0.5, yBot, range[protBlockEnds] + 0.5, yTop,
                 col = chooseColor(frame + 1, fwdProt[[frame]][[seq]][range],
                                   protBlockBegs, protBlockEnds),
                 border=NA)
          }
        }

        ## Plot the genome stops for the current reading frame.
        if (length(stops[[frame]]) > 0) {
          rect(stops[[frame]], yBot, stops[[frame]] + 2, yTop, col = "yellow", border = NA)
        }

        ## If there is conservation evidence, plot the conserved starts.
        if (plotConStarts) {
          ## Determine which of all the possible conserved starts in the current
          ## reading pass the given conservation to coverage ratio threshold.
          conStartIdxs <- which((fwdConStart[range] / fwdCov[range] >= minConCovStart) &
                                  ((range - frame) %% 3) == 0)
          
          ## Plot those conserved starts, picking the shade of gray based on the
          ## strength of conservation (= conservation to coverage ratio).
          if (length(conStartIdxs) > 0) {
            segments(range[conStartIdxs], yBot, range[conStartIdxs], yTop,
                     col = gray(1 - fwdConStart[range[conStartIdxs]] / fwdCov[range[conStartIdxs]],
                                alpha = 0.5))
          }
        }
      } else {
        ## This is a reverse strand frame. Use the corresponding data.
        range <- revRange
        
        ## Determine the top and bottom y-values for current reading frame.
        yBot <- 6 - frame
        yTop <- 7 - frame

        ## If there is proteomics evidence, plot it.
        if (plotProt) {
          ## Find all the proteomic hit blocks in the current reading frame and
          ## then determine the start and end position for each one.
          protBlocks <- revProt[[frame - 3]][[seq]][range] > 0
  
          protBlockBegs <- which(diff(c(0, protBlocks)) == 1)
          protBlockEnds <- which(diff(c(protBlocks, 0)) == -1)
  
          ## Plot the proteomic hit blocks, determining the color & saturation
          ## based on the current reading frame and proteomic hit strengths.
          if (length(protBlockBegs) > 0){
            rect(range[protBlockBegs] - 0.5, yBot, range[protBlockEnds] + 0.5, yTop,
                 col = chooseColor(frame - 2, revProt[[frame - 3]][[seq]][range],
                                   protBlockBegs, protBlockEnds),
                 border=NA)
          }
        }

        ## Plot the genome stops for the current reading frame.
        if (length(stops[[frame]]) > 0){
          rect(stops[[frame]], yBot, stops[[frame]] + 2, yTop, col="yellow", border = NA)
        }

        ## If there is conservation evidence, plot the conserved starts.
        if (plotConStarts) {
          ## Determine which of all the possible conserved starts in the current
          ## reading pass the given conservation to coverage ratio threshold.
          conStartIdxs <- which((revConStart[range] / revCov[range] >= minConCovStart) &
                                  ((range - frame - 3) %% 3) == 0)
          
          ## Plot those conserved starts, picking the shade of gray based on the
          ## strength of conservation (= conservation to coverage ratio).
          if (length(conStartIdxs) > 0) {
            segments(range[conStartIdxs], yBot, range[conStartIdxs], yTop,
                     col = gray(1 - revConStart[range[conStartIdxs]] / revCov[range[conStartIdxs]],
                                alpha = 0.5))
          }
        }
      }
    }
  }
  
  ## ------------------------------------------------------------------------ ##

  ## This function takes in a pair of x-coordinates / genome positions, makes
  ## sure the combination of the two values are valid, and then returns the
  ## range of genome positions between the them.
  checkRange <- function(xCoords) {
    if ((all(xCoords < 0)) || (all(xCoords > genomeLength))) {
      stop("x-axis coordinates are out of bounds.")
    }
    
    if (xCoords[2] < xCoords[1]) {
      intSwap <- xCoords[1]
      xCoords[1] <- xCoords[2]
      xCoords[2] <- intSwap
    }
    
    if (xCoords[1] <= 0) {
      xCoords[1] <- 1L
    }
    
    if (xCoords[2] > genomeLength){
      xCoords[2] <- genomeLength
    }
    
    xCoords <- floor(xCoords)
    
    return(xCoords[1]:xCoords[2])
  }
  
  ## ------------------------------------------------------------------------ ##
  
  ## This variable specifies which sequence of the genome to use. As of now,
  ## AssessORF only works with single-chromosome genomes.
  seq <- 1L
  
  ## ------------------------------------------------------------------------ ##

  ## If not plotting an interactive genome viewer, generate a single plot. If
  ## the user provided an initial range of genomic positions, plot that region
  ## of the genome. Otherwise, plot the whole genome.
  if (!interactivePlot) {
    if (useInitialRange) {
      plotRange(checkRange(c(initialPos1, initialPos2)))
    } else {
      plotRange(seq_len(genomeLength))
    }
    
    ## Reset the graphical device's layout once the single, static plot has been
    ## generated and then exit the function.
    par(mfrow = c(1,1))
    
    return(invisible(mapObj))
  }
  
  ## ------------------------------------------------------------------------ ##

  ## If plotting an interactive genome viewer, start by printing out
  ## instructions on how to interact with the genome viewer.
  cat("How to Interact With the Genome Viewer Via the Locator:\n",
      "Click the graphics window one or more times and then terminate the locator.\n",
      "1 click : Scroll the viewer to the left or the right.\n",
      "2 clicks: Zoom into the horizontal range between the two click points.\n",
      "3 clicks: Zoom out 10-fold.\n",
      "4 clicks: Zoom out completely and view the entire genome.\n",
      "To stop interaction, click zero times then terminate the locator.\n",
      "Depending on the graphical device, terminating the locator can either be done ",
      "by pressing the 'Finish' / 'Stop' button, hitting the 'Esc' key, or right-clicking.\n",
      sep = "")
  
  ## This vector will hold the values returned by the locator. It will be set
  ## such that the whole genome is plotted when going through the interaction
  ## options below for the first time. However, if the user provided an initial
  ## range, the vector is set such that those positions will be plotted when
  ## going through the interaction options for the first time instead.
  locValues <- list(x=integer(4))
  
  if (useInitialRange) {
    locValues$x <- c(initialPos1, initialPos2)
  }

  ## Repeat through the interaction options until the user asks to stop.
  repeat {
    if (length(locValues$x) == 0) {
      ## Stop interaction if the user clicks no times.
      break
    } else if (length(locValues$x) == 1) {
      ## Scroll the viewer left or right with one click. The direction is
      ## determined by where the click is in relation the plot mid-point.
      temp <- list(x = par("usr")[c(1, 2)])
      mid <- (temp$x[1] + temp$x[2])/2
      
      if (locValues$x - mid > 0) {
        ## Move right on the forward strand.
        temp$x <- temp$x + (temp$x[2] - temp$x[1])/2
      } else {
        ## Move left on the forward strand.
        temp$x <- temp$x - (temp$x[2] - temp$x[1])/2
      }
      
      plotRange(checkRange(c(temp$x[1], temp$x[2])))
      
    } else if (length(locValues$x) == 2) {
      ## Plot the given initial range OR with two clicks, zoom to the range of
      ## positions between those clicks.
      plotRange(checkRange(c(locValues$x[1], locValues$x[2])))
      
    } else if (length(locValues$x) == 3) {
      ## Zoom out 10-fold with three clicks.
      temp <- par("usr")[c(1, 2)]
      dtemp <- diff(temp)
      temp <- c(temp[1] - dtemp*4.5,
                temp[2] + dtemp*4.5)
      plotRange(checkRange(temp))
    } else {
      ## Plot the whole genome with four or more clicks. This also triggers with
      ## the first genome viewer plot if no initial range is given.
      plotRange(seq_len(genomeLength))
    }
    
    ## Set up (first iteration only) or reset the locator.
    locValues <- locator()
  }
  
  ## Reset the graphical device's layout once plotting is done and then exit.
  par(mfrow = c(1,1))
  
  return(invisible(mapObj))
}
