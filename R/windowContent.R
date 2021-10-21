library(BSgenome)
library(Biostrings)
library(dplyr)

#' @name windowBaseContent
#' @title **Window Base Content**
#' @description Given an individual sequence, this function splits the provided sequence into windows of a user-specified size and then computes the content in each window for either a base or a vector of individual bases.
#' @param seq (required): Biostring sequence that is to be split into windows, in which the base content is to be calculated for each window.
#' @param windowSize (required): The size of each window, in base pairs.
#' @param bases (required): Character or vector of characters for bases that are to be used for calculating their content proportion in each created window from the sequence. If a vector of multiple characters is provided, this function finds the proportion of the window that is any one of the bases (i.e. not necessarily in a contiguous manner).
#' @param zero.rm (default FALSE): boolean for if the user wishes to remove all proportions or percentages in the returned vector that are equal to zero.
#' @param pct (default FALSE): boolean for if the user wishes to return the output as a percentage (TRUE) or proportion (FALSE, default).
#' @return A vector of proportions or percentages for base content in each window of the user-specified size for a sequence.
#' @export
#'
windowBaseContent <- function(seq, windowSize, bases, zero.rm=FALSE, pct=FALSE){
	inputClass <- class(seq)[1]
	# find start indices for each window, removing the last index so that all windows are equal size
	indices <- seq(1, length(seq), by=windowSize)
	indices <- indices[-length(indices)]
	# create StringSet - check if DNA or RNA and create StringSet accordingly
	if(inputClass == "DNAString"){
		seq.set <- DNAStringSet(seq, start=indices, end=indices + (windowSize-1))
	}
	else if(inputClass == "RNAString"){
		seq.set <- RNAStringSet(seq, start=indices, end=indices + (windowSize-1))
	}
	else{
		stop("Sequence inputted is of class: ", inputClass, "\nPlease provide a DNAString or RNAString for the input sequence.")
	}

	out <- getBaseContent(seq.set, bases, baseOnly=FALSE, pct=FALSE) # baseOnly is forced to be false since the full window size must be considered the denominator regardless of if they are base pair characteres
	if(zero.rm) out <- out[out>0]
	return(out)
}

#' @name windowPatternContent
#' @title **Window Pattern Content**
#' @description Given an individual sequence, this function splits the provided sequence into windows of a user-specified size and then computes the content in each window for a contiguous pattern of bases.
#' @param seq (required): Biostring sequence that is to be split into windows, in which the base content is to be calculated for each window.
#' @param windowSize (required): The size of each window, in base pairs.
#' @param pattern (required): String for a contiguous pattern of bases that are to be used for calculating their observed and expected occurrences in each created window from the sequence.
#' @param na.rm (default FALSE): boolean for if the user wishes to remove all ratios in the returned vector that are NaN.
#' @param log (default FALSE): boolean for if the user wishes to natural log transform the returned vector.
#' @return A vector of observed to expected ratios for base content in each window of the user-specified size for a sequence.
#' @export
#'
windowPatternContent <- function(seq, windowSize, pattern, na.rm=FALSE, log=FALSE){
	inputClass <- class(seq)[1]
	# find start indices for each window, removing the last index so that all windows are equal size
	indices <- seq(1, length(seq), by=windowSize)
	indices <- indices[-length(indices)]
	# create StringSet - check if DNA or RNA and create StringSet accordingly
	if(inputClass == "DNAString"){
		seq.set <- DNAStringSet(seq, start=indices, end=indices + (windowSize-1))
	}
	else if(inputClass == "RNAString"){
		seq.set <- RNAStringSet(seq, start=indices, end=indices + (windowSize-1))
	}
	else{
		stop("Sequence inputted is of class: ", inputClass, "\nPlease provide a DNAString or RNAString for the input sequence.")
	}

	out <- getPatternContent(seq.set, pattern, baseOnly=FALSE, log=log) # baseOnly is forced to be false since the full window size must be considered the denominator regardless of if they are base pair characteres
	if(na.rm) out <- out[!is.nan(out)]
	return(out)
}
