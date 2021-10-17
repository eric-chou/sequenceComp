library(BSgenome)
library(dplyr)

#' @name windowBaseContent
#' @title **Window Base Content**
#' @description Given an individual sequence, this function splits the provided sequence into windows of a user-specified size and then computes the content in eqch window for either a base or a vector of individual bases.
#' @param seq (required): Biostring sequence that is to be split into windows, in which the base content is to be calculated for each window.
#' @param windowSize (required): The size of each window, in base pairs.
#' @param bases (required): Character or vector of characters for bases that are to be used for calculating their content proportion in each created window. If a vector of characters is provided, this function finds the proportion of the sequence that is any one of the bases (i.e. not necessarily in a contiguous manner).
#' @param pct (default FALSE): boolean for if the user wishes to return the output as a percentage (TRUE) or proportion (FALSE, default).
#' @return A vector of proportions or percentages for base content in each window of the user-specified size for a sequence.
windowBaseContent <- function(seq, windowSize, bases, pct=FALSE){
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

	return(getBaseContent(seq.set, bases, baseOnly=FALSE, pct=FALSE)) # baseOnly is forced to be false since the full window size must be considered the denominator regardless of if they are base pair characteres
}

# windowPatternContent <- function(seq, bases,  out=c("histogram", "boxplot", "vector")){

# }
