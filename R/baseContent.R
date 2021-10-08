library(BSgenome)
library(GenomicFeatures)
library(dplyr)

#' @name getBaseContent
#' @title **Get Base Content**
#' @description Computes the content in a sequence for either a base or a vector of individual bases.
#' @param seq (required): Biostring sequence in which the content is to be calculated
#' @param bases (required): Character or vector of characters for bases that are to be used for calculating their content proportion in the sequence provided. If a vector of characters is provided, this function finds the proportion of the sequence that is any one of the bases, not necessarily in a contiguous manner. If one aims to find content of a contiguous pattern, they should use the function `getPatternContent`.
#' @param baseOnly (default TRUE): boolean for if the user wants to calculate the content on using the base pair characters (ex. A, T, C, G for DNA) in the sequence, or from all characters in the sequence.
#' @param pct (default FALSE): boolean for if the user wishes to return the output as a percentage (TRUE) or proportion (FALSE, default).
#' @return Proportion or percentage for the base content of the provided base(s).
getBaseContent <- function(seq, bases, baseOnly=TRUE, pct=FALSE) {
	seq <- alphabetFrequency(seq, baseOnly=baseOnly)
	ntotBases <- ifelse(baseOnly, sum(seq[1:4]), sum(seq))
	baseFreq <- sum(seq[bases] / ntotBases)
	ifelse(pct, return(baseFreq*100), return(baseFreq))
}

#' @name compareBaseContent
#' @title **Compare Base Content**
#' @description Computes the content for each sequence in a list of sequences for either a base or a vector of individual bases, and then outputs them in either ascending or descending order of content.
#' @param seqList (required): List of Biostring sequences, each for which the content is to be calculated. It is recommended that this list be named for readability of output. Otherwise, the output will provide the corresponding index occupied by the sequence in the original list that was provided alongside its base content calculation.
#' @param bases (required): Character or vector of characters for bases that are to be used for calculating their content proportion in the sequence provided. If a vector of characters is provided, this function finds the proportion of the sequence that is any one of the bases, not necessarily in a contiguous manner. If one aims to find content of a contiguous pattern, they should use the function `getPatternContent`.
#' @param baseOnly (default TRUE): boolean for if the user wants to calculate the content on using the base pair characters (ex. A, T, C, G for DNA) in the sequence, or from all characters in the sequence.
#' @param pct (default FALSE): boolean for if the user wishes to return the output as a percentage (TRUE) or proportion (FALSE, default).
#' @param desc (default TRUE): boolean for if the user wishes to return the output in descending order of base content (TRUE, default) or in ascending order (FALSE).
#' @return Named vector of proportions or percentages for base content of each provided sequence, in either descending or ascending order of base content.
compareBaseContent <- function(seqList, bases, baseOnly=TRUE, pct=FALSE, desc=TRUE) {
	if(is.null(names(seqList))){
		warning("The sequence list provided is unnamed, it is recommended to use a named list to label each sequence for more readable output")
	}
	contentVector <- c()
	idxVector <- c()
	for(i in 1:length(seqList)){
		content <- getBaseContent(seqList[[i]], bases=bases, baseOnly=baseOnly, pct=pct)
		contentVector <- c(contentVector, content)
		idx <- ifelse(is.null(names(seqList)), i, names(seqList)[i])
		idxVector <- c(idxVector, idx)
	}

	names(contentVector) <- idxVector
	contentVector <- sort(contentVector, decreasing=desc)
	return(contentVector)
}
