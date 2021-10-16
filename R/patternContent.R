library(BSgenome)
library(GenomicFeatures)
library(dplyr)

#' @name getPatternContent
#' @title **Get Pattern Content**
#' @description Computes the observed content in a sequence for a given pattern of bases, as well as the expected content for that pattern.
#' @param seq (required): Biostring sequence in which the content is to be calculated.
#' @param pattern (required): String for a contiguous pattern of bases that are to be used for calculating their observed and expected occurrences in the sequence provided. 
#' @param baseOnly (default TRUE): boolean for if the user wants to calculate the content on using the base pair characters (ex. A, T, C, G for DNA) in the sequence, or from all characters in the sequence.
#' @param ratio (default FALSE): boolean for if the user wishes to return the output as a named list for observed and expected counts (FALSE, default) or a ratio of observed over expected (TRUE).
#' @return Unlike getBaseContent, this provides not a proportion or percentage, but either a list (default) or the ratio of observed to expected.
getPatternContent <- function(seq, pattern, baseOnly=TRUE, ratio=FALSE){
	# observed count 
	n <- countPattern(pattern, seq)
	freq <- alphabetFrequency(seq, baseOnly=baseOnly)
	ntotBases <- ifelse(baseOnly, sum(freq[1:4]), sum(freq))
	freq <- freq / ntotBases

	# calculate expected - multiply frequencies for each individual base in the given pattern string, then multiply by total bases
	basesVec <- strsplit(pattern, "")[[1]]
	expected <- prod(freq[basesVec]) * ntotBases

	outList <- list(observedCount = n, expectedCount = expected)
	ifelse(ratio, return(n / expected), return(outList))
}

#' @name comparePatternContent
#' @title **Compare Pattern Content**
#' @description Computes the ratio of observed to expected occurrences for a contiguous pattern of baes for each sequence in a list of sequences. These ratios are then outputted in either ascending or descending order of content.
#' @param seqList (required): List or S4 class of Biostring sequences, each for which the observed to expected ratio is to be calculated. It is recommended that this list be named for readability of output. Otherwise, the output will provide the corresponding index occupied by the sequence in the original list that was provided alongside its base content calculation.
#' @param pattern (required): String for a contiguous pattern of bases that are to be used for calculating their observed to expected ratios in each sequence provided.
#' @param baseOnly (default TRUE): boolean for if the user wants to calculate the content on using the base pair characters (ex. A, T, C, G for DNA) in the sequence, or from all characters in the sequence.
#' @param desc (default TRUE): boolean for if the user wishes to return the output in descending order of observed to expected ratio (TRUE, default) or in ascending order (FALSE).
#' @return Named vector of observed to expected ratios of pattern occurrences for each provided sequence, in either descending or ascending order of ratios.
comparePatternContent <- function(seqList, pattern, baseOnly=TRUE, desc=TRUE){
	if(is.null(names(seqList))){
		warning("The sequence list provided is unnamed, it is recommended to use a named list to label each sequence for more readable output")
	}
	ratioVector <- c()
	idxVector <- c()
	for(i in 1:length(seqList)){
		ratio <- getPatternContent(seqList[[i]], pattern=pattern, baseOnly=baseOnly, ratio=TRUE)
		ratioVector <- c(ratioVector, ratio)
		idx <- ifelse(is.null(names(seqList)), i, names(seqList)[i])
		idxVector <- c(idxVector, idx)
	}

	names(ratioVector) <- idxVector
	ratioVector <- sort(ratioVector, decreasing=desc)
	return(ratioVector)
}
