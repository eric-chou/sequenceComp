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

# windowPatternContent <- function(seq, bases, baseOnly=TRUE, out=c("histogram", "boxplot", "vector")){

# }
