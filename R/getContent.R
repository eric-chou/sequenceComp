library(BSgenome)
library(GenomicFeatures)
library(dplyr)

#' @name getBaseContent
#' @title **Get Base Content**
#' @description Computes the content in a sequence for either a base or a vector of individual bases.
#' @param seq (required): Sequence in which the content is to be calculated
#' @param bases (required): Character or vector of characters for bases that are to be used for calculating their content proportion in the sequence provided. If a vector of characters is provided, this function finds the proportion of the sequence that is any one of the bases, not necessarily in a contiguous manner. If one aims to find content of a contiguous pattern, they should use the function `getPatternContent`
#' @param baseOnly (default TRUE): boolean for if the user wants to calculate the content on using the base pair characters (ex. A, T, C, G for DNA) in the sequence, or from all characters in the sequence.
#' @param pct (default FALSE): boolean for if the user wishes to return the output as a percentage (TRUE) or proportion (FALSE, default).
#' @return Proportion or percentage for the base content of the provided base(s)
getBaseContent <- function(seq, bases, baseOnly=TRUE, pct=FALSE) {
	seq <- alphabetFrequency(seq, baseOnly=baseOnly)
	ntotBases <- ifelse(baseOnly, sum(seq[1:4]), sum(seq))
	baseFreq <- sum(seq[bases] / ntotBases)
	ifelse(pct, return(baseFreq*100), return(baseFreq))
}

getPatternContent <- function(seq, pattern, pct=FALSE, baseOnly=TRUE){

}
