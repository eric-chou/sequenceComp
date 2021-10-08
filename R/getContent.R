library(BSgenome)
library(GenomicFeatures)
library(dplyr)

#' @name getBaseContent
#' @title **Get Base Content**
#' @description Computes the content in a sequence for either a base or a vector of individual bases. 
#' @param seq (required): Sequence in which the content is to be calculated
#' @param bases (required): Character or vector of characters for bases that are to be used for calculating their content proportion in the sequence provided. If a vector of characters is provided, this function finds the proportion of the sequence that is any one of the bases, not necessarily in a contiguous manner. If one aims to find content of a contiguous pattern, they should use the function `getPatternContent`
#' @param pct (default FALSE): boolean for if the user wishes to return the output as a percentage (TRUE) or proportion (FALSE, default).
#' @param baseOnly (default TRUE): boolean for if the user wants to calculate the content on using the base pair characters (ex. A, T, C, G for DNA) in the sequence, or from all characters in the sequence.
getBaseContent <- function(seq, bases, pct=FALSE, baseOnly=TRUE) {
	bases = alphabetFrequences(seq, baseOnly=baseOnly)

}

getPatternContent <- function(seq, pattern, pct=FALSE, baseOnly=TRUE){

}
