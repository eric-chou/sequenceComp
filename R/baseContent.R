library(BSgenome)
library(Biostrings)
library(dplyr)

#' @name getBaseContent
#' @title **Get Base Content**
#' @description This function computes the content in a sequence or a StringSet of multiple sequences for either a base or a vector of individual bases.
#' @param seq (required): Biostring sequence, or a StringSet of multiple Biostring sequences for which the content is to be calculated.
#' @param bases (required): Character or vector of characters for bases that are to be used for calculating their content proportion in the sequence(s) provided. If a vector of multiple characters is provided, this function finds the proportion of the sequence that is any one of the bases (i.e. not necessarily in a contiguous manner).
#' @param baseOnly (default TRUE): boolean for if the user wants to calculate the content on using the base pair characters (ex. A, T, C, G for DNA; A, U, C, G for RNA) in the sequence, or from all characters in the sequence.
#' @param pct (default FALSE): boolean for if the user wishes to return the output as a percentage (TRUE) or proportion (FALSE, default).
#' @return Proportion or percentage for the base content of the provided base(s).
#' @export
#'
getBaseContent <- function(seq, bases, baseOnly=TRUE, pct=FALSE) {
	freq <- alphabetFrequency(seq, baseOnly=baseOnly)
	# single sequences will give frequencies as a one dimensional list. StringSet will be matrix
	if(is.matrix(freq)){
		# for StringSet input - multiple rows
		if(baseOnly) rowTotBases <- rowSums(freq[,1:4])
		else rowTotBases <- rowSums(freq)

		baseFreq <- rowSums(freq[,bases]) / rowTotBases

		if(is.null(names(seq))){
			message("The StringSet provided is unnamed. If you wish to directly compare specific sequences, it is recommended that the input is named for readability.")
		}
		names(baseFreq) <- names(seq)
	}
	else{
		# for single sequence input - single dimension list
		ntotBases <- ifelse(baseOnly, sum(freq[1:4]), sum(freq))
		baseFreq <- sum(freq[bases] / ntotBases)
	}
	# return
	ifelse(pct, return(baseFreq*100), return(baseFreq))
}

#' @name compareBaseContent
#' @title **Compare Base Content**
#' @description This function computes the content for each sequence in a genome for either a base or a vector of individual bases. This allows the user to compare between distinct sequences and identify which of the sequences provided has the highest or lowest base content for a given genome. For instance, a user may be interested in finding which chromosome in a genome has the highest CG content.
#' @param seqList (required): BSgenome object or a list of sequences, each for which the content is to be calculated. It is recommended that the elements of this data structure be named for readability of output. Otherwise, the output will provide the corresponding index occupied by the sequence in the original data structure that was provided alongside its base content calculation.
#' @param bases (required): Character or vector of characters for bases that are to be used for calculating their content proportion in the sequence provided. If a vector of characters is provided, this function finds the proportion of the sequence that is any one of the bases (i.e. necessarily in a contiguous manner).
#' @param baseOnly (default TRUE): boolean for if the user wants to calculate the content on using the base pair characters (ex. A, T, C, G for DNA) in the sequence, or from all characters in the sequence.
#' @param pct (default FALSE): boolean for if the user wishes to return the output as a percentage (TRUE) or proportion (FALSE, default).
#' @return Named vector of proportions or percentages for base content of each provided sequence.
#' @export
#'
compareBaseContent <- function(seqList, bases, baseOnly=TRUE, pct=FALSE) {
	if(is.null(names(seqList))){
		message("The sequences in the genome or list provided are unnamed. If you wish to directly compare sequences, it is recommended that the input is named for readability.")
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
	return(contentVector)
}
