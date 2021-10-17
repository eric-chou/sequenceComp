library(BSgenome)
library(GenomicFeatures)
library(dplyr)

#' @name getPatternContent
#' @title **Get Pattern Content**
#' @description This function computes the ratio of observed to expected content for a specific contiguous pattern of bases in a sequence or a StringSet of multiple sequences.
#' @param seq (required): Biostring sequence, or a StringSet of multiple Biostring sequences for which the content is to be calculated.
#' @param pattern (required): String for a contiguous pattern of bases that are to be used for calculating their observed and expected occurrences in the sequence(s) provided. 
#' @param baseOnly (default TRUE): boolean for if the user wants to calculate the content on using the base pair characters (ex. A, T, C, G for DNA; A, U, C, G for RNA) in the sequence, or from all characters in the sequence.
#' @return Unlike getBaseContent, this provides not proportions or percentages for each sequence, but either a ratio of observed to expected counts.
getPatternContent <- function(seq, pattern, baseOnly=TRUE){
	inputClass <- class(seq)[1]
	# convert into StringSet in order to normalize vectorized functionality
	if(inputClass == "DNAString"){
		seq <- DNAStringSet(seq)
	}
	else if(inputClass == "RNAString"){
		seq <- RNAStringSet(seq)
	}
	else if(inputClass %in% c("DNAStringSet", "RNAStringSet")){
		if(is.null(names(seq))){
			message("The StringSet provided is unnamed. If you wish to directly compare specific sequences, it is recommended that the input is named for readability.")
		}
	}
	else{
		stop("Sequence inputted is of class: ", inputClass, "\nPlease provide a DNAString, RNAString, DNAStringSet, or RNAStringSet for the input sequence(s).")
	}
	# observed count 
	n <- vcountPattern(pattern, seq)
	freq <- alphabetFrequency(seq, baseOnly=baseOnly)
	if(nrow(freq) > 1){
		# base R rowSums only works with matrices of more than 1 row
		if(baseOnly) rowTotBases <- rowSums(freq[,1:4])
		else rowTotBases <- rowSums(freq)
	}
	else{
		if(baseOnly) rowTotBases <- sum(freq[1:4])
		else rowTotBases <- sum(freq)
	}
	freq <- freq / rowTotBases

	# calculate expected - multiply frequencies for each individual base in the given pattern string, then multiply by total bases
	basesVec <- strsplit(pattern, "")[[1]]
	if(nrow(freq) > 1) expected <- apply(freq[,basesVec], 1, prod) * rowTotBases # multiplication across rows for a multidimensional matrix requires apply function instead of simple prod()
	else expected <- prod(freq[,basesVec]) * rowTotBases

	out <- n / expected
	names(out) <- names(seq)

	return(out)
}

#' @name genomePatternContent
#' @title **Genome Pattern Content**
#' @description This function omputes the ratio of observed to expected occurrences for a contiguous pattern of bases for each sequence in a genome or a list of sequences. These ratios are then outputted in either ascending or descending order of content.
#' @param seqList (required): BSgenome or a list of sequences, each for which the content is to be calculated. It is recommended that the elements of this data structure be named for readability of output. Otherwise, the output will provide the corresponding index occupied by the sequence in the original datas tructure that was provided alongside its calculated ratio.
#' @param pattern (required): String for a contiguous pattern of bases that are to be used for calculating their observed to expected ratios in each sequence provided.
#' @param baseOnly (default TRUE): boolean for if the user wants to calculate the content on using the base pair characters (ex. A, T, C, G for DNA) in the sequence, or from all characters in the sequence.
#' @return Named vector of observed to expected ratios of pattern occurrences for each provided sequence.
comparePatternContent <- function(seqList, pattern, baseOnly=TRUE){
	if(is.null(names(seqList))){
		message("The sequences in the genome or list provided are unnamed. If you wish to directly compare sequences, it is recommended that the input is named for readability.")
	}
	ratioVector <- c()
	idxVector <- c()
	for(i in 1:length(seqList)){
		ratio <- getPatternContent(seqList[[i]], pattern=pattern, baseOnly=baseOnly)
		ratioVector <- c(ratioVector, ratio)
		idx <- ifelse(is.null(names(seqList)), i, names(seqList)[i])
		idxVector <- c(idxVector, idx)
	}

	names(ratioVector) <- idxVector
	return(ratioVector)
}
