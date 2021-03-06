% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/baseContent.R
\name{getBaseContent}
\alias{getBaseContent}
\title{\strong{Get Base Content}}
\usage{
getBaseContent(seq, bases, baseOnly = TRUE, pct = FALSE)
}
\arguments{
\item{seq}{(required): Biostring sequence, or a StringSet of multiple Biostring sequences for which the content is to be calculated.}

\item{bases}{(required): Character or vector of characters for bases that are to be used for calculating their content proportion in the sequence(s) provided. If a vector of multiple characters is provided, this function finds the proportion of the sequence that is any one of the bases (i.e. not necessarily in a contiguous manner).}

\item{baseOnly}{(default TRUE): boolean for if the user wants to calculate the content on using the base pair characters (ex. A, T, C, G for DNA; A, U, C, G for RNA) in the sequence, or from all characters in the sequence.}

\item{pct}{(default FALSE): boolean for if the user wishes to return the output as a percentage (TRUE) or proportion (FALSE, default).}
}
\value{
Proportion or percentage for the base content of the provided base(s).
}
\description{
This function computes the content in a sequence or a StringSet of multiple sequences for either a base or a vector of individual bases.
}
