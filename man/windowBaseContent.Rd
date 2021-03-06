% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/windowContent.R
\name{windowBaseContent}
\alias{windowBaseContent}
\title{\strong{Window Base Content}}
\usage{
windowBaseContent(seq, windowSize, bases, zero.rm = FALSE, pct = FALSE)
}
\arguments{
\item{seq}{(required): Biostring sequence that is to be split into windows, in which the base content is to be calculated for each window.}

\item{windowSize}{(required): The size of each window, in base pairs.}

\item{bases}{(required): Character or vector of characters for bases that are to be used for calculating their content proportion in each created window from the sequence. If a vector of multiple characters is provided, this function finds the proportion of the window that is any one of the bases (i.e. not necessarily in a contiguous manner).}

\item{zero.rm}{(default FALSE): boolean for if the user wishes to remove all proportions or percentages in the returned vector that are equal to zero.}

\item{pct}{(default FALSE): boolean for if the user wishes to return the output as a percentage (TRUE) or proportion (FALSE, default).}
}
\value{
A vector of proportions or percentages for base content in each window of the user-specified size for a sequence.
}
\description{
Given an individual sequence, this function splits the provided sequence into windows of a user-specified size and then computes the content in each window for either a base or a vector of individual bases.
}
