# Introduction

The goal of the `sequenceComp` R package is to facilitate the characterization of DNA or RNA sequences and comparing their contents. This package works in concert with objects and functions inherited from the `Biostrings` package. The functions implemented by `sequenceComp` expedite and streamline the exploration of DNA or RNA sequences by characterizing their base pair compositions in a flexible, user-specified manner. This package also allows the user to easily compare the compositions of sets of sequences or ranges of a specific length within a single sequence.

# Installation

After loading the `sequenceComp.Rproj` The R package can be loaded into the user's local memory using the `devtools::load_all()` function. 

If the user wishes to continue using `sequenceComp` in future R sessions without having to reload the package into memory, they can install the package using `devtools::install()`. Then, the package will be added to the user's R system library from which the package can be loaded using `library("sequenceComp")`.

# List of functions

After installation, the user can call `?nameOfFunction` within R for more detailed documentation.

### getBaseContent

**Usage:** `getBaseContent(seq, bases, baseOnly = TRUE, pct = FALSE)`

**Description:** This function computes the content in a sequence or a StringSet of multiple sequences for either a base or a vector of individual bases.

### compareBaseContent

**Usage:** `compareBaseContent(seqList, bases, baseOnly = TRUE, pct = FALSE)
`

**Description:** This function computes the content for each sequence in a genome for either a base or a vector of individual bases. This allows the user to compare between distinct sequences and identify which of the sequences provided has the highest or lowest base content for a given genome. For instance, a user may be interested in finding which chromosome in a genome has the highest CG content.

### getPatternContent

**Usage:** `getPatternContent(seq, pattern, baseOnly = TRUE)
`

**Description:** This function computes the ratio of observed to expected content for a specific contiguous pattern of bases in a sequence or a StringSet of multiple sequences.

### comparePatternContent

**Usage:** `comparePatternContent(seqList, pattern, baseOnly = TRUE, log = FALSE)`

**Description:** This function omputes the ratio of observed to expected occurrences for a contiguous pattern of bases for each sequence in a genome or a list of sequences. These ratios are then outputted in either ascending or descending order of content.

### windowBaseContent

**Usage:** `windowBaseContent(seq, windowSize, bases, zero.rm = FALSE, pct = FALSE)`

**Description:** Given an individual sequence, this function splits the provided sequence into windows of a user-specified size and then computes the content in each window for either a base or a vector of individual bases.

### windowPatternContent

**Usage:** `windowPatternContent(seq, windowSize, pattern, na.rm = FALSE, log = FALSE)`

**Description:** Given an individual sequence, this function splits the provided sequence into windows of a user-specified size and then computes the content in each window for a contiguous pattern of bases.
