
library(BSgenome)
available.genomes()
library(BSgenome.Hsapiens.UCSC.hg19)

### get sequence for chromosome 1
Seq=Hsapiens[["chr1"]]
Seq2=Hsapiens[["chr2"]]

bases=alphabetFrequency(Seq2,baseOnly=TRUE)
bases[1:4]
ntotBases=sum(bases[1:4])
baseFreq=bases[1:4]/ntotBases
GCcontent=baseFreq["C"]+baseFreq["G"] # proportion of C or G
GCcontent
# double check
getBaseContent(Seq, c("C", "G"), baseOnly=TRUE, pct=TRUE)
getBaseContent(Seq2, c("C", "G"), baseOnly=TRUE, pct=FALSE)

# unnamed list
seqs <- list(Seq, Seq2) %>% DNAStringSet()
# named list
seqs_named <- list(s1 = Seq, s2 = Seq2) %>% DNAStringSet()

getBaseContent(seqs, c("C", "G"), baseOnly=TRUE, pct=TRUE)
getBaseContent(seqs_named, c("C", "G"), baseOnly=FALSE, pct=FALSE)

compareBaseContent(list(Seq, Seq2), c("C", "G"), baseOnly=TRUE, pct=TRUE)
Hsapiens.GCcontent <- compareBaseContent(Hsapiens, c("C", "G"), baseOnly=TRUE, pct=FALSE) # takes a bit of time but works
hist(Hsapiens.GCcontent) # visualize for each chromosome
sort(Hsapiens.GCcontent, decreasing=TRUE) # sort chromosomes by highest GC content

###########################################################################
############### PATTERN CONTENT
getPatternContent(Seq, "CG", baseOnly=TRUE)
getPatternContent(Seq2, "CG", baseOnly=TRUE)
getPatternContent(seqs, "CG", baseOnly=TRUE)

getPatternContent(Seq, "TG", baseOnly=FALSE)
getPatternContent(Seq2, "TG", baseOnly=FALSE)
getPatternContent(seqs_named, "TG", baseOnly=FALSE)

comparePatternContent(seqs, "CG", baseOnly=TRUE)
comparePatternContent(seqs_named, "CG", baseOnly=TRUE)
comparePatternContent(Hsapiens, "CG", baseOnly=TRUE)
