
library(BSgenome)
available.genomes()
library(BSgenome.Hsapiens.UCSC.hg19)

### get sequence for chromosome 1
Seq=Hsapiens[["chr1"]]
Seq2=Hsapiens[["chr2"]]

bases=alphabetFrequency(Seq,baseOnly=TRUE)
bases[1:4]
ntotBases=sum(bases[1:4])
baseFreq=bases[1:4]/ntotBases
GCcontent=baseFreq["C"]+baseFreq["G"] # proportion of C or G
# double check
getBaseContent(Seq, c("C", "G"), baseOnly=TRUE, pct=TRUE)
getBaseContent(Seq2, c("C", "G"), baseOnly=TRUE, pct=FALSE)

# unnamed list
seqs <- list(Seq, Seq2)
# named list
seqs_named <- list(s1 = Seq, s2 = Seq2)

compareBaseContent(seqs, c("C", "G"), baseOnly=TRUE, pct=TRUE, desc=TRUE)
compareBaseContent(seqs_named, c("C", "G"), baseOnly=TRUE, pct=TRUE, desc=FALSE)
compareBaseContent(Hsapiens, c("C", "G"), baseOnly=TRUE, pct=TRUE, desc=TRUE) # takes a little longer but works

getPatternContent(Seq, "CG", baseOnly=TRUE, ratio=TRUE)
getPatternContent(Seq2, "TG", baseOnly=TRUE, ratio=FALSE)
