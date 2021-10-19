library(sequenceComp)
library(BSgenome)
available.genomes()
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Btaurus.UCSC.bosTau3)

### get sequence for chromosome 1
Seq=Hsapiens[["chr1"]]
Seq2=Hsapiens[["chr2"]]

bases=alphabetFrequency(Seq,baseOnly=TRUE)
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
sort(Hsapiens.GCcontent, decreasing=TRUE) # find highest to lowest GC content in the Hsapiens genome
hist(Hsapiens.GCcontent, main="Histogram of GC content for each chromosome in the Hsapiens genome") # visualize for each chromosome

###########################################################################
############### PATTERN CONTENT
## look at CG dinucleotide content
cg=matchPattern("CG", Seq)
ncg=length(cg) # 2284470 observed CG's
## compute the observed to expected ratio
baseFreq["C"]*baseFreq["G"]*ntotBases # expect 9814098 CG's
ncg/(baseFreq["C"]*baseFreq["G"]*ntotBases) ## this shows CG rarely stay together.

## compare to the observed to expected ratio of TG
tg=matchPattern("TG", Seq)
ntg=length(tg) # 16410688 observed TG's
baseFreq["T"]*baseFreq["G"]*ntotBases # expect 13705214 TG's
ntg/(baseFreq["T"]*baseFreq["G"]*ntotBases) ## this shows TG presented more than expected.

getPatternContent(Seq, "CG", baseOnly=TRUE)
getPatternContent(Seq2, "CG", baseOnly=TRUE)
getPatternContent(seqs, "CG", baseOnly=TRUE)

getPatternContent(Seq, "TG", baseOnly=TRUE)
getPatternContent(Seq2, "TG", baseOnly=TRUE)
getPatternContent(seqs_named, "TG", baseOnly=TRUE)

comparePatternContent(seqs, "CG", baseOnly=TRUE)
comparePatternContent(seqs_named, "CG", baseOnly=TRUE)

Hsapiens.obsExpCG <- comparePatternContent(Hsapiens, "CG", baseOnly=TRUE, log=TRUE)
sort(Hsapiens.obsExpCG, decreasing=TRUE) # sort chromosomes by highest CG dinucleotide content
Hsapiens.obsExpCG %>% hist(main="Histogram of CG dinucleotide obs/expected ratios for each chromosome in the Hsapiens genome")

###########################################################################
############### WINDOW PATTERN CONTENT
ss=seq(1, length(Seq), by=1000)
ss=ss[-length(ss)] ## remove the last one
Seq.set=DNAStringSet(Seq, start=ss, end=ss+999)
ff=alphabetFrequency(Seq.set, baseOnly=TRUE) # frequency of A, C, G, T in each 1000 bp window
pCG=(ff[,"C"]+ff[,"G"])/rowSums(ff) # proportion of C or G in each 1000 bp window
hist(pCG[pCG>0],100)
density(pCG[pCG>0])

pCG1 <- windowBaseContent(Seq, 1000, c("C", "G"), zero.rm=TRUE, pct=FALSE)
pCG1 %>% hist(breaks=100, main="Histogram of GC base content in 1000 bp windows")
density(pCG1) %>% plot()

obsExpCG <- windowPatternContent(Seq, 1000, "CG", na.rm=TRUE)
obsExpCG %>% hist(breaks=100, main="Histogram of CG dinucleotide obs/expected ratios in 1000 bp windows")
density(obsExpCG) %>% plot()

log.obsExpCG <- windowPatternContent(Seq, 1000, "CG", na.rm=TRUE, log=TRUE)
log.obsExpCG %>% hist(breaks=100, main="Histogram of CG dinucleotide logged obs/expected ratios in 1000 bp windows")
density(log.obsExpCG) %>% plot()
