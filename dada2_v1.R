## DADA2 for ITS###
## Input data already trimmed with cutadapt##

### https://benjjneb.github.io/dada2/bigdata.html
### https://github.com/benjjneb/ITS-Workflow/blob/master/ITS_workflow.md

library(dada2); packageVersion("dada2")
# File parsing
pathF <- "/Users/kc178/Documents/Florida_projects/florencia/asilomar_redo/ITS/R1/cut" # CHANGE ME to the directory containing your demultiplexed forward-read fastqs
pathR <- "/Users/kc178/Documents/Florida_projects/florencia/asilomar_redo/ITS/R2/cut" # CHANGE ME ...
filtpathF <- file.path(pathF, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- file.path(pathR, "filtered") # ...
fastqFs <- sort(list.files(pathF, pattern="fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="fastq.gz"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
# Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS
out<-filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
                   rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
                   maxEE=c(2,6), truncQ=11, maxN=0, rm.phix=TRUE, minLen = 50,
                   compress=TRUE, verbose=TRUE, multithread=TRUE)



filtpathF <- "/Users/kc178/Documents/Florida_projects/florencia/asilomar_redo/ITS/R1/cut/filtered" # CHANGE ME to the directory containing your filtered forward fastqs
filtpathR <- "/Users/kc178/Documents/Florida_projects/florencia/asilomar_redo/ITS/R2/cut/filtered" # CHANGE ME ...
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}

rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "/Users/kc178/Documents/Florida_projects/florencia/asilomar_redo/test/seqtab.rds")

head(merger[1])
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))



track <- cbind(out, sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "merged", "nonchim")
rownames(track) <- sample.names
track
write.table(track,"/Users/kc178/Documents/Florida_projects/florencia/asilomar_redo/test/track.txt",quote=F, sep="\t")


tax <- assignTaxonomy(seqtab, "/Users/kc178/Documents/Florida_projects/florencia/asilomar_redo/ITS/sh_general_release_dynamic_01.12.2017.fasta", multithread=TRUE)