#Load libraries

library(dada2)

#-------------------------------------------------------
# Before running
#-------------------------------------------------------
# 1) Demultiplex your samples if necessary

# 2) Save demultiplexed files in folder named "reads" and set working directory to parent directory:

setwd("~/demultiplex_reads")

#-------------------------------------------------------
# Setup
#-------------------------------------------------------

#Paths to demultiplexed reads
reads<-"reads"

# Get the filenames with relative path
# sort to ensure same order of fwd/rev reads
fnFs <- sort(list.files(reads, pattern="-R1.fastq.gz", full.names=TRUE))
fnRs <- sort(list.files(reads, pattern="-R2.fastq.gz", full.names=TRUE))
# Get sample names only (remove path, and everything after the first "-")
# Assuming filenames have format: SAMPLENAME-XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "-R"), `[`, 1)

#list the files (not required)
#list.files(taxpath)

# Get the filenames with relative path
# sort to ensure same order of fwd/rev reads
fnFs <- sort(list.files(reads, pattern="-R1.fastq", full.names=TRUE))
fnRs <- sort(list.files(reads, pattern="-R2.fastq", full.names=TRUE))
# Get sample names only (remove path, and everything after the first "-")
# Assuming filenames have format: SAMPLENAME-XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "-R"), `[`, 1)

#check for duplicated sample names before you move on
message ("###### checking for duplicated sample names")
any(duplicated(sample.names))
# STOP if this is TRUE and check: sample.names[duplicated(sample.names)]

#-------------------------------------------------------
# Check read quality
#-------------------------------------------------------

#This will pick a random subset of 4 samples to look at read quality
#note: you can also plot an aggregate of all fastqs instead using aggregate=TRUE
ids<-round(runif(12,1,length(sample.names)))
pdf("qualprofiles.pdf")
plotQualityProfile(fnFs[ids])
plotQualityProfile(fnRs[ids])
dev.off()

# The distribution of quality scores at each position is shown as a grey-scale heat map,
# with dark colors corresponding to higher frequency.
# The plotted lines show positional summary statistics:
# green is the mean,
# orange is the median,
# and the dashed orange lines are the 25th and 75th quantiles.

#-------------------------------------------------------
# Filter reads based on QC
#-------------------------------------------------------
message ("###### Filtering reads based on QC")

# Make filenames for the filtered fastq files
filtFs <- paste0(reads, "/", sample.names, "-F-filt.fastq.gz")
filtRs <- paste0(reads, "/", sample.names, "-R-filt.fastq.gz")

# This sequencing run (honey bee gut microbiota) had ************12 mer barcodes************

out<-filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                   truncLen=c(295,295),		#optimized based on QC readout                           #225/180 for att_no1 ...220/180 for att_no2
                   truncQ=2, #enter quality score cut off (trunicates after </=2 QC score)
                   trimLeft=c(29,33), #remove 10 bases off the 5' end of R2 sequences only
                   maxN=0,
                   minLen=c(265,261),
                   maxEE=c(2,2),
                   rm.phix=TRUE,
                   compress=TRUE, verbose=TRUE, multithread=TRUE)

write.table(out, file="DADA2_after_filter.txt", sep="\t", col.names=NA, quote=F)

message ("###### Learning error rates - SLOW !!")

# *************** On Mac set multithread=TRUE / On Windows set multithread=FALSE**********
errF <- learnErrors(filtFs, multithread=TRUE, randomize=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE, randomize=TRUE)
#	randomize=TRUE #don't pick the first 1mil for the model, pick a random set

#Plot the error rates
# Do not proceed without a good fit
pdf("err.pdf")
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()

# Dereplication combines all identical sequencing reads into into "unique sequences" with a corresponding "abundance": the number of reads with that unique sequence
# Dereplication substantially reduces computation time by eliminating redundant comparisons.

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#message ("###### Saving your R session...")

#-------------------------------------------------------
# Sample inference, merge paired reads, remove chimeras
#-------------------------------------------------------
message ("###### Inferring the sequence variants in each sample - SLOW!!")

dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool=FALSE, verbose=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool=FALSE, verbose=TRUE)

# overlap the ends of the forward and reverse reads
message ("###### merging the Fwd and Rev reads")
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# make the sequence table, samples are by rows
seqtab <- makeSequenceTable(mergers)

message ("###### summarize the output by sequence length")
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", verbose=TRUE, multithread=TRUE, minFoldParentOverAbundance=2)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

write.table(seqtab.nochim, file="temp_dada2_nochim.txt", sep="\t", col.names=NA, quote=F)
write.table((t(seqtab.nochim)), file="temp_dada2_nochim_transpose.txt", sep="\t", col.names=NA, quote=F)

#-------------------------------------------------------
# Sanity check
#-------------------------------------------------------
message ("###### sanity check - how many reads made it: readsout.txt")
# Check how many reads made it through the pipeline
# This is good to report in your methods/results
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
write.table(track, file="DADA2_readsout.txt", sep="\t", col.names=NA, quote=F)

#-------------------------------------------------------
# Read files back in if starting in a new session
#-------------------------------------------------------

seqtab.nochim.t <-read.table("temp_dada2_nochim_transpose.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
seqtab.nochim <- t(seqtab.nochim.t)

library(dada2)
library(DECIPHER)
library(data.table)

#------------------------------------------------------------------
# Write ASV count table, sequence table, and sequence fasta files
#------------------------------------------------------------------
#:::::::::::::::::::::
# ASV sequence fasta
#:::::::::::::::::::::
#Pull out seqs
dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
#Write FASTA
names(dna) <- paste("SV", seq(from = 0, to = length(dna)-1), sep="_")
writeXStringSet(dna, filepath="ASVseqs1111.fasta", append=FALSE, compress=FALSE, compression_level=NA, format="fasta")

#::::::::::::::::::::::::::::::::::
# ASV count table + sequence table
#::::::::::::::::::::::::::::::::::

seqtab.nochim.t.x <- seqtab.nochim.t
sv.seqs<-rownames(seqtab.nochim.t.x)
sv.num<-paste("SV", seq(from = 0, to = nrow(seqtab.nochim.t.x)-1), sep="_")
rownames(seqtab.nochim.t.x)<-sv.num
seqtab.nochim.t.x$tax.vector <- sv.num
write.table(seqtab.nochim.t.x, file="ASV_counts_table.txt", sep="\t", col.names=NA, quote=F) #Kingdom to Species
write.table(sv.seqs, file="ASV_seqs_table.txt", sep="\t", row.names=sv.num, col.names=F,  quote=F) #Kingdom to Species
