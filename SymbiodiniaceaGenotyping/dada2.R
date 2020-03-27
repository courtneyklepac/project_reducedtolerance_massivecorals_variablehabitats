##Installing external dependencies
##R needs to be version 3.4, Bioconductor 3.5, and dada2 1.4
source("http://bioconductor.org/biocLite.R")
biocLite(suppressUpdates = FALSE)
biocLite("ShortRead", suppressUpdates = FALSE)
biocLite("dada2")
##Install using devtools. Doesn't work. 
# biocLite("devtools")
# library("devtools")
# devtools::install_github("benjjneb/dada2")

##Check package version
packageVersion("dada2")
#If Bioconductor is out of date
biocLite("BiocUpgrade") 
# biocLite()

##Load package
library("dada2")
##documentation
help(package="dada2")

##working with example from website
path <- "~/Documents/Dissertation/Data/dna/MiSeq/fastq2" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

##First, read in fastq file names and perform string manipulation to get lists of forward and reverse fastq files in matched order
# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(path, pattern="_R1_001.fastq_4cut.fastq"))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq_4cut.fastq"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
sample.names <- sapply(strsplit(fnRs, "_"), `[`, 1)
# Specify the full path to the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)
head(fnFs)
##With my data: the string manipulations may have to be modified, especially the extraction of sample names from the file names.

##Look at quality profiles of forward reads
plotQualityProfile(fnFs[1:2]) #looks like 240 for my data as well
plotQualityProfile(fnRs[1:2]) #200, drops off a lot and sooner than example dataset...
##Will want to trim where qualities crash, even though dada2 is robust to low quality sequences

##MUY IMPORTANTE!!!!! READS MUST OVERLAP AFTER TRIMMING IN ORDER TO MERGE LATER. truncLen() needs to be large enough to maintain overlap if primer sets are less-overlapping

##Make a filtered subdirectory. Define the filenames for the filtered
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))
dim(filtFs)
#Filter reads based on filtering parameters: dada2 requires no Ns (maxN=0), truncQ=2, rm.phix=T, maxEE=2 (sets max # of expected errors allowed in read), truncLen is based on quality profile assesment
#ITS amplicon strategies, it is undesirable to truncate reads to a fixed length due to the large amount of length variation at that locus. That is OK, just leave out truncLen. 
#output directory called /Users/corky/Downloads/MiSeq_SOP/filtered
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,200),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)
out <- filterAndTrim(fnFs, filtFs, truncLen=240,
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)
head(out) #reverse didn't work well, fewer # of reads and that's why it didn't work

##DADA2 algorithm depends on a parametric error model (err) and every amplicon dataset has a different set of error rates. The learnErrors method learns the error model from the data, by alternating estimation of the error rates and inference of sample composition until they converge
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
#visualize estimated error rates
plotErrors(errF, nominalQ=TRUE)

##Dereplicate filtered fastqs. Combines all identical reads into unique sequences w/ corresponding abundances (clustering and read counts in one), with a summary of quality info associated with each unique seq.
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

##Infer sequence variants in each sample from dereplicated data.
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
#Inspect the dada-class object returned - real variant # from unique sequences in each sample. [1] is for first sample.
dadaFs[[1]]
dadaRs[[1]]

##Merge paired reads to reduce spurious sequence variants. Depends on F and R reads being in matching order at time of dereplication.
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
#Now have a data.frame for each sample with the merged $sequence, its $abundance, and the indices of merged $forward and $reverse denoised seq. Paired reads that didn't exactly overlap were removed.

##Construct sequence table (which is a matrix of rows=samples and columns=seq variants), higher-res version of OTU table
seqtab <- makeSequenceTable(dadaFs,derepFs)
write.csv(seqtab,"~/Documents/Dissertation/Data/dna/MiSeq/fastqs/filtered/seqtable-dadaF.csv",row.names=T)
## The sequences being tabled vary in length.
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

##Remove chimeras which still remain after dada (removes sub and indel errors)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
#calculate fraction of chimeras, here variants account for only 4% of total seq reads
sum(seqtab.nochim)/sum(seqtab)

##Final check of number of reads that made it through each step in pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab), rowSums(seqtab.nochim))
dim(track)
colnames(track) <- c("input", "filtered", "denoised", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track,"~/Documents/Dissertation/Data/dna/MiSeq/2017_AS-ITS-DADAfsumm.csv",row.names=T)


#If majority of reads failed to merge, revisit truncLen parameter and make sure the truncated reads span amplicon
#If reads failed chimera check, revisit removal of primers

##Assign taxonomy- takes set of sequences and training set of taxon classified seqs, and outputs with at least minBoot bootstrap confidence.
taxa <- assignTaxonomy(seqtab.nochim, "~/Downloads/MiSeq_SOP/rdp_train_set_14.fa.gz", multithread=TRUE)
unname(head(taxa))

##Evaluate accuracy to compare the seq variants inferred by dada2 to the expected composition of the community. Needed?.....
unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop SVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
#DADA2 inferred 20 sample sequences present in the Mock community.
mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
## Of those, 20 were exact matches to the expected reference sequences.

###############END OF DADA2#####################
########PHYLOGENETIC ANALYSIS W PHYLOSEQ########
##dada2 pipeline produces seq table and taxonomy table for further analysis in phyloseq
#Import into phyloseq
install.packages("phyloseq", dependencies=T)
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
# Make a data.frame holding the sample data
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out

# Construct phyloseq object (straightforward from dada2 outputs)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
ps
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 229 taxa and 19 samples ]
## sample_data() Sample Data:       [ 19 samples by 4 sample variables ]
## tax_table()   Taxonomy Table:    [ 229 taxa by 6 taxonomic ranks ]

