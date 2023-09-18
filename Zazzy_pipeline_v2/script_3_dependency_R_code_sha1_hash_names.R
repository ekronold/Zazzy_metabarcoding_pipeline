#### newDADA2 script ####################################
#### originally by LUIS MORGADO 11 DEC 2019  #########################
#### Modified for SWARM pipeline by EIVIND KVERME RONOLD 20 FEB 2023 #########
####
#### for 250 bp FASTERIS
#### General PATH
library(tidyverse)
library(ggplot2)
library(dada2)
library(ShortRead)
library(Biostrings)
library(openxlsx)
library(digest)
#setDadaOpt(BAND_SIZE=16) #set to 16 (default) for 18S, 32 for ITS
path <- getwd()  
list.files(path)
########################## SENSUS - SS analyses #########
##########################
path.SS <- file.path(path, "DADA2_SS")
list.files(path.SS)
fnFs <- sort(list.files(path.SS, pattern = "_R1.fastq", full.names = TRUE))    
fnRs <- sort(list.files(path.SS, pattern = "_R2.fastq", full.names = TRUE))	
#### set the samples names #######
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names

##### look at the primers in the reads ###########

FWD <-"ACTGATCGATCGACGA"    #INSERT_FWD_PRIMER
REV <- "ATGCTAGCATGCATCGATG"    #INSERT_REV_PRIMER

allOrients <- function(primer) { # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

### check for primers in sequences ####
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)))
  return(sum(nhits > 0))
}

### in the demultiplexed files SS ####
Primers.SS<-rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]), 
                  FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[1]]), 
                  REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[1]]), 
                  REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]]))
Primers.SS
tbname_Primers.SS <- file.path(path,"Primers.in.SS.table")
write.table(Primers.SS,tbname_Primers.SS, sep = "\t", col.names = NA,row.names = T, quote = F)

plotQualityProfile(fnFs[1:4])
ggsave("FWD.raw_qualityProfile.SS.pdf", width = 18, height = 14, units = c("in"))
plotQualityProfile(fnRs[1:4])
ggsave("REV.raw_qualityProfile.SS.pdf", width = 18, height = 14, units = c("in"))

### filter the SS reads

filtFs <- file.path(path.SS, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path.SS, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out.SS <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                        maxN=0, maxEE=c(1,1), truncQ=2, minLen = 50, 
                        rm.phix=TRUE, compress=TRUE, multithread=T)

##### the analyses #######
errF <- learnErrors(filtFs, verbose=TRUE, multithread=T)
errR <- learnErrors(filtRs, verbose=TRUE, multithread=T)

plotErrors(errF, nominalQ=TRUE)
ggsave("error_model_FWD.SS.pdf")
plotErrors(errR, nominalQ=TRUE)
ggsave("error_model_REV.SS.pdf")

dadaFs <- dada(filtFs, err=errF, pool = TRUE, multithread=T)
dadaRs <- dada(filtRs, err=errR, pool = TRUE, multithread=T)
names(dadaFs) <- sample.names
names(dadaRs) <- sample.names
mergers.SS <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, trimOverhang=T, verbose=TRUE)
seqtab.SS <- makeSequenceTable(mergers.SS)
seqtab.nochim.SS <- removeBimeraDenovo(seqtab.SS, method="consensus",multithread=TRUE, verbose=TRUE)
#####
1-sum(seqtab.nochim.SS)/sum(seqtab.SS)
stSS <- file.path(path,"seqtab_SS")
stnsSS <- file.path(path,"seqtab.nochim_SS")
saveRDS(seqtab.SS,stSS)
saveRDS(seqtab.nochim.SS,stnsSS)

getN <- function(x) sum(getUniques(x))
track.SS <- cbind(out.SS, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers.SS, getN), rowSums(seqtab.nochim.SS))
colnames(track.SS) <- c("input.SS", "filtered.SS", "denoisedF.SS", "denoisedR.SS", "merged.SS", "nonchim.SS")
rownames(track.SS) <- sample.names
track.SS

########################## Anti-SENSUS - AS analyses #########
##########################
path.AS <- file.path(path, "DADA2_AS")
list.files(path.AS)
fnFs <- sort(list.files(path.AS, pattern = "_R2.fastq", full.names = TRUE)) ### R1 - should be the FORWARD reads - this pattern might change depending on the files outputted from the demultiplexing step
fnRs <- sort(list.files(path.AS, pattern = "_R1.fastq", full.names = TRUE))	### R2 - should be the REVERSE reads - this pattern might change depending on the files outputted from the demultiplexing step
#### set the samples names #######
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names

##### look at the primers in the reads ###########
allOrients <- function(primer) { # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

### check for primers in sequences ####
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)))
  return(sum(nhits > 0))
}

### in the demultiplexed files SS ####
Primers.AS<-rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]), 
                  FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[1]]), 
                  REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[1]]), 
                  REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]]))
Primers.AS
tbname_Primers.AS <- file.path(path,"Primers.in.AS.table")
write.table(Primers.AS,tbname_Primers.AS, sep = "\t", col.names = NA,row.names = T, quote = F)

filtFs <- file.path(path.AS, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path.AS, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out.AS <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                        maxN=0, maxEE=c(1,1), truncQ=2, minLen = 50, 
                        rm.phix=TRUE, compress=TRUE, multithread=T)

##### the analyses #######
errF <- learnErrors(filtFs, verbose=TRUE, multithread=T)
errR <- learnErrors(filtRs, verbose=TRUE, multithread=T)
plotErrors(errF, nominalQ=TRUE)
ggsave("error_model_FWD.AS.pdf")
plotErrors(errR, nominalQ=TRUE)
ggsave("error_model_REV.AS.pdf")
dadaFs <- dada(filtFs, err=errF, pool = "pseudo", multithread=T)
dadaRs <- dada(filtRs, err=errR, pool = "pseudo", multithread=T)
names(dadaFs) <- sample.names
names(dadaRs) <- sample.names
mergers.AS <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, trimOverhang=T, verbose=TRUE)
seqtab.AS <- makeSequenceTable(mergers.AS)
seqtab.nochim.AS <- removeBimeraDenovo(seqtab.AS, method="consensus",multithread=TRUE, verbose=TRUE)

stAS <- file.path(path,"seqtab_AS")
stnsAS <- file.path(path,"seqtab.nochim_AS")
saveRDS(seqtab.AS,stAS)
saveRDS(seqtab.nochim.AS,stnsAS)

getN <- function(x) sum(getUniques(x))
track.AS <- cbind(out.AS, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers.AS, getN), rowSums(seqtab.nochim.AS))
colnames(track.AS) <- c("input.AS", "filtered.AS", "denoisedF.AS", "denoisedR.AS", "merged.AS", "nonchim.AS")
rownames(track.AS) <- sample.names

track.AS
track.SS
write.xlsx(track.AS, "track_reads_AS.xlsx")
write.xlsx(track.SS, "track_reads_SS.xlsx")


#Define a function for combining two or more tables, collapsing samples with similar names:  
sumSequenceTables <- function(table1, table2, ..., orderBy = "abundance") {
  # Combine passed tables into a list
  tables <- list(table1, table2)
  tables <- c(tables, list(...))
  # Validate tables
  if(!(all(sapply(tables, dada2:::is.sequence.table)))) {
    stop("At least two valid sequence tables, and no invalid objects, are expected.")
  }
  sample.names <- rownames(tables[[1]])
  for(i in seq(2, length(tables))) {
    sample.names <- c(sample.names, rownames(tables[[i]]))
  }
  seqs <- unique(c(sapply(tables, colnames), recursive=TRUE))
  sams <- unique(sample.names)
  # Make merged table
  rval <- matrix(0L, nrow=length(sams), ncol=length(seqs))
  rownames(rval) <- sams
  colnames(rval) <- seqs
  for(tab in tables) {
    rval[rownames(tab), colnames(tab)] <- rval[rownames(tab), colnames(tab)] + tab
  }
  # Order columns
  if(!is.null(orderBy)) {
    if(orderBy == "abundance") {
      rval <- rval[,order(colSums(rval), decreasing=TRUE),drop=FALSE]
    } else if(orderBy == "nsamples") {
      rval <- rval[,order(colSums(rval>0), decreasing=TRUE),drop=FALSE]
    }
  }
  rval
}

stAS <- file.path(path,"seqtab_AS")
stnsAS <- file.path(path,"seqtab.nochim_AS")
stSS <- file.path(path,"seqtab_SS")
stnsSS <- file.path(path,"seqtab.nochim_SS")
seqtab.nochim.AS <- readRDS(stnsAS)
seqtab.nochim.SS <- readRDS(stnsSS)
seqtab_AS <- readRDS(stAS)
seqtab_SS <- readRDS(stSS)
Both_sumtable_with_chim <- sumSequenceTables(seqtab_SS,seqtab_AS)
Both_nochim_sumtable <- sumSequenceTables(seqtab.nochim.SS,seqtab.nochim.AS)
stBoth <- file.path(path,"seqtab.with_Chim_Both")
stnsBoth <- file.path(path,"seqtab.nochim_Both")
saveRDS(Both_sumtable_with_chim,stBoth)
saveRDS(Both_nochim_sumtable,stnsBoth)

#### get and save read length
reads.per.seqlen.no.chim <- tapply(colSums(Both_nochim_sumtable), factor(nchar(getSequences(Both_nochim_sumtable))), sum)
df.without.chim <- data.frame(length=as.numeric(names(reads.per.seqlen.no.chim)), count=reads.per.seqlen.no.chim)
write.xlsx(df.without.chim, "No_chimeras_seq_length_table.xlsx")
df.without.chim$count<-as.numeric(df.without.chim$count)
ggplot(data=df.without.chim, aes(x=length, y=count)) + geom_col()
ggsave("ASV_no_chime_length.pdf", width = 18, height = 14, units = c("in"))

#### get all reads counted
track.SS<-as.data.frame(track.SS) %>% 
  rownames_to_column("sample.names")
track.SS[is.na(track.SS)]<-0
track.AS<-as.data.frame(track.AS) %>% 
  rownames_to_column("sample.names")
track.AS[is.na(track.AS)]<-0

track_merged_data<-cbind(rowSums(Both_sumtable_with_chim), rowSums(sign(Both_sumtable_with_chim)), rowSums(Both_nochim_sumtable), rowSums(sign(Both_nochim_sumtable)))
rownames(track_merged_data) <- sample.names
colnames(track_merged_data) <- c("Reads_w_chim","ASVs_w_chim" ,"Reads_nonchim", "ASVs_nonchim")
track_merged_data<-as.data.frame(track_merged_data) %>% 
  rownames_to_column("sample.names")
track_merged_data[is.na(track_merged_data)]<-0

track_both<-track.SS %>% 
  full_join(track.AS) %>%
  full_join(track_merged_data) %>% 
  replace(., is.na(.), 0) %>% 
  mutate(input=(input.SS+input.AS), denoisedF=(denoisedF.SS+denoisedF.AS), 
         denoisedR=(denoisedR.SS+denoisedR.AS)) %>% 
  select(sample.names,input, denoisedF, denoisedR, Reads_w_chim, Reads_nonchim, ASVs_w_chim, ASVs_nonchim)

write.xlsx(track_both,"track_reads_both.xlsx")

#Transpose table, assign names, extract sequences and saving table, for further processing:
trasBoth_nochim_sumtable <- as.data.frame(t(Both_nochim_sumtable))
#Get DNA sequences
sequences <- row.names(trasBoth_nochim_sumtable)
#Assign new rownames
for (i in 1:nrow(trasBoth_nochim_sumtable)){rownames(trasBoth_nochim_sumtable)[i] <- sha1(trasBoth_nochim_sumtable)}
tbname <- file.path(path,"DADA2_noChime.table")
{write.table(trasBoth_nochim_sumtable,tbname,sep="\t",col.names = NA, quote=FALSE)}
#Extract OTUs (sequences)
sinkname <- file.path(path,"DADA2_noChime_raw.otus")
{
  sink(sinkname)
  for (seqX in seq.int(nrow(trasBoth_nochim_sumtable))) {
    header <- paste0(">",rownames(trasBoth_nochim_sumtable)[seqX],"\n")
    cat(header)
    seqq <- paste0(sequences[seqX],"\n")
    cat(seqq)
  }
  sink()
}

#Define function to extract sequences sample-wise
extrSamDADA2 <- function(my_table, table2) {
  out_path <- file.path(path, "DADA2_extracted_samples_no_chim")
  if(!file_test("-d", out_path)) dir.create(out_path)
  for (sampleX in seq(1:dim(my_table)[1])){
    sinkname <- file.path(out_path, paste0(rownames(my_table)[sampleX],".fas"))
    {
      sink(sinkname)
      for (seqX in seq(1:dim(my_table)[2])) {
        if (my_table[sampleX,seqX] > 0) {
          header <- paste0(">",rownames(table2)[seqX],";size=",my_table[sampleX,seqX],";","\n")
          cat(header)
          seqq <- paste0(colnames(my_table)[seqX],"\n")
          cat(seqq)
        }
      }
      sink()
    }
  }
}

#Extract samplewise sequences from the non-chimera table using the above function:
extrSamDADA2(Both_nochim_sumtable, trasBoth_nochim_sumtable)

sink("sessionInfo.txt")
sessionInfo()
sink()
