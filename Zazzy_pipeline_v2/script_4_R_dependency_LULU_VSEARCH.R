#Set the library path to ensure loading correct packages
#Optional, you will generate personal library folders in R by installing packages. To wensure that the corrects packages you will use
#are loaded in the scripts this is advisable to include. The path here will route to the users home folders
.libPaths(c(.libPaths(), "~/R/x86_64-pc-linux-gnu-library/3.6"))

#Install lulu if not already done:
#library(devtools)
#install_github("tobiasgf/lulu")

#Load libraries
library(tidyverse)
library(lulu)
library(openxlsx)

path <- getwd()
OTU_table_path <- list.files(path, pattern="otutable", full.names=TRUE)

OTU_table <- read.csv(OTU_table_path, sep = "\t")

OTU_table$OTUId <- gsub("\\;.*", "", OTU_table$OTUId)

OTU_table1 <- OTU_table %>%
  #Set the amplicon ID as rownames for LULU curation
  column_to_rownames("OTUId")

match_list_path<-list.files(path, pattern = "list.txt", full.names = TRUE)
match_list<-read.csv(match_list_path, sep="\t", header = FALSE)
lulu_curation<-lulu(OTU_table1, match_list)

OTU_table_lulu_curated<-lulu_curation$curated_table %>%
  rownames_to_column("OTUId")


write.xlsx(OTU_table_lulu_curated,"OTU_table_lulu_curated.xlsx")
saveRDS(lulu_curation, "lulu_analyses.R")

#Write out a new centroid file for taxonomic annotation based on just the sequences that are kept after LULU
Seqs <- read.table("Solhomfjell_vsearch.centroids")
Seqs <- as.data.frame(matrix(Seqs$V1, ncol = 2, byrow = TRUE, dimnames = list(NULL, c("OTUId", "Sequence"))))

Seqs$OTUId <- gsub("\\;.*", "", Seqs$OTUId)
Seqs <- Seqs %>% mutate(OTUId = gsub(">", "", OTUId))

Centroids_lulu <- Seqs[which(Seqs$OTUId %in% OTU_table_lulu_curated$OTUId), ]

Centroids_lulu$OTUId <- paste0(">", Centroids_lulu$OTUId)

#Some rearranging necessary to make a fasta file, arranging by X to ensure the correct order
Sequence <- as.data.frame(Centroids_lulu[,2])
names(Sequence) <- "Sequence"
Sequence <- Sequence %>% add_column("X"=seq(1:nrow(Sequence)))

SeqNames <- as.data.frame(Centroids_lulu[,1])
names(SeqNames) <- "Sequence"
SeqNames <- SeqNames %>% add_column("X"=seq(1:nrow(SeqNames)))

FASTA <- bind_rows(SeqNames, Sequence) %>% arrange(X)
#Extract only the first column in order to write out the fasta file
FASTA1 <- data.frame("Sequence" = FASTA[,1])

#Write the new LULU curated centroid file for taxonomic annotation
write.table(FASTA1, "Centroid_lulu_curated.fas", quote = FALSE, col.names = FALSE, row.names = FALSE)

