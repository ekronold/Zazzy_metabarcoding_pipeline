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
OTU_table_path <- list.files(path, pattern="OTU.table", full.names=TRUE)

OTU_table <- read.csv(OTU_table_path, sep = "\t")

#Some cleanup of the raw OTU table
OTU_table2 <- OTU_table %>% 
  #Remove the OTUs tagged as chimeric
  dplyr::filter(chimera == "N") %>%
  #Remove the SWARM stats columns 
  select(-c("OTU", "total", "cloud", "length", "abundance", "chimera", "spread", "sequence", "identity", "taxonomy"))

#Save the SWARM stats to tag on the OTU table after curation
SWARM_stats <- OTU_table %>% 
  #Remove the OTUs tagged as chimeric
  dplyr::filter(chimera == "N") %>%
  #Remove the SWARM stats columns 
  select(c("total", "cloud", "length", "abundance", "chimera", "spread", "sequence", "identity", "taxonomy"))

OTU_table3 <- OTU_table2 %>%
  #Set the amplicon ID as rownames for LULU curation
  column_to_rownames("amplicon")

match_list_path<-list.files(path, pattern = "list.txt", full.names = TRUE)
match_list<-read.csv(match_list_path, sep="\t", header = FALSE)
lulu_curation<-lulu(OTU_table3, match_list)

OTU_table_lulu_curated<-lulu_curation$curated_table %>%
  rownames_to_column("OTUid")

#Reattach the SWARM stats for the sequences that are kept
SWARM_stats <- SWARM_stats %>% rename(OTUid = "identity")
OTU_table4 <- left_join(OTU_table_lulu_curated, SWARM_stats, by = "OTUid")

write.xlsx(OTU_table4,"OTU_table_lulu_curated.xlsx")
saveRDS(lulu_curation, "lulu_analyses.R")

#Write out a new centroid file for taxonomic annotation based on just the sequences that are kept after LULU
SeqNames <- OTU_table4$OTUid
SeqNames <- paste0(">", SeqNames)
SeqNames <- data.frame(cbind("SeqNames"=SeqNames, "X"=seq(1:length(SeqNames))))


Sequence <- OTU_table4$sequence
Sequence <- data.frame(cbind("SeqNames"=as.character(Sequence), "X"=seq(1:length(Sequence))))

#Combine the SeqId and the actual Sequence together in one row, arranging by X to ensure the correct order
FASTA <- bind_rows(SeqNames, Sequence) %>% arrange(X)
#Extract only the first column in order to write out the fasta file
FASTA1 <- data.frame("Sequence" = FASTA[,1])

#Write the new LULU curated centroid file for taxonomic annotation
write.table(FASTA1, "Centroid_lulu_curated.fas", quote = FALSE, col.names = FALSE, row.names = FALSE)

