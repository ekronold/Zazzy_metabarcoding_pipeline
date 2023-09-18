#!/bin/bash

module purge
ml load BLAST+/2.13.0-gompi-2022a

cp [INSERT CENTROID FILE FROM POST CLUSTERING] OTU_centroids

sed -i 's/;.*//' OTU_centroids

makeblastdb -in OTU_centroids -parse_seqids -dbtype nucl

blastn -db OTU_centroids -outfmt '6 qseqid sseqid pident' -out match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query OTU_centroids

#Run the lulu curation:
ml purge
ml load R/4.2.1-foss-2022a

#Select one R script based on clustering method
#R < script_4_R_dependency_LULU_SWARM.R --save
#R < script_4_R_dependency_LULU_VSEARCH.R --save

