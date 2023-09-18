#!/bin/bash


## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error


module purge
module load BLAST+/2.13.0-gompi-2022a

blastn -db [path to database] -query [YOUR CENTROID FILE HERE] -evalue 0.0001 -qcov_hsp_perc 80 -outfmt "6 qseqid sseqid pident qcovs evalue bitscore length mismatch gapopen qstart qend sstart send" -max_target_seqs 5 -num_threads 8 -out blast_hits.txt

