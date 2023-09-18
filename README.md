# Zazzy_metabarcoding_pipeline
Pipeline for analysing ITS and 18S amplicon data used in Oslo Mycology Group

UPDATED PIPELINE FOR OSLO MYCOLOGY GROUP
BY: Eivind Kverme Ronold - 26 June 2023
LAST EDIT: 8 September 2023
Everything you need for the complete pipeline is IN THIS FOLDER

It is best to copy (cp) this entire folder into your own work area since you will be setting your own names on mapping and input/output files
Runs in a unix environment with R, python and perl dependencies

Run the scripts in sequence. Demultiplexing and denoising takes quite a some time to run thorugh

Required modules:
cutadapt v4.2
R v4.2 (R package installation in the R scripts provided)
ITSx v1.1
VSEARCH v2.21
SWARM v3.1.3
BLAST v2.13

Check all the $PATH variables for correct paths to scripts

NOTE: After DADA2 denoising it is possible to run either VSEARCH 97% clustering or SWARM clustering

SWARM is a single-linkage clustering method. This builds clusters based on a centroid and amplicons with 1 bp difference for each step sequentially
-----Details in doi: 10.7717/peerj.593
-----It is possible to increase number of accepted differences by changing 'd' in the base swarm script (default d=1)

VSEARCH cluster is a clustering method based on a set threshold of similarity. It takes the most abundant centroid and lumps amplicons of a certain similarity into an OTU
-----Default threshold is 97%, can be changed in the VSEARCH cluster script
-----Details in doi: 10.7717/peerj.2584


NOTE ON TAXONOMIC ASSIGNMENT:
There are two options for running taxonomic assignment
SINTAX is an algorithm that predicts the taxonomic assignment of a query sequence against a set of reference sequences and outputs a bootstrap value for each taxonomic level
-----Developed for USEARCH, is implemented in VSEARCH with some fewer options than currently available in USEARCH (no important functions are missing)
-----Details of algorithm here: doi: https://doi.org/10.1101/074161

BLAST 
-----Can be run in parallell to compare and control SINTAX output if you have many unknown or unassigned sequences
-----Details on how blast works on the NCBI website


STEP 1 - Demultiplexing raw reads
:
Run slurm_script_1_dem_pairedadapters.sh, this step takes a long time
------Demultiplexing using cutadapt on the raw sequence files
------This script also checks for matching adapters in both ends to account and control for tag-switching
------Set up the batch_file with primer and tag sequences, look at example files in this folder for how to set them up

STEP 2 - Removal of samples with very low read count
:
Run script_2_remove_demult_files_with_very_low_read_number.sh
------This script removes the samples with very low read counts
------IMPORTANT: After running this, check that all your samples in both DADA2_SS and DADA2_AS have both an _R1 and an _R2 .fastq file
------Sometimes one (R1 or R2) go past the threshold and the other dont. This throws an error in the dada2 denoising.
------If you find an unpaired _R1.fastq or _R2.fastq, move it into the samples_with_few_reads subfolder.


STEP 3 - DADA2 amplicon denoising
:
Run slurm_script_2_runDADA2_for_SWARM.sh, this step takes a long time
------Denoised sequences are given unique sha1-hash names.
------These hashes are unique PER SEQUENCE, so the same sequence in different samples will have the same hash
------This is required in order to cluster using SWARM in the next steps and has no impact on how VSEARCH runs


STEP 4 - Clean up sequences with ITSx before clustering ONLY WHEN WORKING WITH ITS AMPLICONS
:
-----Run the script parallell_itsx_fungi.sh
-----New sample fasta files will be generated with only ITS sequence, old sample files are saved in the new folder "uncut_fasta"


######  Now you can choose clustering method for the denoised amplicons

###################
###### SWARM ######
###################

STEP 5.1 - SWARM clustering of DADA2 output
:
-----Run the script SWARM_clust.sh within the DADA2_extracted_samples_no_chim folder
-----Remember to set the name for the output files in the beginning of SWARM_clust.sh (FINAL_FASTA="COOL_PROJECT_NAME.fas")


STEP 5.2 - Build an OTU table based on SWARM output
:
Run the script Build_OTU_table.sh in the DADA2_extracted_samples_no_chim folder
-----File COOL_PROJECT_NAME.OTU.table is the finished OTU Table
-----You can set the name for your output files in the beginning of the SWARM.clust.sh script
-----Go to STEP 6

###################
##### VSEARCH #####
###################

STEP 5 - VSEARCH clustering and building OTU table
:
Run script_5_VSEARCH_cluster.sh inside the DADA2_extracted_samples_no_chim folder
-----Set the names for your output files in the shell script
-----Both a centroid file and a OTU_table is generated with the names you set

##### After either SWARM or VSEARCH

STEP 6 - LULU post-clustering curation of OTU table
:
Run the script script_6_lulu_curation.sh
-----cp the R dependency script (script_4_R_dependency_LULU.R) into DADA2_extracted_samples_no_chim before running the shell script
-----Make sure that the right centroid file is being used (centroid output from either SWARM or VSEARCH)
-----Make sure to run the correct R dependency file (VSEARCH and SWARM OTU_tables need different treatments before the LULU curation)
-----You will have to check the names for both the centroid file and the OTU table as VSEARCH and SWARM by default has slightly different naming conventions
-----stay inside the DADA2_extracted_samples_no_chim when running the shell script
-----Since LULU removes a lot of sequences, this script also outputs a new centroid file for taxonomic annotation
-----Old centroid file is NOT deleted, so if you want to compare pre- and post-lulu that is always possible

STEP 7 - Taxonomic assignment using BLAST or SINTAX
NOTE: Databases needs to be properly formatted, SINTAX has some strict requirements for the FASTA header (explained in the VSEARCH documentation)


###################
#####  BLAST  #####
###################

-----Stay in the DADA2_extracted_samples_no_chim folder
-----Run the appropriate BLAST script (ITS or 18S)
-----Smoothest to copy the appropriate script into the same folder as your centroid file
-----If running scripts as is, the centroid file is: Centroid_lulu_curated.fas, insert this in the blast script
-----If there is very low returns from blast, consider lowering or removing the -qcov_hsp_perc (Set to only accept matches that cover 80% or more of the reference)
-----Output file is "blast_hits.txt", 5 hits per sequence (you can change this if you want, but 5 is recommended)
-----Both BLAST hits and the OTU_table has the same hash ID for the same sequences
-----Currently you will have to filter and combine the OTU table and taxonomy manually yourself if using BLAST

###################
#####  SINTAX #####
###################

-----Stay in the DADA2_extracted_samples_no_chim folder
-----Run the appropriate SINTAX script (ITS or 18S), make sure to path the correct database
-----Smoothest to copy the script into the same folder as your centroid
-----Also copy the R dependency script to combine taxonomy and OTU table
-----Remember to set the name of the centroid file to use and the name of the output file
-----Command --sintax_cutoff will be the minimum bootstrap value for including in the last column of the output.
-----The output will in any case have a column with taxonomic assignment and bootstrap value for each sequence that matched anything in the database at all taxonomic ranks
-----Explanation of output columns:
	1: query ID (hash) 2: Predicted match with bootstrap value per rank 3: strand (always "+" here as we only use the forward strand) 4: Match above cutoff bootstrap value
----Outputs OTU_with_taxonomy_[18S][ITS2].txt depending on your marker, finished OTU-table for downstream analyses in tab-demilited txt format

##### DONE
