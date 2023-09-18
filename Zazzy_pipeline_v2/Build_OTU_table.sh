#!/bin/bash

#Script for building OTU-table with output from SWARM_clust.sh, adapted from Frédéric Mahe at https://github.com/frederic-mahe/swarm/wiki/Fred's-metabarcoding-pipeline
#This script calls a python (.py) script. IMPORTANT: This is not the same as the .py script used in Luis's pipeline
#Make sure the correct python script is being called by setting the path in SCRIPT=".."
#Check that all the names otherwise conform to the output from the SWARM_clust.sh script
#Only thing that might need to be changed is the last line, S[ ]*.fas needs to match the sample names

FASTA="COOL_PROJECT_NAME.fas" #Fasta file name needs to match the final fasta output from SWARM script
SCRIPT="/cluster/projects/nn9338k/eivind/scripts/back_ground_scripts/OTU_contingency_table.py"
STATS="${FASTA/.fas/_1f.stats}"
SWARMS="${FASTA/.fas/_1f.swarms}"
REPRESENTATIVES="${FASTA/.fas/_1f_representatives.fas}"
UCHIME="${FASTA/.fas/_1f_representatives.uchime}"
ASSIGNMENTS="${FASTA/.fas/_dummytax_representatives.results}"
OTU_TABLE="${FASTA/.fas/.OTU.table}"

python \
    "${SCRIPT}" \
    "${REPRESENTATIVES}" \
    "${STATS}" \
    "${SWARMS}" \
    "${UCHIME}" \
    "${ASSIGNMENTS}" \
    S[0-9][0-9][0-9].fas > "${OTU_TABLE}"
