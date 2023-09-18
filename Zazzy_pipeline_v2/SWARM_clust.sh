#!/bin/bash

#VSEARCH=$(which vsearch)
#SWARM=$(which swarm)
TMP_FASTA=$(mktemp --tmpdir=".")
FINAL_FASTA="COOL_PROJECT_NAME.fas"
module purge

# Pool sequences, adapt to sample names your study
cat S[0-9][0-9][0-9].fas > "${TMP_FASTA}"

# Dereplicate (vsearch)
ml VSEARCH/2.21.1-GCC-10.3.0
vsearch      --derep_fulllength "${TMP_FASTA}" \
             --sizein \
             --sizeout \
             --fasta_width 0 \
             --output "${FINAL_FASTA}" > /dev/null

rm -f "${TMP_FASTA}"

# Clustering
THREADS=16
TMP_REPRESENTATIVES=$(mktemp --tmpdir=".")
module purge 
ml swarm/3.1.3-foss-2022a
swarm  \
    -d 1 -f -t ${THREADS} -z \
    -i ${FINAL_FASTA/.fas/_1f.struct} \
    -s ${FINAL_FASTA/.fas/_1f.stats} \
    -w ${TMP_REPRESENTATIVES} \
    -o ${FINAL_FASTA/.fas/_1f.swarms} < ${FINAL_FASTA}

# Sort representatives
module purge
ml VSEARCH/2.21.1-GCC-10.3.0
vsearch      --fasta_width 0 \
             --sortbysize ${TMP_REPRESENTATIVES} \
             --output ${FINAL_FASTA/.fas/_1f_representatives.fas} --sizein --sizeout

# Chimera checking, VSEARCH does this in a different way than DADA2 and tends to catch some chimeras that are overlooked. Useful to run both.
REPRESENTATIVES=${FINAL_FASTA/.fas/_1f_representatives.fas}
UCHIME=${REPRESENTATIVES/.fas/.uchime}
vsearch      --uchime_denovo "${REPRESENTATIVES}" \
             --uchimeout "${UCHIME}"

# Prepare to make table without taxonomic assignment. This can be done separately with PROTAX or BLAST on the centroids and added later
#Assign dummy taxonomy, remember to set the correct file names
grep "^>" ${FINAL_FASTA/.fas/_1f_representatives.fas} | \
    sed 's/^>//
         s/;size=/\t/
         s/;$/\t100.0\tNA\tNA/' > ${FINAL_FASTA/.fas/_dummytax_representatives.results}

rm -f $FINAL_FASTA.uc "${TMP_REPRESENTATIVES}"


