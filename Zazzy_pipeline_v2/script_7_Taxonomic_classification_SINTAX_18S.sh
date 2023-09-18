#!/bin/bash


## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error


module purge
module load VSEARCH/2.22.1-GCC-11.3.0

vsearch --sintax [YOUR CENTROID FILE HERE] --db [path to database formatted for SINTAX] --sintax_cutoff 0.8 --tabbedout Taxonomy_sintax_18S.txt

module purge
module load R/4.2.2-foss-2022b

R < Combine_OTU_with_tax_18S.R --save

