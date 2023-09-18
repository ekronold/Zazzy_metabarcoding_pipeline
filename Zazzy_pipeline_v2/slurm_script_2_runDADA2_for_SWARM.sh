#!/bin/bash

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module purge
module load  R/4.1.2-foss-2021b
R < script_3_dependency_R_code_sha1_hash_names.R --save

