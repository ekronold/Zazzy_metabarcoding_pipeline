#!/bin/bash
mkdir DADA2_AS/samples_with_few_reads_AS | find DADA2_AS/*.fastq.gz -type f -size -5000c |while read file; do mv $file ${file%/*}/samples_with_few_reads_AS; done
mkdir DADA2_SS/samples_with_few_reads_SS | find DADA2_SS/*.fastq.gz -type f -size -5000c |while read file; do mv $file ${file%/*}/samples_with_few_reads_SS; done

