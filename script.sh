#!/bin/bash

cd ~/Assignment1/fastq

mkdir fastqc_result

# Raw data quality check and put uncompressed output files in one directory
fastqc -t 64 -extract -outdir fastqc_result *.fq.gz



