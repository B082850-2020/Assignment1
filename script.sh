#!/bin/bash

cd ~/Assignment1/fastq

mkdir fastqc_result

# Raw data quality check and put uncompressed output files in one directory
fastqc -t 64 -extract -outdir fastqc_result *.fq.gz

# Assess number and quality 
cd ~/Assignment1/fastq/fastqc_result
mk fastqc_output
for dirs, files in fastqc_result:
for name in files:
if (name == "fastqc_data.txt"):
cat fastqc_data.txt |grep -m1 -i filename > fastqc_output
cat fastqc_data.txt |grep -m1 -i total >> fastqc_output
cat fastqc_data.txt |grep -m1 -i quality >> fastqc_output
fi
if (name =="summary.txt"):
cat summary.txt |cut -f 1,2 |awk '{FS="\t"; if ($1 =="FAIL"||$1 =="WARN"){print $0}}' >>fastqc_output
fi
