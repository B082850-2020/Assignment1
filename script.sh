#!/bin/bash

#cd ~/Assignment1/fastq
#mkdir fastqc_result

# Raw data quality check and put uncompressed output files in one directory
#fastqc -t 64 -extract -outdir fastqc_result *.fq.gz

# Find the paths for fastqc output files
cd ~/Assignment1/fastq/fastqc_result
find -type f -name fastqc_data* > output.list
find -type f -name summary* >> output.list
sort -o output.list output.list

# Extract sequence numbers and quality details 
for i in `cat output.list`;
do cat $i|grep -m1 -i filename
cat $i|grep -m1 -i total
cat $i|grep -m1 "Sequence length"
cat $i|grep -m1 -i flag
cat $i|cut -f 1,2 |awk '{FS="\t"; if ($1 =="FAIL"||$1 =="WARN"&&$1 !="PASS"){print $0}}'
echo ------------;
done >> file1 # save the information in file1
cat file1
