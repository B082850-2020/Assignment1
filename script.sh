#!/bin/bash

#cd ~/Assignment1/fastq
mkdir fastqc_result

# Raw data quality check and put uncompressed output files in a directory
fastqc -t 64 -extract -outdir fastqc_result *.fq.gz

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
cat file1 #show results on screen

# Format conversion from .fq to .fa
# Not sure if needed
#cd ~/Assignment1/fastq
#for fq in *.fq.gz
#do /localdisk/home/s1544765/seqtk/seqtk seq -a $fq> $fq\.fa
#done

# Prepare reference genome for alignment
cd /localdisk/data/BPSM/Assignment1/Tbb_genome/
cp Tb927_genome.fasta.gz ~/Assignment1/fastq

# Make indexed reference using bowtie2
cd ~/Assignment1/fastq
mkdir reference_index
gunzip Tb927_genome.fasta.gz
bowtie2-build --threads 64 Tb927_genome.fasta reference_index/

#
bowtie2 --threads 64 -x reference_index/ -1 216_L8_1.fq.gz -2 216_L8_2.fq -S gene_pair1.sam




