#!/bin/bash

# Prepare fastq files
cd /localdisk/data/BPSM/Assignment1/
cp -r fastq ~/Assignment1
cd ~/Assignment1/fastq
mkdir fastqc_result


# Raw data quality check and put uncompressed output files in a directory
fastqc -t 64 -extract -outdir fastqc_result *.fq.gz
echo -e "----------------------------------- \n
fastqc done, files in directory fastqc_result"


# Find the paths for qc output files
cd ~/Assignment1/fastq/fastqc_result
find -type f -name fastqc_data* > output.list	# put all the paths of fastqc_data file in a list
find -type f -name summary* >> output.list	# add paths of summary file to the list
sort -o output.list output.list		# sort the list and overwrite the list 


# Extract sequence numbers and quality details 
for i in `cat output.list`;	# loop though all the paths in the list
do 
	cat $i|grep -m1 -i filename	# print filename
	cat $i|grep -m1 -i total	# print total sequence number
	cat $i|grep -m1 "Sequence length"	# print sequence length 
	cat $i|grep -m1 -i flag		# print poor quality sequence number
	cat $i|cut -f 1,2 |awk '{FS="\t"; if ($1 =="FAIL"||$1 =="WARN"&&$1 !="PASS"){print $0}}'	# find fail or warn in field 1 of summary.txt files 
	echo ------------;
done >> qc_feedback 	# save the information in a file
echo -e "Fastqc feedbacks are as below:"
cat qc_feedback 	#show results on screen


# Prepare reference genome for alignment
cd /localdisk/data/BPSM/Assignment1/Tbb_genome/
cp Tb927_genome.fasta.gz ~/Assignment1/fastq
cd ~/Assignment1/fastq


# Make indexed reference using bowtie2
mkdir reference_index
gunzip Tb927_genome.fasta.gz
echo -e "preparing indexed reference genome..."
bowtie2-build --threads 64 --quiet Tb927_genome.fasta reference_index/		# build index for reference genome


# Prepare sequence lists for pair-ended alignment
find -type f -name "*_1*.gz"|sort > odd.list	# generate 2 list for each one of the gene pairs 
find -type f -name "*_2*.gz"|sort > even.list
gene_pair_number=$(< "odd.list" wc -l) 		# variable as the length of odd.list


# Pair-ended alignment of gene pairs with reference genome
for number in $(seq $gene_pair_number);		# loop commands for the same times as length of odd.list  
do 
	sq1=$(sed -n $number\p odd.list)	# set variables as  each line of path lists
	sq2=$(sed -n $number\p even.list)
	echo -e "------------------------\n start pair-ended alignment of $sq1 and $sq2"
	bowtie2 --threads 64 -x reference_index/ -1 $sq1 -2 $sq2 -S gene_pair$number\.sam	# high-speed alignment, output sam files
	echo -e "Done alighnment of gene pair $number"
	samtools view -b -h -o gene_pair$number\.bam gene_pair$number\.sam	# convert sam files into bam files
	samtools sort gene_pair$number\.bam > gene_pair$number\.srt.bam		# sort bam files
	samtools index gene_pair$number\.srt.bam	# index bam files
	echo -e "SAM->BAM->INDEXED SORT BAM created for gene pair $number"
        echo ------------------------
done


# Prepare the file with gene location information
cd /localdisk/data/BPSM/Assignment1/
cp Tbbgenes.bed ~/Assignment1/fastq
cd ~/Assignment1/fastq


# Prepare a file for mean count
echo -e "Gene\tSlender_216\tSlender_218\tSlender_219\tStumpy_220\tStumpy_221\tStumpy_222" > count_mean.txt	# prepare header


# Generate number of reads
for number in $(seq $gene_pair_number);
do
        bedtools bamtobed -i gene_pair$number\.srt.bam > gene_pair$number\.bed		# conversion of sorted bam to bed 
        bedtools coverage -a Tbbgenes.bed -b gene_pair$number\.bed -mean |cut -f 4,7 > count$number.txt		# mean gene count for each alignment, output in seperate files
	echo -e "mean gene count for gene pair $number is done" 
done


# Make and format output file
paste count*.txt > combine.txt	# combine the columns 
cut -f 1,2,4,6,8,10,12 combine.txt >> count_mean.txt	# keep only one gene name column, append to the file with header
echo -e "------------------------------------- \n
Showing head of the final gene count output file...\n"
cat count_mean.txt| head
