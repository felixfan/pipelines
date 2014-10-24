#!/bin/bash

### typing HLA using exome sequencing data
### run on Ubuntu 14.04 LTS
### @ Felix Yanhui Fan felixfanyh@gmail.com
### @ 22 Oct 2014

### !!! before use it, check and modify all parts with "!!!"

#################### how to use it ##############################
### after run this script, copy all fastq file and run{1,10}.txt to the cluster server to run razers3
### cd $FQDIR
### scp *.fastq yanhui@statgenpro.psychiatry.hku.hk:/home/yanhui/ngs/
### scp run*.txt yanhui@statgenpro.psychiatry.hku.hk:/home/yanhui/ngs/
### on cluster server:
### qsub run1.txt
### qsub run2.txt
### ...
#################################################################

###################### tools installed ##########################
### picard-tools-1.123
### GATK 3.2-2
### samtools 1.1
### bedtools v2.17.0
### razers3 version 3.2 [13859]
### OptiType version 1.0
###################### set directory ############################
### !!! change this part !!! ###
OptiTypeDIR="/home/fan/optiType"         # this is the directory contains OptiType pipeline
DATADIR="/home/fan/ADBWA"                # this directory contains a directory named "bam" contains all bam and bai file, all temporary dir will be created under this dir
############### do not need to change this part #################
REF="$OptiTypeDIR/data/hla_reference_dna.fasta"
BAMDIR="$DATADIR/bam"
TEMPDIR="$DATADIR/tempdir"
FQDIR="$DATADIR/fastq"
mkdir $TEMPDIR
mkdir $FQDIR
######################### BEGIN #################################
### 0 Get individual id (first part of the file name) and number of individuals
ls $BAMDIR/*.bam > fileN.txt
sed -i 's/\./ /' fileN.txt
awk '{print $1}' fileN.txt > fileN2.txt
sed 's/fastq\// /' fileN2.txt > fileN.txt
awk '{print $2}' fileN.txt > fileN3.txt
IFS=$'\n' read -d '' -r -a lines < fileN2.txt   # dir and file names without extension
IFS=$'\n' read -d '' -r -a ids < fileN3.txt     # individual ids
n=$(ls $BAMDIR/*.bam | wc -l)                   # number of individuals
let m=$n-1
#################################################################
### 1 Extract chr6 from BAM
### samtools view -b myfile.bam chr6 > myfile.chr6.bam
### !!! in this analysis "6" instead of "chr6" was used, it depends on the tag in your bam file !!!
for i in $(seq 0 $m)
do
	echo "samtools view -h -b $BAMDIR/${ids[$i]}.bam 6 > $TEMPDIR/${ids[$i]}.chr6.bam" >> batch1ExtractChr6.sh
done

### 2 sort bam by read name
### One can sort the BAM file by query name with samtools
### samtools sort -n aln.bam aln.qsort
for i in $(seq 0 $m)
do
	echo "samtools sort -n $TEMPDIR/${ids[$i]}.chr6.bam $TEMPDIR/${ids[$i]}.chr6.qsort" >> batch2Sort.sh
done

### 3 convert bam to fastq
### bedtools -fq2 Creating two FASTQ files for paired-end sequences
### When using this option, it is required that the BAM file is sorted/grouped by the read name. 
### This keeps the resulting records in the two output FASTQ files in the same order.
for i in $(seq 0 $m)
do
	echo "bedtools bamtofastq -i $TEMPDIR/${ids[$i]}.chr6.qsort.bam -fq $FQDIR/${ids[$i]}.chr6.end1.fastq -fq2 $FQDIR/${ids[$i]}.chr6.end2.fastq" >> batch3BAM2Fastq.sh
done

### 4 filtering reads using razers3
### submit jobs to server
echo '#!/bin/bash' >> qsub.txt
echo '#PBS -l nodes=1' >> qsub.txt               # !!! require one node
echo '#PBS -l walltime=168:00:00' >> qsub.txt    # !!! max time: 168hrs !!! reduce number of individuals for each sumbit, or it will be killed after 168 hrs
echo '#PBS -m abe' >> qsub.txt                   # enable email notify
echo '#PBS -q default' >> qsub.txt               # Queue name "default"
echo '#PBS -N anyName' >> qsub.txt               # name of the job
echo 'cd $PBS_O_WORKDIR' >> qsub.txt             # change to the dir where you submit your job
for i in $(seq 0 $m)
do
	echo "./razers3 -i 90 -m 1 -dr 0 -o ./${ids[$i]}.chr6.end1.fished.sam ./hla_reference_dna.fasta ./${ids[$i]}.chr6.end1.fastq" >> qsub.txt
	echo "./razers3 -i 90 -m 1 -dr 0 -o ./${ids[$i]}.chr6.end2.fished.sam ./hla_reference_dna.fasta ./${ids[$i]}.chr6.end2.fastq" >> qsub.txt
done





