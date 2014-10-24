#!/bin/bash

### typing HLA using exome sequencing data
### run on Ubuntu 14.04 LTS
### @ Felix Yanhui Fan felixfanyh@gmail.com
### @ 22 Oct 2014

######################################################
### picard-tools-1.123
### GATK 3.2-2
### samtools 1.1
### bedtools v2.17.0
### razers3 version 3.2 [13859]
### OptiType version 1.0
######################################################
FQDIR="/home/fan/ADBWA/run/fastq"
OUTDIR="/home/fan/optiType/AD2014"
#####################################################################
### after run razers3, copy all sam file to local ubuntu
### cd $FQDIR
### scp yanhui@statgenpro.psychiatry.hku.hk:/home/yanhui/ngs/*.sam .
#####################################################################
### 5 convert sam to fastq
ls $FQDIR/*.end1.fastq > fileN.txt
sed -i 's/\./ /' fileN.txt
awk '{print $1}' fileN.txt > fileN2.txt
sed 's/fastq\// /' fileN2.txt > fileN.txt
awk '{print $2}' fileN.txt > fileN3.txt
IFS=$'\n' read -d '' -r -a lines < fileN2.txt
IFS=$'\n' read -d '' -r -a ids < fileN3.txt
n=$(ls $FQDIR/*.end1.fastq | wc -l)
let m=$n-1
### convert the sam to bam, then bam to fastq
for i in $(seq 0 $m)
do
	echo "samtools view -Sb ${lines[$i]}.chr6.end1.fished.sam > ${lines[$i]}.end1.temp.bam" >> batch4Sam2fastq.sh
	echo "samtools view -Sb ${lines[$i]}.chr6.end2.fished.sam > ${lines[$i]}.end2.temp.bam" >> batch4Sam2fastq.sh
	echo "bedtools bamtofastq -i ${lines[$i]}.end1.temp.bam -fq ${lines[$i]}.end1.fished.fastq" >> batch4Sam2fastq.sh
	echo "bedtools bamtofastq -i ${lines[$i]}.end2.temp.bam -fq ${lines[$i]}.end2.fished.fastq" >> batch4Sam2fastq.sh
done
# sh batch4Sam2fastq.sh

### 6 HLA typing using OptiType
### copy batch5OptiType.sh to the directory of OptiType and run it from there
for i in $(seq 0 $m)
do
	mkdir $OUTDIR/${ids[$i]}
	echo "python OptiTypePipeline.py -i ${ids[$i]}.end1.fished.fastq ${ids[$i]}.end2.fished.fastq -d -v -o $OUTDIR/${ids[$i]}" >> batch5OptiType.sh
done 
