#!/bin/bash

### typing HLA using exome sequencing data
### run on Ubuntu 14.04 LTS
### @ Felix Yanhui Fan felixfanyh@gmail.com
### @ 28 Oct 2014

### !!! before use it, check and modify all parts with "!!!"

### This pipeline starts from paired-end fastq
### If you only have bam exome data, see OptiTypePipeline part 1-3 to generate fastq files
########## how to use it ############################################
### bash PHLATpipeline.sh
#####################################################################
### phlat-release, need to ask the author to get one copy of this software and data
### bowtie 2.2.4
####################### install python modules ######################
### sudo pip install pysam
#####################################################################
## -tag: name label for the sample associated with the fastq files
## -p: number of threads for running Bowtie2 [default 8]
## -e: url to the home folder of phlat-1.0 package
## -o: url to a directory where results shall be stored
## -pe: flag indicating whether the data shall be treated as paired-end(1) or single-end(0) [default 1]
#####################################################################
phlatdir=/home/fan/phlat-release          # phlat dir
indexdir=/home/fan/phlat-release/b2folder # data with phlat
b2url=/home/fan/bowtie2-2.2.4/bowtie2     # bowtie2 dir
datadir=/home/fan/ADBWA/fastq             # paired-end fastq dir, name must be "fastq" !!!
rsdir=/home/fan/ADExomeHLA/phlat/results  # output dir
#####################################################################
if [ ! -d "$rsdir" ]
then
	mkdir $rsdir
fi
#####################################################################
### 0 Get individual id (first part of the file name) and number of individuals
ls $datadir/*.end1.fastq > fileN.txt
sed -i 's/\./ /' fileN.txt
awk '{print $1}' fileN.txt > fileN2.txt
sed 's/fastq\// /' fileN2.txt > fileN.txt
awk '{print $2}' fileN.txt > fileN3.txt
IFS=$'\n' read -d '' -r -a lines < fileN2.txt   # dir and file names without extension
IFS=$'\n' read -d '' -r -a ids < fileN3.txt     # individual ids
n=$(ls $datadir/*.end1.fastq | wc -l)           # number of individuals
let m=$n-1
### 1 run PHLAT
for i in $(seq 0 $m)
do
	tag=${ids[$i]}       # name label for the sample associated with the fastq files
	fastq1=${tag}".chr6.end1.fastq"           
	fastq2=${tag}".chr6.end2.fastq"
	python2.7 -O ${phlatdir}/dist/PHLAT.py -1 ${datadir}/${fastq1} -2 ${datadir}/${fastq2} -index $indexdir -b2url $b2url -tag $tag -e $phlatdir -o $rsdir -p 1
done


