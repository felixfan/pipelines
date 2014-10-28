#!/bin/bash

### typing HLA using exome sequencing data
### run on Ubuntu 14.04 LTS
### @ Felix Yanhui Fan felixfanyh@gmail.com
### @ 22 Oct 2014

### !!! before use it, check and modify all parts with "!!!"

############### install python modules and other tolls ##########
### https://github.com/felixfan/OptiType
### sudo pip install biopython
### sudo pip install Coopr
### sudo pip install matplotlib
### sudo pip install pandas
### Install hdf5 (http://www.hdfgroup.org/HDF5/)
### Install pyTables (http://www.pytables.org) 
### Install RazerS (http://www.seqan.de/projects/razers/)
### Install glpk: http://ftp.gnu.org/gnu/glpk/
#################### how to use it ##############################
### step 1
### bash OptiTypePipeline.sh
### bash batch1ExtractChr6.sh
### bash batch2Sort.sh
### bash batch3BAM2Fastq.sh
### step 2
### after run step 1, copy all fastq file and run{1,10}.txt to the cluster server to run razers3
### cd $FQDIR
### scp *.fastq yanhui@statgenpro.psychiatry.hku.hk:/home/yanhui/ngs/
### scp run*.txt yanhui@statgenpro.psychiatry.hku.hk:/home/yanhui/ngs/
### on cluster server:
### qsub run1.txt
### qsub run2.txt
### ...
### step 3
### after run razers3, copy all sam file to local ubuntu
### cd $FQDIR
### scp yanhui@statgenpro.psychiatry.hku.hk:/home/yanhui/ngs/*.sam .
### bash batch4Sam2fastq.sh
### copy batch5OptiType.sh to the directory of OptiType and run it from there
### bash batch5OptiType.sh
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
OptiTypeDIR="/home/fan/optiType"     # this is the directory contains OptiType pipeline
BAMDIR="/home/fan/ADBWA/bam"   # this directory must be named "bam" and contains all bam and bai file 
DATADIR="/home/fan/ADBWA"      # all temporary dir will be created under this dir
OUTDIR="/home/fan/optiType/AD2014"  # final HLA typing results
############### do not need to change this part #################
REF="$OptiTypeDIR/data/hla_reference_dna.fasta"        # reference sequence
TEMPDIR="$DATADIR/tempdir"                             # chr6 bam files, sorted chr6 bam files
FQDIR="$DATADIR/fastq"                                 # paired-end fastq files
SAMDIR="$DATADIR/sam"                                  # razers3 output sam files
FFQDIR="$DATADIR/fishedFq"                             # fished fastq files, input of OptiType
#################################################################
if [ ! -d "$TEMPDIR" ]
then
	mkdir $TEMPDIR
fi
###
if [ ! -d "$FQDIR" ]
then
	mkdir $FQDIR
fi
###
if [ ! -d "$SAMDIR" ]
then
	mkdir $SAMDIR
fi
###
if [ ! -d "$FFQDIR" ]
then
	mkdir $FFQDIR
fi
#################################################################
### 1
if [ -e batch1ExtractChr6.sh ]
then
	rm batch1ExtractChr6.sh
fi
### 2
if [ -e batch2Sort.sh ]
then
	rm batch2Sort.sh
fi
### 3
if [ -e batch3BAM2Fastq.sh ]
then
	rm batch3BAM2Fastq.sh
fi
### 4
if [ -e batch4Sam2fastq.sh ]
then
	rm batch4Sam2fastq.sh
fi
### 5
if [ -e batch5OptiType.sh ]
then
	rm batch5OptiType.sh
fi
######################### BEGIN #################################
### 0 Get individual id (first part of the file name) and number of individuals
ls $BAMDIR/*.bam > fileN.txt
sed -i 's/\./ /' fileN.txt
awk '{print $1}' fileN.txt > fileN2.txt
sed 's/bam\// /' fileN2.txt > fileN.txt
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
	echo "samtools view -h -b ${lines[$i]}.bam 6 > $TEMPDIR/${ids[$i]}.chr6.bam" >> batch1ExtractChr6.sh
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
### submit jobs to server:
### qsub qsub.xxxx.txt
### each qsub.xxxx.txt is for one individual
for i in $(seq 0 $m)
do
	echo '#!/bin/bash' >> qsub.${ids[$i]}.txt
	echo '#PBS -l nodes=1' >> qsub.${ids[$i]}.txt             # !!! require one node
	echo '#PBS -l walltime=168:00:00' >> qsub.${ids[$i]}.txt  # !!! max time: 168hrs !!! 
	echo '#PBS -m abe' >> qsub.${ids[$i]}.txt                 # enable email notify
	echo '#PBS -q default' >> qsub.${ids[$i]}.txt             # Queue name "default"
	echo '#PBS -N '${ids[$i]} >> qsub.${ids[$i]}.txt             # name of the job
	echo 'cd $PBS_O_WORKDIR' >> qsub.${ids[$i]}.txt           # change to the dir where you submit your job
	echo "./razers3 -i 90 -m 1 -dr 0 -o ./${ids[$i]}.chr6.end1.fished.sam ./hla_reference_dna.fasta ./${ids[$i]}.chr6.end1.fastq" >> qsub.${ids[$i]}.txt
	echo "./razers3 -i 90 -m 1 -dr 0 -o ./${ids[$i]}.chr6.end2.fished.sam ./hla_reference_dna.fasta ./${ids[$i]}.chr6.end2.fastq" >> qsub.${ids[$i]}.txt	
done

### 5 convert sam to fastq
### convert the sam to bam, then bam to fastq
for i in $(seq 0 $m)
do
	echo "samtools view -Sb $SAMDIR/${ids[$i]}.chr6.end1.fished.sam > $SAMDIR/${ids[$i]}.chr6.end1.temp.bam" >> batch4Sam2fastq.sh
	echo "samtools view -Sb $SAMDIR/${ids[$i]}.chr6.end2.fished.sam > $SAMDIR/${ids[$i]}.chr6.end2.temp.bam" >> batch4Sam2fastq.sh
	echo "bedtools bamtofastq -i $SAMDIR/${ids[$i]}.chr6.end1.temp.bam -fq $FFQDIR/${ids[$i]}.end1.fished.fastq" >> batch4Sam2fastq.sh
	echo "bedtools bamtofastq -i $SAMDIR/${ids[$i]}.chr6.end2.temp.bam -fq $FFQDIR/${ids[$i]}.end2.fished.fastq" >> batch4Sam2fastq.sh
done

### 6 HLA typing using OptiType
### copy batch5OptiType.sh to the directory of OptiType and run it from there
for i in $(seq 0 $m)
do
	if [ ! -d "$OUTDIR/${ids[$i]}" ]
	then
		mkdir $OUTDIR/${ids[$i]}
	fi
	echo "python OptiTypePipeline.py -i $FFQDIR/${ids[$i]}.end1.fished.fastq $FFQDIR/${ids[$i]}.end2.fished.fastq -d -v -o $OUTDIR/${ids[$i]}" >> batch5OptiType.sh
done 

### CLEAR 
rm fileN*.txt


