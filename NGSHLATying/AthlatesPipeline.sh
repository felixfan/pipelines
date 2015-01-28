#!/usr/bin/env bash

### typing HLA using exome sequencing data
### run on Ubuntu 14.04 LTS
### @ Felix Yanhui Fan felixfanyh@gmail.com
### @ 26 Jan 2015

### !!! before use it, check and modify all parts with "!!!"
################### how to use it ###############################
### ./Athlates.sh
### then
### bash run.sh
### or (parameters after -l must can be divided by 51, which is the lines of code to deal with one individual)
### split -l 5100 run.sh runs
### then 
### bash runsaa
### bash runsab
### ...
###################### set directory ############################
FQDIR="/home/felixfan/uwork/ADExomeHLA/results/OptiType/fastq"      # paired-end fastq files, name must be "fastq"
AthlatesDIR="/home/felixfan/ubin/Athlates_2014_04_26"
OUTDIR="/home/felixfan/uwork/ADExomeHLA/results/Athlates/results"
############### do not need to change this part #################
BINDIR="$AthlatesDIR/bin"
REFDIR="$AthlatesDIR/db/ref"
BEDDIR="$AthlatesDIR/db/bed"
MSADIR="$AthlatesDIR/db/msa"
#################################################################
if [ ! -d "$OUTDIR" ]
then
	mkdir $OUTDIR
fi
if [ -e run.sh ]
then
	rm run.sh
fi
#############################PREPARE#############################
cp $MSADIR/DQB_nuc.txt $MSADIR/DQB1_nuc.txt
cp $MSADIR/DRB_nuc.txt $MSADIR/DRB1_nuc.txt
cp $REFDIR/hla.clean.fasta ref.fasta
novoindex ref.nix ref.fasta
######################### BEGIN #################################
### 0 Get individual id (first part of the file name) and number of individuals
ls $FQDIR/*.fastq > fileN.txt
sed -i 's/\./ /' fileN.txt
awk '{print $1}' fileN.txt | sort | uniq > fileN2.txt
sed 's/fastq\// /' fileN2.txt > fileN.txt
awk '{print $2}' fileN.txt > fileN3.txt
IFS=$'\n' read -d '' -r -a lines < fileN2.txt   # dir and file names without extension
IFS=$'\n' read -d '' -r -a ids < fileN3.txt     # individual ids
n=$(ls $FQDIR/*.end1.fastq | wc -l)                           # number of individuals
let m=$n-1
#######################################################################
### novoalign
for i in $(seq 0 $m)
do
	echo "novoalign -d ref.nix -t 10 -o SAM -r All -l 80 -e 100 -i PE 200 140 -F STDFQ -f $FQDIR/${ids[$i]}.chr6.end1.fastq $FQDIR/${ids[$i]}.chr6.end2.fastq > ${ids[$i]}.sam" >> run.sh
	echo "samtools view -bS -h -F 4 ${ids[$i]}.sam > ${ids[$i]}.bam" >> run.sh
	echo "samtools sort ${ids[$i]}.bam ${ids[$i]}.sort" >> run.sh
	for g in {"A","B","C","DQB1","DRB1"}
	do
		echo "samtools view -b -L $BEDDIR/hla.$g.bed ${ids[$i]}.sort.bam > ${ids[$i]}.$g.bam" >> run.sh
		echo "samtools view -b -L $BEDDIR/hla.non-$g.bed ${ids[$i]}.sort.bam > ${ids[$i]}.non.$g.bam" >> run.sh
		echo "samtools view -h -o ${ids[$i]}.$g.sam ${ids[$i]}.$g.bam" >> run.sh
		echo "samtools view -h -o ${ids[$i]}.non.$g.sam ${ids[$i]}.non.$g.bam" >> run.sh
		echo "LC_ALL=C sort -r -k 1,1 -k 3,3 ${ids[$i]}.$g.sam > ${ids[$i]}.$g.sort.sam" >> run.sh
		echo "LC_ALL=C sort -r -k 1,1 -k 3,3 ${ids[$i]}.non.$g.sam > ${ids[$i]}.non.$g.sort.sam" >> run.sh
		echo "samtools view -bS ${ids[$i]}.$g.sort.sam > ${ids[$i]}.$g.sort.bam" >> run.sh
		echo "samtools view -bS ${ids[$i]}.non.$g.sort.sam > ${ids[$i]}.non.$g.sort.bam" >> run.sh
		echo "$BINDIR/typing -bam ${ids[$i]}.$g.sort.bam -exlbam ${ids[$i]}.non.$g.sort.bam -msa $MSADIR/${g}_nuc.txt -o $OUTDIR/${ids[$i]}.$g" >> run.sh
	done
	echo "rm ${ids[$i]}*.bam" >> run.sh
	echo "rm ${ids[$i]}*.sam" >> run.sh
	echo "rm ${ids[$i]}*.bti" >> run.sh
done


### CLEAR 
rm fileN*.txt


