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
REF="/home/fan/optiType/data/hla_reference_dna.fasta"
BAMDIR="/home/fan/ADBWA"
TEMPDIR="/home/fan/ADBWA/tempdir"
FQDIR="/home/fan/ADBWA/fastq"
mkdir $TEMPDIR
mkdir $FQDIR
######################### BEGIN ##############################
### 1 Extract chr6 from BAM
### samtools view -b myfile.bam chr6 > myfile.chr6.bam
ls $BAMDIR/*.bam > fileN.txt
cp fileN.txt fileN2.txt
sed -i 's/\// /g' fileN2.txt
sed -i 's/\./ /' fileN2.txt
awk '{print $(NF-1)}' fileN2.txt > fileN3.txt
paste fileN.txt fileN3.txt > fileN2.txt
# awk -v var="$TEMPDIR" '{print "samtools view -h -b",$1,"chr6 >",var"/"$2".chr6.bam"}' fileN2.txt >batch1ExtractChr6.sh
awk -v var="$TEMPDIR" '{print "samtools view -h -b",$1,"6 >",var"/"$2".chr6.bam"}' fileN2.txt >batch1ExtractChr6.sh
rm fileN*.txt
sh batch1ExtractChr6.sh

### 2 sort bam by read name
### One can sort the BAM file by query name with samtools
### samtools sort -n aln.bam aln.qsort
ls $TEMPDIR/*.bam > fileN.txt
cp fileN.txt fileN2.txt
sed -i 's/\./ /' fileN2.txt
awk '{print $1}' fileN2.txt > fileN3.txt
paste fileN.txt fileN3.txt > fileN2.txt
awk '{print "samtools sort -n",$1,$2".chr6.qsort"}' fileN2.txt >batch2Sort.sh
rm fileN*.txt
sh batch2Sort.sh

### 3 convert bam to fastq
### bedtools -fq2 Creating two FASTQ files for paired-end sequences
### When using this option, it is required that the BAM file is sorted/grouped by the read name. 
### This keeps the resulting records in the two output FASTQ files in the same order.
ls $TEMPDIR/*qsort.bam > fileN.txt
cp fileN.txt fileN2.txt
sed -i 's/\// /g' fileN2.txt
sed -i 's/\./ /' fileN2.txt
awk '{print $(NF-1)}' fileN2.txt > fileN3.txt
paste fileN.txt fileN3.txt > fileN2.txt
awk -v var="$FQDIR" '{print "bedtools bamtofastq -i",$1,"-fq",var"/"$2".chr6.end1.fastq","-fq2",var"/"$2".chr6.end2.fastq"}' fileN2.txt >batch3BAM2Fastq.sh
rm fileN*.txt
sh batch3BAM2Fastq.sh

### 4 filtering reads using razers3

### deal with 10 individuals each time !!!!
### submit jobs to server
### each user can submit 10 jobs
### split the batch file into ten jobs
######################################
### require one node, max time: 100hrs
### qsub0.txt
### #!/bin/bash
### #PBS -l nodes=1
### #PBS -l walltime=100:00:00
### #PBS -m abe
### #PBS -q default
### #PBS -N changeThis
### cd $PBS_O_WORKDIR
######################################
###### get individual ids
ls *.end1.fastq > fileN.txt
sed -i 's/\./ /' fileN.txt
awk '{print $1}' fileN.txt > fileN2.txt
###### Read IDs to array "lines"
IFS=$'\n' read -d '' -r -a lines < fileN2.txt
rm fileN*.txt
###### generate batch files
for i in `seq 1 10`
do
	cp qsub0.txt run$i.txt
	sed -i "s/changeThis/${lines[$i]}/" run$i.txt
	echo "./razers3 -i 90 -m 1 -dr 0 -o ./${lines[$i]}.chr6.end1.fished.sam ./hla_reference_dna.fasta ./${lines[$i]}.chr6.end1.fastq" >> run$i.txt
	echo "./razers3 -i 90 -m 1 -dr 0 -o ./${lines[$i]}.chr6.end2.fished.sam ./hla_reference_dna.fasta ./${lines[$i]}.chr6.end2.fastq" >> run$i.txt
done
###### 

### after run this script, copy all fastq file and run{1,10}.txt to the cluster server to run razers3
### cd $FQDIR
### scp *.fastq yanhui@statgenpro.psychiatry.hku.hk:/home/yanhui/ngs/
### scp run*.txt yanhui@statgenpro.psychiatry.hku.hk:/home/yanhui/ngs/
### on cluster server:
### qsub run1.txt
### qsub run2.txt
### ...
##################################### END OF PART 1 ###################################

