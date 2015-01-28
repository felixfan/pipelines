#!/usr/bin/env bash

### typing HLA using exome sequencing data
### run on Ubuntu 14.04 LTS
### @ Felix Yanhui Fan felixfanyh@gmail.com
### @ 27 Jan 2015

### !!! before use it, check and modify all parts with "!!!"
#################################################################
BAMDIR="/media/felixfan/ad/Mapping/BWA_GATK/bam"
#BAMDIR="/home/felixfan/uwork/ADExomeHLA/results/OptiType/tempdir/"
TOOLDIR="/home/felixfan/ubin/HLA-VBSeq"
OUTDIR="/home/felixfan/uwork/ADExomeHLA/results/HLA-VBSeq/result"
#################################################################
if [ ! -d "$OUTDIR" ]
then
	mkdir $OUTDIR
fi
################################!!!################################
### Get individual id (first part of the file name) and number of individuals
#ls $BAMDIR/*.chr6.bam > fileN.txt
ls $BAMDIR/*.bam > fileN.txt
sed -i 's/\// /g' fileN.txt
sed -i 's/\.bam$/ /g' fileN.txt
#awk '{print $8}' fileN.txt > fileN2.txt
awk '{print $7}' fileN.txt > fileN2.txt
IFS=$'\n' read -d '' -r -a ids < fileN2.txt     # individual ids
# n=$(ls $BAMDIR/*.chr6.bam | wc -l)                   # number of individuals
n=$(ls $BAMDIR/*.bam | wc -l)                   # number of individuals
let m=$n-1
###################################################################
### Build index for reference
bwa index hla_all.fasta
######################################################################

for i in $(seq 0 $m)
do
	ID=${ids[$i]}

	### Build index for bam
	cp $BAMDIR/$ID.bam .
	samtools index $ID.bam

	### Extract a list of read name that were aligned to HLA loci (HLA-A, B, C, DM, DO, DP, DQ, DR, E, F, G, H, J, K, L, P, V, MIC, and TAP)
	samtools view $ID.bam 6:29907037-29915661 6:31319649-31326989 6:31234526-31241863 6:32914391-32922899 6:32900406-32910847 6:32969960-32979389 6:32778540-32786825 6:33030346-33050555 6:33041703-33059473 6:32603183-32613429 6:32707163-32716664 6:32625241-32636466 6:32721875-32733330 6:32405619-32414826 6:32544547-32559613 6:32518778-32554154 6:32483154-32559613 6:30455183-30463982 6:29689117-29699106 6:29792756-29800899 6:29793613-29978954 6:29855105-29979733 6:29892236-29899009 6:30225339-30236728 6:31369356-31385092 6:31460658-31480901 6:29766192-29772202 6:32810986-32823755 6:32779544-32808599 6:29756731-29767588 | gawk '{print $1}' | sort | uniq > $ID.partial.reads.txt

	### Build read name index and search read pairs and their sequences on HLA loci
	java -jar $TOOLDIR/bamNameIndex.jar index $ID.bam --indexFile $ID.bam.idx
	java -jar $TOOLDIR/bamNameIndex.jar search $ID.bam --name $ID.partial.reads.txt --output $ID.partial.sam
	samtools view -Sb $ID.partial.sam > $ID.temp.bam
	bedtools bamtofastq -i $ID.temp.bam -fq $ID.partial.end1.fastq -fq2 $ID.partial.end2.fastq

	### Extract unmapped reads
	samtools view -bh -f 12 $ID.bam > $ID.unmapped.bam
	bedtools bamtofastq -i $ID.unmapped.bam -fq $ID.unmapped.end1.fastq -fq2 $ID.unmapped.end2.fastq

	### Combine reads in FASTQ format
	cat $ID.partial.end1.fastq $ID.unmapped.end1.fastq > $ID.end1.fastq
	cat $ID.partial.end2.fastq $ID.unmapped.end2.fastq > $ID.end2.fastq

	### Alignment by BWA-MEM allowing multiple alignments for each read
	bwa mem -t 8 -P -L 10000 -a hla_all.fasta $ID.end1.fastq $ID.end2.fastq > $ID.sam

	### Estimation of HLA types by HLA-VBSeq 
	java -jar -Xmx32g -Xms32g $TOOLDIR/HLAVBSeq.jar hla_all.fasta $ID.sam $OUTDIR/$ID.result.txt --alpha_zero 0.01 --is_paired

	### parse result
	parse_result.pl Allelelist.txt $OUTDIR/$ID.result.txt | grep "^A\*" | sort -k2 -n -r > $OUTDIR/$ID.HLA_A.txt
	parse_result.pl Allelelist.txt $OUTDIR/$ID.result.txt | grep "^B\*" | sort -k2 -n -r > $OUTDIR/$ID.HLA_B.txt
	parse_result.pl Allelelist.txt $OUTDIR/$ID.result.txt | grep "^C\*" | sort -k2 -n -r > $OUTDIR/$ID.HLA_C.txt
	parse_result.pl Allelelist.txt $OUTDIR/$ID.result.txt | grep "^DRB1\*" | sort -k2 -n -r > $OUTDIR/$ID.HLA_DRB1.txt
	parse_result.pl Allelelist.txt $OUTDIR/$ID.result.txt | grep "^DQB1\*" | sort -k2 -n -r > $OUTDIR/$ID.HLA_DQB1.txt
	parse_result.pl Allelelist.txt $OUTDIR/$ID.result.txt | grep "^DQA1\*" | sort -k2 -n -r > $OUTDIR/$ID.HLA_DQA1.txt

	### clear
	rm ./$ID.*

done

rm fileN*.txt

