#!/usr/bin/env bash

### typing HLA using exome sequencing data
### run on Ubuntu 14.04 LTS
### @ Felix Yanhui Fan felixfanyh@gmail.com
### @ 17 Aug 2016

### From bam to HLA typing
SAMPLE='1094'
INPUT='1094.recal.bam'
TOOLDIR="/home/felixfan/ubin/HLA-VBSeq"
################################################################
### Build index for reference
#bwa index ${TOOLDIR}/hla_all.fasta
################################################################
## EXTRACT chr6
samtools view -b ${INPUT} chr6 > ${SAMPLE}.chr6.bam
################################################################
##INDEX
samtools index ${SAMPLE}.chr6.bam
################################################################
### Extract a list of read name that were aligned to HLA loci (HLA-A, B, C, DM, DO, DP, DQ, DR, E, F, G, H, J, K, L, P, V, MIC, and TAP)
samtools view ${SAMPLE}.chr6.bam chr6:29907037-29915661 chr6:31319649-31326989 chr6:31234526-31241863 chr6:32914391-32922899 chr6:32900406-32910847 chr6:32969960-32979389 chr6:32778540-32786825 chr6:33030346-33050555 chr6:33041703-33059473 chr6:32603183-32613429 chr6:32707163-32716664 chr6:32625241-32636466 chr6:32721875-32733330 chr6:32405619-32414826 chr6:32544547-32559613 chr6:32518778-32554154 chr6:32483154-32559613 chr6:30455183-30463982 chr6:29689117-29699106 chr6:29792756-29800899 chr6:29793613-29978954 chr6:29855105-29979733 chr6:29892236-29899009 chr6:30225339-30236728 chr6:31369356-31385092 chr6:31460658-31480901 chr6:29766192-29772202 chr6:32810986-32823755 chr6:32779544-32808599 chr6:29756731-29767588 | gawk '{print $1}' | sort | uniq > ${SAMPLE}.partial.reads.txt
################################################################
### Build read name index and search read pairs and their sequences on HLA loci
java -Xmx30g -jar $TOOLDIR/bamNameIndex.jar index ${SAMPLE}.chr6.bam --indexFile ${SAMPLE}.chr6.bam.idx
java -Xmx30g -jar $TOOLDIR/bamNameIndex.jar search ${SAMPLE}.chr6.bam --name ${SAMPLE}.partial.reads.txt --output ${SAMPLE}.partial.sam
java -Xmx30g -jar $TOOLDIR/SamToFastq.jar I=${SAMPLE}.partial.sam F=${SAMPLE}.partial.end1.fastq F2=${SAMPLE}.partial.end2.fastq
################################################################
### Extract unmapped reads
samtools view -bh -f 12 ${SAMPLE}.chr6.bam > ${SAMPLE}.unmapped.bam
java -Xmx30g -jar $TOOLDIR/SamToFastq.jar I=${SAMPLE}.unmapped.bam F=${SAMPLE}.unmapped.end1.fastq F2=${SAMPLE}.unmapped.end2.fastq
################################################################
### Combine reads in FASTQ format
cat ${SAMPLE}.partial.end1.fastq ${SAMPLE}.unmapped.end1.fastq > ${SAMPLE}.end1.fastq
cat ${SAMPLE}.partial.end2.fastq ${SAMPLE}.unmapped.end2.fastq > ${SAMPLE}.end2.fastq
################################################################
### Alignment by BWA-MEM allowing multiple alignments for each read
bwa mem -t 8 -P -L 10000 -a $TOOLDIR/hla_all.fasta ${SAMPLE}.end1.fastq ${SAMPLE}.end2.fastq > ${SAMPLE}.hla.sam
################################################################
### Estimation of HLA types by HLA-VBSeq 
java -jar -Xmx30g $TOOLDIR/HLAVBSeq.jar $TOOLDIR/hla_all.fasta ${SAMPLE}.hla.sam ${SAMPLE}.result.txt --alpha_zero 0.01 --is_paired
################################################################
### parse result
perl $TOOLDIR/parse_result.pl $TOOLDIR/Allelelist.txt ${SAMPLE}.result.txt | grep "^A\*" | sort -k2 -n -r > ${SAMPLE}.HLA_A.txt
perl $TOOLDIR/parse_result.pl $TOOLDIR/Allelelist.txt ${SAMPLE}.result.txt | grep "^B\*" | sort -k2 -n -r > ${SAMPLE}.HLA_B.txt
perl $TOOLDIR/parse_result.pl $TOOLDIR/Allelelist.txt ${SAMPLE}.result.txt | grep "^C\*" | sort -k2 -n -r > ${SAMPLE}.HLA_C.txt
perl $TOOLDIR/parse_result.pl $TOOLDIR/Allelelist.txt ${SAMPLE}.result.txt | grep "^DRB1\*" | sort -k2 -n -r > ${SAMPLE}.HLA_DRB1.txt
perl $TOOLDIR/parse_result.pl $TOOLDIR/Allelelist.txt ${SAMPLE}.result.txt | grep "^DQB1\*" | sort -k2 -n -r > ${SAMPLE}.HLA_DQB1.txt
perl $TOOLDIR/parse_result.pl $TOOLDIR/Allelelist.txt ${SAMPLE}.result.txt | grep "^DQA1\*" | sort -k2 -n -r > ${SAMPLE}.HLA_DQA1.txt

