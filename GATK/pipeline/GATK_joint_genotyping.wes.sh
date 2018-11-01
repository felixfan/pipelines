#!/bin/bash

#### pipeline for WGS variant calling - part 2
#### use the best practices for variant calling with GATK
#### software used: GATK
#### How to use:
####          param 1: max memory for java in Gb, e.g. 32
####          param 2: number of threads will be used, e.g. 8
####          param 3: reference genome, e.g. ucsc.hg19.fasta
####          param 4: dbSNP, e.g. dbsnp_138.hg19.vcf
####          param 5: file with gVCF file names, one g.vcf per line
####          param 6: output file name. e.g. out.raw.vcf
#### Example: sh GATK_joint_genotyping.wes.sh 32 8 ucsc.hg19.fasta dbsnp_138.hg19.vcf gvcf.txt out.raw.vcf
#### Note: according to the best practice, if you have more than 200 gVCF files,
####       you need to combine in batches first using CombineGVCFs

echo "@-------------------------------------------------------------@"
echo "|    GATK_joint_genotyping.wes   |   v2.0.0  |  07 Sep 2017   |"
echo "|-------------------------------------------------------------|"
echo "|  (C) 2017 Felix Yanhui Fan, GNU General Public License, v2  |"
echo "|-------------------------------------------------------------|"
echo "|    For documentation, citation & bug-report instructions:   |"
echo "|            http://felixfan.github.io/pipelines              |"
echo "@-------------------------------------------------------------@"

MEM=$1   # MAX MEMORY IN G
NT=$2    # how many core to be used
REF=$3   # reference genome
DBSNP=$4 # dbSNP
GVCFS=$5 # file of gvcfs, one gvcf per line
OUT=$6   # output name

comms="java -Xmx${MEM}g -jar GenomeAnalysisTK.jar -nt ${NT} -T GenotypeGVCFs -R ${REF} --dbsnp ${DBSNP} -o ${OUT}"

while read -r line
do
    comms="$comms -V $line"
done < "$GVCFS"

echo $comms

$comms
