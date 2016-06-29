#!/bin/bash

#### pipeline for WGS variant calling - part 2
#### use the best practices for variant calling with GATK
#### software used: GATK
#### How to use:
####          param 1: number of threads will be used, e.g. 8
####          param 2: reference genome, e.g. ucsc.hg19.fasta
####          param 3: prefix of output
####          param 4: the first gVCF file name
####          param n(n>4): the second to the last gVCF file name
#### Example: sh GATK_joint_genotyping.sh 8 ucsc.hg19.fasta pht 1.g.vcf 2.g.vcf 3.g.vcf
#### Note: according to the best practice, if you have more than 200 gVCF files, 
####       you need to combine in batches first using CombineGVCFs

echo "@-------------------------------------------------------------@"
echo "|     GATK_joint_genotyping   |   v1.0.0    |   29 Jun 2016   |"
echo "|-------------------------------------------------------------|"
echo "|  (C) 2016 Felix Yanhui Fan, GNU General Public License, v2  |"
echo "|-------------------------------------------------------------|"
echo "|    For documentation, citation & bug-report instructions:   |"
echo "|            http://felixfan.github.io/pipelines              |"
echo "@-------------------------------------------------------------@"

comms="java -Xmx32g -jar GenomeAnalysisTK.jar -T GenotypeGVCFs"

n=0
for gvcf in $@
do
    n=$(($n+1))
    if [ $n -eq 1 ]; then
        comms="$comms -nt $gvcf"
    elif [ $n -eq 2 ]; then
        comms="$comms -R $gvcf"
    elif [ $n -eq 3 ]; then
        comms="$comms -o ${gvcf}.raw.vcf"
    else
        comms="$comms -V $gvcf"
    fi
done

$comms
