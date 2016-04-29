################################################################################
####                  pipeline for ngs analysis: part 1                     ####
####                  1 FROM FASTQ TO ANALYSIS-READY BAM                    ####
####                  2 Identify potential variants in each sample          ####
####                  !!!! CAN NOT RUN WITH 'nohup'!!!!                     ####
####                  !!!! For PAIR-END FASTQ ONLY !!!!                     ####
####                  last modified on 8 Oct 2015                           ####
####                  @ Felix Yanhui Fan, felixfanyh@gmail.com              ####
################################################################################
#### usage: sh wgsPipeline1.sh 1 '@RG\tID:AIS7588\tSM:7588\tPL:ILLUMINA\tLB:lib7588\tPU:H52M2CCXX:1:none' 7588_R1.fastq.gz 7588_R2.fastq.gz 7588
CORE=$1
### e.g. 8
RG=$2
### e.g. '@RG\tID:AIS7588\tSM:7588\tPL:ILLUMINA\tLB:lib7588\tPU:H52M2CCXX:1:none'
### e.g. '@RG\tID:AIS7590\tSM:7590\tPL:ILLUMINA\tLB:lib7590\tPU:H52M2CCXX:2:none'
FQ1=$3
### e.g. 7588_R1.fastq.gz
FQ2=$4
### e.g. 7588_R2.fastq.gz
OUT=$5
### prefix of output e.g. 7588
################################################################
GOLD="/home/felixfan/ulib/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
KG="/home/felixfan/ulib/hg19/1000G_phase1.indels.hg19.sites.vcf"
DBSNP="/home/felixfan/ulib/hg19/dbsnp_138.hg19.vcf"
REF="/home/felixfan/ulib/hg19/ucsc.hg19.fasta"
### 1 Mapping the data to the reference—Generate a SAM file containing aligned reads by bwa
bwa mem -M -t $CORE -R $RG $REF $FQ1 $FQ2 > $OUT.align.sam
rm $FQ1
rm $FQ2
### 2  sorting and marking duplicates
java -jar SortSam.jar INPUT=$OUT.align.sam OUTPUT=$OUT.sorted.bam SORT_ORDER=coordinate
rm $OUT.align.sam
java -jar MarkDuplicates.jar INPUT=$OUT.sorted.bam OUTPUT=$OUT.dedup.bam METRICS_FILE=$OUT.metrics.txt
rm $OUT.sorted.bam
java -jar BuildBamIndex.jar INPUT=$OUT.dedup.bam
### 3 Local realignment around indels
java -jar GenomeAnalysisTK.jar -nt $CORE -T RealignerTargetCreator -R $REF -I $OUT.dedup.bam -known $GOLD -known $KG -o $OUT.intervals
java -jar GenomeAnalysisTK.jar -T IndelRealigner -R $REF -targetIntervals $OUT.intervals -I $OUT.dedup.bam -known $GOLD -known $KG -o $OUT.realigned.bam
rm $OUT.dedup.bam
### 4 Base quality score recalibration (BQSR)
java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R $REF -I $OUT.realigned.bam -knownSites $DBSNP -knownSites $GOLD -knownSites $KG -o $OUT.recal.grp
java -jar GenomeAnalysisTK.jar -T PrintReads -R $REF -I $OUT.realigned.bam -BQSR $OUT.recal.grp -o $OUT.recal.bam
rm $OUT.realigned.bam
### 5 Identify	potential variants in each sample (HaplotypeCaller)
### output ends with .g.vcf
java -jar -Xmx20g GenomeAnalysisTK.jar -nct $CORE -T HaplotypeCaller -R $REF -I $OUT.recal.bam -o $OUT.g.vcf -ERC GVCF
#########################END######################################################
