################################################################################
####                  pipeline for ngs analysis: part 1                     ####
####                  1 FROM FASTQ TO ANALYSIS-READY BAM                    ####
####                  2 Identify potential variants in each sample          ####
####                  !!!! CAN NOT RUN WITH 'nohup'!!!!                     ####
####                  !!!! For PAIR-END FASTQ ONLY !!!!                     ####
####                  last modified on 8 Oct 2015                           ####
####                  @ Felix Yanhui Fan, felixfanyh@gmail.com              ####
################################################################################
#### usage: sh wgsPipeline1.sh 1 '@RG\tID:AIS7588\tSM:7588\tPL:ILLUMINA\tLB:lib7588\tPU:H52M2CCXX:1:none' 7588_R1.fastq.gz 7588_R2.fastq.gz 7588 True True True
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
REDUCE=$6
### data compression with ReduceReads, e.g. True
RMTEMP=$7
### remove generated file in each last step to save disk sapce, e.g. False
### keep fastq input & *dedup.bam, remove all others 
DPRM=$8
### remove generated file in each last step to save disk sapce, e.g. False
### remove fastq input & remove *dedup.bam
################################################################
GOLD="Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
KG="1000G_phase1.indels.hg19.sites.vcf"
DBSNP="dbsnp_138.hg19.vcf"
REF="ucsc.hg19.fasta"
### 1 Mapping the data to the reference—Generate a SAM file containing aligned reads by bwa
bwa mem -M -t $CORE -R $RG $REF $FQ1 $FQ2 > $OUT.align.sam
if [ $DPRM == "True" ]; then
    rm $FQ1
    rm $FQ2
fi
### 2 Converting to BAM
samtools view -bS $OUT.align.sam -o $OUT.align.bam
if [ $RMTEMP == "True" ]; then
    rm $OUT.align.sam
fi
### 3  sorting and marking duplicates
java -jar SortSam.jar INPUT=$OUT.align.bam OUTPUT=$OUT.sorted.bam SORT_ORDER=coordinate
java -jar MarkDuplicates.jar INPUT=$OUT.sorted.bam OUTPUT=$OUT.dedup.bam METRICS_FILE=$OUT.metrics.txt
java -jar BuildBamIndex.jar INPUT=$OUT.dedup.bam
if [ $RMTEMP == "True" ]; then
    rm $OUT.align.bam
    rm $OUT.sorted.bam
fi
### 4 Local realignment around indels
java -jar GenomeAnalysisTK.jar -nt $CORE -T RealignerTargetCreator -R $REF -I $OUT.dedup.bam -known $GOLD -known $KG -o $OUT.intervals
java -jar GenomeAnalysisTK.jar -T IndelRealigner -R $REF -targetIntervals $OUT.intervals -I $OUT.dedup.bam -known $GOLD -known $KG -o $OUT.realigned.bam
### 5 Base quality score recalibration (BQSR)
java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R $REF -I $OUT.realigned.bam -knownSites $DBSNP -knownSites $GOLD -knownSites $KG -o $OUT.recal.grp
java -jar GenomeAnalysisTK.jar -T PrintReads -R $REF -I $OUT.realigned.bam -BQSR $OUT.recal.grp -o $OUT.recal.bam
java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R $REF -I $OUT.realigned.bam -knownSites $DBSNP -knownSites $GOLD -knownSites $KG -BQSR $OUT.recal.grp -o $OUT.post.recal.grp
java -jar GenomeAnalysisTK.jar -T AnalyzeCovariates -R $REF -before $OUT.recal.grp -after $OUT.post.recal.grp -plots $OUT.recal.plots.pdf
if [ $RMTEMP == "True" ]; then
    rm $OUT.realigned.bam
fi
### 6 Identify	potential variants in each sample (HaplotypeCaller)
### output ends with .g.vcf
java -jar -Xmx32g GenomeAnalysisTK.jar -nct 8 -T HaplotypeCaller -R $REF -I $OUT.recal.bam -o $OUT.g.vcf -ERC GVCF
#########################END######################################################
