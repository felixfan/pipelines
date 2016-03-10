### 1 download reference data
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/dbsnp_138.hg19.vcf.gz
gunzip ucsc.hg19.fasta.gz
### 2 Generate the BWA index
bwa index ucsc.hg19.fasta
### about 50 min
### 3 Generate fasta file index by samtools
samtools faidx ucsc.hg19.fasta
### about 1 min
### 4 Generate the sequence dictionary by Picard
java -jar CreateSequenceDictionary.jar REFERENCE=ucsc.hg19.fasta OUTPUT=ucsc.hg19.dict
### about 1 min
### 5 Mapping the data to the reference—Generate a SAM file containing aligned reads by bwa
# -M Mark shorter split hits as secondary (for Picard compatibility).
# -t Number of threads [1]
# -R Complete read group header line. ’\t’ can be used in STR and will be converted to a TAB in the output SAM.
# a="bwa mem -M -t 8 -R "
# c="'@RG\tID:AIS"
# d="\tSM:"
# e="\tPL:ILLUMINA\tLB:lib"
# f="\tPU:H52M2CCXX:"
# g=":none' ucsc.hg19.fasta"
# n=1
# for i in 7588 7590 7591 7592 7701 7702 7703 7704B 7706 7708 7713
# do
#     echo $a$c$i$d$i$e$i$f$n$g $i"_R1.fastq.gz" $i"_R2.fastq.gz >" $i"_aligned_reads.sam"
#     n=$(($n+1))
# done
bwa mem -M -t 8 -R '@RG\tID:AIS7588\tSM:7588\tPL:ILLUMINA\tLB:lib7588\tPU:H52M2CCXX:1:none' ucsc.hg19.fasta 7588_R1.fastq.gz 7588_R2.fastq.gz > 7588_aligned_reads.sam
bwa mem -M -t 8 -R '@RG\tID:AIS7590\tSM:7590\tPL:ILLUMINA\tLB:lib7590\tPU:H52M2CCXX:2:none' ucsc.hg19.fasta 7590_R1.fastq.gz 7590_R2.fastq.gz > 7590_aligned_reads.sam
bwa mem -M -t 8 -R '@RG\tID:AIS7591\tSM:7591\tPL:ILLUMINA\tLB:lib7591\tPU:H52M2CCXX:3:none' ucsc.hg19.fasta 7591_R1.fastq.gz 7591_R2.fastq.gz > 7591_aligned_reads.sam
bwa mem -M -t 8 -R '@RG\tID:AIS7592\tSM:7592\tPL:ILLUMINA\tLB:lib7592\tPU:H52M2CCXX:4:none' ucsc.hg19.fasta 7592_R1.fastq.gz 7592_R2.fastq.gz > 7592_aligned_reads.sam
bwa mem -M -t 8 -R '@RG\tID:AIS7701\tSM:7701\tPL:ILLUMINA\tLB:lib7701\tPU:H52M2CCXX:5:none' ucsc.hg19.fasta 7701_R1.fastq.gz 7701_R2.fastq.gz > 7701_aligned_reads.sam
bwa mem -M -t 8 -R '@RG\tID:AIS7702\tSM:7702\tPL:ILLUMINA\tLB:lib7702\tPU:H52M2CCXX:6:none' ucsc.hg19.fasta 7702_R1.fastq.gz 7702_R2.fastq.gz > 7702_aligned_reads.sam
bwa mem -M -t 8 -R '@RG\tID:AIS7703\tSM:7703\tPL:ILLUMINA\tLB:lib7703\tPU:H52M2CCXX:7:none' ucsc.hg19.fasta 7703_R1.fastq.gz 7703_R2.fastq.gz > 7703_aligned_reads.sam
bwa mem -M -t 8 -R '@RG\tID:AIS7704B\tSM:7704B\tPL:ILLUMINA\tLB:lib7704B\tPU:H52M2CCXX:8:none' ucsc.hg19.fasta 7704B_R1.fastq.gz 7704B_R2.fastq.gz > 7704B_aligned_reads.sam
bwa mem -M -t 8 -R '@RG\tID:AIS7706\tSM:7706\tPL:ILLUMINA\tLB:lib7706\tPU:H52M2CCXX:9:none' ucsc.hg19.fasta 7706_R1.fastq.gz 7706_R2.fastq.gz > 7706_aligned_reads.sam
bwa mem -M -t 8 -R '@RG\tID:AIS7708\tSM:7708\tPL:ILLUMINA\tLB:lib7708\tPU:H52M2CCXX:10:none' ucsc.hg19.fasta 7708_R1.fastq.gz 7708_R2.fastq.gz > 7708_aligned_reads.sam
bwa mem -M -t 8 -R '@RG\tID:AIS7713\tSM:7713\tPL:ILLUMINA\tLB:lib7713\tPU:H52M2CCXX:11:none' ucsc.hg19.fasta 7713_R1.fastq.gz 7713_R2.fastq.gz > 7713_aligned_reads.sam
### 6 Converting to BAM
# for i in 7588 7590 7591 7592 7701 7702 7703 7704B 7706 7708 7713
# do
# 	echo "samtools view -bS "$i"_aligned_reads.sam -o "$i"_aligned_reads.bam"
# done
samtools view -bS 7588_aligned_reads.sam -o 7588_aligned_reads.bam
samtools view -bS 7590_aligned_reads.sam -o 7590_aligned_reads.bam
samtools view -bS 7591_aligned_reads.sam -o 7591_aligned_reads.bam
samtools view -bS 7592_aligned_reads.sam -o 7592_aligned_reads.bam
samtools view -bS 7701_aligned_reads.sam -o 7701_aligned_reads.bam
samtools view -bS 7702_aligned_reads.sam -o 7702_aligned_reads.bam
samtools view -bS 7703_aligned_reads.sam -o 7703_aligned_reads.bam
samtools view -bS 7704B_aligned_reads.sam -o 7704B_aligned_reads.bam
samtools view -bS 7706_aligned_reads.sam -o 7706_aligned_reads.bam
samtools view -bS 7708_aligned_reads.sam -o 7708_aligned_reads.bam
samtools view -bS 7713_aligned_reads.sam -o 7713_aligned_reads.bam
### 7  sorting and marking duplicates
# for i in 7588 7590 7591 7592 7701 7702 7703 7704B 7706 7708 7713
# do
#     echo "java -jar SortSam.jar INPUT="$i"_aligned_reads.bam OUTPUT="$i"_sorted_reads.bam SORT_ORDER=coordinate"
#     echo "java -jar MarkDuplicates.jar INPUT="$i"_sorted_reads.bam OUTPUT="$i"_dedup_reads.bam METRICS_FILE="$i"_metrics.txt"
#     echo "java -jar BuildBamIndex.jar INPUT="$i"_dedup_reads.bam" 
# done
java -jar SortSam.jar INPUT=7588_aligned_reads.bam OUTPUT=7588_sorted_reads.bam SORT_ORDER=coordinate
java -jar MarkDuplicates.jar INPUT=7588_sorted_reads.bam OUTPUT=7588_dedup_reads.bam METRICS_FILE=7588_metrics.txt
java -jar BuildBamIndex.jar INPUT=7588_dedup_reads.bam
java -jar SortSam.jar INPUT=7590_aligned_reads.bam OUTPUT=7590_sorted_reads.bam SORT_ORDER=coordinate
java -jar MarkDuplicates.jar INPUT=7590_sorted_reads.bam OUTPUT=7590_dedup_reads.bam METRICS_FILE=7590_metrics.txt
java -jar BuildBamIndex.jar INPUT=7590_dedup_reads.bam
java -jar SortSam.jar INPUT=7591_aligned_reads.bam OUTPUT=7591_sorted_reads.bam SORT_ORDER=coordinate
java -jar MarkDuplicates.jar INPUT=7591_sorted_reads.bam OUTPUT=7591_dedup_reads.bam METRICS_FILE=7591_metrics.txt
java -jar BuildBamIndex.jar INPUT=7591_dedup_reads.bam
java -jar SortSam.jar INPUT=7592_aligned_reads.bam OUTPUT=7592_sorted_reads.bam SORT_ORDER=coordinate
java -jar MarkDuplicates.jar INPUT=7592_sorted_reads.bam OUTPUT=7592_dedup_reads.bam METRICS_FILE=7592_metrics.txt
java -jar BuildBamIndex.jar INPUT=7592_dedup_reads.bam
java -jar SortSam.jar INPUT=7701_aligned_reads.bam OUTPUT=7701_sorted_reads.bam SORT_ORDER=coordinate
java -jar MarkDuplicates.jar INPUT=7701_sorted_reads.bam OUTPUT=7701_dedup_reads.bam METRICS_FILE=7701_metrics.txt
java -jar BuildBamIndex.jar INPUT=7701_dedup_reads.bam
java -jar SortSam.jar INPUT=7702_aligned_reads.bam OUTPUT=7702_sorted_reads.bam SORT_ORDER=coordinate
java -jar MarkDuplicates.jar INPUT=7702_sorted_reads.bam OUTPUT=7702_dedup_reads.bam METRICS_FILE=7702_metrics.txt
java -jar BuildBamIndex.jar INPUT=7702_dedup_reads.bam
java -jar SortSam.jar INPUT=7703_aligned_reads.bam OUTPUT=7703_sorted_reads.bam SORT_ORDER=coordinate
java -jar MarkDuplicates.jar INPUT=7703_sorted_reads.bam OUTPUT=7703_dedup_reads.bam METRICS_FILE=7703_metrics.txt
java -jar BuildBamIndex.jar INPUT=7703_dedup_reads.bam
java -jar SortSam.jar INPUT=7704B_aligned_reads.bam OUTPUT=7704B_sorted_reads.bam SORT_ORDER=coordinate
java -jar MarkDuplicates.jar INPUT=7704B_sorted_reads.bam OUTPUT=7704B_dedup_reads.bam METRICS_FILE=7704B_metrics.txt
java -jar BuildBamIndex.jar INPUT=7704B_dedup_reads.bam
java -jar SortSam.jar INPUT=7706_aligned_reads.bam OUTPUT=7706_sorted_reads.bam SORT_ORDER=coordinate
java -jar MarkDuplicates.jar INPUT=7706_sorted_reads.bam OUTPUT=7706_dedup_reads.bam METRICS_FILE=7706_metrics.txt
java -jar BuildBamIndex.jar INPUT=7706_dedup_reads.bam
java -jar SortSam.jar INPUT=7708_aligned_reads.bam OUTPUT=7708_sorted_reads.bam SORT_ORDER=coordinate
java -jar MarkDuplicates.jar INPUT=7708_sorted_reads.bam OUTPUT=7708_dedup_reads.bam METRICS_FILE=7708_metrics.txt
java -jar BuildBamIndex.jar INPUT=7708_dedup_reads.bam
java -jar SortSam.jar INPUT=7713_aligned_reads.bam OUTPUT=7713_sorted_reads.bam SORT_ORDER=coordinate
java -jar MarkDuplicates.jar INPUT=7713_sorted_reads.bam OUTPUT=7713_dedup_reads.bam METRICS_FILE=7713_metrics.txt
java -jar BuildBamIndex.jar INPUT=7713_dedup_reads.bam
### 8 Local realignment around indels
for i in 7588 7590 7591 7592 7701 7702 7703 7704B 7706 7708 7713
do
    echo "java -jar GenomeAnalysisTK.jar -nt 8 -T RealignerTargetCreator -R ucsc.hg19.fasta -I "$i"_dedup_reads.bam -known Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known 1000G_phase1.indels.hg19.sites.vcf -o "$i"_target_intervals.list"
    echo "java -jar GenomeAnalysisTK.jar -T IndelRealigner -R ucsc.hg19.fasta -I "$i"_dedup_reads.bam -known Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known 1000G_phase1.indels.hg19.sites.vcf -o "$i"_realigned_reads.bam -targetIntervals "$i"_target_intervals.list"  
done
java -jar GenomeAnalysisTK.jar -nt 8 -T RealignerTargetCreator -R ucsc.hg19.fasta -I 7588_dedup_reads.bam -known Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known 1000G_phase1.indels.hg19.sites.vcf -o 7588_target_intervals.list
java -jar GenomeAnalysisTK.jar -T IndelRealigner -R ucsc.hg19.fasta -I 7588_dedup_reads.bam -known Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known 1000G_phase1.indels.hg19.sites.vcf -o 7588_realigned_reads.bam -targetIntervals 7588_target_intervals.list
java -jar GenomeAnalysisTK.jar -nt 8 -T RealignerTargetCreator -R ucsc.hg19.fasta -I 7590_dedup_reads.bam -known Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known 1000G_phase1.indels.hg19.sites.vcf -o 7590_target_intervals.list
java -jar GenomeAnalysisTK.jar -T IndelRealigner -R ucsc.hg19.fasta -I 7590_dedup_reads.bam -known Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known 1000G_phase1.indels.hg19.sites.vcf -o 7590_realigned_reads.bam -targetIntervals 7590_target_intervals.list
java -jar GenomeAnalysisTK.jar -nt 8 -T RealignerTargetCreator -R ucsc.hg19.fasta -I 7591_dedup_reads.bam -known Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known 1000G_phase1.indels.hg19.sites.vcf -o 7591_target_intervals.list
java -jar GenomeAnalysisTK.jar -T IndelRealigner -R ucsc.hg19.fasta -I 7591_dedup_reads.bam -known Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known 1000G_phase1.indels.hg19.sites.vcf -o 7591_realigned_reads.bam -targetIntervals 7591_target_intervals.list
java -jar GenomeAnalysisTK.jar -nt 8 -T RealignerTargetCreator -R ucsc.hg19.fasta -I 7592_dedup_reads.bam -known Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known 1000G_phase1.indels.hg19.sites.vcf -o 7592_target_intervals.list
java -jar GenomeAnalysisTK.jar -T IndelRealigner -R ucsc.hg19.fasta -I 7592_dedup_reads.bam -known Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known 1000G_phase1.indels.hg19.sites.vcf -o 7592_realigned_reads.bam -targetIntervals 7592_target_intervals.list
java -jar GenomeAnalysisTK.jar -nt 8 -T RealignerTargetCreator -R ucsc.hg19.fasta -I 7701_dedup_reads.bam -known Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known 1000G_phase1.indels.hg19.sites.vcf -o 7701_target_intervals.list
java -jar GenomeAnalysisTK.jar -T IndelRealigner -R ucsc.hg19.fasta -I 7701_dedup_reads.bam -known Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known 1000G_phase1.indels.hg19.sites.vcf -o 7701_realigned_reads.bam -targetIntervals 7701_target_intervals.list
java -jar GenomeAnalysisTK.jar -nt 8 -T RealignerTargetCreator -R ucsc.hg19.fasta -I 7702_dedup_reads.bam -known Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known 1000G_phase1.indels.hg19.sites.vcf -o 7702_target_intervals.list
java -jar GenomeAnalysisTK.jar -T IndelRealigner -R ucsc.hg19.fasta -I 7702_dedup_reads.bam -known Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known 1000G_phase1.indels.hg19.sites.vcf -o 7702_realigned_reads.bam -targetIntervals 7702_target_intervals.list
java -jar GenomeAnalysisTK.jar -nt 8 -T RealignerTargetCreator -R ucsc.hg19.fasta -I 7703_dedup_reads.bam -known Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known 1000G_phase1.indels.hg19.sites.vcf -o 7703_target_intervals.list
java -jar GenomeAnalysisTK.jar -T IndelRealigner -R ucsc.hg19.fasta -I 7703_dedup_reads.bam -known Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known 1000G_phase1.indels.hg19.sites.vcf -o 7703_realigned_reads.bam -targetIntervals 7703_target_intervals.list
java -jar GenomeAnalysisTK.jar -nt 8 -T RealignerTargetCreator -R ucsc.hg19.fasta -I 7704B_dedup_reads.bam -known Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known 1000G_phase1.indels.hg19.sites.vcf -o 7704B_target_intervals.list
java -jar GenomeAnalysisTK.jar -T IndelRealigner -R ucsc.hg19.fasta -I 7704B_dedup_reads.bam -known Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known 1000G_phase1.indels.hg19.sites.vcf -o 7704B_realigned_reads.bam -targetIntervals 7704B_target_intervals.list
java -jar GenomeAnalysisTK.jar -nt 8 -T RealignerTargetCreator -R ucsc.hg19.fasta -I 7706_dedup_reads.bam -known Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known 1000G_phase1.indels.hg19.sites.vcf -o 7706_target_intervals.list
java -jar GenomeAnalysisTK.jar -T IndelRealigner -R ucsc.hg19.fasta -I 7706_dedup_reads.bam -known Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known 1000G_phase1.indels.hg19.sites.vcf -o 7706_realigned_reads.bam -targetIntervals 7706_target_intervals.list
java -jar GenomeAnalysisTK.jar -nt 8 -T RealignerTargetCreator -R ucsc.hg19.fasta -I 7708_dedup_reads.bam -known Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known 1000G_phase1.indels.hg19.sites.vcf -o 7708_target_intervals.list
java -jar GenomeAnalysisTK.jar -T IndelRealigner -R ucsc.hg19.fasta -I 7708_dedup_reads.bam -known Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known 1000G_phase1.indels.hg19.sites.vcf -o 7708_realigned_reads.bam -targetIntervals 7708_target_intervals.list
java -jar GenomeAnalysisTK.jar -nt 8 -T RealignerTargetCreator -R ucsc.hg19.fasta -I 7713_dedup_reads.bam -known Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known 1000G_phase1.indels.hg19.sites.vcf -o 7713_target_intervals.list
java -jar GenomeAnalysisTK.jar -T IndelRealigner -R ucsc.hg19.fasta -I 7713_dedup_reads.bam -known Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known 1000G_phase1.indels.hg19.sites.vcf -o 7713_realigned_reads.bam -targetIntervals 7713_target_intervals.list
### 9 Base quality score recalibration (BQSR)
for i in 7588 7590 7591 7592 7701 7702 7703 7704B 7706 7708 7713
do
	echo "java -jar GenomeAnalysisTK.jar -nct 8 -T BaseRecalibrator -R ucsc.hg19.fasta -I "$i"_realigned_reads.bam -knownSites dbsnp_138.hg19.vcf -knownSites Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites 1000G_phase1.indels.hg19.sites.vcf -o "$i"_recal_data.grp"
	echo "java -jar GenomeAnalysisTK.jar -nct 8 -T PrintReads -R ucsc.hg19.fasta -I "$i"_realigned_reads.bam -BQSR "$i"_recal_data.grp -o "$i"_recal_reads.bam"
	echo "java -jar GenomeAnalysisTK.jar -nct 8 -T BaseRecalibrator -R ucsc.hg19.fasta -I "$i"_realigned_reads.bam -knownSites dbsnp_138.hg19.vcf -knownSites Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites 1000G_phase1.indels.hg19.sites.vcf -BQSR "$i"_recal_data.grp -o "$i"_post_recal_data.grp"
	echo "java -jar GenomeAnalysisTK.jar -T AnalyzeCovariates -R ucsc.hg19.fasta -before "$i"_recal_data.grp -after "$i"_post_recal_data.grp -plots "$i"_recal_plots.pdf"
done
java -jar GenomeAnalysisTK.jar -nct 8 -T BaseRecalibrator -R ucsc.hg19.fasta -I 7588_realigned_reads.bam -knownSites dbsnp_138.hg19.vcf -knownSites Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites 1000G_phase1.indels.hg19.sites.vcf -o 7588_recal_data.grp
java -jar GenomeAnalysisTK.jar -nct 8 -T PrintReads -R ucsc.hg19.fasta -I 7588_realigned_reads.bam -BQSR 7588_recal_data.grp -o 7588_recal_reads.bam
java -jar GenomeAnalysisTK.jar -nct 8 -T BaseRecalibrator -R ucsc.hg19.fasta -I 7588_realigned_reads.bam -knownSites dbsnp_138.hg19.vcf -knownSites Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites 1000G_phase1.indels.hg19.sites.vcf -BQSR 7588_recal_data.grp -o 7588_post_recal_data.grp
java -jar GenomeAnalysisTK.jar -T AnalyzeCovariates -R ucsc.hg19.fasta -before 7588_recal_data.grp -after 7588_post_recal_data.grp -plots 7588_recal_plots.pdf
java -jar GenomeAnalysisTK.jar -nct 8 -T BaseRecalibrator -R ucsc.hg19.fasta -I 7590_realigned_reads.bam -knownSites dbsnp_138.hg19.vcf -knownSites Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites 1000G_phase1.indels.hg19.sites.vcf -o 7590_recal_data.grp
java -jar GenomeAnalysisTK.jar -nct 8 -T PrintReads -R ucsc.hg19.fasta -I 7590_realigned_reads.bam -BQSR 7590_recal_data.grp -o 7590_recal_reads.bam
java -jar GenomeAnalysisTK.jar -nct 8 -T BaseRecalibrator -R ucsc.hg19.fasta -I 7590_realigned_reads.bam -knownSites dbsnp_138.hg19.vcf -knownSites Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites 1000G_phase1.indels.hg19.sites.vcf -BQSR 7590_recal_data.grp -o 7590_post_recal_data.grp
java -jar GenomeAnalysisTK.jar -T AnalyzeCovariates -R ucsc.hg19.fasta -before 7590_recal_data.grp -after 7590_post_recal_data.grp -plots 7590_recal_plots.pdf
java -jar GenomeAnalysisTK.jar -nct 8 -T BaseRecalibrator -R ucsc.hg19.fasta -I 7591_realigned_reads.bam -knownSites dbsnp_138.hg19.vcf -knownSites Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites 1000G_phase1.indels.hg19.sites.vcf -o 7591_recal_data.grp
java -jar GenomeAnalysisTK.jar -nct 8 -T PrintReads -R ucsc.hg19.fasta -I 7591_realigned_reads.bam -BQSR 7591_recal_data.grp -o 7591_recal_reads.bam
java -jar GenomeAnalysisTK.jar -nct 8 -T BaseRecalibrator -R ucsc.hg19.fasta -I 7591_realigned_reads.bam -knownSites dbsnp_138.hg19.vcf -knownSites Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites 1000G_phase1.indels.hg19.sites.vcf -BQSR 7591_recal_data.grp -o 7591_post_recal_data.grp
java -jar GenomeAnalysisTK.jar -T AnalyzeCovariates -R ucsc.hg19.fasta -before 7591_recal_data.grp -after 7591_post_recal_data.grp -plots 7591_recal_plots.pdf
java -jar GenomeAnalysisTK.jar -nct 8 -T BaseRecalibrator -R ucsc.hg19.fasta -I 7592_realigned_reads.bam -knownSites dbsnp_138.hg19.vcf -knownSites Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites 1000G_phase1.indels.hg19.sites.vcf -o 7592_recal_data.grp
java -jar GenomeAnalysisTK.jar -nct 8 -T PrintReads -R ucsc.hg19.fasta -I 7592_realigned_reads.bam -BQSR 7592_recal_data.grp -o 7592_recal_reads.bam
java -jar GenomeAnalysisTK.jar -nct 8 -T BaseRecalibrator -R ucsc.hg19.fasta -I 7592_realigned_reads.bam -knownSites dbsnp_138.hg19.vcf -knownSites Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites 1000G_phase1.indels.hg19.sites.vcf -BQSR 7592_recal_data.grp -o 7592_post_recal_data.grp
java -jar GenomeAnalysisTK.jar -T AnalyzeCovariates -R ucsc.hg19.fasta -before 7592_recal_data.grp -after 7592_post_recal_data.grp -plots 7592_recal_plots.pdf
java -jar GenomeAnalysisTK.jar -nct 8 -T BaseRecalibrator -R ucsc.hg19.fasta -I 7701_realigned_reads.bam -knownSites dbsnp_138.hg19.vcf -knownSites Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites 1000G_phase1.indels.hg19.sites.vcf -o 7701_recal_data.grp
java -jar GenomeAnalysisTK.jar -nct 8 -T PrintReads -R ucsc.hg19.fasta -I 7701_realigned_reads.bam -BQSR 7701_recal_data.grp -o 7701_recal_reads.bam
java -jar GenomeAnalysisTK.jar -nct 8 -T BaseRecalibrator -R ucsc.hg19.fasta -I 7701_realigned_reads.bam -knownSites dbsnp_138.hg19.vcf -knownSites Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites 1000G_phase1.indels.hg19.sites.vcf -BQSR 7701_recal_data.grp -o 7701_post_recal_data.grp
java -jar GenomeAnalysisTK.jar -T AnalyzeCovariates -R ucsc.hg19.fasta -before 7701_recal_data.grp -after 7701_post_recal_data.grp -plots 7701_recal_plots.pdf
java -jar GenomeAnalysisTK.jar -nct 8 -T BaseRecalibrator -R ucsc.hg19.fasta -I 7702_realigned_reads.bam -knownSites dbsnp_138.hg19.vcf -knownSites Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites 1000G_phase1.indels.hg19.sites.vcf -o 7702_recal_data.grp
java -jar GenomeAnalysisTK.jar -nct 8 -T PrintReads -R ucsc.hg19.fasta -I 7702_realigned_reads.bam -BQSR 7702_recal_data.grp -o 7702_recal_reads.bam
java -jar GenomeAnalysisTK.jar -nct 8 -T BaseRecalibrator -R ucsc.hg19.fasta -I 7702_realigned_reads.bam -knownSites dbsnp_138.hg19.vcf -knownSites Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites 1000G_phase1.indels.hg19.sites.vcf -BQSR 7702_recal_data.grp -o 7702_post_recal_data.grp
java -jar GenomeAnalysisTK.jar -T AnalyzeCovariates -R ucsc.hg19.fasta -before 7702_recal_data.grp -after 7702_post_recal_data.grp -plots 7702_recal_plots.pdf
java -jar GenomeAnalysisTK.jar -nct 8 -T BaseRecalibrator -R ucsc.hg19.fasta -I 7703_realigned_reads.bam -knownSites dbsnp_138.hg19.vcf -knownSites Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites 1000G_phase1.indels.hg19.sites.vcf -o 7703_recal_data.grp
java -jar GenomeAnalysisTK.jar -nct 8 -T PrintReads -R ucsc.hg19.fasta -I 7703_realigned_reads.bam -BQSR 7703_recal_data.grp -o 7703_recal_reads.bam
java -jar GenomeAnalysisTK.jar -nct 8 -T BaseRecalibrator -R ucsc.hg19.fasta -I 7703_realigned_reads.bam -knownSites dbsnp_138.hg19.vcf -knownSites Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites 1000G_phase1.indels.hg19.sites.vcf -BQSR 7703_recal_data.grp -o 7703_post_recal_data.grp
java -jar GenomeAnalysisTK.jar -T AnalyzeCovariates -R ucsc.hg19.fasta -before 7703_recal_data.grp -after 7703_post_recal_data.grp -plots 7703_recal_plots.pdf
java -jar GenomeAnalysisTK.jar -nct 8 -T BaseRecalibrator -R ucsc.hg19.fasta -I 7704B_realigned_reads.bam -knownSites dbsnp_138.hg19.vcf -knownSites Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites 1000G_phase1.indels.hg19.sites.vcf -o 7704B_recal_data.grp
java -jar GenomeAnalysisTK.jar -nct 8 -T PrintReads -R ucsc.hg19.fasta -I 7704B_realigned_reads.bam -BQSR 7704B_recal_data.grp -o 7704B_recal_reads.bam
java -jar GenomeAnalysisTK.jar -nct 8 -T BaseRecalibrator -R ucsc.hg19.fasta -I 7704B_realigned_reads.bam -knownSites dbsnp_138.hg19.vcf -knownSites Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites 1000G_phase1.indels.hg19.sites.vcf -BQSR 7704B_recal_data.grp -o 7704B_post_recal_data.grp
java -jar GenomeAnalysisTK.jar -T AnalyzeCovariates -R ucsc.hg19.fasta -before 7704B_recal_data.grp -after 7704B_post_recal_data.grp -plots 7704B_recal_plots.pdf
java -jar GenomeAnalysisTK.jar -nct 8 -T BaseRecalibrator -R ucsc.hg19.fasta -I 7706_realigned_reads.bam -knownSites dbsnp_138.hg19.vcf -knownSites Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites 1000G_phase1.indels.hg19.sites.vcf -o 7706_recal_data.grp
java -jar GenomeAnalysisTK.jar -nct 8 -T PrintReads -R ucsc.hg19.fasta -I 7706_realigned_reads.bam -BQSR 7706_recal_data.grp -o 7706_recal_reads.bam
java -jar GenomeAnalysisTK.jar -nct 8 -T BaseRecalibrator -R ucsc.hg19.fasta -I 7706_realigned_reads.bam -knownSites dbsnp_138.hg19.vcf -knownSites Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites 1000G_phase1.indels.hg19.sites.vcf -BQSR 7706_recal_data.grp -o 7706_post_recal_data.grp
java -jar GenomeAnalysisTK.jar -T AnalyzeCovariates -R ucsc.hg19.fasta -before 7706_recal_data.grp -after 7706_post_recal_data.grp -plots 7706_recal_plots.pdf
java -jar GenomeAnalysisTK.jar -nct 8 -T BaseRecalibrator -R ucsc.hg19.fasta -I 7708_realigned_reads.bam -knownSites dbsnp_138.hg19.vcf -knownSites Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites 1000G_phase1.indels.hg19.sites.vcf -o 7708_recal_data.grp
java -jar GenomeAnalysisTK.jar -nct 8 -T PrintReads -R ucsc.hg19.fasta -I 7708_realigned_reads.bam -BQSR 7708_recal_data.grp -o 7708_recal_reads.bam
java -jar GenomeAnalysisTK.jar -nct 8 -T BaseRecalibrator -R ucsc.hg19.fasta -I 7708_realigned_reads.bam -knownSites dbsnp_138.hg19.vcf -knownSites Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites 1000G_phase1.indels.hg19.sites.vcf -BQSR 7708_recal_data.grp -o 7708_post_recal_data.grp
java -jar GenomeAnalysisTK.jar -T AnalyzeCovariates -R ucsc.hg19.fasta -before 7708_recal_data.grp -after 7708_post_recal_data.grp -plots 7708_recal_plots.pdf
java -jar GenomeAnalysisTK.jar -nct 8 -T BaseRecalibrator -R ucsc.hg19.fasta -I 7713_realigned_reads.bam -knownSites dbsnp_138.hg19.vcf -knownSites Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites 1000G_phase1.indels.hg19.sites.vcf -o 7713_recal_data.grp
java -jar GenomeAnalysisTK.jar -nct 8 -T PrintReads -R ucsc.hg19.fasta -I 7713_realigned_reads.bam -BQSR 7713_recal_data.grp -o 7713_recal_reads.bam
java -jar GenomeAnalysisTK.jar -nct 8 -T BaseRecalibrator -R ucsc.hg19.fasta -I 7713_realigned_reads.bam -knownSites dbsnp_138.hg19.vcf -knownSites Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites 1000G_phase1.indels.hg19.sites.vcf -BQSR 7713_recal_data.grp -o 7713_post_recal_data.grp
java -jar GenomeAnalysisTK.jar -T AnalyzeCovariates -R ucsc.hg19.fasta -before 7713_recal_data.grp -after 7713_post_recal_data.grp -plots 7713_recal_plots.pdf
### 10 Identify	potential variants in each sample (HaplotypeCaller)
### output ends with .g.vcf
for i in 7588 7590 7591 7592 7701 7702 7703 7704B 7706 7708 7713
do
    echo "java -Xmx32g -jar GenomeAnalysisTK.jar -nct 8 -T HaplotypeCaller -R ucsc.hg19.fasta -I "$i"_recal_reads.bam -o "$i".g.vcf -ERC GVCF"
done
java -Xmx32g -jar GenomeAnalysisTK.jar -nct 8 -T HaplotypeCaller -R ucsc.hg19.fasta -I 7588_recal_reads.bam -o 7588.g.vcf -ERC GVCF
java -Xmx32g -jar GenomeAnalysisTK.jar -nct 8 -T HaplotypeCaller -R ucsc.hg19.fasta -I 7590_recal_reads.bam -o 7590.g.vcf -ERC GVCF
java -Xmx32g -jar GenomeAnalysisTK.jar -nct 8 -T HaplotypeCaller -R ucsc.hg19.fasta -I 7591_recal_reads.bam -o 7591.g.vcf -ERC GVCF
java -Xmx32g -jar GenomeAnalysisTK.jar -nct 8 -T HaplotypeCaller -R ucsc.hg19.fasta -I 7592_recal_reads.bam -o 7592.g.vcf -ERC GVCF
java -Xmx32g -jar GenomeAnalysisTK.jar -nct 8 -T HaplotypeCaller -R ucsc.hg19.fasta -I 7701_recal_reads.bam -o 7701.g.vcf -ERC GVCF
java -Xmx32g -jar GenomeAnalysisTK.jar -nct 8 -T HaplotypeCaller -R ucsc.hg19.fasta -I 7702_recal_reads.bam -o 7702.g.vcf -ERC GVCF
java -Xmx32g -jar GenomeAnalysisTK.jar -nct 8 -T HaplotypeCaller -R ucsc.hg19.fasta -I 7703_recal_reads.bam -o 7703.g.vcf -ERC GVCF
java -Xmx32g -jar GenomeAnalysisTK.jar -nct 8 -T HaplotypeCaller -R ucsc.hg19.fasta -I 7704B_recal_reads.bam -o 7704B.g.vcf -ERC GVCF
java -Xmx32g -jar GenomeAnalysisTK.jar -nct 8 -T HaplotypeCaller -R ucsc.hg19.fasta -I 7706_recal_reads.bam -o 7706.g.vcf -ERC GVCF
java -Xmx32g -jar GenomeAnalysisTK.jar -nct 8 -T HaplotypeCaller -R ucsc.hg19.fasta -I 7708_recal_reads.bam -o 7708.g.vcf -ERC GVCF
java -Xmx32g -jar GenomeAnalysisTK.jar -nct 8 -T HaplotypeCaller -R ucsc.hg19.fasta -I 7713_recal_reads.bam -o 7713.g.vcf -ERC GVCF
### 11 Perform joint genotyping on the cohort (GenotypeGVCFs)
java -Xmx32g -jar GenomeAnalysisTK.jar -nt 8 -T GenotypeGVCFs -R ucsc.hg19.fasta -V 7588.g.vcf -V 7590.g.vcf -V 7591.g.vcf -V 7592.g.vcf -V 7701.g.vcf -V 7702.g.vcf -V 7703.g.vcf -V 7704B.g.vcf -V 7706.g.vcf -V 7708.g.vcf -V 7713.g.vcf -o ais.raw.vcf
### 12 Recalibrate variant quality scores for SNPs
java -Xmx32g -jar GenomeAnalysisTK.jar -nt 8 -T VariantRecalibrator -R ucsc.hg19.fasta -input ais.raw.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg19.sites.vcf -resource:omni,known=false,training=true,truth=true,prior=12.0 1000G_omni2.5.hg19.sites.vcf -resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.indels.hg19.sites.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_138.hg19.vcf -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile recalibrate_SNP.recal -tranchesFile recalibrate_SNP.tranches -rscriptFile recalibrate_SNP_plots.R
java -Xmx32g -jar GenomeAnalysisTK.jar -nt 8 -T ApplyRecalibration -R ucsc.hg19.fasta -input ais.raw.vcf -mode SNP --ts_filter_level 99.0 -recalFile recalibrate_SNP.recal -tranchesFile recalibrate_SNP.tranches -o ais_recalibrated_snps_raw_indels.vcf
### 13 Recalibrate variant quality scores for Indels
java -Xmx32g -jar GenomeAnalysisTK.jar -nt 8 -T VariantRecalibrator -R ucsc.hg19.fasta -input ais_recalibrated_snps_raw_indels.vcf -resource:mills,known=true,training=true,truth=true,prior=12.0 Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -an DP -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 -recalFile recalibrate_INDEL.recal -tranchesFile recalibrate_INDEL.tranches -rscriptFile recalibrate_INDEL_plots.R
java -Xmx32g -jar GenomeAnalysisTK.jar -nt 8 -T ApplyRecalibration -R ucsc.hg19.fasta -input ais_recalibrated_snps_raw_indels.vcf -mode INDEL --ts_filter_level 99.0 -recalFile recalibrate_INDEL.recal -tranchesFile recalibrate_INDEL.tranches -o ais_recalibrated_variants.vcf

