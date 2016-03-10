################################################################################
####                  pipeline for ngs analysis: part 3                     ####
####                  4 Recalibrate variant quality scores                  ####
####                  last modified on 8 Oct 2015                           ####
####                  @ Felix Yanhui Fan, felixfanyh@gmail.com              ####
################################################################################
### usage: sh wgsPipeline3.sh 8 32 ais.raw.vcf ais
### e.g. uses 8 core & 32g memory
### !!!!!! add path !!!!
CORE=$1
MEM=$2
IN=$3
out=$4
REF="/home/felixfan/ulib/hg19/ucsc.hg19.fasta"
GOLD="/home/felixfan/ulib/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
KG="/home/felixfan/ulib/hg19/1000G_phase1.indels.hg19.sites.vcf"
DBSNP="/home/felixfan/ulib/hg19/dbsnp_138.hg19.vcf"
HAPMAP="/home/felixfan/ulib/hg19/hapmap_3.3.hg19.sites.vcf"
OMNI="/home/felixfan/ulib/hg19/1000G_omni2.5.hg19.sites.vcf"
### 12 Recalibrate variant quality scores for SNPs
mymem="-Xmx"$MEM"g"
java $mymem -jar GenomeAnalysisTK.jar -nt $CORE -T VariantRecalibrator -R $REF -input $IN -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP -resource:omni,known=false,training=true,truth=true,prior=12.0 $OMNI -resource:1000G,known=false,training=true,truth=false,prior=10.0 $KG -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile recalibrate_SNP.recal -tranchesFile recalibrate_SNP.tranches -rscriptFile recalibrate_SNP_plots.R
java $mymem -jar GenomeAnalysisTK.jar -nt $CORE -T ApplyRecalibration -R $REF -input $IN -mode SNP --ts_filter_level 99.0 -recalFile recalibrate_SNP.recal -tranchesFile recalibrate_SNP.tranches -o $out.recalibrated_snps_raw_indels.vcf
### 13 Recalibrate variant quality scores for Indels
java $mymem -jar GenomeAnalysisTK.jar -nt $CORE -T VariantRecalibrator -R $REF -input $out.recalibrated_snps_raw_indels.vcf -resource:mills,known=true,training=true,truth=true,prior=12.0 $GOLD -an DP -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 -recalFile recalibrate_INDEL.recal -tranchesFile recalibrate_INDEL.tranches -rscriptFile recalibrate_INDEL_plots.R
java $mymem -jar GenomeAnalysisTK.jar -nt $CORE -T ApplyRecalibration -R $REF -input $out.recalibrated_snps_raw_indels.vcf -mode INDEL --ts_filter_level 99.0 -recalFile recalibrate_INDEL.recal -tranchesFile recalibrate_INDEL.tranches -o $out.recalibrated_variants.vcf

