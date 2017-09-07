#!/bin/bash

#### pipeline for WGS variant calling - part 3
#### use the best practices for variant calling with GATK
#### software used: GATK
#### how to use:
#### example:
#### sh GATK_variants_recal.wes.sh --vcf=input.raw.vcf --ref=ucsc.hg19.fasta --kg=1000G_phase1.indels.hg19.sites.vcf --gold=Mills_and_1000G_gold_standard.indels.hg19.sites.vcf --dbsnp=dbsnp_138.hg19.vcf --hapmap=hapmap_3.3.hg19.sites.vcf --omni=1000G_omni2.5.hg19.sites.vcf --nt=8 --out=output

echo "@-------------------------------------------------------------@"
echo "|    GATK_variants_recal.wes   |   v2.0.0   |   07 Sep 2017   |"
echo "|-------------------------------------------------------------|"
echo "|  (C) 2017 Felix Yanhui Fan, GNU General Public License, v2  |"
echo "|-------------------------------------------------------------|"
echo "|    For documentation, citation & bug-report instructions:   |"
echo "|            http://felixfan.github.io/pipelines              |"
echo "@-------------------------------------------------------------@"
for i in "$@"
do
case $i in
    --vcf=*)
    VCF="${i#*=}"
    shift # past argument=value
    ;;
    --ref=*)
    REF="${i#*=}"
    shift # past argument=value
    ;;
    --kg=*)
    KG="${i#*=}"
    shift # past argument=value
    ;;
    --gold=*)
    GOLD="${i#*=}"
    shift # past argument=value
    ;;
    --dbsnp=*)
    DBSNP="${i#*=}"
    shift # past argument=value
    ;;
    --hapmap=*)
    HAPMAP="${i#*=}"
    shift # past argument=value
    ;;
    --omni=*)
    OMNI="${i#*=}"
    shift # past argument=value
    ;;
    --out=*)
    OUT="${i#*=}"
    shift # past argument=value
    ;;
    --nt=*)
    CORE="${i#*=}"
    shift # past argument=value
    ;;
    *)
            # unknown option
    ;;
esac
done
################################################################################
#### check arguments
if [ ! $VCF ]; then
    echo 'missing INPUT VCF file!'
    exit 1
fi
if [ ! $REF ]; then
    echo 'missing the reference genome fasta file!'
    exit 1
fi
if [ ! $KG ]; then
    echo 'missing the 1000G indels file'
    exit 1
fi
if [ ! $GOLD ]; then
    echo 'missing the gold standard indels file'
    exit 1
fi
if [ ! $DBSNP ]; then
    echo 'missing the dbSNPs file'
    exit 1
fi
if [ ! $HAPMAP ]; then
    echo 'missing the HapMap VCF file'
    exit 1
fi
if [ ! $OMNI ]; then
    echo 'missing the OMNI VCF file'
    exit 1
fi
if [ ! $OUT ]; then
    OUT='output'
fi
if [ ! $CORE ]; then
    CORE=2
fi
#### effected options
echo -e "\n\tOptions in effect:"
echo -e "\t--vcf=$VCF"
echo -e "\t--ref=$REF"
echo -e "\t--kg=$KG"
echo -e "\t--gold=$GOLD"
echo -e "\t--dbsnp=$DBSNP"
echo -e "\t--hapmap=$HAPMAP"
echo -e "\t---omni=$OMNI"
echo -e "\t--nt=$CORE"
echo -e "\t--out=$OUT"
echo

# exit

#######################################################################################
#### step 1 Recalibrate variant quality scores for SNPs
java -Xmx32g -jar GenomeAnalysisTK.jar -nt $CORE -T VariantRecalibrator -R $REF -input $VCF -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP -resource:omni,known=false,training=true,truth=true,prior=12.0 $OMNI -resource:1000G,known=false,training=true,truth=false,prior=10.0 $KG -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile ${OUT}.SNP.recal -tranchesFile ${OUT}.SNP.tranches -rscriptFile ${OUT}.SNP.plots.R
java -Xmx32g -jar GenomeAnalysisTK.jar -nt $CORE -T ApplyRecalibration -R $REF -input $VCF -mode SNP --ts_filter_level 99.0 -recalFile ${OUT}.SNP.recal -tranchesFile ${OUT}.SNP.tranches -o ${OUT}.recal.SNPs.vcf

### step 2 Recalibrate variant quality scores for Indels
java -Xmx32g -jar GenomeAnalysisTK.jar -nt $CORE -T VariantRecalibrator -R $REF -input ${OUT}.recal.SNPs.vcf -resource:mills,known=true,training=true,truth=true,prior=12.0 $GOLD -an DP -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 -recalFile ${OUT}.INDELs.recal -tranchesFile ${OUT}.INDELs.tranches -rscriptFile ${OUT}.INDELs.plots.R
java -Xmx32g -jar GenomeAnalysisTK.jar -nt $CORE -T ApplyRecalibration -R $REF -input ${OUT}.recal.SNPs.vcf -mode INDEL --ts_filter_level 99.0 -recalFile ${OUT}.INDELs.recal -tranchesFile ${OUT}.INDELs.tranches -o ${OUT}.recal.variants.vcf
