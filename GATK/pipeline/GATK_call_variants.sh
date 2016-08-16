#!/bin/bash

#### pipeline for WGS variant calling - part 1
#### use the best practices for variant calling with GATK
#### software used: bwa, samtools, picard, GATK
#### total seven steps:
#### step 1 Download reference files [--download]
#### step 2 create the index and dictionary of the reference genome [--indexing]
#### step 3 mapping with bwa mem [--mapping]
#### step 4  sorting and marking duplicates [--sort]
#### step 5 Local realignment around indels [--realign]
#### step 6 Base quality score recalibration (BQSR) [--bqsr] [--plot]
#### step 7 Identify  potential variants in each sample (HaplotypeCaller) [--haplotypeCaller]
#### how to run:
#### example 1: run steps 1-7 with plot
#### sh GATK_call_variants.sh --full --fq1=588_R1.fastq.gz --fq2=7588_R2.fastq.gz --sample=7588 --ref=ucsc.hg19.fasta --kg=1000G_phase1.indels.hg19.sites.vcf --gold=Mills_and_1000G_gold_standard.indels.hg19.sites.vcf --dbsnp=dbsnp_138.hg19.vcf
#### example 2: run steps 3 - 7 with plot
#### sh GATK_call_variants.sh --default --fq1=588_R1.fastq.gz --fq2=7588_R2.fastq.gz --sample=7588 --ref=ucsc.hg19.fasta --kg=1000G_phase1.indels.hg19.sites.vcf --gold=Mills_and_1000G_gold_standard.indels.hg19.sites.vcf --dbsnp=dbsnp_138.hg19.vcf
#### example 3: run steps 3 - 7 without plot
#### sh GATK_call_variants.sh --fast --fq1=588_R1.fastq.gz --fq2=7588_R2.fastq.gz --sample=7588 --ref=ucsc.hg19.fasta --kg=1000G_phase1.indels.hg19.sites.vcf --gold=Mills_and_1000G_gold_standard.indels.hg19.sites.vcf --dbsnp=dbsnp_138.hg19.vcf
#### example 4: only run steps 1 - 2
#### sh GATK_call_variants.sh --download --indexing
#### Default setting:
#### By default, 2 threads are used (--nt=2). This can be changed, e.g. --nt=8
#### By default, all temporary files will be deleted after each step when --full, --default, or --fast was used. Only BAM file of BQSR was kept when --plot was used, else, both BAM files of BQSR and realignment were kept, so you can make the plot later.  
#### By default, RG is "@RG\tID:${ID}\tSM:${SM}\tPL:${PL}\tLB:${LB}", where ID="${SAMPLE}.id", SM=${SAMPLE}, PL="ILLUMINA", LB="${SAMPLE}.lib"
echo "@-------------------------------------------------------------@"
echo "|     GATK_call_variants    |     v1.0.0    |   29 Jun 2016   |"
echo "|-------------------------------------------------------------|"
echo "|  (C) 2016 Felix Yanhui Fan, GNU General Public License, v2  |"
echo "|-------------------------------------------------------------|"
echo "|    For documentation, citation & bug-report instructions:   |"
echo "|            http://felixfan.github.io/pipelines              |"
echo "@-------------------------------------------------------------@"
for i in "$@"
do
case $i in
    --download)
    DOWNLOAD=true
    shift # past argument with no value
    ;;
    --indexing)
    INDEX=true
    shift # past argument with no value
    ;;
    --mapping)
    MAP=true
    shift # past argument with no value
    ;;
    --sort)
    SORT=true
    shift # past argument with no value
    ;;
    --realign)
    REALIGN=true
    shift # past argument with no value
    ;;
    --bqsr)
    BQSR=true
    shift # past argument with no value
    ;;
    --plot)
    PLOT=true
    shift # past argument with no value
    ;;
    --haplotypeCaller)
    HTC=true
    shift # past argument with no value
    ;;
    --rmtmp)
    DEL=true
    shift # past argument with no value
    ;;
    --full)
    FULL=true
    shift # past argument with no value
    ;;
    --fast)
    FAST=true
    shift # past argument with no value
    ;;
    --default)
    DEFAULT=true
    shift # past argument with no value
    ;;
    --fq1=*)
    FQ1="${i#*=}"
    shift # past argument=value
    ;;
    --fq2=*)
    FQ2="${i#*=}"
    shift # past argument=value
    ;;
    --sample=*)
    SAMPLE="${i#*=}"
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
#### check which steps to run
flag=false
if [ $DEFAULT ]; then
    DOWNLOAD=false
    INDEX=false
    MAP=true
    SORT=true
    REALIGN=true
    BQSR=true
    PLOT=true
    HTC=true
    DEL=true
    flag=true
elif [ $FULL ]; then
    DOWNLOAD=true
    INDEX=true
    MAP=true
    SORT=true
    REALIGN=true
    BQSR=true
    PLOT=true
    HTC=true
    DEL=true
    flag=true
elif [ $FAST ]; then
    DOWNLOAD=false
    INDEX=false
    MAP=true
    SORT=true
    REALIGN=true
    BQSR=true
    PLOT=false
    HTC=true
    DEL=true
    flag=true
elif [ $DOWNLOAD ]; then
    DOWNLOAD=true
elif [ $INDEX ]; then
    INDEX=true
elif [ $MAP ]; then
    MAP=true
elif [ $SORT ]; then
    SORT=true
elif [ $REALIGN ]; then
    REALIGN=true
elif [ $BQSR ]; then
    BQSR=true
elif [ $PLOT ]; then
    PLOT=true
elif [ $HTC ]; then
    HTC=true
else
    echo 'nothing to do, exit now'
    exit 2
fi
if [ $flag != true ]; then
    if [ ! $DOWNLOAD ]; then
        DOWNLOAD=false
    fi
    if [ ! $INDEX ]; then
        INDEX=false
    fi
    if [ ! $MAP ]; then
        MAP=false
    fi
    if [ ! $SORT ]; then
        SORT=false
    fi
    if [ ! $REALIGN ]; then
        REALIGN=false
    fi
    if [ ! $BQSR ]; then
        BQSR=false
    fi
    if [ ! $PLOT ]; then
        PLOT=false
    fi
    if [ ! $HTC ]; then
        HTC=false
    fi
    if [ ! $DEL ]; then
        DEL=false
    fi
fi
#### check input files
if [ $MAP = true ] || [ $SORT = true ] || [ $REALIGN = true ] || [ $BQSR = true ] || [ $PLOT = true ] || [ $HTC = true ]; then
    if [ ! $FQ1 ]; then
        echo 'missing the first fastq file!'
        exit 1
    fi
    if [ ! $FQ2 ]; then
        echo 'missing the second fastq file!'
        exit 1
    fi
    if [ ! $SAMPLE ]; then
        echo 'missing the sample ID!'
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
fi
if [ ! $CORE ]; then
    CORE=2
fi
#### effected options
echo "\n\tOptions in effect:"
if [ $DOWNLOAD = true ]; then
    echo "\t--download"
fi
if [ $INDEX = true ]; then
    echo "\t--indexing"
fi
if [ $MAP = true ]; then
    echo "\t--mapping"
fi
if [ $SORT = true ]; then
    echo "\t--sort"
fi
if [ $REALIGN = true ]; then
    echo "\t--realign"
fi
if [ $BQSR = true ]; then
    echo "\t--bqsr"
fi
if [ $PLOT = true ]; then
    echo "\t--plot"
fi
if [ $HTC = true ]; then
    echo "\t--haplotypeCaller"
fi
if [ $DEL = true ]; then
    echo "\t--rmtmp"
fi
echo "\t--nt=$CORE"
if [ $MAP = true ] || [ $SORT = true ] || [ $REALIGN = true ] || [ $BQSR = true ] || [ $PLOT = true ] || [ $HTC = true ]; then
    echo "\t--fq1=$FQ1"
    echo "\t--fq2=$FQ2"
    echo "\t--sample=$SAMPLE"
    echo "\t--ref=$REF"
    echo "\t--kg=$KG"
    echo "\t--gold=$GOLD"
    echo "\t--dbsnp=$DBSNP"
fi
echo
####
#for key in DOWNLOAD INDEX MAP SORT REALIGN BQSR PLOT HTC DEL FQ1 FQ2 SAMPLE REF KG GOLD DBSNP CORE
#do
#    echo "\t${key} = ${!key}"
#done
#### RG
ID="${SAMPLE}.id"
SM=${SAMPLE}
PL="ILLUMINA"
LB="${SAMPLE}.lib"
RG="@RG\tID:${ID}\tSM:${SM}\tPL:${PL}\tLB:${LB}"
#### out
OUT=$SAMPLE

# exit

################################################################################
#### step 1 Download reference files
#### download bundle data from gatk and index
#### down load data from GATK ftp://ftp.broadinstitute.org/bundle/2.8/hg19/
if [ $DOWNLOAD == 'true' ]; then
    wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta.gz
    wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/dbsnp_138.hg19.vcf.gz
    wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
    wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/1000G_phase1.indels.hg19.sites.vcf.gz
    wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/1000G_omni2.5.hg19.sites.vcf.gz
    wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/hapmap_3.3.hg19.sites.vcf.gz
    gunzip ucsc.hg19.fasta.gz
    gunzip dbsnp_138.hg19.vcf.gz
    gunzip Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
    gunzip 1000G_phase1.indels.hg19.sites.vcf.gz
    gunzip 1000G_omni2.5.hg19.sites.vcf.gz
    gunzip hapmap_3.3.hg19.sites.vcf.gz
fi
################################################################################
#### step 2 create the index and dictionary of the reference genome
if [ $INDEX == 'true' ]; then
    bwa index ucsc.hg19.fasta
    samtools faidx ucsc.hg19.fasta
    #java -jar CreateSequenceDictionary.jar REFERENCE=ucsc.hg19.fasta OUTPUT=ucsc.hg19.dict
    java -jar picard.jar CreateSequenceDictionary R=ucsc.hg19.fasta O=ucsc.hg19.dict # picard 2.6
fi
################################################################################
#### step 3 mapping with bwa mem
#### paired-end alignment
#### ID:<unique id> LB:<library name> SM:<sample name> PL:<platform name>
if [ $MAP == 'true' ]; then
    bwa mem -M -t $CORE -R $RG $REF $FQ1 $FQ2 > $OUT.align.sam
    if [ $DEL == 'true' ]; then
        rm $FQ1 $FQ2
    fi
fi
################################################################################
#### step 4  sorting and marking duplicates
#### sort sam and output as bam
#### marking duplicates and index the bam
if [ $SORT == 'true' ]; then
    #java -jar SortSam.jar INPUT=$OUT.align.sam OUTPUT=$OUT.sorted.bam SORT_ORDER=coordinate
    java -jar picard.jar SortSam I=$OUT.align.sam O=$OUT.sorted.bam SORT_ORDER=coordinate
    if [ $DEL == 'true' ]; then
        rm $OUT.align.sam
    fi
    #java -jar MarkDuplicates.jar INPUT=$OUT.sorted.bam OUTPUT=$OUT.dedup.bam METRICS_FILE=$OUT.metrics.txt
    java -jar picard.jar MarkDuplicates I=$OUT.sorted.bam O=$OUT.dedup.bam M=$OUT.metrics.txt
    if [ $DEL == 'true' ]; then
        rm $OUT.sorted.bam
    fi
    #java -jar BuildBamIndex.jar INPUT=$OUT.dedup.bam
    java -jar picard.jar BuildBamIndex I=$OUT.dedup.bam
fi
################################################################################
#### step 5 Local realignment around indels
#### create a target list of intervals to be realigned
#### perform realignment
if [ $REALIGN == 'true' ]; then
    java -jar GenomeAnalysisTK.jar -nt $CORE -T RealignerTargetCreator -R $REF -I $OUT.dedup.bam -known $GOLD -known $KG -o $OUT.intervals
    java -jar GenomeAnalysisTK.jar -T IndelRealigner -R $REF -targetIntervals $OUT.intervals -I $OUT.dedup.bam -known $GOLD -known $KG -o $OUT.realigned.bam
    if [ $DEL == 'true' ]; then
        rm $OUT.dedup.bam
    fi
fi
################################################################################
#### step 6 Base quality score recalibration (BQSR)
#### analyze patterns of covariation and builds recalibration model
#### apply the recalibration
if [ $BQSR == 'true' ]; then
    java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R $REF -I $OUT.realigned.bam -knownSites $DBSNP -knownSites $GOLD -knownSites $KG -o $OUT.recal.grp
    java -jar GenomeAnalysisTK.jar -T PrintReads -R $REF -I $OUT.realigned.bam -BQSR $OUT.recal.grp -o $OUT.recal.bam
fi
#### second pass evaluateds what the data looks like after recalibration
#### makes plots based on before/after recalibration tables
if [ $PLOT == 'true' ]; then
    java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R $REF -I $OUT.realigned.bam -knownSites $DBSNP -knownSites $GOLD -knownSites $KG -BQSR $OUT.recal.grp -o $OUT.post.recal.grp
    java -jar GenomeAnalysisTK.jar -T AnalyzeCovariates -R $REF -before $OUT.recal.grp -after $OUT.post.recal.grp -plots $OUT.recal.plots.pdf
    if [ $DEL == 'true' ]; then
        rm $OUT.realigned.bam
    fi
fi
################################################################################
#### step 7 Identify  potential variants in each sample (HaplotypeCaller)
#### output records for all sites in gVCF format
if [ $HTC == 'true' ]; then
    java -jar -Xmx32g GenomeAnalysisTK.jar -nct $CORE -T HaplotypeCaller -R $REF -I $OUT.recal.bam -o $OUT.g.vcf -ERC GVCF
fi
################################################################################
#### The END ####
