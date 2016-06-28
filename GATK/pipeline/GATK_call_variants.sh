#!/bin/bash

for i in "$@"
do
case $i in
    -d=*|--download=*)
    DOWLOAD="${i#*=}"
    shift # past argument=value
    ;;
    -i=*|--indexing=*)
    INDEX="${i#*=}"
    shift # past argument=value
    ;;
    -m=*|--mapping=*)
    MAP="${i#*=}"
    shift # past argument=value
    ;;
    -s=*|--sort=*)
    SORT="${i#*=}"
    shift # past argument=value
    ;;
    -r=*|--realign=*)
    REALIGN="${i#*=}"
    shift # past argument=value
    ;;
    -q=*|--bqsr=*)
    BQSR="${i#*=}"
    shift # past argument=value
    ;;
    -p=*|--plot=*)
    PLOT="${i#*=}"
    shift # past argument=value
    ;;
    -h=*|--haplotypeCaller=*)
    HTC="${i#*=}"
    shift # past argument=value
    ;;
    -z=*|--rmtmp=*)
    DEL="${i#*=}"
    shift # past argument=value
    ;;
    -a=*|--fq1=*)
    FQ1="${i#*=}"
    shift # past argument=value
    ;;
    -b=*|--fq2=*)
    FQ2="${i#*=}"
    shift # past argument=value
    ;;
    -g=*|--ref=*)
    REF="${i#*=}"
    shift # past argument=value
    ;;
    -k=*|--kg=*)
    KG="${i#*=}"
    shift # past argument=value
    ;;
    -x=*|--gold=*)
    GOLD="${i#*=}"
    shift # past argument=value
    ;;
    -y=*|--dbsnp=*)
    DBSNP="${i#*=}"
    shift # past argument=value
    ;;
    -c=*|--nt=*)
    CORE="${i#*=}"
    shift # past argument=value
    ;;
    -o=*|--out=*)
    OUT="${i#*=}"
    shift # past argument=value
    ;;
    --default)
    DEFAULT=YES
    shift # past argument with no value
    ;;
    *)
            # unknown option
    ;;
esac
done

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
    java -jar CreateSequenceDictionary.jar REFERENCE=ucsc.hg19.fasta OUTPUT=ucsc.hg19.dict
fi

################################################################################
#### step 3 mapping with bwa mem
#### paired-end alignment
#### ID:<unique id> LB:<library name> SM:<sample name> PL:<platform name>
if [ $MAP == 'true' ]; then
    bwa mem -M -t $CORE $REF $FQ1 $FQ2 > $OUT.align.sam
    if [ $DEL == 'true' ]; then
        rm $FQ1 $FQ2
    fi
fi

################################################################################
#### step 4  sorting and marking duplicates
#### sort sam and output as bam
#### marking duplicates and index the bam
if [ $SORT == 'true' ]; then
    java -jar SortSam.jar INPUT=$OUT.align.sam OUTPUT=$OUT.sorted.bam SORT_ORDER=coordinate
    if [ $DEL == 'true' ]; then
        rm $OUT.align.sam
    fi
    java -jar MarkDuplicates.jar INPUT=$OUT.sorted.bam OUTPUT=$OUT.dedup.bam METRICS_FILE=$OUT.metrics.txt
    if [ $DEL == 'true' ]; then
        rm $OUT.sorted.bam
    fi
    java -jar BuildBamIndex.jar INPUT=$OUT.dedup.bam
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
fi

if [ $DEL == 'true' ]; then
    rm $OUT.realigned.bam
fi

################################################################################
#### step 7 Identify  potential variants in each sample (HaplotypeCaller)
#### output records for all sites in gVCF format
if [ $HTC == 'true' ]; then
    java -jar -Xmx32g GenomeAnalysisTK.jar -nct $CORE -T HaplotypeCaller -R $REF -I $OUT.recal.bam -o $OUT.g.vcf -ERC GVCF
fi

################################################################################
#### The END ####
