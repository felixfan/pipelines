################################################################################
####                  pipeline for ngs analysis: part 0                     ####
####                  1 Download files                                      ####
####                  2 Identify potential variants in each sample          ####
####                  !!!! CAN NOT RUN WITH 'nohup'!!!!                     ####
####                  !!!! For PAIR-END FASTQ ONLY !!!!                     ####
####                  last modified on 8 Oct 2015                           ####
####                  @ Felix Yanhui Fan, felixfanyh@gmail.com              ####
################################################################################
### download bundle data from gatk and index
#### down load data from GATK ftp://ftp.broadinstitute.org/bundle/2.8/hg19/
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
bwa index ucsc.hg19.fasta
samtools faidx ucsc.hg19.fasta
java -jar CreateSequenceDictionary.jar REFERENCE=ucsc.hg19.fasta OUTPUT=ucsc.hg19.dict

