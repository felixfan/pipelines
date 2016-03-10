################################################################################
####                  pipeline for ngs analysis: part 2                     ####
####                  3 Perform joint genotyping on the cohort              ####
####                  last modified on 6 Oct 2015                           ####
####                  @ Felix Yanhui Fan, felixfanyh@gmail.com              ####
################################################################################
### usage: sh wgsPipeline2.sh 8 32 ais
### e.g. uses 8 core & 32g memory
CORE=$1
MEM=$2
OUT=$3
REF="ucsc.hg19.fasta"
### Perform joint genotyping on the cohort (GenotypeGVCFs)
comm="java -Xmx"$MEM"g -jar GenomeAnalysisTK.jar -nt "$CORE" -T GenotypeGVCFs -R "$REF
### !!!! ADD ALL .g.vcf !!!!
### !!!!!!!!!!! list of all file name (individal id) !!!!!!!!!
for i in 7588 7590 7591 7592 7701 7702 7703 7704B 7706 7708 7713
do
    comm=$comm" -V "$i".g.vcf"
done
comm=$comm" -o "$OUT".raw.vcf"
echo $comm > tempwgs.sh
chmod +x tempwgs.sh
sh tempwgs.sh
rm tempwgs.sh

