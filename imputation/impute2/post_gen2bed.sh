### input is impute2 output genotypes
### include only SNPs
### convert gen to ped
### convert missing 'N' to '0'
### update SNP names
### output bed/bim/fam

### GTOOL (http://www.well.ox.ac.uk/~cfreeman/software/gwas/gtool.html)

### by Yanhui Fan (felixfanyh@gmail.com)
### revised on 29 April 2016

### usage: sh post_gen2bed.sh -g=out.phased.impute2 -s=user.sample -c=1 -t=0.9 -p=phenotype -x=sex -o=out.chr1

for i in "$@"
do
case $i in
    -g=*|--genotype=*)
    GENOTYPE="${i#*=}"
    shift # past argument=value
    ;;
    -s=*|--sample=*)
    SAMPLE="${i#*=}"
    shift # past argument=value
    ;;
    -c=*|--chr=*)
    CHR="${i#*=}"
    shift # past argument=value
    ;;
    -t=*|--threshold=*)
    THRESHOLD="${i#*=}"
    shift # past argument=value
    ;;
    -p=*|--phenotype=*)
    PHENOTYPE="${i#*=}"
    shift # past argument=value
    ;;
    -x=*|--sex=*)
    SEX="${i#*=}"
    shift # past argument=value
    ;;
    -o=*|--out=*)
    OUTPUT="${i#*=}"
    shift # past argument=value
    ;;
    *)
            # unknown option
    ;;
esac
done

### extract SNPs only
awk '$4~/^[ACGT]$/ && $5~/^[ACGT]$/ {print $2}' > tmp.snps
gtool -S --g ${GENOTYPE} --s ${SAMPLE} --inclusion tmp.snps --og tmp.gens

### convert gen to ped
gtool -G --g tmp.gens --s ${SAMPLE} --ped tmp.ped --map tmp.map --phenotype ${PHENOTYPE} --sex ${SEX} --chr ${CHR} --threshold ${THRESHOLD} --snp
rm tmp.*

### gtool code missing as 'N',  => plink code missing as '0'
python recodeMissing.py tmp.ped tmp2.ped
mv tmp2.ped tmp.ped

### update ID  e.g. rs123:2277777:C:T => rs123; 1:222333:G:A => 1:222333
python updateName.py tmp2.ped tmp.name.txt
plink --file tmp --update-map tmp.name.txt --update-name --make-bed --out ${OUT}
