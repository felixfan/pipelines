### input is genome-wide genotype (input.bed/input.bim/input.fam)
### split into chromsome, only include autosomes
### convert ped to gen
### output is output.chr1.gen, output.chr2.sample, ...

### GTOOL (http://www.well.ox.ac.uk/~cfreeman/software/gwas/gtool.html)
### PLINK (http://pngu.mgh.harvard.edu/~purcell/plink/)

### by Yanhui Fan (felixfanyh@gmail.com)
### revised on 29 April 2016

### usage: sh pre_ped2gen.sh -i=input -o=output

for i in "$@"
do
case $i in
    -i=*|--input=*)
    INPUT="${i#*=}"
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

for chr in $(seq 1 22)
do
    plink --bfile ${INPUT} --chr ${chr} --recode --out tmp
    gtool -P --ped tmp.ped --map tmp.map --binary_phenotype --og ${OUTPUT}.chr${chr}.gen --os ${OUTPUT}.chr${chr}.sample
    rm tmp*
done
