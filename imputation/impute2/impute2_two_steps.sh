#!/bin/bash

#### impute2 pipeline
#### imputation with one phased reference panel (pre-phasing)
#### autosome: chr1-chr22
#### step1: pre-phasing
#### step2: imputation into pre-phased haplotypes

#### !!! IMPORTANT !!! ####
#### STRAND ALIGNMENT OPTIONS (https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#strand_options)
#### Here I used THE -align_by_maf_g option to activate the program's internal strand alignment procedure for the -g file.
#### If you know the strand information, use the -strand_g option to provide a strand file to the program.

#### suggested Ne=20000

#### input: reference(map/legend/haplotype), genotypes, region (e.g. from 1bp to 5000000 bp), Ne (e.g. 20000)
#### output: out.phased.impute2, ...

#### IMPUTE2 (https://mathgen.stats.ox.ac.uk/impute/impute_v2.html)

#### by Yanhui Fan felixfanyh@gmail.com
#### last revised on 29 April 2016

#### usage: sh impute2_two_steps.sh -m=ref.map -h=ref.haps -l=ref.legend -g=user.gens -f=1 -t=5000000 -n=20000 -o=out

for i in "$@"
do
case $i in
    -m=*|--map=*)
    MAP="${i#*=}"
    shift # past argument=value
    ;;
    -h=*|--haplotype=*)
    HAPLOTYPE="${i#*=}"
    shift # past argument=value
    ;;
    -l=*|--legend=*)
    LEGEND="${i#*=}"
    shift # past argument=value
    ;;
    -g=*|--genotype=*)
    GENOTYPE="${i#*=}"
    shift # past argument=value
    ;;
    -f=*|--from=*)
    FROM="${i#*=}"
    shift # past argument=value
    ;;
    -t=*|--to=*)
    TO="${i#*=}"
    shift # past argument=value
    ;;
    -n=*|--Ne=*)
    NE="${i#*=}"
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

### STEP 1
COMM1="impute2 -prephase_g -m ${MAP} -g ${GENOTYPE} -int ${FROM} ${TO} -Ne ${NE} -o ${OUT}.prephasing.impute2"
echo "step 1: pre-phasing"
echo $COMM1
$COMM1

### STEP 2
COMM2="impute2 -use_prephased_g -m ${MAP} -h ${HAPLOTYPE} -l ${LEGEND} -known_haps_g ${OUT}.prephasing.impute2_haps -int ${FROM} ${TO} -Ne ${NE} -o ${OUT}.phased.impute2 -phase -align_by_maf_g"
echo "step2: imputation"
echo $COMM2
$COMM2

