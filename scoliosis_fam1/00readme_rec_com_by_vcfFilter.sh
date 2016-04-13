mv ../*/*.vcf .
### 2 Ref genotype filter
###### 7591 is hom REF
python ~/ubin/vcfFilter.py -vcf 1505KHX-0015_7588_SNP_INDEL.vcf -gtp hom-ref -ind 7588 -mv rm -o 7588.vcf
###### 7590 is het
python ~/ubin/vcfFilter.py -vcf 1505KHX-0015_7590_SNP_INDEL.vcf -gtp het -ind 7590 -mv rm -o 7590.vcf
###### 7703 is het
python ~/ubin/vcfFilter.py -vcf 1505KHX-0015_7703_SNP_INDEL.vcf -gtp het -ind 7703 -mv rm -o 7703.vcf
###### 7588 is hom REF or missing
python ~/ubin/vcfFilter.py -vcf 1505KHX-0015_7588_SNP_INDEL.vcf -gtp hom-ref -ind 7588 -mv keep -o 7588.vcf
###### 7592 is hom REF or missing
python ~/ubin/vcfFilter.py -vcf 1505KHX-0015_7592_SNP_INDEL.vcf -gtp hom-ref -ind 7592 -mv keep -o 7592.vcf
###### 7702 is hom REF or missing
python ~/ubin/vcfFilter.py -vcf 1505KHX-0015_7702_SNP_INDEL.vcf -gtp hom-ref -ind 7702 -mv keep -o 7702.vcf
###### 7704B is hom REF or missing
python ~/ubin/vcfFilter.py -vcf 1505KHX-0015_7704B_SNP_INDEL.vcf -gtp hom-ref -ind 7704B -mv keep -o 7704B.vcf
###### 7713 is hom REF or missing
python ~/ubin/vcfFilter.py -vcf 1505KHX-0015_7713_SNP_INDEL.vcf -gtp hom-ref -ind 7713 -mv keep -o 7713.vcf
###### 7701 is not hom REF or missing
python ~/ubin/vcfFilter.py -vcf 1505KHX-0015_7701_SNP_INDEL.vcf -gtp not-hom-ref -ind 7701 -mv keep -o 7701.vcf
###### 7706 is not hom REF or missing
python ~/ubin/vcfFilter.py -vcf 1505KHX-0015_7706_SNP_INDEL.vcf -gtp not-hom-ref -ind 7706 -mv keep -o 7706.vcf
###### 7708 is not hom REF or missing
python ~/ubin/vcfFilter.py -vcf 1505KHX-0015_7708_SNP_INDEL.vcf -gtp not-hom-ref -ind 7708 -mv keep -o 7708.vcf
###### COMMON = 0 !!!!
### 4 Alt genotype filter
###### 7591 is hom Alt
python ~/ubin/vcfFilter.py -vcf 1505KHX-0015_7591_SNP_INDEL.vcf -gtp hom-alt -ind 7591 -mv rm -o 7591.vcf
###### 7590 is het
python ~/ubin/vcfFilter.py -vcf 1505KHX-0015_7590_SNP_INDEL.vcf -gtp het -ind 7590 -mv rm -o 7590.vcf
###### 7703 is het
python ~/ubin/vcfFilter.py -vcf 1505KHX-0015_7703_SNP_INDEL.vcf -gtp het -ind 7703 -mv rm -o 7703.vcf
###### 7588 is hom Alt or missing
python ~/ubin/vcfFilter.py -vcf 1505KHX-0015_7588_SNP_INDEL.vcf -gtp hom-alt -ind 7588 -mv keep -o 7588.vcf
###### 7592 is hom Alt or missing
python ~/ubin/vcfFilter.py -vcf 1505KHX-0015_7592_SNP_INDEL.vcf -gtp hom-alt -ind 7592 -mv keep -o 7592.vcf
###### 7702 is hom Alt or missing
python ~/ubin/vcfFilter.py -vcf 1505KHX-0015_7702_SNP_INDEL.vcf -gtp hom-alt -ind 7702 -mv keep -o 7702.vcf
###### 7704B is hom Alt or missing
python ~/ubin/vcfFilter.py -vcf 1505KHX-0015_7704B_SNP_INDEL.vcf -gtp hom-alt -ind 7704B -mv keep -o 7704B.vcf
###### 7713 is hom Alt or missing
python ~/ubin/vcfFilter.py -vcf 1505KHX-0015_7713_SNP_INDEL.vcf -gtp hom-alt -ind 7713 -mv keep -o 7713.vcf
###### 7701 is not hom Alt or missing
python ~/ubin/vcfFilter.py -vcf 1505KHX-0015_7701_SNP_INDEL.vcf -gtp not-hom-alt -ind 7701 -mv keep -o 7701.vcf
###### 7706 is not hom Alt or missing
python ~/ubin/vcfFilter.py -vcf 1505KHX-0015_7706_SNP_INDEL.vcf -gtp not-hom-alt -ind 7706 -mv keep -o 7706.vcf
###### 7708 is not hom Alt or missing
python ~/ubin/vcfFilter.py -vcf 1505KHX-0015_7708_SNP_INDEL.vcf -gtp not-hom-alt -ind 7708 -mv keep -o 7708.vcf
###### common
for n in 7588 7590 7591 7592 7701 7702 7703 7704B 7706 7708 7713
do
    grep -v '#' $n.vcf | awk '{print $1,$2,$2}' | sort -u > $n.txt
done
comm -12 7588.txt 7590.txt > t1
comm -12 7591.txt 7592.txt > t2
comm -12 7701.txt 7702.txt > t3
comm -12 7703.txt 7704B.txt > t4
comm -12 7706.txt 7708.txt > t5
comm -12 t1 t2 > t12
comm -12 t3 t4 > t34
comm -12 t5 7713.txt > t56
comm -12 t12 t34 > t1234
comm -12 t1234 t56 > t123456
awk '$1!="chrX" && $1!="chrY" && $1!="chrM"' t123456 > t123456auto
echo -e "chr\tstart\tend" > t.txt
cat t.txt t123456auto > alt.comm.txt
rm 7*.txt
rm t*
###### extract common gtp
for n in 7588 7590 7591 7592 7701 7702 7703 7704B 7706 7708 7713
do
    vcftools --vcf "1505KHX-0015_"$n"_SNP_INDEL.vcf" --bed alt.comm.txt --recode --recode-INFO-all --out $n.alt
done
###### annotation
grep -v '#' 7591.alt.recode.vcf | awk '{print $1,$2,$2,$4,$5}'  > alt.annovar.input

