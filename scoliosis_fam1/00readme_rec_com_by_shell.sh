mv ../*/*.vcf .
### 2 Ref genotype filter
###### 7591 is hom REF
grep -v '#' 1505KHX-0015_7591_SNP_INDEL.vcf | awk '$10~"0/0"' | awk '{print $1,$2,$2}' | sort -u > 7591.txt
###### 7590 is het
grep -v '#' 1505KHX-0015_7590_SNP_INDEL.vcf | awk '$10!~"0/0" && $10~"0/"' | awk '{print $1,$2,$2}' | sort -u > 7590.txt
###### 7703 is het
grep -v '#' 1505KHX-0015_7703_SNP_INDEL.vcf | awk '$10!~"0/0" && $10~"0/"' | awk '{print $1,$2,$2}' | sort -u > 7703.txt
###### 7588 is hom REF or missing
grep -v '#' 1505KHX-0015_7588_SNP_INDEL.vcf | awk '$10~"0/0" || $10~/\./' | awk '{print $1,$2,$2}' | sort -u > 7588.txt
###### 7592 is hom REF or missing
grep -v '#' 1505KHX-0015_7592_SNP_INDEL.vcf | awk '$10~"0/0" || $10~/\./' | awk '{print $1,$2,$2}' | sort -u > 7592.txt
###### 7702 is hom REF or missing
grep -v '#' 1505KHX-0015_7702_SNP_INDEL.vcf | awk '$10~"0/0" || $10~/\./' | awk '{print $1,$2,$2}' | sort -u > 7702.txt
###### 7704B is hom REF or missing
grep -v '#' 1505KHX-0015_7704B_SNP_INDEL.vcf | awk '$10~"0/0" || $10~/\./' | awk '{print $1,$2,$2}' | sort -u > 7704B.txt
###### 7713 is hom REF or missing
grep -v '#' 1505KHX-0015_7713_SNP_INDEL.vcf | awk '$10~"0/0" || $10~/\./' | awk '{print $1,$2,$2}' | sort -u > 7713.txt
###### 7701 is not hom REF or missing
grep -v '#' 1505KHX-0015_7701_SNP_INDEL.vcf | awk '$10!~"0/0" || $10~/\./' | awk '{print $1,$2,$2}' | sort -u > 7701.txt
###### 7706 is not hom REF or missing
grep -v '#' 1505KHX-0015_7706_SNP_INDEL.vcf | awk '$10!~"0/0" || $10~/\./' | awk '{print $1,$2,$2}' | sort -u > 7706.txt
###### 7708 is not hom REF or missing
grep -v '#' 1505KHX-0015_7708_SNP_INDEL.vcf | awk '$10!~"0/0" || $10~/\./' | awk '{print $1,$2,$2}' | sort -u > 7708.txt
###### COMMON = 0 !!!!
### 4 Alt genotype filter
###### 7591 is hom Alt
grep -v '#' 1505KHX-0015_7591_SNP_INDEL.vcf | awk '$10!~"0/0" && $10!~"0/"' | awk '{print $1,$2,$2}' | sort -u > 7591.txt
###### 7590 is het
grep -v '#' 1505KHX-0015_7590_SNP_INDEL.vcf | awk '$10!~"0/0" && $10~"0/"' | awk '{print $1,$2,$2}' | sort -u > 7590.txt
###### 7703 is het
grep -v '#' 1505KHX-0015_7703_SNP_INDEL.vcf | awk '$10!~"0/0" && $10~"0/"' | awk '{print $1,$2,$2}' | sort -u > 7703.txt
###### 7588 is hom Alt or missing
grep -v '#' 1505KHX-0015_7588_SNP_INDEL.vcf | awk '$10!~"0/0" && $10!~"0/" || $10~/\./' | awk '{print $1,$2,$2}' | sort -u > 7588.txt
###### 7592 is hom Alt or missing
grep -v '#' 1505KHX-0015_7592_SNP_INDEL.vcf | awk '$10!~"0/0" && $10!~"0/" || $10~/\./' | awk '{print $1,$2,$2}' | sort -u > 7592.txt
###### 7702 is hom Alt or missing
grep -v '#' 1505KHX-0015_7702_SNP_INDEL.vcf | awk '$10!~"0/0" && $10!~"0/" || $10~/\./' | awk '{print $1,$2,$2}' | sort -u > 7702.txt
###### 7704B is hom Alt or missing
grep -v '#' 1505KHX-0015_7704B_SNP_INDEL.vcf | awk '$10!~"0/0" && $10!~"0/" || $10~/\./' | awk '{print $1,$2,$2}' | sort -u > 7704B.txt
###### 7713 is hom Alt or missing
grep -v '#' 1505KHX-0015_7713_SNP_INDEL.vcf | awk '$10!~"0/0" && $10!~"0/" || $10~/\./' | awk '{print $1,$2,$2}' | sort -u > 7713.txt
###### 7701 is not hom Alt or missing
grep -v '#' 1505KHX-0015_7701_SNP_INDEL.vcf | awk '$10~"0/0" || $10~"0/" || $10~/\./' | awk '{print $1,$2,$2}' | sort -u > 7701.txt
###### 7706 is not hom Alt or missing
grep -v '#' 1505KHX-0015_7706_SNP_INDEL.vcf | awk '$10~"0/0" || $10~"0/" || $10~/\./' | awk '{print $1,$2,$2}' | sort -u > 7706.txt
###### 7708 is not hom Alt or missing
grep -v '#' 1505KHX-0015_7708_SNP_INDEL.vcf | awk '$10~"0/0" || $10~"0/" || $10~/\./' | awk '{print $1,$2,$2}' | sort -u > 7708.txt
###### common
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

