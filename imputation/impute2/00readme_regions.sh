### step 1 convert ped to gen
sh pre_ped2gen.sh -i=sco.hg19 -o=sco
sh pre_ped2gen.sh -i=ddd.hg19 -o=ddd
### step 2 imputation
map1="/home/felixfan/ulib/impute2ref/genetic_map_chr1_combined_b37.txt"
legend1="/home/felixfan/ulib/impute2ref/1000GP_Phase3_chr1.legend.gz"
hap1="/home/felixfan/ulib/impute2ref/1000GP_Phase3_chr1.hap.gz"
gen1="/home/felixfan/uwork/ais_gwas_2016/imputation/impute2/test/sco.chr1.gen"
gen1d="/home/felixfan/uwork/ais_gwas_2016/imputation/impute2/test/ddd.chr1.gen"
sh impute2_two_steps.sh -m=${map1} -g=${gen1} -f=224000000 -t=229000000 -n=20000 -o=sco.chr1.r1 -h=${hap1} -l=${legend1}
sh impute2_two_steps.sh -m=${map1} -g=${gen1d} -f=224000000 -t=229000000 -n=20000 -o=ddd.chr1.r1 -h=${hap1} -l=${legend1}
map2="/home/felixfan/ulib/impute2ref/genetic_map_chr2_combined_b37.txt"
legend2="/home/felixfan/ulib/impute2ref/1000GP_Phase3_chr2.legend.gz"
hap2="/home/felixfan/ulib/impute2ref/1000GP_Phase3_chr2.hap.gz"
gen2="/home/felixfan/uwork/ais_gwas_2016/imputation/impute2/test/sco.chr2.gen"
gen2d="/home/felixfan/uwork/ais_gwas_2016/imputation/impute2/test/ddd.chr2.gen"
sh impute2_two_steps.sh -m=${map2} -g=${gen2} -f=40000000 -t=45000000 -n=20000 -o=sco.chr2.r1 -h=${hap2} -l=${legend2}
sh impute2_two_steps.sh -m=${map2} -g=${gen2d} -f=40000000 -t=45000000 -n=20000 -o=ddd.chr2.r1 -h=${hap2} -l=${legend2}
map9="/home/felixfan/ulib/impute2ref/genetic_map_chr9_combined_b37.txt"
legend9="/home/felixfan/ulib/impute2ref/1000GP_Phase3_chr9.legend.gz"
hap9="/home/felixfan/ulib/impute2ref/1000GP_Phase3_chr9.hap.gz"
gen9="/home/felixfan/uwork/ais_gwas_2016/imputation/impute2/test/sco.chr9.gen"
gen9d="/home/felixfan/uwork/ais_gwas_2016/imputation/impute2/test/ddd.chr9.gen"
sh impute2_two_steps.sh -m=${map9} -g=${gen9d} -f=14000000 -t=19000000 -n=20000 -o=ddd.chr9.r1 -h=${hap9} -l=${legend9}
sh impute2_two_steps.sh -m=${map9} -g=${gen9} -f=14000000 -t=19000000 -n=20000 -o=sco.chr9.r1 -h=${hap9} -l=${legend9}
map20="/home/felixfan/ulib/impute2ref/genetic_map_chr20_combined_b37.txt"
legend20="/home/felixfan/ulib/impute2ref/1000GP_Phase3_chr20.legend.gz"
hap20="/home/felixfan/ulib/impute2ref/1000GP_Phase3_chr20.hap.gz"
gen20="/home/felixfan/uwork/ais_gwas_2016/imputation/impute2/test/sco.chr20.gen"
gen20d="/home/felixfan/uwork/ais_gwas_2016/imputation/impute2/test/ddd.chr20.gen"
sh impute2_two_steps.sh -m=${map20} -g=${gen20d} -f=42000000 -t=48000000 -n=20000 -o=ddd.chr20.r1 -h=${hap20} -l=${legend20}
sh impute2_two_steps.sh -m=${map20} -g=${gen20} -f=42000000 -t=48000000 -n=20000 -o=sco.chr20.r1 -h=${hap20} -l=${legend20}
#### step 3 convert gen to ped (SNPs only)
sh post_gen2bed.sh -g=sco.chr1.r1.phased.impute2 -s=sco.chr1.sample -c=1 -t=0.9 -p=phenotype -x=sex -o=sco.chr1
sh post_gen2bed.sh -g=sco.chr2.r1.phased.impute2 -s=sco.chr2.sample -c=2 -t=0.9 -p=phenotype -x=sex -o=sco.chr2
sh post_gen2bed.sh -g=sco.chr9.r1.phased.impute2 -s=sco.chr9.sample -c=9 -t=0.9 -p=phenotype -x=sex -o=sco.chr9
sh post_gen2bed.sh -g=sco.chr20.r1.phased.impute2 -s=sco.chr20.sample -c=20 -t=0.9 -p=phenotype -x=sex -o=sco.chr20
sh post_gen2bed.sh -g=ddd.chr1.r1.phased.impute2 -s=ddd.chr1.sample -c=1 -t=0.9 -p=phenotype -x=sex -o=ddd.chr1
sh post_gen2bed.sh -g=ddd.chr2.r1.phased.impute2 -s=ddd.chr2.sample -c=2 -t=0.9 -p=phenotype -x=sex -o=ddd.chr2
sh post_gen2bed.sh -g=ddd.chr9.r1.phased.impute2 -s=ddd.chr9.sample -c=9 -t=0.9 -p=phenotype -x=sex -o=ddd.chr9
sh post_gen2bed.sh -g=ddd.chr20.r1.phased.impute2 -s=ddd.chr20.sample -c=20 -t=0.9 -p=phenotype -x=sex -o=ddd.chr20
#### step 4 merge
echo -e "sco.chr2.r1.bed\tsco.chr2.r1.bim\tsco.chr2.r1.fam"> files.txt
echo -e "sco.chr9.r1.bed\tsco.chr9.r1.bim\tsco.chr9.r1.fam">> files.txt
echo -e "sco.chr20.r1.bed\tsco.chr20.r1.bim\tsco.chr20.r1.fam">> files.txt
echo -e "ddd.chr1.r1.bed\tddd.chr1.r1.bim\tddd.chr1.r1.fam">> files.txt
echo -e "ddd.chr2.r1.bed\tddd.chr2.r1.bim\tddd.chr2.r1.fam"> files.txt
echo -e "ddd.chr9.r1.bed\tddd.chr9.r1.bim\tddd.chr9.r1.fam">> files.txt
echo -e "ddd.chr20.r1.bed\tddd.chr20.r1.bim\tddd.chr20.r1.fam">> files.txt
plink --bfile sco.chr1.r1 --merge-list files.txt --make-bed --out ais_impute_regions
#### step 5 update id
cat ais_impute_regions.fam | tr "_" " " | awk '{print $1,$1}' > new.fam
awk '{print $1,$2}' ais_impute_regions.fam > old.fam
paste old.fam new.fam > update.id.txt
plink --bfile ais_impute_regions --update-ids update.id.txt --make-bed --out ais_impute_regions_upID
#### step 6 extract inds
plink --bfile ais_impute_regions_upID --keep ais_inds2.txt --make-bed --out ais_impute2_final
plink --bfile ais_impute2_final --filter-females --make-bed --out ais_impute2_final_female
### step 6 assoc
plink --bfile ais_impute2_final --assoc --ci 0.95 --out ais_impute2_final_assoc
plink --bfile ais_impute2_final --model --out ais_impute2_final_model
plink --bfile ais_impute2_final_female --model --out ais_impute2_final_female_model
plink --bfile ais_impute2_final_female --assoc --ci 0.95 --out ais_impute2_final_female_assoc
### extract SNP infor
###==================================================================
grep rs10511625 sco_ddd_auto_qc2_update_female_model.model
grep rs4961736 sco_ddd_auto_qc2_update_female_model.model
grep rs1888207 sco_ddd_auto_qc2_update_female_model.model
grep rs10738446 sco_ddd_auto_qc2_update_female_model.model
grep rs10962542 sco_ddd_auto_qc2_update_female_model.model
grep rs10810593 sco_ddd_auto_qc2_update_female_model.model
grep rs1475711 sco_ddd_auto_qc2_update_female_model.model
grep rs59404763 sco_ddd_auto_qc2_update_female_model.model
grep rs10196262 sco_ddd_auto_qc2_update_female_model.model
grep rs1475712 sco_ddd_auto_qc2_update_female_model.model
grep rs3904778 sco_ddd_auto_qc2_update_female_model.model
grep rs10738445 sco_ddd_auto_qc2_update_female_model.model
###==================================================================
grep rs10511625 sco_ddd_auto_qc2_update_female_assoc.assoc
grep rs4961736 sco_ddd_auto_qc2_update_female_assoc.assoc
grep rs1888207 sco_ddd_auto_qc2_update_female_assoc.assoc
grep rs10738446 sco_ddd_auto_qc2_update_female_assoc.assoc
grep rs10962542 sco_ddd_auto_qc2_update_female_assoc.assoc
grep rs10810593 sco_ddd_auto_qc2_update_female_assoc.assoc
grep rs1475711 sco_ddd_auto_qc2_update_female_assoc.assoc
grep rs59404763 sco_ddd_auto_qc2_update_female_assoc.assoc
grep rs10196262 sco_ddd_auto_qc2_update_female_assoc.assoc
grep rs1475712 sco_ddd_auto_qc2_update_female_assoc.assoc
grep rs3904778 sco_ddd_auto_qc2_update_female_assoc.assoc
grep rs10738445 sco_ddd_auto_qc2_update_female_assoc.assoc
###==================================================================
grep rs10511625 sco_ddd_auto_qc2_update_model.model
grep rs4961736 sco_ddd_auto_qc2_update_model.model
grep rs1888207 sco_ddd_auto_qc2_update_model.model
grep rs10738446 sco_ddd_auto_qc2_update_model.model
grep rs10962542 sco_ddd_auto_qc2_update_model.model
grep rs10810593 sco_ddd_auto_qc2_update_model.model
grep rs1475711 sco_ddd_auto_qc2_update_model.model
grep rs59404763 sco_ddd_auto_qc2_update_model.model
grep rs10196262 sco_ddd_auto_qc2_update_model.model
grep rs1475712 sco_ddd_auto_qc2_update_model.model
grep rs3904778 sco_ddd_auto_qc2_update_model.model
grep rs10738445 sco_ddd_auto_qc2_update_model.model
###==================================================================
grep rs10511625 sco_ddd_assoc2.assoc
grep rs4961736 sco_ddd_assoc2.assoc
grep rs1888207 sco_ddd_assoc2.assoc
grep rs10738446 sco_ddd_assoc2.assoc
grep rs10962542 sco_ddd_assoc2.assoc
grep rs10810593 sco_ddd_assoc2.assoc
grep rs1475711 sco_ddd_assoc2.assoc
grep rs59404763 sco_ddd_assoc2.assoc
grep rs10196262 sco_ddd_assoc2.assoc
grep rs1475712 sco_ddd_assoc2.assoc
grep rs3904778 sco_ddd_assoc2.assoc
grep rs10738445 sco_ddd_assoc2.assoc
###==================================================================
grep rs10511625 ais_impute2_final_model.model
grep rs4961736 ais_impute2_final_model.model
grep rs1888207 ais_impute2_final_model.model
grep rs10738446 ais_impute2_final_model.model
grep rs10962542 ais_impute2_final_model.model
grep rs10810593 ais_impute2_final_model.model
grep rs1475711 ais_impute2_final_model.model
grep rs59404763 ais_impute2_final_model.model
grep rs10196262 ais_impute2_final_model.model
grep rs1475712 ais_impute2_final_model.model
grep rs3904778 ais_impute2_final_model.model
grep rs10738445 ais_impute2_final_model.model
###==================================================================
grep rs10511625 ais_impute2_final_assoc.assoc
grep rs4961736 ais_impute2_final_assoc.assoc
grep rs1888207 ais_impute2_final_assoc.assoc
grep rs10738446 ais_impute2_final_assoc.assoc
grep rs10962542 ais_impute2_final_assoc.assoc
grep rs10810593 ais_impute2_final_assoc.assoc
grep rs1475711 ais_impute2_final_assoc.assoc
grep rs59404763 ais_impute2_final_assoc.assoc
grep rs10196262 ais_impute2_final_assoc.assoc
grep rs1475712 ais_impute2_final_assoc.assoc
grep rs3904778 ais_impute2_final_assoc.assoc
grep rs10738445 ais_impute2_final_assoc.assoc
###==================================================================
grep rs10511625 ais_impute2_final_female_model.model
grep rs4961736 ais_impute2_final_female_model.model
grep rs1888207 ais_impute2_final_female_model.model
grep rs10738446 ais_impute2_final_female_model.model
grep rs10962542 ais_impute2_final_female_model.model
grep rs10810593 ais_impute2_final_female_model.model
grep rs1475711 ais_impute2_final_female_model.model
grep rs59404763 ais_impute2_final_female_model.model
grep rs10196262 ais_impute2_final_female_model.model
grep rs1475712 ais_impute2_final_female_model.model
grep rs3904778 ais_impute2_final_female_model.model
grep rs10738445 ais_impute2_final_female_model.model
###==================================================================
grep rs10511625 ais_impute2_final_female_assoc.assoc
grep rs4961736 ais_impute2_final_female_assoc.assoc
grep rs1888207 ais_impute2_final_female_assoc.assoc
grep rs10738446 ais_impute2_final_female_assoc.assoc
grep rs10962542 ais_impute2_final_female_assoc.assoc
grep rs10810593 ais_impute2_final_female_assoc.assoc
grep rs1475711 ais_impute2_final_female_assoc.assoc
grep rs59404763 ais_impute2_final_female_assoc.assoc
grep rs10196262 ais_impute2_final_female_assoc.assoc
grep rs1475712 ais_impute2_final_female_assoc.assoc
grep rs3904778 ais_impute2_final_female_assoc.assoc
grep rs10738445 ais_impute2_final_female_assoc.assoc

