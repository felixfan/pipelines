### read all individuals passed QC
mapfile -t inds < /home/felixfan/uwork/ADExomeHLA2016/infor/adExonQCinds.txt
### read all individuals have processed
ls /home/felixfan/uwork/ADExomeHLA2016/HLA-VBSeq_results/*.result.txt > tmp999.txt
sed -i 's/\// /g' tmp999.txt
sed -i 's/\./ /g' tmp999.txt
awk '{print $6}' tmp999.txt > tmp666.txt
mapfile -t finished < tmp666.txt
rm tmp666.txt tmp999.txt
### passed QC & not processed
for i in "${inds[@]}"
do
    if [[ ! ${finished[*]} =~ $i ]]; then
        fq1=$i"_1.fastq.gz"
        fq2=$i"_2.fastq.gz"
        if [ -f "$fq1" ] && [ -f "$fq2" ] && [ ! $i = '1094' ]; then
            bash HLA_VBSeq_Fastq.sh $i $fq1 $fq2
            # echo $i
        fi
    fi
done
