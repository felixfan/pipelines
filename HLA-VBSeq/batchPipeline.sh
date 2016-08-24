### read all individuals passed QC
mapfile -t inds < /home/felixfan/uwork/ADExomeHLA2016/infor/adExonQCinds.txt
###
for i in "${inds[@]}"
do
    fq1=$i"_1.fastq.gz"
    fq2=$i"_2.fastq.gz"
    if [ -f "$fq1" ] && [ -f "$fq2" ] && [ ! $i = '1094' ]; then
        bash HLA_VBSeq_Fastq.sh $i $fq1 $fq2
    fi
done
