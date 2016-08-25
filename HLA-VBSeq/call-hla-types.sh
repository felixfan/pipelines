ls *.result.txt > filelist.txt
sed -i 's/\./ /g' filelist.txt 
awk '{print $1}' filelist.txt > iid.txt
mapfile -t iids < iid.txt
touch hla_types.txt
for iid in ${iids[@]}
do
    python HLA-VBSeq-call-types.py $iid 4 >> hla_types.txt
done
