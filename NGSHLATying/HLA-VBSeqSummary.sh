## usage: bash summary.sh

if [ -e HLAtype.txt ]
then
    rm HLAtype.txt
fi

#####################################################################
ls ./result/*.result.txt > file.txt
sed -i 's/\// /g' file.txt
sed -i 's/\./ /g' file.txt
awk '{print $2}' file.txt > file2.txt
IFS=$'\n' read -d '' -r -a lines < file2.txt
n=$(ls ./result/*.result.txt | wc -l)    
let m=$n-1
rm file*.txt
#####################################################################
for i in $(seq 0 $m)
do
	echo ${lines[$i]} > temp

	for g in {'HLA_A','HLA_B','HLA_C','HLA_DQA1','HLA_DQB1','HLA_DRB1'}
	do
    	awk 'NR==1{print $1}' ./result/${lines[$i]}.$g.txt > $g.1
    	awk 'NR==2{print $1}' ./result/${lines[$i]}.$g.txt > $g.2
    done

    paste temp HLA_A.1 HLA_A.2 HLA_B.1 HLA_B.2 HLA_C.1 HLA_C.2 HLA_DQA1.1 HLA_DQA1.2 HLA_DQB1.1 HLA_DQB1.2 HLA_DRB1.1 HLA_DRB1.2 >> HLAtype.txt
    rm temp HLA_A.1 HLA_A.2 HLA_B.1 HLA_B.2 HLA_C.1 HLA_C.2 HLA_DQA1.1 HLA_DQA1.2 HLA_DQB1.1 HLA_DQB1.2 HLA_DRB1.1 HLA_DRB1.2
done
