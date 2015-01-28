## usage: bash summary.sh

if [ -e HLAtype.txt ]
then
    rm HLAtype.txt
fi

ls ADHLATyping2014/*/*/*.tsv > file.txt
cp file.txt file2.txt
sed -i 's/\// /' file2.txt
sed -i 's/\// /' file2.txt
awk '{print $2}' file2.txt > file3.txt

IFS=$'\n' read -d '' -r -a lines < file.txt
IFS=$'\n' read -d '' -r -a ids < file3.txt

n=$(ls ADHLATyping2014/*/*/*.tsv | wc -l)    
let m=$n-1
for i in $(seq 0 $m)
do
    cp ${lines[$i]} temp
    sed -i '1d' temp
    echo ${ids[$i]} > n
    awk '{print $2,$3,$4,$5,$6,$7}' temp > g
    paste n g >> HLAtype.txt
    rm temp n g
done

rm file*.txt

