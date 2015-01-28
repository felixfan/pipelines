## usage: bash summary.sh

if [ -e HLAtype.txt ]
then
    rm HLAtype.txt
fi

ls ./results/*.sum > file.txt
IFS=$'\n' read -d '' -r -a lines < file.txt
cp file.txt file2.txt
sed -i 's/\// /g' file2.txt
sed -i 's/\_/ /g' file2.txt
awk '{print $3}' file2.txt > file3.txt
IFS=$'\n' read -d '' -r -a ids < file3.txt
n=$(ls ./results/*.sum | wc -l)    
let m=$n-1
for i in $(seq 0 $m)
do
    cp ${lines[$i]} temp
    sed -i '1d' temp
    echo ${ids[$i]} > n
    grep "HLA_A" temp | awk '{print $2}' > a1
    grep "HLA_B" temp | awk '{print $2}' > b1
    grep "HLA_C" temp | awk '{print $2}' > c1
    grep "HLA_DQA1" temp | awk '{print $2}' > dqa1
    grep "HLA_DQB1" temp | awk '{print $2}' > dqb1
    grep "HLA_DRB1" temp | awk '{print $2}' > drb1
    grep "HLA_A" temp | awk '{print $3}' > a2
    grep "HLA_B" temp | awk '{print $3}' > b2
    grep "HLA_C" temp | awk '{print $3}' > c2
    grep "HLA_DQA1" temp | awk '{print $3}' > dqa2
    grep "HLA_DQB1" temp | awk '{print $3}' > dqb2
    grep "HLA_DRB1" temp | awk '{print $3}' > drb2
    # check empty or not
    if [ ! -s a1 ]
    then
        echo "NA" > a1
    fi
    if [ ! -s a2 ]
    then
        echo "NA" > a2
    fi
    if [ ! -s b1 ]
    then
        echo "NA" > b1
    fi
    if [ ! -s b2 ]
    then
        echo "NA" > b2
    fi
    if [ ! -s c1 ]
    then
        echo "NA" > c1
    fi
    if [ ! -s c2 ]
    then
        echo "NA" > c2
    fi
    if [ ! -s dqa1 ]
    then
        echo "NA" > dqa1
    fi
    if [ ! -s dqa2 ]
    then
        echo "NA" > dqa2
    fi
    if [ ! -s dqb1 ]
    then
        echo "NA" > dqb1
    fi
    if [ ! -s dqb2 ]
    then
        echo "NA" > dqb2
    fi
    if [ ! -s drb1 ]
    then
        echo "NA" > drb1
    fi
    if [ ! -s drb2 ]
    then
        echo "NA" > drb2
    fi
    # paste file
    paste n a1 a2 b1 b2 c1 c2 dqa1 dqa2 dqb1 dqb2 drb1 drb2 >> HLAtype.txt
    rm n a1 a2 b1 b2 c1 c2 dqa1 dqa2 dqb1 dqb2 drb1 drb2 temp
done

rm file*.txt

# check ambigous
awk '$2 ~/,/ {$2="NA"} 1' HLAtype.txt > temp.txt
awk '$3 ~/,/ {$3="NA"} 1' temp.txt > temp1.txt
awk '$4 ~/,/ {$4="NA"} 1' temp1.txt > temp.txt
awk '$5 ~/,/ {$5="NA"} 1' temp.txt > temp1.txt
awk '$6 ~/,/ {$6="NA"} 1' temp1.txt > temp.txt
awk '$7 ~/,/ {$7="NA"} 1' temp.txt > temp1.txt
awk '$8 ~/,/ {$8="NA"} 1' temp1.txt > temp.txt
awk '$9 ~/,/ {$9="NA"} 1' temp.txt > temp1.txt
awk '$10 ~/,/ {$10="NA"} 1' temp1.txt > temp.txt
awk '$11 ~/,/ {$11="NA"} 1' temp.txt > temp1.txt
awk '$12 ~/,/ {$12="NA"} 1' temp1.txt > temp.txt
awk '$13 ~/,/ {$13="NA"} 1' temp.txt > HLAtype.txt

rm temp*.txt
