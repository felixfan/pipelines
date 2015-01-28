## usage: bash summary.sh

if [ -e HLAtype.txt ]
then
    rm HLAtype.txt
fi

#####################################################################
ls ./results/*.A.typing.txt > file.txt
sed -i 's/\// /g' file.txt
sed -i 's/\./ /g' file.txt
awk '{print $2}' file.txt > file2.txt
IFS=$'\n' read -d '' -r -a lines < file2.txt
n=$(ls ./results/*.A.typing.txt | wc -l)    
let m=$n-1
rm file*.txt
#####################################################################
for i in $(seq 0 $m)
do
    echo ${lines[$i]} > N
    ###########################################
    if [ -e results/${lines[$i]}.A.typing.txt ]
    then
        tail -n 1 results/${lines[$i]}.A.typing.txt | awk '{print $1,$2}' > A
    else
        echo "NA NA" > A
    fi
    
    if [ -e results/${lines[$i]}.B.typing.txt ]
    then
        tail -n 1 results/${lines[$i]}.B.typing.txt | awk '{print $1,$2}' > B
    else
        echo "NA NA" > B
    fi
    
    if [ -e results/${lines[$i]}.C.typing.txt ]
    then
        tail -n 1 results/${lines[$i]}.C.typing.txt | awk '{print $1,$2}' > C
    else
        echo "NA NA" > C
    fi
    
    if [ -e results/${lines[$i]}.DQB1.typing.txt ]
    then
        tail -n 1 results/${lines[$i]}.DQB1.typing.txt | awk '{print $1,$2}' > DQB
    else
        echo "NA NA" > DQB
    fi
    
    if [ -e results/${lines[$i]}.DRB1.typing.txt ]
    then
        tail -n 1 results/${lines[$i]}.DRB1.typing.txt | awk '{print $1,$2}' > DRB
    else
        echo "NA NA" > DRB
    fi
    ################rm blank/or only white sapce line#########################
    sed -i '/^[[:space:]]*$/d' A
    sed -i '/^[[:space:]]*$/d' B
    sed -i '/^[[:space:]]*$/d' C
    sed -i '/^[[:space:]]*$/d' DQB
    sed -i '/^[[:space:]]*$/d' DRB
    ### empty file
    if [ ! -s A ]
    then
        echo "NA NA" > A
    fi
    
    if [ ! -s B ]
    then
        echo "NA NA" > B
    fi
    
    if [ ! -s C ]
    then
        echo "NA NA" > C
    fi
    
    if [ ! -s DQB ]
    then
        echo "NA NA" > DQB
    fi
    
    if [ ! -s DRB ]
    then
        echo "NA NA" > DRB
    fi
    ####################################
    paste N A B C DQB DRB >> HLAtype.txt
    # rm N A B C DQB DRB
done
