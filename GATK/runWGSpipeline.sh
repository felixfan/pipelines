echo sh wgsPipeline0.sh
n=1
for i in 7588 7590 7591 7592 7701 7702 7703 7704B 7706 7708 7713
do
    echo "sh wgsPipeline1.sh 8" "'@RG\tID:AIS"$i"\tSM:"$i"\tPL:ILLUMINA\tLB:lib"$i"\tPU:H52M2CCXX:"$n":none'" $i"_R1.fastq.gz" $i"_R2.fastq.gz" $i "True True False False"
    n=$(($n+1))
done
echo sh wgsPipeline2.sh 8 32 ais
echo sh wgsPipeline3.sh 8 32 ais.raw.vcf ais

