#antes de executar o script pela primeira vez, eh necessario
#executar esse comando: chmod +x script_coleta.sh

#!/bin/bash
#remove todas as subpastas e o arquivo masa.csv
rm -rf */
rm masa.csv
n=4
k=1
n1=$((x=$n, y=1, x-y))
ARQUIVO="alignment.00.txt"
for i in $(seq 1 $n1)
do
    ni=$((x=$i, y=1, x+y))
    for j in $(seq $ni $n)
    do    
        ../masa-serial --work-dir=$k seq$i.fasta seq$j.fasta
        k=$((x=$k, y=1, x+y))
    done
	cd $k
    coleta=$(awk '/Total Score/ {print $3}' $ARQUIVO)
    coleta+=","    
    coleta+=$(awk '/Matches/ {print $2}' $ARQUIVO)
    coleta+=","
    coleta+=$(awk '/Mismatches/ {print $2}' $ARQUIVO)
    coleta+=","
    coleta+=$(awk '/Gap Openings/ {print $3}' $ARQUIVO)
    coleta+=","
    coleta+=$(awk '/Gap Extentions/ {print $3}' $ARQUIVO)
    echo $coleta >> ../masa.csv
    cd ..
	rm -rf */
done