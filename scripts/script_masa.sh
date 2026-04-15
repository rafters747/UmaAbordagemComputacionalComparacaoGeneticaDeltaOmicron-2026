#antes de executar o script pela primeira vez, eh necessario
#executar esse comando: chmod +x script_masa.sh

#para executar o script basta colocar no prompt ./script_masa.sh

#remove todas as subpastas
rm -rf */
#!/bin/bash
n=4
k=1
n1=$((x=$n, y=1, x-y))
for i in $(seq 1 $n1)
do
    echo $i    
    ni=$((x=$i, y=1, x+y))
    for j in $(seq $ni $n)
    do    
        ../masa-serial --work-dir=$k seq$i.fasta seq$j.fasta
        k=$((x=$k, y=1, x+y))
    done
done
