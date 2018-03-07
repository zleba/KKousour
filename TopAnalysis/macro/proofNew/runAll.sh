#!/bin/bash
N=10

add=$(pwd)

for per in {B..H}
do
    res=$(./GetEntries.py /nfs/dust/cms/user/zlebcr/JEC/ntuplesNewFormat/merged/jets${per}.root |& tail -1 |grep -o "[0-9]*$")
    ((N =  res / 3000000 + 1))
    for ((i=0; i < N; ++i))
    do
        echo ./run.sh $N $i $per
        qsub -e $add/out/${per}_${N}_${i}.err  -o $add/out/${per}_${N}_${i}.out -N  ${per}_${N}_${i} run.sh $N $i $per
    done
done
