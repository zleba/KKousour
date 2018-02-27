
for run in B  C D E F G H
do
    echo $run
    root -l -b -q "match.C(\"/nfs/dust/cms/user/zlebcr/JEC/ntuplesJCunc/merged/jets${run}.root\", \"AugErrAltM\",   -1)" &
    #root -l -b -q "match.C(\"/nfs/dust/cms/user/zlebcr/JEC/ntuplesJCuncSpring/merged/jets${run}.root\", \"SpringAltM\",   -1)" &
done
wait
