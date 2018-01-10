
for run in B C D E F G H
do
    echo $run
    root -l -b -q "match.C(\"/nfs/dust/cms/user/zlebcr/JEC/ntuplesAug/merged/jets${run}.root\", -1)" &
done
