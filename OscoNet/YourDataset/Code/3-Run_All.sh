#!/bin/bash
root=../OsconetInput
outdir=../OsconetOutput
mkdir -p $outdir
log=run_all.log
experiment=YourDataset # change this to reflect the 'case' name set in script 1
for fn in ls $root/*$experiment*; do
 echo "file $fn"
 python collaboration.py $fn 128
 mv $root/*.out.csv $outdir/ 
 mv $root/*.psi_ng.csv $outdir/
done
