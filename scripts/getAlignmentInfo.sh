#!/bin/bash

## -------------------------------------------------------
## Example: getAlignmentInfo.sh path/to/MAFFT
## -------------------------------------------------------

mywd=$1

echo -e OG'\t'nSp'\t'len'\t'nGap'\t'pGap'\t'nX'\t'pX >${mywd}/alignInfo.txt

for file in ${mywd}/*.phy; do
  OG=`basename ${file%".phy"}`
  nSp=`cat ${file} | sed -n '1p' | cut -d ' ' -f1`
  len=`cat ${file} | sed -n '1p' | cut -d ' ' -f2`
  nAA=$(( $nSp * $len ))
  nGap=`cut -d ' ' -f2 ${file} | grep -o "-" | wc -w`
  nX=`cut -d ' ' -f2 ${file} | grep -oi "X" | wc -w`
  pGap=`echo "scale=4; $nGap/$nAA" | bc -l`
  pX=`echo "scale=4; $nX/$nAA" | bc -l`
  echo -e $OG'\t'$nSp'\t'$len'\t'$nGap'\t'$pGap'\t'$nX'\t'$pX >>${mywd}/alignInfo.txt
done


## -------------------------------------------------------
