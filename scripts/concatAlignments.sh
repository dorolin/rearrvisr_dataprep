#!/bin/bash

## -------------------------------------------------------
## Example: concatAlignments.sh path/to/MAFFT
## -------------------------------------------------------

mywd=$1

declare -a mysp=("ANA" "ERE" "GRI" "MEL" "MOJ" "PER" "PSE" "SEC" "SIM" "VIR" "WIL" "YAK")

touch tmp.phy
for sp in "${mysp[@]}"
do
  echo "$sp " >> tmp.phy
done

mylen=0

while read file; do
  len=`cat ${mywd}/${file} | sed -n '1p' | cut -d ' ' -f2`
  ids=`tail -n +2 ${mywd}/${file} | cut -d ' ' -f1 | cut -d '_' -f1 | sed 's/\n/ /g'`
  tmpsp=($ids)

  diff=$(diff <(printf "%s\n" "${mysp[@]}") <(printf "%s\n" "${tmpsp[@]}"))

  if [[ -z "$diff" ]]; then
      paste -d '' <(cat tmp.phy) <(cat ${mywd}/${file} | tail -n +2 | cut -d ' ' -f2) >tmpout
      mv -f tmpout tmp.phy
      mylen=$(( $mylen + $len ))
  else
    echo "${file} does not contain the correct species"
  fi
done < ${mywd}/filteredAlignments.txt

echo "${#mysp[@]} $mylen" >${mywd}/concatOGAlignment.phy
cat tmp.phy >>${mywd}/concatOGAlignment.phy
rm -f tmp.phy

