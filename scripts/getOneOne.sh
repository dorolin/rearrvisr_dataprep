#!/bin/bash

## -------------------------------------------------------
## extract 1:1 orthologs for all species pairs and all
## orthologous groups defined in 'fullOGs.txt', which was
## generated with 'getOGs.R'

## Example: getOneOne.sh path/to/OMA/Output
## -------------------------------------------------------


mywd=$1

rm -f ${mywd}/fullOneOneOGfiles.txt
touch ${mywd}/fullOneOneOGfiles.txt

while read OG
do

  isoneone=0

  for file in ${mywd}/PairwiseOrthologs/*.txt
  do

    class=`grep -P "\t$OG$" $file | cut -f 5`

    if [ "$class" != "1:1" ]; then

      isoneone=1

    fi

  done

  if [ $isoneone -eq 0 ]; then

    echo "OG$OG.fa" >> ${mywd}/fullOneOneOGfiles.txt

  fi

done < ${mywd}/fullOGs.txt

