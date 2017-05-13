#!/bin/bash
for file in *mixture.txt ; do
 # GME_imi_0.80_bim_.20_mixture.txt
 CIFFile=$(echo $file | sed 's/\.txt/\.cif/g')
 topol=$(echo $file | sed 's/_/ /g' | awk '{print $1}')
 n1=$(echo $file | sed 's/_/ /g' | awk '{print $2}')
 n2=$(echo $file | sed 's/_/ /g' | awk '{print $4}')
 l1t=$(echo $file | sed 's/_/ /g' | awk '{print $3}')
 l2t=$(echo $file | sed 's/_/ /g' | awk '{print $5}')
 l1=$(grep 'molar fractions:' $file | tail -n1 | awk '{print $7}')
 l2=$(grep 'molar fractions:' $file | tail -n1 | awk '{print $8}')
 # echo $l1t $l2t $l1 $l2
 mv $file ${topol}_${n1}_${l1}_${n2}_${l2}_mixture.txt
 mv $CIFFile ${topol}_${n1}_${l1}_${n2}_${l2}_mixture.cif
done
