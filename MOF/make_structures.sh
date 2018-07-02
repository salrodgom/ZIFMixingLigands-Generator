#!/bin/bash
declare -A ligands
make
ligands[1]="li1"
ligands[2]="li2"
ligands[3]="li3"
ligands[4]="li4"
ligands[5]="li5"
ligands[6]="li6"
ligands[7]="li7"
ligands[8]="li8"
ligands[9]="li9"
ligands[10]="l10"
ligands[11]="l11"
ligands[12]="l12"
for topology in "I10" "CP2" "LP2" ; do
 for (( i=1; i<=12; i++ )) ; do
   linker_1=${ligands[${i}]}
   echo "$topology ${linker_1}"
   #./generator -t $topology -l 1 ${linker_1} 1.0 > ${topology}_${linker_1}_pure.txt
   #cp ${topology}_${linker_1}_mixture.cif ${topology}_${linker_1}_pure.cif
 done
 for (( i=1; i<=12; i++ )) ; do
  linker_1=${ligands[${i}]}
  #for (( j=i+1; j<=12; j++ )) ; do
   j=2
   linker_2=${ligands[${j}]}
   if [ "$linker_1" != "$linker_2" ] ; then 
    for molar_fraction in $(seq 0.10 0.10 0.90 ) ; do
     echo "$topology ${linker_1} ${molar_fraction} ${linker_2} $(bc <<<"scale=2;1.0-${molar_fraction}")"
 #   #./generator -t $topology -l 2 ${linker_1} ${molar_fraction} ${linker_2} $(bc <<<"scale=2;1.0-${molar_fraction}") > ${topology}_${linker_1}_${molar_fraction}_${linker_2}_$(bc <<<"scale=2;1.0-${molar_fraction}")_mixture.txt
 #   # mv ${topology}_${linker_1}${linker_2}_mixture.cif ${topology}_${linker_1}_${molar_fraction}_${linker_2}_$(bc <<<"scale=2;1.0-${molar_fraction}")_mixture.cif
     done
   fi
  #done
 done
done
