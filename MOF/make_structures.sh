#!/bin/bash
function look_for_jam {
 n_max_proc=16
 n_proc=$(ps aux | grep "generator" | sed "/grep/d" | wc -l | awk '{print $1}')
 while [ $(echo "${n_proc} >= ${n_max_proc}" | bc -lq ) == 1 ] ; do
  sleep 30
  n_proc=$(ps aux | grep "generator" | sed "/grep/d" | wc -l | awk '{print $1}')
 done
 sleep 0.1
}
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
# main code:
make
for topology in "I10" "CP2" "LP2" ; do
 macro_folder=Pure
 if [ ! -d $macro_folder ] ; then mkdir $macro_folder ; fi
 for (( i=1; i<=12; i++ )) ; do
   linker_1=${ligands[${i}]}
   folder=${topology}_${linker_1}_pure
   if [ ! -d $macro_folder/$folder ] ; then mkdir $macro_folder/$folder ; fi
   cd ${macro_folder}/$folder
    look_for_jam
    echo "[send] $folder"
    ln -s ../../generator
    ln -s ../../zif_${topology}_cif_gin_all
    nohup ./generator -t $topology -l 1 ${linker_1} 1.0 > ${folder}.txt &
   cd ../..
 done
 macro_folder=Binary
 if [ ! -d $macro_folder ] ; then mkdir $macro_folder ; fi
 for (( i=1; i<=12; i++ )) ; do
  linker_1=${ligands[${i}]}
   j=2
   linker_2=${ligands[${j}]}
   if [ "$linker_1" != "$linker_2" ] ; then 
    for molar_fraction in $(seq 0.10 0.10 0.90 ) ; do
     folder=${topology}_${linker_1}_${molar_fraction}_${linker_2}_$(bc <<<"scale=2;1.0-${molar_fraction}")_binary
     if [ ! -d $macro_folder/$folder ] ; then mkdir $macro_folder/$folder ; fi
     cd ${macro_folder}/$folder
      look_for_jam
      echo "[send] $folder"
      ln -s ../../generator
      ln -s ../../zif_${topology}_cif_gin_all
      nohup ./generator -t $topology -l 2 ${linker_1} ${molar_fraction} ${linker_2} $(bc <<<"scale=2;1.0-${molar_fraction}") > ${folder}.txt &
     cd ../..
    done
   fi
 done
done
