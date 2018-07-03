#!/bin/bash
# main code:
for topology in SOD SZ7 ; do
 macro_folder=Pure
 for linker_1 in imi mim bim ; do
   folder=${topology}_${linker_1}_pure
   cd ${macro_folder}/$folder
    CIFFIle=$(ls *.cif)
    for CIFFile in *.cif ; do
     mv $CIFFile ${folder}.cif
    done
   cd ../..
 done
 macro_folder=Binary
 if [ ! -d $macro_folder ] ; then mkdir $macro_folder ; fi
 for linker_1 in imi mim ; do
  for linker_2 in mim bim ; do
   if [ "${linker_1}" != "${linker_2}" ] ; then
    for molar_fraction in $(seq 0.10 0.10 0.90 ) ; do
     folder=${topology}_${linker_1}_${molar_fraction}_${linker_2}_$(bc <<<"scale=2;1.0-${molar_fraction}")_binary
     cd ${macro_folder}/$folder
      for CIFFile in *.cif ; do
       echo $CIFFile
       mv $CIFFile ${folder}.cif
      done
     cd ../..
    done
   fi
  done
 done
 macro_folder=Ternary
 if [ ! -d $macro_folder ] ; then mkdir $macro_folder ; fi
 link_1="imi"
 link_2="mim"
 link_3="bim"
 for i in $(seq 0.2 0.2 0.8 ) ; do
  for j in $(seq 0.2 0.2 0.8 ) ; do
   for k in $(seq 0.2 0.2 0.8 ) ; do
    sum=$(bc <<<"scale=5;$i + $j + $k")
    iu=$(bc <<<"scale=4;$i / $sum")
    ju=$(bc <<<"scale=4;$j / $sum")
    ku=$(bc <<<"scale=4;$k / $sum")
    folder=${topology}_${link_1}_${iu}_${link_2}_${ju}_${link_3}_${ku}_ternary
    if [ ! -d $macro_folder/$folder ] ; then mkdir $macro_folder/$folder ; fi
    cd ${macro_folder}/$folder
     for CIFFile in *.cif ; do
      mv $CIFFile ${folder}.cif
     done
    cd ../..
   done
  done
 done
done
