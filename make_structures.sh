#!/bin/bash
make
for topology in CAN GME AFI ; do
 for linker_1 in imi mim bim; do
   ./generator -t $topology -l 1 ${linker_1} 1.0 > ${topology}_${linker_1}_pure.txt
   cp ${topology}_${linker_1}_mixture.cif ${topology}_${linker_1}_pure.cif
 done
 for linker_1 in imi ; do
  for linker_2 in mim bim ; do
   for molar_fraction in $(seq 0.10 0.10 0.90 ) ; do
    ./generator -t $topology -l 2 ${linker_1} ${molar_fraction} ${linker_2} $(bc <<<"scale=2;1.0-${molar_fraction}") > ${topology}_${linker_1}_${molar_fraction}_${linker_2}_$(bc <<<"scale=2;1.0-${molar_fraction}")_mixture.txt
    cp ${topology}_${linker_1}${linker_2}_mixture.cif ${topology}_${linker_1}_${molar_fraction}_${linker_2}_$(bc <<<"scale=2;1.0-${molar_fraction}")_mixture.cif
   done
  done
 done
done
