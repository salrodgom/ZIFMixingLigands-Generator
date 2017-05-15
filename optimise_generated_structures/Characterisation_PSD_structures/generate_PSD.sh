#!/bin/bash
for CIFFile in ../CIFFilesFinished/*.cif ; do
 name=$(echo $CIFFile | sed 's/CIFFilesFinished\///g' | sed 's/\.cif//g' | sed 's/\.\.\///g')
 echo $name
 if [ ! -f ${name} ] ; then
  mkdir $name
  cd $name
   sed "s/STRUCTURE/$name/g" ../fff_raspa/INPUT > simulation.input
   cp ../fff_raspa/*.def .
   cp ../input_raspa.slurm .
   cp ../../CIFFilesFinished/${name}.cif . 
   sbatch -J ${name} input_raspa.slurm 
  cd ..
 fi
 sleep 0.1
done
