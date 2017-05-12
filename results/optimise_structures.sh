#!/bin/bash
# variables:
nCPU=16
# functions:
function count_used_CPUs {
 n_used=0
 for process in "simulate" "lmp_mpi" "gulp" ; do
  n=$(ps aux | grep ${process} | sed '/grep/d' | wc -l | awk '{print $1}')
  n_used=$((${n}+${n_used}))
 done
}
function get_seed {
 #choose a seed for coding the IO properly:
 seed=$(od --read-bytes=3 --address-radix=n --format=uL /dev/urandom | tr --delete " ")
 while [ $(echo "$seed > 900000000" | bc -l) == 1 ] || [ $(echo "$seed < 0" | bc -l) == 1 ] ; do
  seed=$(od --read-bytes=3 --address-radix=n --format=uL /dev/urandom | tr --delete " ")
  sleep 0.5
 done
}
function go_lammps {
 count_used_CPUs
 while [ $(echo "${n_used} >= ${nCPU}" | bc -l) == 1 ] ; do
  sleep 30
  count_used_CPUs
 done
 lammps
 sleep 1
}
function lammps {
 nohup mpirun --np $((${nCPU}-2)) lmp_mpi -in in.lmp -sf opt &
}
function wait_for_lammps {
 lammps_end=$( grep "Total wall time:" logs/log.main | wc -l | awk '{print $1}')
 while [ $(echo "${lammps_end} < 1" | bc -l) == 1 ] ; do
  sleep 30
  lammps_end=$( grep "Total wall time:" logs/log.main | wc -l | awk '{print $1}')
 done
}
function make_binaries {
 gfortran -O3 lammpstrj2pdb.f90 -o lammpstrj2pdb
 gfortran -O3 pdb2cif.f90 -o pdb2cif
 gfortran -O3 cif2lammps.f90 -o cif2lammps
}
function lammps_raspa {
 ./lammpstrj2pdb < movs/opti.lammpstrj
 n_lines=$(wc -l out.pdb |awk '{print $1}')
 line=$(sed -n '/MODEL/{=;p}' out.pdb | sed '{N;s/\n/ /}' | tail -n1 | awk '{print $1}')
 end=$((n_lines - line + 1))
 tail -n$end out.pdb > input.pdb
 ./pdb2cif
 mv p1.cif ${folder}_optimised.cif
}
# main section
make_binaries
for CIFFile in *.cif ; do
 get_seed
 folder=$(echo $CIFFile | sed 's/\.cif//g')
 if [ ! -d $folder ] ; then
  mkdir $folder
  cp forcefield.lib in.lmp $CIFFile $folder
  cp lammpstrj2pdb pdb2cif cif2lammps $folder/.
  cd $folder
   ./cif2lammps -c $CIFFile -wq -S > out_cif2lammps_${folder}.txt
   if [ $(echo "$(wc -l out_cif2lammps_*.txt | awk '{print $1}') > 50" | bc -l) == 1 ] ; then
    mv ${folder}.data ${folder}.lmp
    sed -i "s/RANDOMSEED/${seed}/g"   in.lmp
    sed -i "s/DATAFILE/${folder}.lmp/g"   in.lmp
    elements=$(cat atom_types_for_dump.txt | sed 's/[0-9]//g' | sed 's/  / /g')
    sed -i "s/ELEMENTS/${elements}/g" in.lmp
    go_lammps
    wait_for_lammps
    lammps_raspa
   else
    echo "cif2lammps ERROR"
   fi
  cd ..
 fi
done
