#!/bin/bash

# This example file creates the directory structure and MESA input files for the SMSP section of the set. The initial donor mass (Msol) and log10(orbital period/days) values are set by mvals and pvals in the code.

mvals=($(seq 0.95 0.05 4.0))
pvals=($(seq -0.6 0.02 1.64))

for i in ${mvals[@]}
do
   mkdir "m_$i"
   echo "Mass $i"
   for j in ${pvals[@]}
   do
      m1="$(printf "%8.4f" "$i")"

      # You can change initial BH mass (Msol) by replacing "5.00" in this line.
      m2="$(printf "%8.4f" "5.00")"

      p="$(printf "%8.4f" "$(echo "scale=6;e($j*l(10.0))" | bc -l)")"
      mkdir "m_$i/p_$j"
      echo "Period $j"
      cp ref_inlist inlist
      cp ref_inlist1 inlist1
      cp ref_profile_columns.list profile_columns.list
      sed -i -e "s/init_m1/$m1/g" inlist
      sed -i -e "s/init_m2/$m2/g" inlist
      sed -i -e "s/init_p/$p/g" inlist
      mv inlist inlist1 profile_columns.list ./m_$i/p_$j/
      
      # Add any files that you want copied into all simulation directories here.
      cp rn ./m_$i/p_$j/
   done
done  
