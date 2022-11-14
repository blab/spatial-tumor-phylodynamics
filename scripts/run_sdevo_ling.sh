#!/bin/sh


#Ling et al
## Run in ling-application directory
## XMLs are created with write_ling_state_clocks.R
rsync -a ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/ling-application/xmls mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/ling-application

seed=1009
for file in xmls/hcc-wes_unidir_state_fixed_cr_rep[1-3].xml
do
  echo $file
  seed=$(($seed + 1))
  echo $seed
  sbatch submit_sdevo.sh $file $seed
done

seed=1929
for file in xmls/hcc-wes_unidir_state_fixed_strict_cr_rep[1-3].xml
do
  echo $file
  seed=$(($seed + 1))
  echo $seed
  sbatch submit_sdevo.sh $file $seed
done


#java -jar SDevo.jar xmls/hcc-wes_unidir_state_rep1.xml
#java -jar SDevo.jar -resume -seed 11 xmls/hcc-wes_unidir_state_rep2.xml
#java -jar SDevo.jar -resume -seed 210 xmls/hcc-wes_unidir_state_rep3.xml
#java -jar SDevo.jar -seed 44 xmls/hcc-wes_unidir_state_strict_clock_rep3.xml
#java -jar SDevo.jar -seed 33 xmls/hcc-wes_unidir_state_strict_clock_rep2.xml
#java -jar SDevo.jar xmls/hcc-wes_unidir_state_strict_clock_rep1.xml
