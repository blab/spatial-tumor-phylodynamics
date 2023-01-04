#!/bin/sh
##script to find combine tree files
rsync -a ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/li-application/xmls mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/li-application
rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/li-application/*.log ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/li-application/out
rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/li-application/*.trees ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/li-application/out
#Run SDevo
seed=10
for file in xmls/T[1-2]_wgs_oristates_unidir_[0-9]_strict_rep[0-9].xml
do
  echo $file
  seed=$(($seed + 1))
  echo $seed
  sbatch submit_sdevo.sh $file $seed
done

seed=100
for file in xmls/T[1-2]_wgs_oristates_unidir_[0-9]_state_rep[0-9]_rt1l13.xml
do
  echo $file
  seed=$(($seed + 1))
  echo $seed
  sbatch submit_sdevo.sh $file $seed
done

seed=210
for file in xmls/T[1-2]_wgs_oristates_unidir_[0-9]_state_rep[0-9].xml
do
  echo $file
  seed=$(($seed + 1))
  echo $seed
  sbatch submit_sdevo.sh $file $seed
done

seed=301
for file in xmls/T[1-2]_wgs_newstates_unidir_[0-9]_state_rep[0-9].xml
do
  echo $file
  seed=$(($seed + 1))
  echo $seed
  sbatch submit_sdevo.sh $file $seed
done

#Tumor 1
## New states / unidirectional
/Applications/BEAST\ 2.6.2/bin/LogCombiner -log T1_wgs_newstates_unidir_1_state_rep[0-9].HCCtumor.typed.node.trees \
-o T1_wgs_newstates_unidir_state_comb.HCCtumor.typed.node.trees \
-b 20

#/Applications/BEAST\ 2.6.2/bin/LogCombiner -log T1_wgs_newstates_unidir_[0-9]_state_rep[0-9].log \
#-o T1_wgs_newstates_unidir_state_comb.log \
#-b 10

## Old states / unidirectional
/Applications/BEAST\ 2.6.2/bin/LogCombiner -log T1_wgs_oristates_unidir_1_state_rep[0-9].HCCtumor.typed.node.trees \
-o T1_wgs_oristates_unidir_state_comb.HCCtumor.typed.node.trees \
-b 20

#/Applications/BEAST\ 2.6.2/bin/LogCombiner -log T1_wgs_oristates_unidir_[0-9]_state_rep[0-9].log \
#-o T1_wgs_oristates_unidir_state_comb.log \
#-b 10

#Tumor 2
## New states / unidirectional
/Applications/BEAST\ 2.6.2/bin/LogCombiner -log T2_wgs_newstates_unidir_1_state_rep[0-9].HCCtumor.typed.node.trees \
-o T2_wgs_newstates_unidir_state_comb.HCCtumor.typed.node.trees \
-b 20

#/Applications/BEAST\ 2.6.2/bin/LogCombiner -log T2_wgs_newstates_unidir_[0-9]_state_rep[0-9].log \
#-o T2_wgs_newstates_unidir_state_comb.log \
#-b 10

## Old states / unidirectional
/Applications/BEAST\ 2.6.2/bin/LogCombiner -log T2_wgs_oristates_unidir_1_state_rep[023].HCCtumor.typed.node.trees \
-o T2_wgs_oristates_unidir_state_comb.HCCtumor.typed.node.trees \
-b 20

#/Applications/BEAST\ 2.6.2/bin/LogCombiner -log T2_wgs_oristates_unidir_[0-9]_state_rep[023].log \
#-o T2_wgs_oristates_unidir_state_comb.log \
#-b 10

#Strict clock

#Tumor 1
/Applications/BEAST\ 2.6.2/bin/LogCombiner -log T1_wgs_oristates_unidir_1_strict_rep[0-9].HCCtumor.typed.node.trees \
-o T1_wgs_oristates_unidir_strict_comb.HCCtumor.typed.node.trees \
-b 20

#/Applications/BEAST\ 2.6.2/bin/LogCombiner -log T1_wgs_oristates_unidir_[0-9]_strict_rep[0-9].log \
#-o T1_wgs_oristates_unidir_strict_comb.log \
#-b 10

#Tumor 2
/Applications/BEAST\ 2.6.2/bin/LogCombiner -log T2_wgs_oristates_unidir_1_strict_rep[0-9].HCCtumor.typed.node.trees \
-o T2_wgs_oristates_unidir_strict_comb.HCCtumor.typed.node.trees \
-b 20

#/Applications/BEAST\ 2.6.2/bin/LogCombiner -log T2_wgs_oristates_unidir_[0-9]_strict_rep[0-9].log \
#-o T2_wgs_oristates_unidir_strict_comb.log \
#-b 10


# Remove single tip for tumor one, only state clock and original states
/Applications/BEAST\ 2.6.2/bin/LogCombiner -log T1red_wgs_oristates_unidir_1_state_rep[0-9].HCCtumor.typed.node.trees \
-o T1red_wgs_oristates_unidir_state_comb.HCCtumor.typed.node.trees \
-b 20

#/Applications/BEAST\ 2.6.2/bin/LogCombiner -log T1_wgs_oristates_unidir_[0-9]_state_rep[0-9]_rt1l13.log \
#-o T1_wgs_oristates_unidir_state_comb_rt1l13.log \
#-b 10


for trees_file in *_comb.HCCtumor.typed.node.trees
do
  MODEL=$(basename $trees_file .typed.node.trees)
  /Applications/BEAST\ 2.6.2/bin/treeannotator -burnin 0 -heights median $trees_file ${MODEL}_mcc.tree
done
