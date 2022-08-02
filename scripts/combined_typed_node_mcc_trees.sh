#!/bin/sh
##script to find combine tree files

#Tumor 1
## New states / unidirectional
/Applications/BEAST\ 2.6.2/bin/LogCombiner -log T1_wgs_newstates_unidir_state_rep0.HCCtumor.typed.node.trees \
-log T1_wgs_newstates_unidir_state_rep1.HCCtumor.typed.node.trees \
-log T1_wgs_newstates_unidir_state_rep2.HCCtumor.typed.node.trees \
-o T1_wgs_newstates_unidir_state_comb.HCCtumor.typed.node.trees \
-b 10
## Old states / unidirectional
/Applications/BEAST\ 2.6.2/bin/LogCombiner -log T1_wgs_oristates_unidir_state_rep0.HCCtumor.typed.node.trees \
-log T1_wgs_oristates_unidir_state_rep1.HCCtumor.typed.node.trees \
-log T1_wgs_oristates_unidir_state_rep2.HCCtumor.typed.node.trees \
-o T1_wgs_oristates_unidir_state_comb.HCCtumor.typed.node.trees \
-b 10

#Tumor 2
## New states / unidirectional
/Applications/BEAST\ 2.6.2/bin/LogCombiner -log T2_wgs_newstates_unidir_state_rep0.HCCtumor.typed.node.trees \
-log T2_wgs_newstates_unidir_state_rep1.HCCtumor.typed.node.trees \
-log T2_wgs_newstates_unidir_state_rep2.HCCtumor.typed.node.trees \
-o T2_wgs_newstates_unidir_state_comb.HCCtumor.typed.node.trees \
-b 10
## Old states / unidirectional
/Applications/BEAST\ 2.6.2/bin/LogCombiner -log T1_wgs_oristates_unidir_state_rep0_r.HCCtumor.typed.node.trees \
-log T2_wgs_oristates_unidir_state_rep1.HCCtumor.typed.node.trees \
-log T2_wgs_oristates_unidir_state_rep2.HCCtumor.typed.node.trees \
-o T2_wgs_oristates_unidir_state_comb.HCCtumor.typed.node.trees \
-b 10

#Strict clock, old stats
/Applications/BEAST\ 2.6.2/bin/LogCombiner -log T1_wgs_oristates_unidir_strict_rep0.HCCtumor.typed.node.trees \
-log T1_wgs_oristates_unidir_strict_rep1.HCCtumor.typed.node.trees \
-log T1_wgs_oristates_unidir_strict_rep2.HCCtumor.typed.node.trees \
-o T1_wgs_oristates_unidir_strict_comb.HCCtumor.typed.node.trees \
-b 10

/Applications/BEAST\ 2.6.2/bin/LogCombiner -log T2_wgs_oristates_unidir_strict_rep0.HCCtumor.typed.node.trees \
-log T2_wgs_oristates_unidir_strict_rep1.HCCtumor.typed.node.trees \
-log T2_wgs_oristates_unidir_strict_rep2.HCCtumor.typed.node.trees \
-o T2_wgs_oristates_unidir_strict_comb.HCCtumor.typed.node.trees \
-b 10

# Remove single tip
/Applications/BEAST\ 2.6.2/bin/LogCombiner -log T1_wgs_oristates_unidir_state_rep0_rt1l13.HCCtumor.typed.node.trees \
-log T1_wgs_oristates_unidir_state_rep1_rt1l13.HCCtumor.typed.node.trees \
-log T1_wgs_oristates_unidir_state_rep2_rt1l13.HCCtumor.typed.node.trees \
-o T1_wgs_oristates_unidir_state_comb_rt1l13.HCCtumor.typed.node.trees\
-b 10

for trees_file in *_comb.HCCtumor.typed.node.trees
do
  MODEL=$(basename $trees_file .typed.node.trees)
  /Applications/BEAST\ 2.6.2/bin/treeannotator -burnin 0 -heights median $trees_file ${MODEL}_mcc.tree
done

/Applications/BEAST\ 2.6.2/bin/treeannotator -burnin 0 -heights median T1_wgs_oristates_unidir_state_comb_rt1l13.HCCtumor.typed.node.trees T1_wgs_oristates_unidir_state_comb_rt1l13_mcc.tree
