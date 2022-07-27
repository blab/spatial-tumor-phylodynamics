#!/bin/sh
##script to find MCC tree for typed.node.trees in single directory
for trees_file in *.typed.node.trees
do
  MODEL=$(basename $trees_file .typed.node.trees)
  /Applications/BEAST\ 2.6.2/bin/treeannotator -burnin 10 -heights median $trees_file ${MODEL}_mcc.tree
done
