
#Li et al tumors
fasttree -nt data/li_t1_wgs.fasta > li_t1_tree_raw.nwk
#augur tree \
#  --alignment data/li_t1_wgs.fasta \
#  --method fasttree \
#  --output li_t1_tree_raw.nwk

#augur tree \
#      --alignment data/li_t1_wgs.fasta \
#      --method iqtree \
#      --output li_t1_tree_raw.nwk \
#      --tree-builder-args '-bb 2000' \
#      --override-default-args

#Root based on trees made by Li et al
augur refine \
    --alignment data/li_t1_wgs.fasta \
    --tree li_t1_tree_raw.nwk \
    --metadata data/t1_metadata.tsv \
    --output-tree li_t1_tree.nwk \
    --root blood \
    --output-node-data t1_branch_lengths.json

augur ancestral \
  --tree li_t1_tree.nwk \
  --alignment data/li_t1_wgs.fasta \
  --output-node-data t1_nt_muts.json

augur export v2 \
  --tree li_t1_tree.nwk \
  --metadata data/t1_metadata.tsv \
  --node-data t1_branch_lengths.json \
              t1_nt_muts.json \
  --color-by-metadata edge \
  --output t1_tumor.json


#augur tree \
#      --alignment data/li_t2_wgs.fa \
#      --method fasttree \
#      --output li_t2_tree_raw.nwk

#augur tree \
#            --alignment data/li_t2_wgs.fasta \
#            --method iqtree \
#            --output li_t2_tree_raw.nwk \
#            --tree-builder-args '-bb 2000 -czb -m GTR+ASC+I+R' \
#            --override-default-args

#iqtree -s data/li_t2_wgs.fasta -czb -m JC+ASC+I+R --bcon 50 -redo
#iqtree -s data/li_t2_wgs.fasta -czb -m GTR+ASC+I+R+T -redo --root-seq data/li_t2_wgs.fasta,blood
#cat data/li_t2_wgs.fasta.treefile > li_t2_tree_raw.nwk
fasttree -nt data/li_t2_wgs.fasta > li_t2_tree_raw.nwk
augur refine \
          --alignment data/li_t2_wgs.fasta\
          --tree li_t2_tree_raw.nwk \
          --metadata data/t2_metadata.tsv \
          --output-tree li_t2_tree.nwk \
          --root blood \
          --output-node-data t2_branch_lengths.json

augur ancestral \
      --tree li_t2_tree.nwk \
      --alignment data/li_t2_wgs.fasta \
      --output-node-data t2_nt_muts.json

augur export v2 \
    --tree li_t2_tree.nwk \
    --metadata data/t2_metadata.tsv \
    --node-data t2_branch_lengths.json \
                t2_nt_muts.json \
    --color-by-metadata edge \
    --output t2_tumor.json
