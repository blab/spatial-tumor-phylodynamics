#Ling et al tumors

augur tree \
  --alignment ../data/ling-punch-sequences.fa \
  --method fasttree \
  --output ling_tree_raw.nwk

augur refine \
  --alignment ../data/ling-punch-sequences.fa \
  --tree ling_tree_raw.nwk \
  --metadata ../data/punch_metadata.tsv \
  --output-tree ling_tree.nwk \
  --root Normal \
  --output-node-data branch_lengths.json

augur ancestral \
  --tree ling_tree.nwk \
  --alignment ../data/ling-punch-sequences.fa \
  --output-node-data nt_muts.json


augur export v2 \
  --tree ling_tree.nwk \
  --metadata ../data/punch_metadata.tsv \
  --node-data branch_lengths.json \
              nt_muts.json \
  --color-by-metadata edge \
  --output tumor.json

#Li et al tumors

augur tree \
  --alignment data/li_t1_wgs.fa \
  --method fasttree \
  --output li_t1_tree_raw.nwk

#Root based on trees made by Li et al
augur refine \
    --alignment data/li_t1_wgs.fa \
    --tree li_t1_tree_raw.nwk \
    --metadata data/t1_metadata.tsv \
    --output-tree li_t1_tree.nwk \
    --root blood \
    --output-node-data t1_branch_lengths.json

augur ancestral \
  --tree li_t1_tree.nwk \
  --alignment data/li_t1_wgs.fa \
  --output-node-data t1_nt_muts.json

augur export v2 \
  --tree li_t1_tree.nwk \
  --metadata data/t1_metadata.tsv \
  --node-data t1_branch_lengths.json \
              t1_nt_muts.json \
  --color-by-metadata edge \
  --output t1_tumor.json


augur tree \
      --alignment data/li_t2_wgs.fa \
      --method fasttree \
      --output li_t2_tree_raw.nwk
augur refine \
          --alignment data/li_t2_wgs.fa\
          --tree li_t2_tree_raw.nwk \
          --metadata data/t2_metadata.tsv \
          --output-tree li_t2_tree.nwk \
          --root blood \
          --output-node-data t2_branch_lengths.json

augur ancestral \
      --tree li_t2_tree.nwk \
      --alignment data/li_t2_wgs.fa \
      --output-node-data t2_nt_muts.json

augur export v2 \
    --tree li_t2_tree.nwk \
    --metadata data/t2_metadata.tsv \
    --node-data t2_branch_lengths.json \
                t2_nt_muts.json \
    --color-by-metadata edge \
    --output t2_tumor.json
