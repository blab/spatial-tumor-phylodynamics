augur tree \
  --alignment ../data/ling-punch-sequences.fa \
  --method fasttree \
  --output ling_tree_raw.nwk

augur refine \
    --alignment ../data/ling-punch-sequences.fa \
    --tree ling_tree_raw.nwk \
    --metadata ../data/punch_metadata.tsv \
    --output-tree ling_tree.nwk \
    --keep-root
    --root Normal\
    --output-node-data branch_lengths.json

  augur refine \
        --alignment ../data/ling-punch-sequences.fa \
        --tree ling_tree_raw.nwk \
        --metadata ../data/punch_metadata.tsv \
        --output-tree ling_tree.nwk \
        --keep-root \
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


data/punch_metadata.tsv
