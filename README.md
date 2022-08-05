# State-dependent evolutionary models reveal modes of solid tumor growth

Maya A. Lewinsohn<sup>1,2</sup>, Trevor Bedford<sup>1,2,3</sup>, Nicola F. Müller<sup>2</sup>, Alison F. Feder<sup>1</sup>

<sup>1</sup>Department of Genome Sciences, University of Washington, Seattle, WA, USA;
<sup>2</sup>Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Research Center, Seattle, WA, USA;
<sup>3</sup>Howard Hughes Medical Institute, Seattle, WA, USA

## Abstract
Spatial properties of tumor growth have profound implications for cancer progression, therapeutic resistance and
metastasis. Yet, how spatial position governs tumor cell division remains difficult to evaluate in clinical tumors.
Here, we demonstrate that elevated cellular growth rates on the tumor periphery leave characteristic patterns in
the genomes of cells sampled from different parts of a tumor, which become evident when they are used to
construct a tumor phylogenetic tree. Namely, rapidly-dividing peripheral lineages branch more extensively and
acquire more mutations than slower-dividing lineages in the tumor center. We develop a Bayesian
state-dependent evolutionary phylodynamic model (SDevo) that quantifies these patterns to infer the differential
cell division rates between peripheral and central cells jointly from the branching and mutational patterns of
single-time point, multi-region sequencing data. We validate this approach on simulated tumors by
demonstrating its ability to accurately infer spatially-varying birth rates under a range of growth conditions and
sampling strategies. We then show that SDevo outperforms state-of-the-art, non-cancer multi-state
phylodynamic methods which ignore differential mutational acquisition. Finally, we apply SDevo to
multi-region sequencing data from clinical hepatocellular carcinomas and find evidence that cells on the tumor
edge divide 2-4x faster than those in the center. As multi-region and single-cell sequencing increase in
resolution and availability, we anticipate that SDevo will be useful in interrogating spatial restrictions on tumor
growth and could be extended to model non-spatial factors that influence tumor progression, including hypoxia
and immune infiltration.

## Analyses and figures

### Eden simulation studies
Scripts to generate simulated tumors under boundary-driven and unrestricted growth in a 2D lattice can be found in `simulated_data/spatial_tumor_simulation.ipynb`.
Simulated tumors are recorded .csv files recording all cells in tumor simulation. For pushing simulation, cell positions in lattice are
recorded in locs.csv at each time slice.

`scripts/reconstruct_simulated_trees.R` contains code to convert simulated cells records into S4 tree objects containing edge/center state information.
Install local R package _tumortree_ for necessary functions.

### Physicell simulation studies
Simulation data is generated from https://github.com/federlab/PhysiCellTrees and should be put `phyiscell/simulation_data` for downstream analyses.

### Run SDEvo on simulated trees
Script to generate input XML files from simulated Eden simulation can be found here: `scripts/set_up_simulated_tumors_xmls.R`
Script to generate input XML files from PhysiCell outputs can be found here: `scripts/write_state_clocks_xml_from_physicell.R`

### Run strict clock model on simulated trees
Make XML files with strict clock by `running create_strict_clock.R`, where input is state clocks XML file generated above.

### HCC Tumor analysis with SDevo
Scripts to generate input DNA sequences HCC sequence data can be found in `scripts/process_li_data.R` and `li-application/`. Data will be available pending approval of original authors. Temporarily, a template xml is included for reference.

### Figures
Local R package _tumortree_ is needed for most figures.
