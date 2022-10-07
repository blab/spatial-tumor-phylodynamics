rsync -a ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/scripts/run_eden_simulation*.py mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/simulation_data
rsync -a ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/scripts/*.py mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/scripts
rsync -a ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/scripts/* mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/scripts

for file in run_eden_simulations[0-9].py
do
  sbatch submit_sim_batch.sh $file
done

for file in run_eden_simulations2[4-7].py
do
  sbatch submit_sim_batch.sh $file
done

rsync -a ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/scripts/reconstruct_sim_trees_large*.Rscript mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/scripts


for file in set_up_simulated_tumors_xmls_large[0-9].Rscript
do
  sbatch submit_Rscript.sh $file
done
