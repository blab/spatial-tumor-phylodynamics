#!/bin/bash

#combine posterior summary statistics (means, HPD intervals)
## Run in eden/stats/posteriors

#Add headers
head -n 1 death_rate_validation_pop_1000_dr_0.00_n_100_state_clock_estimate_dr_posterior_summary.tsv > posterior_summary_all.tsv

for file in *_posterior_summary.tsv
do
  #add summary line from each file
  tail -n 1 $file >> posterior_summary_all.tsv
done

#Get posterior summary
rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/stats/posteriors/posterior_summary_all.tsv /Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/eden/stats/posteriors

#Get example logs and trees
rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/death_rate_validation_pop_1000_dr_0.29_n_100_state_clock_estimate_dr.log /Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/eden/logs
rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/death_rate_validation_pop_1000_dr_0.29_n_100_state_clock_estimate_dr.typed.node.trees /Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/eden/trees
