#!/bin/bash

#combine posterior summary statistics (means, HPD intervals)
## Run in eden/stats/posteriors


#Add headers
head -n 1 bulk_sampling_death_rate_validation_pop_1000_dr_0.00_posterior_summary.tsv > bulk_posterior_summary_all.tsv

for file in *_posterior_summary.tsv
do
  #add summary line from each file
  tail -n 1 $file >> bulk_posterior_summary_all.tsv
done


# BULK

## Run in eden/bulk/stats/posteriors
head -n 1 death_rate_validation_pop_1000_dr_0.00_n_100_state_clock_estimate_dr_posterior_summary.tsv > posterior_summary_all.tsv

for file in *_posterior_summary.tsv
do
  #add summary line from each file
  tail -n 1 $file >> posterior_summary_all.tsv
done

# PHYSICELL runs

# 2D neutral BDG
head -n 1 sampconfig_m0_w1_d0.1_t1_mg1_mm1_l2e+08_i1_s482_posterior_summary.tsv > 2D_neut_bdg_posterior_summary_all.tsv

for file in *_posterior_summary.tsv
do
  #add summary line from each file
  tail -n 1 $file >> 2D_neut_bdg_posterior_summary_all.tsv
done

# 2D neutral BDG -- motility
head -n 1 sampconfig_m0_w1_d0.2_t1_mg1_mm1_l2e+08_i1_s4522_diversified_m1_n100_posterior_summary.tsv > 2D_neut_bdg_motility_posterior_summary_all.tsv

for file in *_posterior_summary.tsv
do
  #add summary line from each file
  tail -n 1 $file >> 2D_neut_bdg_motility_posterior_summary_all.tsv
done

# 2D neutral BDG -- sigmoid
head -n 1 sampconfig_m0_w1_d0.2_t1_mg1_mm1_l2e+08_i1_s15745_diversified_m1_n100_posterior_summary.tsv  > 2D_neut_bdg_sigmoid_posterior_summary_all.tsv

for file in *_posterior_summary.tsv
do
  #add summary line from each file
  tail -n 1 $file >> 2D_neut_bdg_sigmoid_posterior_summary_all.tsv
done

# 2D sel bdg
head -n 1 sampconfig_m0_w1.1_d0.1_t1_mg0.99_mm1_l2e+08_i1_s20363_diversified_m1_n100_posterior_summary.tsv > 2D_sel_bdg_posterior_summary_all.tsv

for file in *_posterior_summary.tsv
do
  #add summary line from each file
  tail -n 1 $file >> 2D_sel_bdg_posterior_summary_all.tsv
done

# 3D netural BDG

head -n 1 sampconfig_m0_w1_d0.1_t1_mg1_mm1_l2e+08_i1_s25599_posterior_summary.tsv > 3D_neut_bdg_posterior_summary_all.tsv
for file in *_posterior_summary.tsv
do
  #add summary line from each file
  tail -n 1 $file >> 3D_neut_bdg_posterior_summary_all.tsv
done

#3D sel bdg
head -n 1 sampconfig_m0_w1.1_d0.2_t1_mg0.99_mm1_l2e+08_i1_s35083_diversified_m1_n100_posterior_summary.tsv > 3D_sel_bdg_posterior_summary_all.tsv
for file in *_posterior_summary.tsv
do
  #add summary line from each file
  tail -n 1 $file >> 3D_sel_bdg_posterior_summary_all.tsv
done
#Get posterior summary

#Eden simulations
rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/stats/posteriors/posterior_summary_all.tsv /Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/eden/stats/posteriors
rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/bulk/stats/posteriors/bulk_posterior_summary_all.tsv /Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/eden/stats/posteriors

#Physicell simulations
rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/physicell/simulation_data/2D_neut_bdg/stats/posteriors/2D_neut_bdg_posterior_summary_all.tsv /Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/physicell/stats/posteriors
rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/physicell/simulation_data/2D_neut_bdg_sigmoid/stats/posteriors/2D_neut_bdg_sigmoid_posterior_summary_all.tsv /Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/physicell/stats/posteriors
rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/physicell/simulation_data/2D_neut_bdg_motility/stats/posteriors/2D_neut_bdg_motility_posterior_summary_all.tsv /Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/physicell/stats/posteriors

rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/physicell/simulation_data/2D_sel_bdg/stats/posteriors/2D_sel_bdg_posterior_summary_all.tsv /Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/physicell/stats/posteriors
rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/physicell/simulation_data/3D_sel/stats/posteriors/3D_sel_bdg_posterior_summary_all.tsv /Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/physicell/stats/posteriors
rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/physicell/simulation_data/3D_neut_bdg/stats/posteriors/3D_neut_bdg_posterior_summary_all.tsv /Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/physicell/stats/posteriors

#Get example logs and trees
rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/death_rate_validation_pop_1000_dr_0.29_n_100_state_clock_estimate_dr.log /Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/eden/logs
rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/death_rate_validation_pop_1000_dr_0.29_n_100_state_clock_estimate_dr.typed.node.trees /Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/eden/trees
