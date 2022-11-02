#for dr in `seq -f "%f" 0.275 0.005 0.43`
#do
#  echo $dr
#  sbatch submit_indiv_eden_sim.sh $dr
#done

for dr in `seq -f "%f" 0.000 0.005 0.100`
do
  echo $dr
  sbatch submit_indiv_eden_sim_lowmem.sh $dr
done


for dr in `seq -f "%f" 0.405 0.005 0.43`
do
  echo $dr
  sbatch submit_indiv_eden_sim_lowmem.sh $dr
done

rsync -a /fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/simulation_data2/cells_death_rate_validation_pop_10000_mu_1_dr_0.[0-3]* /fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/simulation_data
rsync -a /fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/xml2/death_rate_validation_pop_10000_mu_1_dr_0.2[8-9][0-9]_n_50_state_clock_estimate_dr.xml /fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/xml

#submit test tumor
#sbatch submit_indiv_eden_sim_lowmem_smallwq.sh 0.252
#sbatch submit_indiv_eden_sim_lowmem_smallwq.sh 0.411

for file in eden/xml2/death_rate_validation2_pop_10000_mu_1_dr_0.[2-3][0-9][0-9]_n_50_state_clock_estimate_dr.xml
do
  echo $file
  sbatch scripts/submit_beast_run2.sh $file
done

for file in eden/xml2/death_rate_validation2_pop_10000_mu_1_dr_0.[0-3][0-9][0-9]_n_50_state_clock_estimate_dr_strict_clock.xml
do
  echo $file
  sbatch scripts/submit_beast_run2.sh $file
done

for file in eden/xml2/death_rate_validation2_pop_10000_mu_1_dr_0.205_n_100_state_clock_estimate_dr.xml
do
  echo $file
  sbatch scripts/submit_beast_run2.sh $file
done

for file in eden/xml/death_rate_validation_pop_10000_mu_1_dr_0.[0-9][0-9][0-9]_n_100_state_clock_estimate_dr_strict_clock.xml
do
  echo $file
  sbatch scripts/submit_beast_run.sh $file
done

for file in eden/xml/death_rate_validation_pop_10000_mu_1_dr_0.[0-9][0-9][0-9]_n_100_state_clock_estimate_dr_random_sampling.xml
do
  echo $file
  sbatch scripts/submit_beast_run.sh $file
done

for file in eden/xml/death_rate_validation_pop_10000_mu_1_dr_0.[0-9][0-9][0-9]_n_100_state_clock_estimate_dr_random_sampling_strict_clock.xml
do
  echo $file
  sbatch scripts/submit_beast_run.sh $file
done

#Special runs
sbatch scripts/submit_beast_run.sh eden/xml/death_rate_validation_pop_10000_mu_1_dr_0.195_n_100_state_clock_estimate_dr.xml
sbatch scripts/submit_beast_run.sh eden/xml2/death_rate_validation_pop_10000_mu_1_dr_0.095_n_50_state_clock_estimate_dr.xml
sbatch scripts/submit_beast_run.sh eden/xml/death_rate_validation_pop_10000_mu_1_dr_0.195_n_50_state_clock_estimate_dr.xml
sbatch scripts/submit_beast_run.sh eden/xml/death_rate_validation_pop_10000_mu_1_dr_0.295_n_50_state_clock_estimate_dr.xml

for file in eden/xml/death_rate_validation_pop_10000_mu_1_dr_0.[0-9][0-9][0-9]_n_100_state_clock_estimate_dr.xml
do
  sbatch scripts/submit_beast_run.sh $file
done

for file in eden/xml/death_rate_validation_pop_10000_mu_1_dr_0.[0-9][0-9][0-9]_n_100_state_clock_estimate_dr_random_sampling.xml
do
  sbatch scripts/submit_beast_run.sh $file
done

for file in eden/xml/death_rate_validation_pop_10000_mu_1_dr_0.[0-9][0-9][0-9]_n_50_state_clock_estimate_dr.xml
do
  sbatch scripts/submit_beast_run.sh $file
done

for file in eden/xml/death_rate_validation_pop_10000_mu_1_dr_0.[0-9][0-9][0-9]_n_50_state_clock_estimate_dr_strict_clock.xml
do
  sbatch scripts/submit_beast_run.sh $file
done

for file in eden/xml/*strict_clock.xml
do
  sbatch scripts/submit_beast_run.sh $file
done

for file in eden/xml/bulk_sampling_diversifieddeath_rate_validation_pop_10000_mu_1_dr_0.0*
do
  echo $file
  sbatch scripts/submit_beast_run2.sh $file
done

#Script to run to desired ESS
for file in eden/xml2/death_rate_validation2_pop_10000_mu_1_dr_0.[0-9][0-9][0-9]_n_50_state_clock_estimate_dr_strict_clock.xml
do
  echo $file
  sbatch scripts/run_beast_to_ess2.sh $file
done

for file in eden/xml2/death_rate_validation2_pop_10000_mu_1_dr_0.[0-9][0-9][0-9]_n_50_state_clock_estimate_dr.xml
do
  echo $file
  sbatch scripts/run_beast_to_ess2.sh $file
done

for file in eden/xml2/death_rate_validation2_pop_10000_mu_1_dr_0.[0-9][0-9][0-9]_n_100_state_clock_estimate_dr.xml
do
  echo $file
  sbatch scripts/run_beast_to_ess2.sh $file
done

for file in eden/xml/death_rate_validation2_pop_10000_mu_1_dr_0.[0-9][0-9][0-9]_n_100_state_clock_estimate_dr_strict_clock.xml
do
  echo $file
  sbatch scripts/run_beast_to_ess2.sh $file
done

for file in eden/xml2/death_rate_validation2_pop_10000_mu_1_dr_0.[0-9][0-9][0-9]_n_50_state_clock_estimate_dr.xml
do
  echo $file
  sbatch scripts/run_get_posteriors.sh $file
done

for file in eden/xml2/death_rate_validation2_pop_10000_mu_1_dr_0.[0-9][0-9][0-9]_n_50_state_clock_estimate_dr_strict_clock.xml
do
  echo $file
  sbatch scripts/run_get_posteriors.sh $file
done

#sh scripts/run_beast_to_esssq.sh eden/xml/death_rate_validation_pop_10000_mu_1_dr_0.260_n_50_state_clock_estimate_dr.xml
#physicell

#submitted
#test: sh submit_single_batch.sh physicell/simulation_data/2D_neut_bdg_motility/diversified_100/to_beast_format/xml_files/sampconfig_m2_w1_d0.2_t1_mg1_mm1_l2e+08_i8_s34287_diversified_m1_n100.xml
for file in physicell/simulation_data/2D_neut_bdg_motility/diversified_100/to_beast_format/xml_files/*
do
  sbatch scripts/submit_beast_run2.sh $file
  #qsub -cwd submit_single_batch.sh $file
done

for file in physicell/simulation_data/2D_neut_bdg_sigmoid/diversified_100/to_beast_format/xml_files/*
do
  sbatch scripts/submit_beast_run2.sh $file
  #qsub -cwd -l mem_free=2G submit_single_batch.sh $file
done

for file in physicell/simulation_data/2D_sel_bdg_deleterious/diversified_100/to_beast_format/xml_files/*
do
  sbatch scripts/submit_beast_run2.sh $file
  #qsub -cwd -l mem_free=2G submit_single_batch.sh $file
done

#TODO: Run xml files (transfer first to gs cluster)
for file in physicell/simulation_data/3D_sel/diversified_100/to_beast_format/xml_files/*
do
  sbatch scripts/submit_beast_run.sh $file
  #qsub -cwd -l mem_free=2G submit_single_batch.sh $file

done

for file in physicell/simulation_data/2D_neut_bdg_usecarryingcapacity/diversified_100/to_beast_format/xml_files/*
do
  sbatch scripts/submit_beast_run2.sh $file
  #qsub -cwd -l mem_free=2G submit_single_batch.sh $file
done

for file in physicell/simulation_data/2D_sel_bdg/diversified_100/to_beast_format/xml_files/*
do
  sbatch scripts/submit_beast_run2.sh $file
  #qsub -cwd -l mem_free=2G submit_single_batch.sh $file
done


for file in *_state_clock_estimate_dr.xml
do
  sh create_strict_clock_xml.sh $file
done

for file in *_random_sampling.xml
do
  sh create_strict_clock_xml.sh $file
done
sbatch sub_xmls_indiv.sh ../eden/simulation_data/cells_death_rate_validation_pop_10000_mu_1_dr_0.275.csv 50
sbatch sub_xmls_indiv.sh ../eden/simulation_data/cells_death_rate_validation_pop_10000_mu_1_dr_0.195.csv 100
sbatch sub_xmls_indiv.sh ../eden/simulation_data/cells_death_rate_validation_pop_10000_mu_1_dr_0.195.csv 50
sh sub_xmls_indiv.sh ../eden/simulation_data/cells_death_rate_validation_pop_10000_mu_1_dr_0.195.csv 50
sh sub_xmls_indiv.sh ../eden/simulation_data/cells_death_rate_validation_pop_1000_mu_1_dr_0.411.csv

#redo 0.275
for file in ../eden/simulation_data/cells_death_rate_validation_pop_10000_mu_1_dr_[0-3][0-9][0-9].csv
do
  echo $file
  sbatch sub_xmls_indiv.sh $file 50
done

for file in ../eden/simulation_data/cells_death_rate_validation_pop_10000_mu_1_dr_0.[0-3][0-9][0-9].csv
do
  echo $file
  sbatch sub_xmls_indiv.sh $file 100
done

for file in ../eden/simulation_data/cells_death_rate_validation_pop_10000_mu_1_dr_0.400.csv
do
  echo $file
  sbatch sub_xmls_indiv_random.sh $file 100
done


rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/simulation_data/cells_death_rate_validation_pop_1000_mu_1_dr_0.252.csv ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/eden/simulation_data
#test run for new simulations low mem
#sbatch sub_xmls_indiv.sh ../eden/simulation_data/cells_death_rate_validation_pop_10000_mu_1_dr_0.295.csv 50

#Reconstruct simtrees
for file in ../eden/simulation_data/cells_death_rate_validation_pop_10000_mu_1_dr_0.[0-4][0-9][0-9].csv
do
  echo $file
  sbatch tree_recon_submit.sh $file
done

#sbatch tree_recon_submit.sh ../eden/simulation_data/cells_death_rate_validation_pop_10000_mu_1_dr_0.150.csv
#sbatch tree_recon_submit.sh ../eden/simulation_data/cells_death_rate_validation_pop_10000_mu_1_dr_0.300.csv
for file in ../eden/simulation_data/cells_pushing_pop_10000_mu_1_dr_0.[0-9][0-9][0-9].csv
do
  echo $file
  sbatch tree_recon_submit.sh $file
done


#To test script
#sh sub_xmls_indiv.sh ../eden/simulation_data/cells_death_rate_validation_pop_10000_mu_1_dr_0.320.csv 20

# To create all state clock XMLS
## Sampling should be reproducible for each individual run
## Will not overwrite existing XML files unless set overwrite=TRUE

for file in ../eden/simulation_data/cells_death_rate_validation_pop_10000_mu_1_dr_0.400.csv
do
  for sample_size in 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
  do

    echo $file
    echo $sample_size
    sbatch sub_xmls_indiv.sh $file $sample_size

  done
done

for file in ../eden/simulation_data/cells_death_rate_validation_pop_10000_mu_1_dr_0.400.csv
do
  for sample_size in 20 30 40 60 70 80 90
  do

    echo $file
    echo $sample_size
    sbatch sub_xmls_indiv.sh $file $sample_size

  done
done



for file in ../eden/simulation_data/cells_death_rate_validation_pop_10000_mu_1_dr_0.[0-9][05][0-9].csv
do
  echo $file
  #sbatch sub_xmls_indiv_random.sh $file 100
done

#test
#sbatch bulk_sample_diversified_submit.sh ../eden/simulation_data/cells_death_rate_validation_pop_10000_mu_1_dr_0.005.csv
#Rscript bulk_sample_simulated_tumor_large_indiv_diversified.Rscript ../eden/simulation_data/cells_death_rate_validation_pop_10000_mu_1_dr_0.235.csv
#Rscript bulk_sample_simulated_tumor_large_indiv_diversified2.Rscript ../eden/simulation_data/cells_death_rate_validation_pop_1000_dr_0.10.csv
#Rscript bulk_sample_simulated_tumor_large_indiv_diversified.Rscript ../eden/simulation_data/cells_death_rate_validation_pop_1000_dr_0.10.csv
for file in ../eden/simulation_data/cells_death_rate_validation_pop_10000_mu_1_dr_0.1[0-9][0-9].csv
do
  echo $file
  sbatch bulk_sample_diversified_submit.sh $file
done

#sbatch rates_submit.sh ../eden/simulation_data/cells_death_rate_validation_pop_10000_mu_1_dr_0.320.csv
#sh rates_submit.sh ../eden/simulation_data/cells_death_rate_validation_pop_10000_mu_1_dr_0.005.csv

#TODO: replace all 0.2 runs
for file in ../eden/simulation_data2/cells_death_rate_validation_pop_10000_mu_1_dr_0.[0-9][0-9][0-9].csv
do
  echo $file
  sbatch rates_submit.sh $file
done

#Calculalte distance from edge to number of divisions
for file in ../eden/simulation_data/cells_death_rate_validation_pop_10000_mu_1_dr_0.05_i_[0-9].csv
do
  echo $file
  sbatch edge_dist_v_div_submit.sh $file
done

for file in ../eden/simulation_data/cells_pushing_pop_10000_mu_1_dr_0.050_i_[0-9].csv
do
  echo $file
  sbatch edge_dist_v_div_submit.sh $file
done

#Combine files -- boundary-driven growth
##Run in eden/stats
head -n1 edge_dist_vs_divisions_mu_1_dr_0.05_i_0_boundary_driven.tsv > edge_dist_vs_divisions_mu_1_dr_0.050_boundary_driven.tsv

#Add all rows starting after row 2 (should be only 2 rows per file, but could combine multiple simulations)

for file in edge_dist_vs_divisions_mu_1_dr_0.05_i_[0-9]_boundary_driven.tsv
do
    echo $file
    tail -n+2 $file >> edge_dist_vs_divisions_mu_1_dr_0.050_boundary_driven.tsv
done

#Combine files -- unrestricted
head -n1 edge_dist_vs_divisions_mu_1_dr_0.050_i_9_unrestricted.tsv > edge_dist_vs_divisions_mu_1_dr_0.050_unrestricted.tsv

#Add all rows starting after row 2 (should be only 2 rows per file, but could combine multiple simulations)
for file in edge_dist_vs_divisions_mu_1_dr_0.050_i_[0-9]_unrestricted.tsv
do
    echo $file
    tail -n+2 $file >> edge_dist_vs_divisions_mu_1_dr_0.050_unrestricted.tsv
done

rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/stats/edge_dist_vs_divisions_mu_1_dr_0.050* ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/eden/stats

# Rsync to gs cluster to run eden beast runs
## Run from eden/xml on rhino
rsync -a *.xml lewinsom@nexus.gs.washington.edu:/net/gs/vol1/home/lewinsom/StateClocks/eden/xml_files

#Rsync physicell simulations to gs cluster from rhino
## RUn from physicell/simulation_data
rsync -a * lewinsom@nexus.gs.washington.edu:/net/gs/vol1/home/lewinsom/StateClocks/physicell/simulation_data

#Run from home directory with Sdevo.jrar
rsync -a SDevo.jar lewinsom@nexus.gs.washington.edu:/net/gs/vol1/home/lewinsom/StateClocks
#Submit all runs to gs cluster
ssh lewinsom@nexus.gs.washington.edu
module load modules modules-init modules-gs
module load gsits-util/1.0
ssh lewinsom@grid.gs.washington.edu
cd /net/gs/vol1/home/lewinsom/StateClocks

rsync -a lewinsom@nexus.gs.washington.edu:/net/gs/vol1/home/lewinsom/StateClocks/physicell/logs  ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/physicell
for file in eden/xml_files/*.xml
do
  echo $file
  qsub -cwd -l mem_free=4G submit_single_batch.sh $file
done

#Run in eden/stats to combine all simulations into single file
#Add header
head -n1 validation_growth_and_death_rates_weighted_large2_0.000.tsv > validation_growth_and_death_rates_weighted_large2.tsv

#Add all rows after header (starting at line 2)
for file in validation_growth_and_death_rates_weighted_large2_0.[0-9][0-9][0-9].tsv
do
    tail -n+2 $file >> validation_growth_and_death_rates_weighted_large2.tsv
done

for file in ../eden/simulation_data/*_mu_1_i_[0-9]_dr_0.050.csv
do
  echo $file
  sbatch bl_stats_submit.sh $file
done

for file in ../eden/simulation_data/*_mu_1_i_[0-9]_dr_0.050.csv
do
  echo $file
  sbatch bl_stats_submit.sh $file
done

for file in ../eden/simtrees/*pop_10000_mu_1_dr_0.[0-9][0-9][0-9].rds
do
  echo $file
  sbatch corr_edge_stats_submit.sh  $file
done

for file in ../eden/simtrees/cells_death_rate_validation_pop_10000_mu_1_dr_0.26[0-9].rds
do
  echo $file
  sbatch corr_edge_stats_submit.sh  $file
done







#Run in eden/stats to combine all iterations into single file

#Add header
head -n1 terminal_bl_clock_rate_stats_full_large_boundary_driven_0.tsv > terminal_bl_clock_rate_stats_full_large.tsv

#Add all rows starting after row 2 (should be only 2 rows per file, but could combine multiple simulations)
for file in terminal_bl_clock_rate_stats_full_large_*
do
    echo $file
    tail -n+2 $file >> terminal_bl_clock_rate_stats_full_large.tsv
done

#Do same with time spend on edge and center / terminal branch length stats

head -n1 edge_center_stats_large_boundary_driven_dr_0.000.tsv > edge_center_stats_large.tsv

#Add all rows starting after row 2 (should be only 2 rows per file, but could combine multiple simulations)
for file in edge_center_stats_large_*_dr_0.[0-9][0-9][0-9].tsv
do
    echo $file
    tail -n+2 $file >> edge_center_stats_large.tsv
done

head -n1 death_rate_validation2_pop_10000_mu_1_dr_0.000_n_50_state_clock_estimate_dr_posterior_summary.tsv > posterior_summary.tsv

#Add all rows starting after row 2 (should be only 2 rows per file, but could combine multiple simulations)
for file in death_rate_validation2_pop_10000_mu_1_dr_0.[0-9][0-9][0-9]_n_50_*_posterior_summary.tsv
do
    echo $file
    tail -n+2 $file >> posterior_summary.tsv
done

death_rate_validation_pop_10000_mu_1_dr_0.015_n_50_state_clock_estimate_dr_strict_clock_posterior_summary.tsv

#Example tumor for figure 3
## Filter alive cells

sbatch filter_submit.sh ../eden/simulation_data/cells_death_rate_validation_pop_10000_mu_1_dr_0.195.csv
for file in ../eden/simulation_data/cells_death_rate_validation_pop_10000_mu_1_dr_0.[0-9][0-9][0-9].csv
do
  echo $file
  sbatch filter_submit.sh $file
done
#Then get log, trees and alive cells for plotting\\
rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/logs/* ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/eden/logs

rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/logs/death_rate_validation_pop_10000_mu_1_dr_0.195_n_50_state_clock_estimate_dr.log ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/eden/logs

rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/logs/death_rate_validation_pop_10000_mu_1_dr_0.195_n_100_state_clock_estimate_dr.log ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/eden/logs
rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/trees2/death_rate_validation2_pop_10000_mu_1_dr_0.205_n_50_state_clock_estimate_dr.typed.node.trees ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/eden/trees
rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/simulation_data/alive_cells_death_rate_validation_pop_10000_mu_1_dr_0.195.csv ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/eden/simulation_data
rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/simtrees/cells_death_rate_validation_pop_10000_mu_1_dr_0.195.rds ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/eden/simtrees


#Find MCC tree
#Ling et al
java -jar SDevo.jar xmls/hcc-wes_unidir_state_rep1.xml
java -jar SDevo.jar -seed 11 xmls/hcc-wes_unidir_state_rep2.xml
java -jar SDevo.jar -seed 22 xmls/hcc-wes_unidir_state_rep3.xml
java -jar SDevo.jar -seed 44 xmls/hcc-wes_unidir_state_strict_clock_rep3.xml
java -jar SDevo.jar -seed 33 xmls/hcc-wes_unidir_state_strict_clock_rep2.xml
java -jar SDevo.jar xmls/hcc-wes_unidir_state_strict_clock_rep1.xml
/Applications/BEAST\ 2.6.2/bin/treeannotator -heights median -burnin 10 hcc-wes.HCCtumor.typed.node.trees hcc-wes.HCCtumor.typed.node.mcc.tree

/Applications/BEAST\ 2.6.2/bin/treeannotator -heights median -burnin 10 death_rate_validation2_pop_10000_mu_1_dr_0.205_n_50_state_clock_estimate_dr.typed.node.trees death_rate_validation_pop_10000_mu_1_dr_0.205_n_50_state_clock_estimate_dr.typed.node.mcc.tree
rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/stats/* ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/eden/stats

rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/stats/* ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/eden/stats

#rsyncing
rsync -a /fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/xml2/* /fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/xml
rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/scripts/*.sh ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/scripts
rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/tumortree/** ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/tumortree

rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/simtrees/cells_death_rate_validation_pop_10000_mu_1_dr_0.050.rds ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/eden/simtrees
rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/simtrees/cells_pushing_pop_10000_mu_1_dr_0.050.rds ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/eden/simtrees

rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/simulation_data/cells_pushing_pop_10000_mu_1_dr_0.050.csv ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/eden/simulation_data
rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/simulation_data/cells_pushing_pop_10000_mu_1_dr_0.050.csv ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/eden/simulation_data

rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/physicell/logs/*.log ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/physicell/logs
rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/logs/*.log ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/eden/logs
