for dr in `seq -f "%f" 0.275 0.005 0.43`
do
  echo $dr
  sbatch submit_indiv_eden_sim.sh $dr
done

for file in eden/xml/death_rate_validation_pop_10000_mu_1_dr_0.[0-9][0-9][0-9]_n_50_state_clock_estimate_dr.xml
do
  sbatch scripts/submit_beast_run.sh $file
done

for file in eden/xml/death_rate_validation_pop_10000_mu_1_dr_0.[0-9][0-9][0-9]_n_100_state_clock_estimate_dr.xml
do
  sbatch scripts/submit_beast_run.sh $file
done

for file in eden/xml/death_rate_validation_pop_10000_mu_1_dr_0.[0-9][0-9][0-9]_n_100_state_clock_estimate_dr_random_sampling.xml
do
  sbatch scripts/submit_beast_run.sh $file
done

for file in eden/xml/death_rate_validation_pop_10000_mu_1_dr_0.[0-9][0-9][0-9]_n_30_state_clock_estimate_dr.xml
do
  sbatch scripts/submit_beast_run.sh $file
done

for file in eden/xml/*strict_clock.xml
do
  sbatch scripts/submit_beast_run.sh $file
done

for file in eden/xml/bulk*.xml
do
  sbatch scripts/submit_beast_run2.sh $file
done

#physicell

#submitted
for file in physicell/simulation_data/2D_neut_bdg_motility/diversified_100/to_beast_format/xml_files/*
do
  sbatch scripts/submit_beast_run.sh $file
done

for file in physicell/simulation_data/2D_neut_bdg_sigmoid/diversified_100/to_beast_format/xml_files/*
do
  sbatch scripts/submit_beast_run.sh $file
done

for file in physicell/simulation_data/2D_sel_bdg_deleterious/diversified_100/to_beast_format/xml_files/*
do
  sbatch scripts/submit_beast_run.sh $file
done

for file in physicell/simulation_data/3D_sel/diversified_100/to_beast_format/xml_files/*
do
  sbatch scripts/submit_beast_run.sh $file
done

for file in physicell/simulation_data/2D_neut_bdg_usecarryingcapacity/diversified_100/to_beast_format/xml_files/*
do
  sbatch scripts/submit_beast_run.sh $file
done

for file in physicell/simulation_data/2D_sel_bdg/diversified_100/to_beast_format/xml_files/*
do
  sbatch scripts/submit_beast_run.sh $file
done


for file in *estimate_dr.xml
do
  sh create_strict_clock_xml.sh $file
done

for file in *_random_sampling.xml
do
  sh create_strict_clock_xml.sh $file
done

for file in ../eden/simulation_data/cells_death_rate_validation_pop_10000_mu_1_dr_0.[0-9][0-9][0-9].csv
do
  echo $file
  #sbatch sub_xmls_indiv.sh $file 40
done

for file in ../eden/simulation_data/cells_death_rate_validation_pop_10000_mu_1_dr_0.[0-9][0-9][0-9].csv
do
  for sample_size in 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 30 40 50 60 70 80 90 100
  do

    echo $file
    echo $sample_size
    sbatch sub_xmls_indiv.sh $file $sample_size

  done
done

for file in ../eden/simulation_data/cells_death_rate_validation_pop_10000_mu_1_dr_0.[0-9][0-9][0-9].csv
do
  for sample_size in 5 6 7 8 9
  do

    echo $file
    echo $sample_size
    sbatch sub_xmls_indiv.sh $file $sample_size

  done
done


for file in ../eden/simulation_data/cells_death_rate_validation_pop_10000_mu_1_dr_0.[0-9][0-9][0-9].csv
do
  echo $file
  sbatch sub_xmls_indiv_random.sh $file 100
done

for file in ../eden/simulation_data/cells_death_rate_validation_pop_10000_mu_1_dr_0.[0-9][0-9][0-9].csv
do
  echo $file
  sbatch bulk_sample_submit.sh $file
done

#rsyncing
rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/simtrees/cells_death_rate_validation_pop_10000_mu_1_dr_0.050.rds ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/eden/simtrees
rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/simtrees/cells_pushing_pop_10000_mu_1_dr_0.050.rds ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/eden/simtrees

rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/simulation_data/cells_pushing_pop_10000_mu_1_dr_0.050.csv ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/eden/simulation_data
rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/simulation_data/cells_pushing_pop_10000_mu_1_dr_0.050.csv ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/eden/simulation_data

rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/physicell/logs/*.log ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/physicell/logs
rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/eden/logs/*.log ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/eden/logs
