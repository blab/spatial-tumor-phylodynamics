#Revisons

rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/samp*.log ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/physicell/logs
#rsync -a ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/SDevo.jar mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics


for file in write_state_clocks_xml_from_physicell[0-9].Rscript
do
 echo $file
 sbatch submit_Rscript.sh $file
done

#submitted
#test: sh submit_single_batch.sh physicell/simulation_data/2D_neut_bdg_motility/diversified_100/to_beast_format/xml_files/sampconfig_m2_w1_d0.2_t1_mg1_mm1_l2e+08_i8_s34287_diversified_m1_n100.xml
cat SDevo.jar > physicell/simulation_data/2D_neut_bdg_motility/diversified_100/to_beast_format/SDevo.jar
cat scripts/submit_beast_run2.sh > physicell/simulation_data/2D_neut_bdg_motility/diversified_100/to_beast_format/submit_beast_run2.sh

cd physicell/simulation_data/2D_neut_bdg_motility/diversified_100/to_beast_format
for file in xml_files/*
do
  #sbatch submit_beast_run2.sh $file
  sbatch ../../../../../scripts/run_physicell_beast_to_ess.sh $file
  #qsub -cwd submit_single_batch.sh $file
done

#rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/physicell/simulation_data/2D_neut_bdg_motility/diversified_100/to_beast_format/*.log ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/physicell/logs/2D_neut_bdg_motility


cat SDevo.jar > physicell/simulation_data/2D_neut_bdg_sigmoid/diversified_100/to_beast_format/SDevo.jar
cat scripts/submit_beast_run2.sh > physicell/simulation_data/2D_neut_bdg_sigmoid/diversified_100/to_beast_format/submit_beast_run2.sh

cd physicell/simulation_data/2D_neut_bdg_sigmoid/diversified_100/to_beast_format


for file in xml_files/*
do
  sbatch submit_beast_run2.sh $file
  #qsub -cwd -l mem_free=2G submit_single_batch.sh $file
done

cat SDevo.jar > physicell/simulation_data/2D_sel_bdg_deleterious/diversified_100/to_beast_format/SDevo.jar
cat scripts/submit_beast_run2.sh > physicell/simulation_data/2D_sel_bdg_deleterious/diversified_100/to_beast_format/submit_beast_run2.sh

cd physicell/simulation_data/2D_sel_bdg_deleterious/diversified_100/to_beast_format


for file in xml_files/*
do
  sbatch scripts/submit_beast_run2.sh $file
  #qsub -cwd -l mem_free=2G submit_single_batch.sh $file
done

cat SDevo.jar > physicell/simulation_data/3D_sel/diversified_100/to_beast_format/SDevo.jar
cat scripts/submit_beast_run2.sh > physicell/simulation_data/3D_sel/diversified_100/to_beast_format/submit_beast_run2.sh

cd physicell/simulation_data/3D_sel/diversified_100/to_beast_format


#TODO: Run xml files (transfer first to gs cluster)
for file in xml_files/*
do
  sbatch submit_beast_run2.sh $file
  #qsub -cwd -l mem_free=2G submit_single_batch.sh $file

done

cat SDevo.jar > physicell/simulation_data/2D_neut_bdg_usecarryingcapacity/diversified_100/to_beast_format/SDevo.jar
cat scripts/submit_beast_run2.sh > physicell/simulation_data/2D_neut_bdg_usecarryingcapacity/diversified_100/to_beast_format/submit_beast_run2.sh

cd physicell/simulation_data/2D_neut_bdg_usecarryingcapacity/diversified_100/to_beast_format


for file in xml_files/*
do
  sbatch scripts/submit_beast_run2.sh $file
  #qsub -cwd -l mem_free=2G submit_single_batch.sh $file
done

cat SDevo.jar > physicell/simulation_data/2D_sel_bdg/diversified_100/to_beast_format/SDevo.jar
cat scripts/submit_beast_run2.sh >physicell/simulation_data/2D_sel_bdg/diversified_100/to_beast_format/submit_beast_run2.sh

cd physicell/simulation_data/2D_sel_bdg/diversified_100/to_beast_format


for file in xml_files/*
do
  sbatch submit_beast_run2.sh $file
  #qsub -cwd -l mem_free=2G submit_single_batch.sh $file
done


#To get all posteriors
for file in physicell/simulation_data/*/diversified_100/to_beast_format/*.log
do
  echo $file
  sbatch scripts/run_get_posteriors_physicell.sh $file
done

#Summarize posteriors

for file in physicell/simulation_data/2D_neut_bdg_motility/diversified_100/to_beast_format/*.log
do

  echo $file
  sbatch scripts/run_get_posteriors_physicell.sh $file
done

for file in physicell/simulation_data/2D_neut_bdg_sigmoid/diversified_100/to_beast_format/*.log
do
  sbatch scripts/run_get_posteriors_physicell.sh $file
done

for file in physicell/simulation_data/2D_sel_bdg_deleterious/diversified_100/to_beast_format/*.log
do
  sbatch scripts/run_get_posteriors_physicell.sh $file
done

for file in physicell/simulation_data/3D_sel/diversified_100/to_beast_format/*.log
do
  sbatch scripts/run_get_posteriors_physicell.sh $file
done

for file in physicell/simulation_data/2D_neut_bdg_usecarryingcapacity/diversified_100/to_beast_format/*.log
do
  sbatch scripts/run_get_posteriors_physicell.sh $file
done

for file in physicell/simulation_data/2D_neut_bdg_usecarryingcapacity/diversified_100/to_beast_format/*.log
do
  sbatch scripts/run_get_posteriors_physicell.sh $file
done

head -n1 sampconfig_m0.5_w1_d0.2_t1_mg1_mm1_l2e+08_i1_s3755_diversified_m1_n100_posterior_summary.tsv > posterior_summary.tsv

#Add all rows starting after row 2 (should be only 2 rows per file, but could combine multiple simulations)
for file in *_posterior_summary.tsv
do
    echo $file
    tail -n+2 $file >> posterior_summary.tsv
done

rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/physicell/stats/posteriors ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/physicell/stats

#Get true growth rate differences
for file in physicell/simulation_data/*/diversified_100/to_beast_format/xml_files/*.xml
do
  echo $file
  sbatch scripts/get_true_rates_physicell_submit.sh  $file
done

#Rscript scripts/get_true_rates_physicell.Rscript physicell/simulation_data/2D_neut_bdg_motility/diversified_100/to_beast_format/xml_files/sampconfig_m0.5_w1_d0.2_t1_mg1_mm1_l2e+08_i1_s3755_diversified_m1_n100.xml

head -n1 sampconfig_m0.5_w1_d0.2_t1_mg1_mm1_l2e+08_i1_s3755_diversified_m1_n100_rates.tsv > true_birth_rates_extended_physicell.tsv

#Add all rows starting after row 2 (should be only 2 rows per file, but could combine multiple simulations)
for file in *_rates.tsv
do
    echo $file
    tail -n+2 $file >> true_birth_rates_extended_physicell.tsv
done

rsync -a mlewinso@rhino:/fh/fast/bedford_t/users/mlewinsohn/tumors_sims/spatial-tumor-phylodynamics/physicell/stats/true_birth_rates ~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/physicell/stats


for N in 1000 2000 3000 4000 5000
do
  echo $N
  for itr in 1 2 3
    do
      echo $itr
      sbatch submit_indiv_eden_sim_N.sh 0.100 $N $itr
    done
done
