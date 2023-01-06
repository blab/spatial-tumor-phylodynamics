ssh lewinsom@nexus.gs.washington.edu
ssh lewinsom@grid.gs.washington.edu
cd lewinsom@nexus.gs.washington.edu:/net/gs/vol1/home/lewinsom/StateClocks
rsync -a /Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/bulk/xml_files/* lewinsom@nexus.gs.washington.edu:/net/gs/vol1/home/lewinsom/StateClocks/eden/bulk/xml
rsync -a /Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/SDevo.jar lewinsom@nexus.gs.washington.edu:/net/gs/vol1/home/lewinsom/StateClocks

rsync -a lewinsom@nexus.gs.washington.edu:/net/gs/vol1/home/lewinsom/StateClocks/eden/bulk/xml/* eden/bulk/xml
for file in xml/*
do
  echo $file
  qsub -cwd submit_single_batch.sh $file
done

rsync -a xml lewinsom@nexus.gs.washington.edu:/net/gs/vol1/home/lewinsom/StateClocks/physicell/simulation_data/2D_neut_bdg_motility
rsync -a xml lewinsom@nexus.gs.washington.edu:/net/gs/vol1/home/lewinsom/StateClocks/physicell/simulation_data/2D_neut_bdg_sigmoid
rsync -a xml lewinsom@nexus.gs.washington.edu:/net/gs/vol1/home/lewinsom/StateClocks/physicell/simulation_data/2D_sel_bdg
rsync -a xml lewinsom@nexus.gs.washington.edu:/net/gs/vol1/home/lewinsom/StateClocks/physicell/simulation_data/3D_sel
rsync -a xml lewinsom@nexus.gs.washington.edu:/net/gs/vol1/home/lewinsom/StateClocks/physicell/simulation_data/2D_neut_bdg

for file in xml_to_run/*
do
  echo $file
  qsub ../resume_single_batch.sh $file
done


rsync -a lewinsom@nexus.gs.washington.edu:/net/gs/vol1/home/lewinsom/StateClocks/eden/bulk_sampling_death_rate_validation_pop_1000_dr_*_strict_clock.xml /Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/bulk/logs

for file in death_rate_validation_pop_1000_dr_0.[0-9][0-9]_n_100_state_clock_estimate_dr_random_sampling.xml
do
    echo $file
    sh create_strict_clock.sh $file
done

for file in death_rate_validation_pop_1000_dr_0.[0-9][0-9]_n_100_state_clock_estimate_dr_random_sampling_strict_clock.xml
do
    echo $file
    sh convert_to_new_version.sh $file
done

for file in death_rate_validation_pop_1000_dr_0.[0-9][0-9]_n_1[1-9]_state_clock_estimate_dr.xml
do
    echo $file
    sh create_strict_clock.sh $file
done

for file in death_rate_validation_pop_1000_dr_0.[0-9][0-9]_n_1[1-9]_state_clock_estimate_dr.xml
do
    echo $file
    sh convert_to_new_version.sh $file
done

for file in death_rate_validation_pop_1000_dr_0.[0-9][0-9]_n_1[1-9]_state_clock_estimate_dr_strict_clock.xml
do
    echo $file
    sh convert_to_new_version.sh $file
done

for file in death_rate_validation_pop_1000_dr_0.[0-9][0-9]_n_[1-9]_state_clock_estimate_dr.xml
do
    echo $file
    sh create_strict_clock.sh $file
done

for file in death_rate_validation_pop_1000_dr_0.[0-9][0-9]_n_[1-9]_state_clock_estimate_dr_strict_clock.xml
do
    echo $file
    sh convert_to_new_version.sh $file
done

rsync -a lewinsom@nexus.gs.washington.edu:/net/gs/vol1/home/lewinsom/StateClocks/physicell/simulation_data/2D_neut_bdg_motility/xml .
rsync -a lewinsom@nexus.gs.washington.edu:/net/gs/vol1/home/lewinsom/StateClocks/physicell/simulation_data/2D_neut_bdg_motility/*.log .
rsync -a lewinsom@nexus.gs.washington.edu:/net/gs/vol1/home/lewinsom/StateClocks/physicell/simulation_data/2D_neut_bdg/*.log .
rsync -a lewinsom@nexus.gs.washington.edu:/net/gs/vol1/home/lewinsom/StateClocks/physicell/simulation_data/2D_sel_bdg/*.log .
rsync -a lewinsom@nexus.gs.washington.edu:/net/gs/vol1/home/lewinsom/StateClocks/physicell/simulation_data/3D_sel/*.log .
rsync -a lewinsom@nexus.gs.washington.edu:/net/gs/vol1/home/lewinsom/StateClocks/physicell/simulation_data/2D_neut_bdg_sigmoid/*.log .
rsync -a xml_to_run lewinsom@nexus.gs.washington.edu:/net/gs/vol1/home/lewinsom/StateClocks/physicell/simulation_data/3D_sel
rsync -a xml_to_run lewinsom@nexus.gs.washington.edu:/net/gs/vol1/home/lewinsom/StateClocks/physicell/simulation_data/2D_neut_bdg
rsync -a xml_to_run lewinsom@nexus.gs.washington.edu:/net/gs/vol1/home/lewinsom/StateClocks/physicell/simulation_data/2D_sel_bdg
rsync -a xml_to_run lewinsom@nexus.gs.washington.edu:/net/gs/vol1/home/lewinsom/StateClocks/physicell/simulation_data/3D_neut_bdg
rsync -a xml_to_run lewinsom@nexus.gs.washington.edu:/net/gs/vol1/home/lewinsom/StateClocks/physicell/simulation_data/2D_neut_bdg_sigmoid
rsync -a ../2D_neut_bdg_motility/diversified_100/updated_growth_rate_differences diversified_100
2D_neut_bdg_SGE
