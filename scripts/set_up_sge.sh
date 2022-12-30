ssh lewinsom@nexus.gs.washington.edu
ssh lewinsom@grid.gs.washington.edu

rsync -a /Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/bulk/xml_files/* lewinsom@nexus.gs.washington.edu:/net/gs/vol1/home/lewinsom/StateClocks/eden/bulk/xml
rsync -a /Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/SDevo.jar lewinsom@nexus.gs.washington.edu:/net/gs/vol1/home/lewinsom/StateClocks

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

for file in xml/*
do
  echo $file
  qsub -cwd ../submit_single_batch.sh $file
done


rsync -a lewinsom@nexus.gs.washington.edu:/net/gs/vol1/home/lewinsom/StateClocks/eden/bulk_sampling_death_rate_validation_pop_1000_dr_*_strict_clock.xml /Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/bulk/logs
