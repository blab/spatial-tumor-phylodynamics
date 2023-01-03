#!/bin/bash

# For all: run at least twice -- as needed to run until ESS threshold (200) and then will output posterior summaries
## Run in eden directory

## Run all death rates for N=50 (Figure 3)

## State-dependent clock
for file in xml/death_rate_validation_pop_1000_dr_0.[0-9][0-9]_n_50_state_clock_estimate_dr.xml
do
  sbatch ../scripts/run_beast_to_ess2.sh $file
done

## Strict clock
#Done!
for file in xml/death_rate_validation_pop_1000_dr_0.[0-9][0-9]_n_50_state_clock_estimate_dr_strict_clock.xml
do
  sbatch ../scripts/run_beast_to_ess2.sh $file
done

# Run all death rates for N=100 (Figure S4)

## State-depedent clock
for file in xml/death_rate_validation_pop_1000_dr_0.[0-9][0-9]_n_100_state_clock_estimate_dr.xml
do
  sbatch ../scripts/run_beast_to_ess2.sh $file
done


## Strict clock
for file in xml/death_rate_validation_pop_1000_dr_0.[0-9][0-9]_n_100_state_clock_estimate_dr_strict_clock.xml
do
  sbatch ../scripts/run_beast_to_ess2.sh $file
done

## State-depedent clock -- random sampling
for file in xml/death_rate_validation_pop_1000_dr_0.[0-9][0-9]_n_100_state_clock_estimate_dr_random_sampling.xml
do
  sbatch ../scripts/run_beast_to_ess2.sh $file
done

## Strict clock -- random sampling

for file in xml/death_rate_validation_pop_1000_dr_0.[0-9][0-9]_n_100_state_clock_estimate_dr_random_sampling_strict_clock.xml
do
  sbatch ../scripts/run_beast_to_ess2.sh $file
done

## Run all death rates for N=60, 70, 80, 90 strict clock (Figure S3)
for file in xml/death_rate_validation_pop_1000_dr_0.[0-9][0-9]_n_[6789]0_state_clock_estimate_dr_strict_clock.xml
do
  sbatch ../scripts/run_beast_to_ess2.sh $file
done


## Run subset death rates for all Ns

### State-dependent, N < 10
for file in xml/death_rate_validation_pop_1000_dr_0.[0-9][05]_n_[0-9]_state_clock_estimate_dr.xml
do
  sbatch ../scripts/run_beast_to_ess2.sh $file
done


### Strict clock, N < 10
for file in xml/death_rate_validation_pop_1000_dr_0.[0-9][05]_n_[0-9]_state_clock_estimate_dr_strict_clock.xml
do
  sbatch ../scripts/run_beast_to_ess2.sh $file
done


### State-dependent, N >= 10

for file in xml/death_rate_validation_pop_1000_dr_0.[0-9][05]_n_[1-4][0-9]_state_clock_estimate_dr.xml
do
  sbatch ../scripts/run_beast_to_ess2.sh $file
done


for file in xml/death_rate_validation_pop_1000_dr_0.[0-9][05]_n_[6-8][0-9]_state_clock_estimate_dr.xml
do
  sbatch ../scripts/run_beast_to_ess2.sh $file
done



for file in xml/death_rate_validation_pop_1000_dr_0.[0-9][05]_n_90_state_clock_estimate_dr.xml
do
  sbatch ../scripts/run_beast_to_ess2.sh $file
done

for file in xml/death_rate_validation_pop_1000_dr_0.[0-1][05]_n_40_state_clock_estimate_dr.xml
do
  sbatch ../scripts/run_beast_to_ess2.sh $file
done

for file in xml/death_rate_validation_pop_1000_dr_0.[0-1][05]_n_60_state_clock_estimate_dr.xml
do
  sbatch ../scripts/run_beast_to_ess2.sh $file
done

for file in xml/death_rate_validation_pop_1000_dr_0.[0-1][05]_n_70_state_clock_estimate_dr.xml
do
  sbatch ../scripts/run_beast_to_ess2.sh $file
done

for file in xml/death_rate_validation_pop_1000_dr_0.[0-1][05]_n_80_state_clock_estimate_dr.xml
do
  sbatch ../scripts/run_beast_to_ess2.sh $file
done

for file in xml/death_rate_validation_pop_1000_dr_0.[0-1][05]_n_90_state_clock_estimate_dr.xml
do
  sbatch ../scripts/run_beast_to_ess2.sh $file
done

#To run
sbatch ../scripts/run_beast_to_ess2.sh xml/death_rate_validation_pop_1000_dr_0.20_n_17_state_clock_estimate_dr.xml

### Strict N >= 10

for file in xml/death_rate_validation_pop_1000_dr_0.[0-9][05]_n_[1-4][0-9]_state_clock_estimate_dr_strict_clock.xml
do
  sbatch ../scripts/run_beast_to_ess2.sh $file
done

## Run in eden/bulk directory
##### BULK SAMPLING ##############
cd bulk
for file in xml/bulk_sampling_death_rate_validation_pop_1000_dr_0.[0-8][0-9].xml
do
  sbatch ../../scripts/run_beast_to_ess_bulk.sh $file
done

## Physicell runs ####

for file in xml/sampconfig_m0_w1*
do
  sbatch ../../../scripts/run_beast_to_ess_physicell.sh $file
done

for file in xml/sampconfig_m0_w1_d0.6_*
do
  sbatch ../../../scripts/run_beast_to_ess_physicell.sh $file
done

#sbatch ../../../scripts/run_beast_to_ess_physicell.sh xml/sampconfig_m0_w1.1_d0.1_t1_mg0.99_mm1_l2e+08_i8_s50357_diversified_m1_n100.xml