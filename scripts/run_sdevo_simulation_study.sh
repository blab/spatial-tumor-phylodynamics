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


for file in xml/death_rate_validation_pop_1000_dr_0.[0-9][05]_n_[6-9][0-9]_state_clock_estimate_dr.xml
do
  sbatch ../scripts/run_beast_to_ess2.sh $file
done


### Strict N >= 10

for file in xml/death_rate_validation_pop_1000_dr_0.[0-9][05]_n_[1-4][0-9]_state_clock_estimate_dr_strict_clock.xml
do
  sbatch ../scripts/run_beast_to_ess2.sh $file
done
