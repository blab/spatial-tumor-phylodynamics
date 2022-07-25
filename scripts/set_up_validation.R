#### set_up_validation.R #### 

#Script to generate XML files from simulation output

library(tumortree)
library(tidyverse)

#estimate clock rate

#sim_cells_file_test <- sim_cells_files[1]
dr_vec <- format(round(seq(0,0.85,by=0.01),2), nsmall = 2)
#sample_sizes <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100))
sample_sizes <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)

#sample_sizes <- c(50)
sim_cells_files <- paste0("outputs/raw_simulation_results/validation/",
                         list.files(path = "outputs/raw_simulation_results/validation",
                                    pattern = "cells_death_rate_validation_pop_1000_dr_0.[0-9][05].*"), sep = "")

sim_cells_files <- sim_cells_files[! str_detect(sim_cells_files, "sampled")]
sim_cells_files <- sim_cells_files[! str_detect(sim_cells_files,"_i_")]
xml_template <- "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/xml_files/state_clocks_template.xml"
xml_template_single_dr <- "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/xml_files/state_clocks_single_dr_template.xml"
#sampled_cells_files <- paste("outputs/raw_simulation_results/validation/sampled_cells_death_rate_validation_pop_1000_dr_", dr_vec, "_n_", sample_size, "_diversified.csv", sep = "")



set.seed(31241)
overwrite = FALSE
for (sample_size in sample_sizes) {
    for (i in 1:length(sim_cells_files)) {
    
    sim_cells_file <- sim_cells_files[i]
        

        dr_extract <- regmatches(basename(sim_cells_file),
                                 gregexpr("(?<=dr_)[[:digit:]]+\\.[[:digit:]][[:digit:]]", basename(sim_cells_file), perl = TRUE))[[1]]
        print(dr_extract)
        xml_file <- paste("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/xml_files/death_rate_validation_pop_1000_dr_", dr_extract,
                           "_n_", sample_size, "_state_clock_estimate_dr.xml", sep = "")
        
        xml_file_single_dr <- paste("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/xml_files/death_rate_validation_pop_1000_dr_", dr_extract,
                          "_n_", sample_size, "_state_clock_estimate_single_dr.xml", sep = "")
        
        sampled_cells_file <- paste("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/raw_simulation_results/validation/sampled_cells_death_rate_validation_pop_1000_dr_", dr_extract, "_n_", sample_size, "_diversified.csv", sep = "")
        
        if (file.exists(xml_file) & (! overwrite)) {
            
            message("XML file exists, skipping....")
            
        } else {
            
            write_state_clocks_xml(sim_cells_file = sim_cells_file,
                           sampled_cells_file = sampled_cells_file,
                           xml_file = xml_file,
                           template_file =  xml_template,
                           diversified = TRUE,
                           resample = FALSE,
                           n_samples = sample_size)
            
            write_state_clocks_xml(sim_cells_file = sim_cells_file,
                                   sampled_cells_file = sampled_cells_file,
                                   xml_file = xml_file_single_dr,
                                   template_file =  xml_template_single_dr,
                                   diversified = TRUE,
                                   resample = FALSE,
                                   n_samples = sample_size)
        }
    }
}

#random sampling
set.seed(21221)
overwrite = TRUE

#sample_sizes <- c(50)
sim_cells_files <- paste0("outputs/raw_simulation_results/validation/",
                          list.files(path = "outputs/raw_simulation_results/validation",
                                     pattern = "cells_death_rate_validation_pop_1000_dr_0.[0-9][0-9].*"), sep = "")

sim_cells_files <- sim_cells_files[! str_detect(sim_cells_files, "sampled")]
sim_cells_files <- sim_cells_files[! str_detect(sim_cells_files,"_i_")]
xml_template <- "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/xml_files/state_clocks_template.xml"
#sampled_cells_files <- paste("outputs/raw_simulation_results/validation/sampled_cells_death_rate_validation_pop_1000_dr_", dr_vec, "_n_", sample_size, "_diversified.csv", sep = "")




for (i in 1:length(sim_cells_files)) {
        
        sim_cells_file <- sim_cells_files[i]
        sample_size <- 100
        
        dr_extract <- regmatches(basename(sim_cells_file),
                                 gregexpr("(?<=dr_)[[:digit:]]+\\.[[:digit:]][[:digit:]]", basename(sim_cells_file), perl = TRUE))[[1]]
        print(dr_extract)
        xml_file <- paste("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/xml_files/death_rate_validation_pop_1000_dr_", dr_extract,
                          "_n_", sample_size, "_state_clock_estimate_dr_random_sampling.xml", sep = "")
        
        sampled_cells_file <- paste("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/raw_simulation_results/validation/sampled_cells_death_rate_validation_pop_1000_dr_", dr_extract, "_n_", sample_size, "_random_sampling.csv", sep = "")
        
        if (file.exists(xml_file) & (! overwrite)) {
            
            message("XML file exists, skipping....")
            
        } else {
            
            write_state_clocks_xml(sim_cells_file = sim_cells_file,
                                   sampled_cells_file = sampled_cells_file,
                                   xml_file = xml_file,
                                   template_file =  xml_template,
                                   diversified = FALSE,
                                   resample = TRUE,
                                   n_samples = sample_size)
            

    }
}


#extended power analysis for N <20 
set.seed(5774)
overwrite = FALSE
sim_cells_files <- paste0("outputs/raw_simulation_results/validation/",
                          list.files(path = "outputs/raw_simulation_results/validation",
                                     pattern = "cells_death_rate_validation_pop_1000_dr_0.[0-9][05].*"), sep = "")

sim_cells_files <- sim_cells_files[! str_detect(sim_cells_files, "sampled")]
sim_cells_files <- sim_cells_files[! str_detect(sim_cells_files,"_i_")]
sample_sizes <- 5:19
for (sample_size in sample_sizes) {
    for (i in 1:length(sim_cells_files)) {
        
        sim_cells_file <- sim_cells_files[i]
        
        
        dr_extract <- regmatches(basename(sim_cells_file),
                                 gregexpr("(?<=dr_)[[:digit:]]+\\.[[:digit:]][[:digit:]]", basename(sim_cells_file), perl = TRUE))[[1]]
        print(dr_extract)
        xml_file <- paste("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/xml_files/death_rate_validation_pop_1000_dr_", dr_extract,
                          "_n_", sample_size, "_state_clock_estimate_dr.xml", sep = "")
        
        # xml_file_single_dr <- paste("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/xml_files/death_rate_validation_pop_1000_dr_", dr_extract,
        #                             "_n_", sample_size, "_state_clock_estimate_single_dr.xml", sep = "")
        
        sampled_cells_file <- paste("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/raw_simulation_results/validation/sampled_cells_death_rate_validation_pop_1000_dr_", dr_extract, "_n_", sample_size, "_diversified.csv", sep = "")
        
        if (file.exists(xml_file) & (! overwrite)) {
            
            message("XML file exists, skipping....")
            
        } else {
            
            write_state_clocks_xml(sim_cells_file = sim_cells_file,
                                   sampled_cells_file = sampled_cells_file,
                                   xml_file = xml_file,
                                   template_file =  xml_template,
                                   diversified = TRUE,
                                   resample = FALSE,
                                   n_samples = sample_size)
            
            # write_state_clocks_xml(sim_cells_file = sim_cells_file,
            #                        sampled_cells_file = sampled_cells_file,
            #                        xml_file = xml_file_single_dr,
            #                        template_file =  xml_template_single_dr,
            #                        diversified = TRUE,
            #                        resample = FALSE,
            #                        n_samples = sample_size)
        }
    }
}


