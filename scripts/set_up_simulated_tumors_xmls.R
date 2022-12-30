#### set_up_simulated_tumor_xmls.R #### 

#Script to generate XML files from tumor simulations

library(tumortree)
library(tidyverse)




#dr_vec <- format(round(seq(0,0.85,by=0.01),2), nsmall = 2)

sample_sizes <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)

#Get file names for simulation data
## These data are outputs of simulation_data/spatial_tumor_simulation.ipynb

sim_cells_files <- list.files(path = "../simulation_data",
                                    pattern = "cells_death_rate_validation_pop_1000_dr_0.[0-9][05].*",
                              full.names = TRUE)

# Filter to only get range of simulations accross dr=0-0.87 and not other simulation studies for not
sim_cells_files <- sim_cells_files[! str_detect(sim_cells_files, "sampled")]
sim_cells_files <- sim_cells_files[! str_detect(sim_cells_files,"_i_")]

xml_template <- "analysis/xml/state-clocks-template.xml"

set.seed(31241)

#set to TRUE if want to regenerated XML files
overwrite = FALSE
for (sample_size in sample_sizes) {
    for (i in 1:length(sim_cells_files)) {
    
    sim_cells_file <- sim_cells_files[i]
        

        dr_extract <- regmatches(basename(sim_cells_file),
                                 gregexpr("(?<=dr_)[[:digit:]]+\\.[[:digit:]][[:digit:]]", basename(sim_cells_file), perl = TRUE))[[1]]
        print(dr_extract)
        xml_file <- paste("../analysis/xml/death_rate_validation_pop_1000_dr_", dr_extract,
                           "_n_", sample_size, "_state_clock_estimate_dr.xml", sep = "")
        
        #save sampled cells as record
        sampled_cells_file <- paste("../simulation_data/sampled_cells_death_rate_validation_pop_1000_dr_",
                                    dr_extract, "_n_",
                                    sample_size,
                                    "_diversified.csv", sep = "")
        
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
            
        }
    }
}

# Also generate files for random sampling
set.seed(21221)

#Set overwrite to TRUE if want to regenerate XMLS that already exist
overwrite = FALSE

for (i in 1:length(sim_cells_files)) {
        
        sim_cells_file <- sim_cells_files[i]
        
        # All random sampling at 100 for now
        sample_size <- 100
        
        dr_extract <- regmatches(basename(sim_cells_file),
                                 gregexpr("(?<=dr_)[[:digit:]]+\\.[[:digit:]][[:digit:]]", basename(sim_cells_file), perl = TRUE))[[1]]
        print(dr_extract)
        xml_file <- paste("../analysis/xml/death_rate_validation_pop_1000_dr_", dr_extract,
                          "_n_", sample_size, "_state_clock_estimate_dr_random_sampling.xml", sep = "")
        
        sampled_cells_file <- paste("../simulation_data/validation/sampled_cells_death_rate_validation_pop_1000_dr_",
                                    dr_extract,
                                    "_n_",
                                    sample_size,
                                    "_random_sampling.csv",
                                    sep = "")
        
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

sample_sizes <- 5:19
for (sample_size in sample_sizes) {
    for (i in 1:length(sim_cells_files)) {
        
        sim_cells_file <- sim_cells_files[i]
        
        
        dr_extract <- regmatches(basename(sim_cells_file),
                                 gregexpr("(?<=dr_)[[:digit:]]+\\.[[:digit:]][[:digit:]]", basename(sim_cells_file), perl = TRUE))[[1]]
        print(dr_extract)
        xml_file <- paste("../analysis/xml/death_rate_validation_pop_1000_dr_",
                          dr_extract,
                          "_n_",
                          sample_size,
                          "_state_clock_estimate_dr.xml", sep = "")
        
        # xml_file_single_dr <- paste("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/validation/xml_files/death_rate_validation_pop_1000_dr_", dr_extract,
        #                             "_n_", sample_size, "_state_clock_estimate_single_dr.xml", sep = "")
        
        sampled_cells_file <- paste("../simulation_data/sampled_cells_death_rate_validation_pop_1000_dr_",
                                    dr_extract, "_n_", sample_size,
                                    "_diversified.csv", sep = "")
        
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
            
        }
    }
}


