#De-identify SNVs for publication
## Note: This can not be run with public repository
## Run process_li_wgs_data.R to generate XMLs from de-identified SNVs

library(tumortree)
library(tidyverse)

set.seed(812)

#Local repo
setwd("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li")

##Whole genome estimated to be ~3Gb
###41% GC content

###DATA PROCESSING ####

#tumor 1
t1_wgs_snvs_files <- list.files(path = "data/All_SNV_MBE2022Li", pattern = "t1", full.names = TRUE)

#tumor 2
t2_wgs_snvs_files <- list.files(path = "data/All_SNV_MBE2022Li", pattern = "t2", full.names = TRUE)

header_info <- read.table(file = "data/All_SNV_MBE2022Li/head.csv")

#note: files are named w/ .csv extension but are actually tab-delimited
subclone_labels <- read.csv("data/T1_subclone_labels.csv") %>% 
    dplyr::mutate(Punch = tolower(Punch))

#read in raw WGS variants for each punch

##tumor 1
t1_snv_list <- list()
t1_all_sites_list <- list()

for (f1 in t1_wgs_snvs_files) {
    
    sample_extract <- gsub(".csv", "", basename(f1))
    punch_extract <- gsub("t1", "", sample_extract)
    slice_extract <- gsub("[[:digit:]]","", punch_extract)
    punch_num_extract <- gsub("[a-z]","", punch_extract)
    
    # punch_snvs <- read.table(text = gsub("\t", ",", readLines(f1)), sep = ",")
    punch_snvs <- read.table(f1, header = FALSE, sep = "\t")
    colnames(punch_snvs) <- header_info
    punch_snvs <- punch_snvs %>%
        add_column(sample = sample_extract, punch = punch_extract, slice = slice_extract, punch_num = punch_num_extract)
    
    t1_snv_list[[sample_extract]] <- punch_snvs
    
    punch_sites <- punch_snvs %>%
        dplyr::filter(judgement == "KEEP") %>%
        dplyr::select(contig, position, ref_allele)
    
    t1_all_sites_list[[sample_extract]] <- punch_sites
    
}



#tumor 2
t2_snv_list <- list()
t2_all_sites_list <- list()

for (f2 in t2_wgs_snvs_files) {
    
    sample_extract <- gsub(".csv", "", basename(f2))
    punch_extract <- gsub("t2", "", sample_extract)
    slice_extract <- gsub("[[:digit:]]","", punch_extract)
    punch_num_extract <- gsub("[a-z]","", punch_extract)
    
    # punch_snvs <- read.table(text = gsub("\t", ",", readLines(f1)), sep = ",")
    punch_snvs <- read.table(f2, header = FALSE, sep = "\t")
    colnames(punch_snvs) <- header_info
    punch_snvs <- punch_snvs %>%
        add_column(sample = sample_extract, punch = punch_extract, slice = slice_extract, punch_num = punch_num_extract)
    
    t2_snv_list[[sample_extract]] <- punch_snvs
    
    punch_sites <- punch_snvs %>%
        dplyr::filter(judgement == "KEEP") %>%
        dplyr::select(contig, position, ref_allele)
    
    t2_all_sites_list[[sample_extract]] <- punch_sites
    
}

t1_all_sites_df <- t1_all_sites_list %>%
    bind_rows %>%
    distinct %>% 
    arrange(contig, position)

#Add de-identified variant id
t1_all_sites_df$var_id <- paste0("VAR", 1:nrow(t1_all_sites_df))

t2_all_sites_df <- t2_all_sites_list %>%
    bind_rows %>%
    distinct %>% 
    arrange(contig, position)

t2_all_sites_df$var_id <- paste0("VAR", 1:nrow(t2_all_sites_df))

#Amplicon variants
T1_data_amp <- read.csv("data/T1_snvs.csv", header = TRUE) %>% 
    dplyr::select(-T1F47.1)
T2_data_amp <- read.csv("data/T2_snvs.csv", header = TRUE)

min_var_allele_freq <- 0.05 #threshold to call presense of variant allele

#pivot data, this takes some wrangling because there are multiple columns per punch
T1_data_amp_long <- T1_data_amp %>%
    tidyr::pivot_longer(cols = contains("T1"), names_to = c("Punch"),values_to = "VAF")  %>%
    dplyr::mutate("basecall" = ifelse(is.na(VAF), "N", ifelse(VAF>= min_var_allele_freq, Alt, Ref))) %>%
    dplyr::filter(! is.na(Position)) %>% 
    dplyr::mutate(Punch = tolower(Punch))  %>% 
    dplyr::rename(contig = "Chromosome", position =  "Position", ref_allele = "Ref", alt_allele = "Alt", t_alt_vaf = "VAF") %>% 
    dplyr::filter(Punch != "t1f47.1") %>% 
    dplyr::left_join(., t1_all_sites_df, by = c("contig", "position", "ref_allele"))

T2_data_amp_long <- T2_data_amp %>%
    tidyr::pivot_longer(cols = contains("T2"), names_to = c("Punch"), values_to = "VAF") %>%
    dplyr::mutate("basecall" = ifelse(is.na(VAF), "N", ifelse(VAF>= min_var_allele_freq, Alt, Ref))) %>%
    dplyr::filter(! is.na(Position)) %>%
    dplyr::mutate(Punch = tolower(Punch)) %>%
    dplyr::rename(contig = "Chromosome", position =  "Position", ref_allele = "Ref", alt_allele = "Alt", t_alt_vaf = "VAF") %>% 
    dplyr::left_join(., t2_all_sites_df, by = c("contig", "position", "ref_allele"))


t1_all_sites_deident <- t1_all_sites_df %>% 
    dplyr::select(-contig, -position)

write_csv()
t2_all_sites_deident <- t2_all_sites_df %>% 
    dplyr::select(-contig, -position)

T1_data_amp_long_deident <- T1_data_amp_long %>% 
    dplyr::select(-contig, -position)

T2_data_amp_long_deident <- T2_data_amp_long %>% 
    dplyr::select(-contig, -position)

#Save to make publically available

write.csv(t1_all_sites_deident, "/Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/li-application/li_hcc_data/t1_all_sites_deident.csv")
write.csv(t2_all_sites_deident, "/Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/li-application/li_hcc_data/t2_all_sites_deident.csv")
                    
write.csv(T1_data_amp_long_deident, "/Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/li-application/li_hcc_data/T1_data_amp_long_deident.csv")
write.csv(T2_data_amp_long_deident, "/Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/li-application/li_hcc_data/T2_data_amp_long_deident.csv")


