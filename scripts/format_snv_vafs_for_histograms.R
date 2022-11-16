#format_snv_vafs_for_histograms.R
library(tumortree)
library(tidyverse)

#for now this needs the raw data which is not public
setwd("/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li")

##Whole genome estimated to be ~3Gb
###41% GC content

###DATA PROCESSING ####

#Published edge/center labels from Li et al (Table S8)
t1_li_edge_labels <- data.frame("sample" = c("t1z5", "t1l13", "t1z1", "t1f24", "t1z3", "t1f23", "t1f11", "t1l8", "t1f14", "t1f9", "t1l1", "t1l10", "t1l3", "t1f2", "t1l6", "t1f4"),
                                "edgeP" = c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))


t2_li_edge_labels <- data.frame("sample" = c("t2z11", "t2z13", "t2z1", "t2f9", "t2f5", "t2f2", "t2f13", "t2z9", "t2z6"),
                                "edgeP" = c(1, 1, 1, 0, 0, 0, 0, 0, 0))

# tumor_purities <- read_csv("../li-application/li_hcc_data/li_table_s2_purities.csv") %>% 
#     dplyr::mutate(Sample = tolower(Sample))
#colnames(tumor_purities) <- c("sample", "mapping_rate", "depth", "n_snvs", "purity", "ploidy")

#Extracted from last rows of Tables S5 and S6 in Li et al Supplement
#T1
# t1_truncal_mutations <- read_csv("../li-application/li_hcc_data/t1_truncal_mutations.csv", col_names = FALSE) 
# colnames(t1_truncal_mutations) <- c("contig", "position")
# t1_truncal_mutations <- t1_truncal_mutations %>% 
#     add_column("mark" = "truncal")
#T2
# t2_truncal_mutations <- read_csv("../li-application/li_hcc_data/t2_truncal_mutations.csv", col_names = FALSE)
# colnames(t2_truncal_mutations) <- c("contig", "position")
# t2_truncal_mutations <- t2_truncal_mutations %>% 
#     add_column("mark" = "truncal")
#tumor 1
#note: files are named w/ .csv extension but are actually tab-delimited

#Local directories
t1_wgs_snvs_files <- list.files(path = "data/All_SNV_MBE2022Li", pattern = "t1", full.names = TRUE)

#tumor 2
t2_wgs_snvs_files <- list.files(path = "data/All_SNV_MBE2022Li", pattern = "t2", full.names = TRUE)

header_info <- read.table(file = "data/All_SNV_MBE2022Li/head.csv")

#Uncomment when data is added 
#t1_wgs_snvs_files <- list.files(path = "../li-appliation/li_hcc_data/All_SNV_MBE2022Li", pattern = "t1", full.names = TRUE)

#tumor 2
#t2_wgs_snvs_files <- list.files(path = "../li-appliation/li_hcc_data/All_SNV_MBE2022Li", pattern = "t2", full.names = TRUE)

#header_info <- read.table(file = "../li-appliation/li_hcc_data/data/All_SNV_MBE2022Li/head.csv")

read in raw WGS variants for each punch

##tumor 1
t1_snv_list <- list()
t1_all_sites_list <- list()

for (f1 in t1_wgs_snvs_files) {
    
    sample_extract <- gsub(".tsv", "", gsub(".csv", "", basename(f1)))
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
    
    sample_extract <- gsub(".tsv", "", gsub(".csv", "", basename(f2)))
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
    mutate("mut_id" = paste0(contig, position, sep = "_")) %>% 
    distinct

t2_all_sites_df <- t2_all_sites_list %>%
    bind_rows %>%
    mutate("mut_id" = paste0(contig, position, sep = "_")) %>% 
    distinct

#Read in masked sites
source("~/Documents/PhD/Bedford_lab/spatial-tumor-phylodynamics/li-application/hp.R")

t1_all_sites_df$masked <- "no"
t1_all_sites_df$masked[t1_muts] <- "yes"

t1_all_sites_df_masked <- t1_all_sites_df %>% 
    dplyr::filter(masked == "yes") 

t2_all_sites_df$masked <- "no"
t2_all_sites_df$masked[t2_muts] <- "yes"

t2_all_sites_df_masked <- t2_all_sites_df %>% 
    dplyr::filter(masked == "yes") 
# #function to make vafs for wgs
make_vafs <- function(punch_snvs, all_sites, vaf_cutoff = 0.1) {
    
    punch_snvs_vafs <- punch_snvs %>%
        dplyr::select(contig, position, ref_allele, alt_allele, t_ref_count, t_alt_count) %>%
        left_join(all_sites, .) %>%
        dplyr::mutate(t_alt_vaf = t_alt_count /(t_ref_count + t_alt_count)) %>%
        dplyr::mutate(basecall = ifelse(is.na(t_alt_vaf), ref_allele, ifelse(t_alt_vaf <= vaf_cutoff, ref_allele, alt_allele))) %>%
        dplyr::mutate(variant = ifelse(basecall == alt_allele, 1, 0)) %>%
        tibble::add_column(sample = punch_snvs$sample[1] )
    
    
    return(punch_snvs_vafs)
    
}

vaf_cutoff = 0.05
t1_vaf_list <- purrr::map(t1_snv_list,
                          function(x) make_vafs(punch_snvs = x,
                                                all_sites = t1_all_sites_df,
                                                vaf_cutoff = vaf_cutoff))
t2_vaf_list <- purrr::map(t2_snv_list,
                          function(x) make_vafs(punch_snvs = x,
                                                all_sites = t2_all_sites_df,
                                                vaf_cutoff = vaf_cutoff))

t1_vaf_df <- t1_vaf_list %>%
    bind_rows %>% 
    left_join(.,t1_li_edge_labels)
     dplyr::filter(mark == "truncal")

t2_vaf_df <- t2_vaf_list %>%
    bind_rows %>% 
    left_join(., t2_li_edge_labels) 

t1_vaf_df_truncal <- t1_vaf_df %>% 
    group_by(contig, position) %>% 
    summarize(n_samples = sum(variant, na.rm = TRUE)) %>% 
    dplyr::filter(n_samples == length(t1_snv_list)) %>% 
    add_column("mark" = "truncal")

t2_vaf_df_truncal <- t2_vaf_df %>% 
    group_by(contig, position) %>% 
    summarize(n_samples = sum(variant, na.rm = TRUE)) %>% 
    dplyr::filter(n_samples == length(t2_snv_list)) %>% 
    add_column("mark" = "truncal")

t1_vaf_df <- t1_vaf_df %>% 
    left_join(t1_vaf_df_truncal)

t2_vaf_df <- t2_vaf_df %>% 
    left_join(t2_vaf_df_truncal)


#tumor 1
## write data to regenerate plot
t1_vaf_freqs <- t1_vaf_df %>% 
    dplyr::filter(t_alt_vaf >= 0.05) %>% 
    dplyr::select(contig, position, edgeP, mark, sample, t_alt_vaf) %>% 
    mutate("mut_id" = paste0(contig, position, sep = "_")) %>% 
    mutate('masked' = mut_id %in% t1_all_sites_df_masked$mut_id)



t2_vaf_freqs <- t2_vaf_df %>% 
    dplyr::filter(t_alt_vaf >= 0.05) %>% 
    dplyr::select(contig, position, edgeP, mark, sample, t_alt_vaf) %>% 
    mutate("mut_id" = paste0(contig, position, sep = "_")) %>% 
    mutate('masked' = mut_id %in% t2_all_sites_df_masked$mut_id)

write_tsv(t1_vaf_freqs,
          file = "../li-application/li_hcc_data/t1_vaf_freqs.tsv")

write_tsv(t2_vaf_freqs,
          file = "../li-application/li_hcc_data/t2_vaf_freqs.tsv")

