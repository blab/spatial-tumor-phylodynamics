#process_li_wgs_data_with_amplicon.R

library(tumortree)
library(tidyverse)
set.seed(812)
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

t2_all_sites_df <- t2_all_sites_list %>%
    bind_rows %>%
    distinct %>% 
    arrange(contig, position)


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
    dplyr::filter(Punch != "t1f47.1")

T2_data_amp_long <- T2_data_amp %>%
    tidyr::pivot_longer(cols = contains("T2"), names_to = c("Punch"), values_to = "VAF") %>%
    dplyr::mutate("basecall" = ifelse(is.na(VAF), "N", ifelse(VAF>= min_var_allele_freq, Alt, Ref))) %>%
    dplyr::filter(! is.na(Position)) %>%
    dplyr::mutate(Punch = tolower(Punch)) %>%
    dplyr::rename(contig = "Chromosome", position =  "Position", ref_allele = "Ref", alt_allele = "Alt", t_alt_vaf = "VAF")


#function to make vafs for wgs
make_vafs <- function(punch_snvs, punch_name, all_sites, vaf_cutoff = 0.05) {

    punch_snvs_vafs <- punch_snvs %>%
        dplyr::select(contig, position, ref_allele, alt_allele, t_ref_count, t_alt_count) %>%
        left_join(all_sites, .) %>%
        dplyr::mutate(t_alt_vaf = t_alt_count /(t_ref_count + t_alt_count)) %>%
        dplyr::mutate(basecall = ifelse(is.na(t_alt_vaf), ref_allele, ifelse(t_alt_vaf <= vaf_cutoff, ref_allele, alt_allele))) %>%
        dplyr::mutate(variant = ifelse(basecall == alt_allele, 1, 0)) %>%
        tibble::add_column(sample = punch_snvs$sample[1], Punch = punch_name)


    return(punch_snvs_vafs)

}


###### MAKE SEQUENCES #####
T1_all_amp_punches <- unique(T1_data_amp_long$Punch)
T2_all_amp_punches <- unique(T2_data_amp_long$Punch)

get_wgs_sequence_from_vaf <- function(punch, wgs_snvs_vaf_df, include_invariable = TRUE, amplicon_only = FALSE)  {
    
    punch_snvs <- wgs_snvs_vaf_df %>% 
        dplyr::filter(Punch == punch)
    
    if (include_invariable) {
        wgs_sequence <- paste(c(punch_snvs$basecall,"N", "A", "G", "C", "T"), collapse = "")
    
    } else if (amplicon_only) {
        
        wgs_sequence <- paste(punch_snvs$basecall, collapse = "")
    
    } else {
        wgs_sequence <- paste(c(punch_snvs$basecall,"N"), collapse = "")
    }

    return(wgs_sequence)
}

get_wgs_seq_from_amp_panel <- function(punch, amp_long_df, wgs_all_sites, include_invariable = TRUE, amplicon_only = FALSE) {
    
    
    punch_all_sites <- amp_long_df %>% 
        dplyr::filter(Punch == punch) %>% 
        dplyr::left_join(wgs_all_sites, ., by = c("contig", "position", "ref_allele")) %>% 
        dplyr::mutate("basecall" = ifelse(is.na(basecall), "N", basecall), "Punch" = punch)
    
    #sanity check 
    #print(sum(punch_all_sites$basecall != "N"))
    if (include_invariable) {
        
        amp_sequence <- paste(c(punch_all_sites$basecall,"N", "N", "N", "N", "N"), collapse = "")
        
    } else if (amplicon_only) {
        
        amp_sequence <- paste(punch_all_sites$basecall, collapse = "")
    } else {
        
        amp_sequence <- paste(c(punch_all_sites$basecall,"N"), collapse = "")
    }
    
    return(amp_sequence)
}

get_sequence <- function(punch, wgs_snvs_vaf_df, amp_long_df, wgs_all_sites, include_invariable = TRUE, amplicon_only = FALSE) {
    
    if (punch %in% unique(wgs_snvs_vaf_df$Punch) & (! amplicon_only)) {
        
        sequence <- get_wgs_sequence_from_vaf(punch, wgs_snvs_vaf_df,
                                              include_invariable = include_invariable,
                                              amplicon_only = amplicon_only) 
        
    } else {
        
        sequence <- get_wgs_seq_from_amp_panel(punch, amp_long_df,
                                               wgs_all_sites,
                                               include_invariable = include_invariable,
                                               amplicon_only = amplicon_only)
    }
    
    return(sequence)
}


generate_subsampled_sequences <- function(all_amp_punches, wgs_snv_list, amp_long_df, wgs_all_sites, vaf_cutoff = 0.1,
                                          n_sites = 10000,
                                          include_invariable = TRUE,
                                          amplicon_only = FALSE, 
                                          sites_file=NULL) {
    
    if (amplicon_only) {
        wgs_subset_sites <- amp_long_df %>% 
            
            dplyr::distinct(contig, position, ref_allele) %>% 
            arrange(contig,position)
        
        if(! is.null(sites_file)) {
          write_csv(wgs_subset_sites, sites_file)
        }
        
        unsampled_sites <- 0
        
    #} else if (! is.na(n_amplicon_sampled_sites)) {
        
       #  #first break sites in amplified or wgs-only 
       #  
       #  
       #  ##make mut_id
       #  
       #  amp_sites_mut_id_df <- amp_long_df %>% 
       #      tidyr::unite("mut_id", c(contig, position), remove = FALSE) %>% 
       #      dplyr::select(mut_id, contig, position, ref_allele) %>% 
       #      dplyr::distinct()
       #      
       #  
       #  wgs_all_mut_id_sites <- wgs_all_sites %>% 
       #      tidyr::unite("mut_id", c(contig, position), remove = FALSE)
       #  
       #  amp_sites_mut_id_df <- amp_sites_mut_id_df %>% 
       #      dplyr::filter(mut_id %in% wgs_all_mut_id_sites$mut_id)
       # 
       #  
       #  ## get ids for sites in amplicon panel and not
       #  amp_mut_ids <- amp_sites_mut_id_df$mut_id
       #  
       #  
       #  wgs_not_amp_sites <-  wgs_all_mut_id_sites %>% 
       #      dplyr::filter(! mut_id %in% amp_mut_ids)
       #  
       #  wgs_not_amp_ids <- wgs_not_amp_sites$mut_id
       #  
       #  #then sample proportional to bias
       #  if (n_amplicon_sampled_sites > length(amp_mut_ids)) {
       #      
       #      warning("Number of desired amplicon sites is greater than the total, using all amplicon sites")
       #      
       #      n_amplicon_sampled_sites <- length(amp_mut_ids)
       #      
       #  } 
       #  
       #  #sample desired number of amplicon sites
       #  amp_sampled_ids  <- sample(amp_mut_ids, size = n_amplicon_sampled_sites, replace = FALSE)
       #  
       # # base::print(n_sites)
       #  n_wgs_sampled <- n_sites - n_amplicon_sampled_sites
       #  
       #  #base::print(n_wgs_sampled)
       #  #sample remaining sites from wgs
       #  wgs_sampled_ids <- sample(wgs_not_amp_ids, size = n_wgs_sampled, replace = FALSE)
       #  
       #  #base::print(length(wgs_sampled_ids))
       #  all_sampled_mut_ids <- c(amp_sampled_ids, wgs_sampled_ids)
       #  
       #  #base::print(length(all_sampled_mut_ids))
       #  
       #  wgs_subset_sites <- wgs_all_mut_id_sites %>% 
       #      dplyr::filter(mut_id %in% all_sampled_mut_ids) %>% 
       #      dplyr::select(-mut_id)
       #  
       #  base::print(nrow(wgs_subset_sites))
       #  
       #  unsampled_sites <- nrow(wgs_all_sites) - nrow(wgs_subset_sites)
        
      # base::print(unsampled_sites)
        
    } else {
        wgs_subset_sites <- wgs_all_sites[sort(sample(1:nrow(wgs_all_sites), size = n_sites, replace = FALSE)), ]
        wgs_subset_sites <- wgs_subset_sites %>% 
          
          dplyr::distinct(contig, position, ref_allele) %>% 
          arrange(contig,position)
        
        print(nrow(wgs_subset_sites))
        
        if(! is.null(sites_file)) {
          
          write_csv(wgs_subset_sites, sites_file)
          
        }
        unsampled_sites <- nrow(wgs_all_sites) - nrow(wgs_subset_sites)

    }
    
    
    vaf_list <- purrr::map2(wgs_snv_list, names(wgs_snv_list),
                              function(x, n) make_vafs(punch_snvs = x,
                                                       punch_name = n,
                                                    all_sites = wgs_subset_sites,
                                                    vaf_cutoff = vaf_cutoff))
    
    
    vaf_df <- vaf_list %>%
        bind_rows
    

    
    sequences <- purrr::map(all_amp_punches, function(p) get_sequence(punch = p, 
                                                                            wgs_snvs_vaf_df = vaf_df,
                                                                            amp_long_df = amp_long_df,
                                                                            wgs_all_sites = wgs_subset_sites,
                                                                      include_invariable = include_invariable,
                                                                      amplicon_only = amplicon_only))
    
    

    
    return(list("sequences" = sequences, "n_unsampled_sites" = unsampled_sites))
    
}

# n_sites = 10000
# t1_sequences <- generate_subsampled_sequences(all_amp_punches = T1_all_amp_punches,
#                                              wgs_snv_list = t1_snv_list,
#                                              amp_long_df =  T1_data_amp_long,
#                                              wgs_all_sites = t1_all_sites_df,
#                                              n_sites = n_sites,
#                                              n_amplicon_sampled_sites = 300)
# 
# 
# t1_sequences_all <- generate_subsampled_sequences(all_amp_punches = T1_all_amp_punches,
#                                               wgs_snv_list = t1_snv_list,
#                                               amp_long_df =  T1_data_amp_long,
#                                               wgs_all_sites = t1_all_sites_df,
#                                               n_sites = nrow(t1_all_sites_df),
#                                               sites_file =  "data/t1_all_sites.csv",
#                                               include_invariable = TRUE)

t1_sequences_amp_sites_only_wgs_punches <- generate_subsampled_sequences(all_amp_punches = names(t1_snv_list),
                                                  wgs_snv_list = t1_snv_list,
                                                  amp_long_df =  T1_data_amp_long,
                                                  wgs_all_sites = t1_all_sites_df,
                                                  n_sites = nrow(t1_all_sites_df),
                                                  amplicon_only = TRUE,
                                                  include_invariable = FALSE,
                                                  sites_file = "data/t1_amp_sites.csv")

t1_sequences_amp_sites_only <- generate_subsampled_sequences(all_amp_punches = T1_all_amp_punches,
                                                             wgs_snv_list = t1_snv_list,
                                                             amp_long_df =  T1_data_amp_long,
                                                             wgs_all_sites = t1_all_sites_df,
                                                             n_sites = nrow(t1_all_sites_df),
                                                             amplicon_only = TRUE,
                                                             include_invariable = FALSE)


# t2_sequences_all <- generate_subsampled_sequences(all_amp_punches = T2_all_amp_punches,
#                                                   wgs_snv_list = t2_snv_list,
#                                                   amp_long_df =  T2_data_amp_long,
#                                                   wgs_all_sites = t2_all_sites_df,
#                                                   n_sites = nrow(t2_all_sites_df),
#                                                   sites_file =  "data/t2_all_sites.csv", 
#                                                   include_invariable = TRUE)

t2_sequences_amp_sites_only_wgs_punches <- generate_subsampled_sequences(all_amp_punches = names(t2_snv_list),
                                                                         wgs_snv_list = t2_snv_list,
                                                                         amp_long_df =  T2_data_amp_long,
                                                                         wgs_all_sites = t2_all_sites_df,
                                                                         n_sites = nrow(t2_all_sites_df),
                                                                         amplicon_only = TRUE,
                                                                         include_invariable = FALSE,
                                                                         sites_file = "data/t2_amp_sites.csv")

t2_sequences_amp_sites_only <- generate_subsampled_sequences(all_amp_punches = T2_all_amp_punches,
                                                                         wgs_snv_list = t2_snv_list,
                                                                         amp_long_df =  T2_data_amp_long,
                                                                         wgs_all_sites = t2_all_sites_df,
                                                                         n_sites = nrow(t2_all_sites_df),
                                                                         amplicon_only = TRUE,
                                                                         include_invariable = FALSE)



# #write fasta file
# fasta_file <- "/Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/data/li_wgs_amp.fa"
# #
# first_line = TRUE
# for (i in 1:length(T1_all_amp_punches)) {
# 
#     punch <- T1_all_amp_punches[i]
# 
#     write(paste0(">", punch, collapse = ""), file=fasta_file, append=!first_line)
# 
#     first_line = FALSE
# 
#     write(t1_sequences_all$sequences[[i]], file=fasta_file, append=!first_line)
# 
# 
# }


### write fasta file with all punches from wgs tree with only amplicon sites
# 
# t1_amp_sites <- T1_data_amp_long %>% 
#   
#   dplyr::distinct(contig, position, ref_allele) %>% 
#   arrange(contig,position)

t1_fasta_file_wgs_punches_amp_sites <- "/Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/data/li_t1_wgs_punches_amp_sites.fa"
first_line = TRUE
for (i in 1:length(names(t1_snv_list))) {
  
  punch <- names(t1_snv_list)[i]
  
  write(paste0(">", punch, collapse = ""), file=t1_fasta_file_wgs_punches_amp_sites, append=!first_line)
  
  first_line = FALSE
  
  write(t1_sequences_amp_sites_only_wgs_punches$sequences[[i]], file=t1_fasta_file_wgs_punches_amp_sites, append=!first_line)

}
write(">blood", file=t1_fasta_file_wgs_punches_amp_sites, append=!first_line)
write(paste0(t1_amp_sites$ref_allele, collapse=""), file=t1_fasta_file_wgs_punches_amp_sites, append=!first_line)


# t2_amp_sites <- T2_data_amp_long %>% 
#   
#   dplyr::distinct(contig, position, ref_allele) %>% 
#   arrange(contig,position)

t2_fasta_file_wgs_punches_amp_sites <- "/Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/data/li_t2_wgs_punches_amp_sites.fa"
first_line = TRUE
for (i in 1:length(names(t2_snv_list))) {
  
  punch <- names(t2_snv_list)[i]
  
  write(paste0(">", punch, collapse = ""), file=t2_fasta_file_wgs_punches_amp_sites, append=!first_line)
  
  first_line = FALSE
  
  write(t2_sequences_amp_sites_only_wgs_punches$sequences[[i]], file=t2_fasta_file_wgs_punches_amp_sites, append=!first_line)
  
}

# add blood
write(">blood", file=t2_fasta_file_wgs_punches_amp_sites, append=!first_line)
write(paste0(t2_amp_sites$ref_allele, collapse=""), file=t2_fasta_file_wgs_punches_amp_sites, append=!first_line)






#sanity checks
#are all punches included?
# length(t1_sequence) == length(T1_all_amp_punches)
# 
# #are all sequences equal to desired length?
# all(unlist(purrr::map(t1_sequence, nchar)) == n_sites + 4)
# 
# #count how many Ns are in each sequence
# nNs <- unlist(purrr::map(t1_sequence, function(s) str_count(s, pattern = "N")))
# sum(nNs > 1000) == length(T1_all_amp_punches) - length(t1_snv_list)
# sum(nNs < 1000) == length(t1_snv_list)
# 
# t1_sequences_df <- purrr::map2(t1_sequences$sequences, T1_all_amp_punches, function(s, p) data.frame("mut_id" = paste("mut_", 1:nchar(s), sep = ""),
#                                                                     "basepair"= unlist(str_split(s, pattern = "")), 
#                                                                     "Punch" = p)) %>% 
#     bind_rows %>% 
#     left_join(., subclone_labels, by = "Punch") %>% 
#     mutate("Subclone" = as.factor(Subclone))
# 
# 
# t1_sequences_df_summary <- t1_sequences_df%>% 
#     group_by(mut_id) %>% 
#     summarise("Ns" = sum(basepair== "N")) %>% 
#     filter(Ns < 100)
# 
# t1_punches_df_summary <- t1_sequences_df%>% 
#     group_by(Punch) %>% 
#     summarise("Ns" = sum(basepair == "N")) %>% 
#     filter(Ns < 4000)
# 
# amplicon_panel <- t1_sequences_df_summary$mut_id
# wgs_punches <- t1_punches_df_summary$Punch

#mut_subset <- unique(t1_sequences_df$mut_id)[1:500]

# t1_sequences_df %>% 
#     filter(mut_id %in% amplicon_panel) %>% 
#     ggplot(., aes(x= as.factor(mut_id), y = Punch, fill = basepair)) +geom_tile() + scale_fill_brewer(type = "qual")
# 
# t1_sequences_df %>% 
#     filter(Punch %in%  wgs_punches, mut_id %in% mut_subset) %>% 
#     ggplot(., aes(x= as.factor(mut_id), y = Punch, fill = basepair)) +geom_tile() + scale_fill_brewer(type = "qual")
# 
# t1_sequences_df %>% 
#     filter( mut_id %in% mut_subset) %>% 
#     ggplot(., aes(x= as.factor(mut_id), y = Punch, fill = basepair)) +geom_tile() + scale_fill_brewer(type = "qual")

# hist(nNs)
# t1_sequences <- purrr::map(T1_all_amp_punches, function(p) get_sequence(punch = p, 
#                                                                         wgs_snvs_vaf_df = t1_vaf_df,
#                                                                         amp_long_df = T1_data_amp_long,
#                                                                         wgs_all_sites = t1_all_sites_df))
# 
# t2_sequences <- purrr::map(T2_all_amp_punches, function(p) get_sequence(punch = p, 
#                                                                         wgs_snvs_vaf_df = t2_vaf_df,
#                                                                         amp_long_df = T2_data_amp_long,
#                                                                         wgs_all_sites = t2_all_sites_df))


###### ASSIGN STATES #####
#based on Li et al supplement 
# t1_li_edge_labels <- data.frame("punch" = c("t1z5", "t1l13", "t1z1", "t1f24", "t1z3", "t1f23", "t1f11", "t1l8", "t1f14", "t1f9", "t1l1", "t1l10", "t1l3", "t1f2", "t1l6", "t1f4"),
#                                  "edge" = c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
# 
# 
# t2_li_edge_labels <- data.frame("punch" = c("t2z11", "t2z13", "t2z1", "t2f9", "t2f5", "t2f2", "t2f13", "t2z9", "t2z6"),
# "edge" = c(1, 1, 1, 0, 0, 0, 0, 0, 0))


#Spatial data from Figure 2 tumor schematics converted to coordinates by WebDigiziter2 
## These data contain both punches and boundary of each tumor slice
## Within each slice, punch coordinates are ordered by punch number
##(units are in mm, but coordinates are arbitrary)
T1_spatial_data <- read.csv("data/T1_spatial_data.csv", header = TRUE)
T2_spatial_data <- read.csv("data/T2_spatial_data.csv", header = TRUE)


# T1F_slices <- c("T1F2","T1F4", "T1F9", "T1F11", "T1F14", "T1F23", "T1F24", "T1F25", "T1F27", "T1F28", "T1F29", "T1F30", "T1F31", "T1F32", "T1F33", "T1F34", "T1F35", "T1F36",
# "T1F37", "T1F38", "T1F39", "T1F40", "T1F41", "T1F42", "T1F43", "T1F44", "T1F45", "T1F46", "T1F47", "T1F48")

create_punch_df <- function(slice, T_spatial_data) {
    #print(slice)
    #punches <- colnames(T_data)[grep(pattern = slice, x = colnames(T_data))]
    # if(slice == "T1F") {
    #     ordered_punches <- T1F_slices
    # } else {
    #     ordered_punches <- punches[order(as.integer(gsub(pattern = slice, replacement = "", fixed = TRUE, x = punches)))]
    # }
    slice_data <- T_spatial_data %>%
        dplyr::select(contains(slice) & !contains("boundary")) %>%
        na.omit %>%
        tibble::add_column("Slice" = slice)
    
    colnames(slice_data) <- c("Punch", "X", "Y", "Slice")
    return(slice_data)
    
}

separate_boundary_points <- function(slice, T_spatial_data) {
    #print(slice)
    
    slice_data <- T_spatial_data %>%
        dplyr::select(contains(slice) & contains("boundary")) %>%
        na.omit
    
    #order points in circle
    
    curr_point <- slice_data[1,1:2]
    ordered_slice_data <- curr_point
    remaining_points <- slice_data[2:nrow(slice_data),1:2]
    
    for (i in 2:nrow(slice_data)) {
        new_row <- which.min(raster::pointDistance(p1 = curr_point, p2 = remaining_points[,1:2], lonlat = FALSE))[[1]]
        curr_point <- remaining_points[new_row ,1:2]
        ordered_slice_data <- ordered_slice_data %>%
            add_row(tibble::tibble_row(curr_point))
        remaining_points <- remaining_points[-new_row, 1:2]
    }
    
    ordered_slice_data <- ordered_slice_data %>%
        tibble::add_column("Slice" = slice)
    
    colnames(ordered_slice_data) <- c("X", "Y","Slice")
    return(ordered_slice_data)
    
}

calculate_circumference_from_boundary_points <- function(slice, T_spatial_data) {
    #print(slice)
    
    slice_data <- T_spatial_data %>%
        dplyr::select(contains(slice) & contains("boundary")) %>%
        na.omit
    
    #order points in circle
    
    curr_point <- slice_data[1,1:2]
    ordered_slice_data <- curr_point
    remaining_points <- slice_data[2:nrow(slice_data),1:2]
    
    circ <- 0
    
    for (i in 2:nrow(slice_data)) {
        new_row <- which.min(raster::pointDistance(p1 = curr_point, p2 = remaining_points[,1:2], lonlat = FALSE))[[1]]
        between_point_dist <- min(raster::pointDistance(p1 = curr_point, p2 = remaining_points[,1:2], lonlat = FALSE))[[1]]
        curr_point <- remaining_points[new_row ,1:2]
        
        remaining_points <- remaining_points[-new_row, 1:2]
        circ <- circ + between_point_dist
    }
    
    return(circ)
    
}



T1_slices <- unique(gsub("\\..*", replacement = "", colnames(T1_spatial_data)[grep(pattern = "T1", x = colnames(T1_spatial_data))]))
T1_slices <- T1_slices[! grepl("_boundary", T1_slices)]

T2_slices <- unique(gsub("\\..*", replacement = "", colnames(T2_spatial_data)[grep(pattern = "T2", x = colnames(T2_spatial_data))]))
T2_slices <- T2_slices[! grepl("_boundary", T2_slices)]

T1_punch_coordinates <- purrr::map(T1_slices, function(slice) create_punch_df(slice = slice, T_spatial_data = T1_spatial_data)) %>%
    bind_rows

T2_punch_coordinates <- purrr::map(T2_slices, function(slice) create_punch_df(slice = slice, T_spatial_data = T2_spatial_data)) %>%
  bind_rows



T1_boundary_coordinates <- purrr::map(T1_slices, function(slice) separate_boundary_points(slice = slice, T_spatial_data = T1_spatial_data)) %>%
    bind_rows

T2_boundary_coordinates <- purrr::map(T2_slices, function(slice) separate_boundary_points(slice = slice, T_spatial_data = T2_spatial_data)) %>%
  bind_rows

height_of_slice <- 0.2 #mm (200uM)
hypothetical_slices <- c(LETTERS, "AA", "AB")

T1_boundary_coordinates$z <- map_dbl(sub("T1", "",
                                      T1_boundary_coordinates$Slice), function(l) (length(hypothetical_slices) - which(l == hypothetical_slices)) * height_of_slice + height_of_slice)


# T1_boundary_circumference <- purrr::map_dbl(T1_slices, function(slice) calculate_circumference_from_boundary_points(slice = slice, T_spatial_data = T1_spatial_data))
# names(T1_boundary_circumference) <- T1_slices

T2_boundary_coordinates$z <- map_dbl(sub("T2", "",
                                         T2_boundary_coordinates$Slice),
                                     function(l) (length(hypothetical_slices) - which(l == hypothetical_slices)) * height_of_slice + height_of_slice)


# T2_boundary_circumference <- purrr::map_dbl(T2_slices, function(slice) calculate_circumference_from_boundary_points(slice = slice, T_spatial_data = T2_spatial_data))
# names(T2_boundary_circumference) <- T2_slices



find_closest_dist_to_boundary <- function(x,y, boundary_points) {
    
    return(min(raster::pointDistance(p1 = c(x,y), p2 = boundary_points, lonlat = FALSE)))
}

T1_punch_coordinates$dist_to_boundary <- map2_dbl(T1_punch_coordinates$X, T1_punch_coordinates$Y,
                                                  function(x,y) find_closest_dist_to_boundary(x,y, boundary_points = T1_boundary_coordinates[,1:2]))

T2_punch_coordinates$dist_to_boundary <- map2_dbl(T2_punch_coordinates$X, T2_punch_coordinates$Y,
                                                   function(x,y) find_closest_dist_to_boundary(x,y, boundary_points = T2_boundary_coordinates[,1:2]))
T1_edge_cutoff <- 2 #mm
T2_edge_cutoff <- 1.5 #mm
# T1_punch_coordinates <- T1_punch_coordinates %>%
#     mutate("edge" = dist_to_boundary <= edge_cutoff | Slice == "T1AB" | Slice == "T1Z") %>% 
#     mutate(Punch = tolower(Punch)) %>% 
#     dplyr::filter(Punch %in% T1_all_amp_punches)
T1_edge_slices <- unique(T1_boundary_coordinates$Slice[T1_boundary_coordinates$z < T1_edge_cutoff])
T2_edge_slices <- unique(T2_boundary_coordinates$Slice[T2_boundary_coordinates$z < T2_edge_cutoff])


T1_punch_coordinates <- T1_punch_coordinates %>%
    dplyr::mutate("edge" = dist_to_boundary <= T1_edge_cutoff | Slice %in% T1_edge_slices)# %>%
    #dplyr::mutate(Punch = tolower(Punch)) %>%
    #dplyr::filter(Punch %in% T1_all_amp_punches)
edge_center_colors <- c("edge" = "#89352F", "center" = "#A2D2E2")

T2_punch_coordinates <- T2_punch_coordinates %>%
  dplyr::mutate("edge" = dist_to_boundary < T2_edge_cutoff| Slice %in% T2_edge_slices)# %>%

t1_map_plot <- ggplot(T1_punch_coordinates, aes(X,-Y)) + geom_polygon(data = T1_boundary_coordinates, fill = "lightgrey", alpha = 0.5,  color = "black", aes(group = Slice)) +
    geom_point(size = 2, shape = 21, color = "black",  aes(fill = ifelse(edge, "edge", "center"))) +
    #scale_fill_manual(values = sim_colors ) +
    # labs(fill = "") +
    facet_wrap(~factor(Slice, levels = T1_slices), scales = "free") +
    theme_void() + theme(legend.position = "none") +
    scale_fill_manual(values = edge_center_colors)

t2_map_plot <- ggplot(T2_punch_coordinates, aes(X,-Y)) + geom_polygon(data = T2_boundary_coordinates, fill = "lightgrey", alpha = 0.5,  color = "black", aes(group = Slice)) +
  geom_point(size = 2, shape = 21, color = "black",  aes(fill = ifelse(edge, "edge", "center"))) +
  #scale_fill_manual(values = sim_colors ) +
  # labs(fill = "") +
  facet_wrap(~factor(Slice, levels = T2_slices), scales = "free") +
  theme_void() + theme(legend.position = "none") +
  scale_fill_manual(values = edge_center_colors)

    #geom_label(aes(label = Punch)) 
t1_map_plot + geom_label_repel(aes(label=Punch))
t2_map_plot+ geom_label(aes(label=Punch))

#Write spatial info for plotting
write_csv(T1_punch_coordinates, file = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/data/T1_punch_coordinates.csv")
write_csv(T2_punch_coordinates, file = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/data/T2_punch_coordinates.csv")

write_csv(T1_boundary_coordinates, file = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/data/T1_boundary_coordinates.csv")
write_csv(T2_boundary_coordinates, file = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/data/T2_boundary_coordinates.csv")
ggsave("/Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial_tumor_growth_simulation/manuscript/figures/li_t1_spatial_map.png",
       t1_map_plot)

ggsave("/Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial_tumor_growth_simulation/manuscript/figures/li_t2_spatial_map.png",
       t2_map_plot)
######### GET CLOCK RATE #############

# expected_origin_time <- 600 #days
# 
# t1_n_mutations <- data.frame(sample = names(t1_snv_list),
#                              n_mutations = unlist(purrr::map(t1_snv_list, nrow)))
# 
# t2_n_mutations <- data.frame(sample = names(t2_snv_list),
#                              n_mutations = unlist(purrr::map(t2_snv_list, nrow)))
# 
# t1_clock_rate <- mean(t1_n_mutations$n_mutations / expected_origin_time / nrow(t1_all_sites_df))
# t2_clock_rate <- mean(t2_n_mutations$n_mutations / expected_origin_time / nrow(t2_all_sites_df))
# 
# T1_clock_rate_string <- as.character(round(t1_clock_rate, 5))
# T2_clock_rate_string <- as.character(round(t2_clock_rate, 5))
# 
# ###### COALESCENT MODEL POP SIZ ##
# 
# t1_pop_size <- length(t1_sequences) / 0.0001
# t2_pop_size <- length(t2_sequences) / 0.0001
######## WRITE ALIGNMENT AND TAXON STRINGS ############

write_li_xml <- function(wgs_all_sites, amp_long_df, wgs_snv_list, all_amp_punches, punch_coordinates,
                         n_sites = 10000, boundary_coordinates, edge_cutoff = 2, n_sequenced_sites = 3E9,
                         state_clocks_template_xml = "state_clocks_template.xml", 
                         state_clocks_xml = "state_clocks.xml", include_invariable = TRUE,
                         amplicon_only = FALSE,  sites_file = NULL, amplicon_weights_file = NULL) {
    

    sequences <- generate_subsampled_sequences(all_amp_punches = all_amp_punches, 
                                                  wgs_snv_list = wgs_snv_list, 
                                                  amp_long_df =  amp_long_df, 
                                                  wgs_all_sites = wgs_all_sites, 
                                                  n_sites = n_sites,
                                               include_invariable = include_invariable,
                                               amplicon_only = amplicon_only,
                                               sites_file = sites_file)
    
    alignment_string <- ""

    for (i in 1:length(sequences$sequences)) {
        alignment_string <- paste0(alignment_string, paste0('\t', '<sequence id="seq_',
                                                            all_amp_punches[i], '" spec="Sequence" taxon="',
                                                            all_amp_punches[i], '" totalcount="4" value="',
                                                                sequences$sequences[i], '"/>', sep = ""),
                                   sep = "\n")
    }
    
    #taxon states
    
    #distance to boundary is only calculated by x-y coordinates, so also call based on Z slice
    edge_slices <- unique(boundary_coordinates$Slice[boundary_coordinates$z < edge_cutoff])
    punch_coordinates <- punch_coordinates %>%
        dplyr::mutate("edge" = dist_to_boundary < edge_cutoff | (Slice %in% edge_slices)) %>% 
        dplyr::mutate(Punch = tolower(Punch)) %>% 
        dplyr::filter(Punch %in% all_amp_punches)
    
    g <- ggplot(punch_coordinates, aes(X,-Y)) + geom_polygon(data = boundary_coordinates, fill = "lightgrey", alpha = 0.5,  color = "black", aes(group = Slice)) +
        geom_point(size = 2, shape = 21, color = "black",  aes(fill = ifelse(edge, "edge", "center"))) +
        theme_void() + theme(legend.position = "none") + ggtitle(paste0("Edge cutoff = ", edge_cutoff, "mm", sep = ""))
    
    base::print(g)
    
    taxon_traits <- c()
    
    for (i in 1:nrow(punch_coordinates)) {
        
        taxon_traits <- c(taxon_traits, paste0(tolower(punch_coordinates$Punch)[i], "=loc", as.integer(punch_coordinates$edge[i]), sep = ""))
    }
    
    taxon_traits_string <- paste(taxon_traits, collapse = ",")
    
    #weights string
    
    if (include_invariable) {
      #assume 60% GC content
      #assume human genome ~3Gb with 70% being a high confidence set
      ##https://www.pnas.org/doi/10.1073/pnas.1613365113
      n_wgs_sites_gc <- format(n_sequenced_sites * 0.7 * 0.4 * 0.5, scientific = FALSE)
      n_wgs_sites_at <- format(n_sequenced_sites * 0.7 * 0.6 * 0.5, scientific = FALSE)
      
      if (! is.null(amplicon_weights_file)) {
        #amplicon_weights_file <- "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/t1_amp_site_weights.csv"
        #amplicon weights from solve_branch_weights_li.R
        amplicon_weights <- read_csv(amplicon_weights_file) %>% 
          dplyr::select(snv_mut_id, amp_wgs_ratio, norm_ratio, branch_weight) 
        
        #just to double check that sites are going to be exactly the same
        #sites_file = "data/t1_all_sites_file.csv"
        sites <- read_csv(sites_file) %>% 
          dplyr::mutate(snv_mut_id = paste0(contig, "_", position)) %>% 
          dplyr::left_join(., amplicon_weights, by = "snv_mut_id")
        
        #non-amplicon sites have weight of 1
        sites$branch_weight[which(is.na(sites$branch_weight))] <- 1
        

        weights_string <- paste(c(sites$branch_weight, sequences$n_unsampled_sites, n_wgs_sites_at,
                                  n_wgs_sites_gc, n_wgs_sites_gc, n_wgs_sites_at), collapse = ",")
      } else {

        weights_string <- paste(c(rep("1", n_sites), sequences$n_unsampled_sites, n_wgs_sites_at,
                                  n_wgs_sites_gc, n_wgs_sites_gc, n_wgs_sites_at), collapse = ",")
        
      }
        
        
    } else if (amplicon_only) {
        
        weights_string <- NULL

    
    } else {
        
        weights_string <- paste(c(rep("1", n_sites), sequences$n_unsampled_sites), collapse = ",")
    }
    
    write_alignment_to_state_xml_from_template(state_clocks_template_xml = state_clocks_template_xml,
                                               alignment_string = alignment_string,
                                               taxon_traits_string = taxon_traits_string,
                                               state_clocks_xml = state_clocks_xml,
                                               weights_string = weights_string)
    
    
}

######## WRITE XML FILE FROM TEMPLATE ############


write_alignment_to_state_xml_from_template <- function(state_clocks_template_xml,
                                                       alignment_string,
                                                       taxon_traits_string,
                                                       weights_string, 
                                                       state_clocks_xml) {

    con = file(state_clocks_template_xml, "r")

    first_line = TRUE

    while(length(x <- readLines(con, n = 1)) > 0) {



        if (grepl("insert_sequence", x, fixed = TRUE)) { #alignment insertion

            write(alignment_string, file=state_clocks_xml, append=!first_line)

        } else if (grepl("insert_trait", x, fixed = TRUE)) { #taxon insertion

            write(gsub("insert_trait", taxon_traits_string , x), file=state_clocks_xml, append = !first_line)

        } else if (grepl("insert_weights", x, fixed = TRUE)) { #taxon insertion

            write(gsub("insert_weights", weights_string , x), file=state_clocks_xml, append = !first_line)

        } else {

            write(x, file=state_clocks_xml, append=!first_line)
        }

        first_line = FALSE


    }
    close(con)

}
#Files for templates and xml files
# state_clocks_template_xml <- "xml_files/T1_state_clocks_wgs_amplicon_template.xml"
# 
# li_T1_state_clocks_xml <- "xml_files/T1_state_clocks_wgs_amplicon.xml"


# t1_sequences <- generate_subsampled_sequences(all_amp_punches = T1_all_amp_punches, 
#                                               wgs_snv_list = t1_snv_list, 
#                                               amp_long_df =  T1_data_amp_long, 
#                                               wgs_all_sites = t1_all_sites_df, 
#                                               n_sites = n_sites)

# set.seed(2812)
# 
# edge_cutoff_vec <- c(1,2,3,4, 5)
# for (c in edge_cutoff_vec) {
#     
#     state_clocks_xml_file <- paste0("xml_files/T1_state_clocks_wgs_amplicon_c", c, ".xml", sep = "")
#     
#     write_li_xml( wgs_all_sites =t1_all_sites_df,
#               amp_long_df =T1_data_amp_long, 
#               wgs_snv_list = t1_snv_list,
#               all_amp_punches =T1_all_amp_punches,
#               punch_coordinates = T1_punch_coordinates,
#               boundary_coordinates = T1_boundary_coordinates,
#               n_sites = 10000, edge_cutoff = c,
#               n_sequenced_sites = 3E6,
#               state_clocks_template_xml = "xml_files/T1_state_clocks_test_temp.xml", 
#               state_clocks_xml = state_clocks_xml_file)
# }

# set.seed(33452)
# write_li_xml( wgs_all_sites =t1_all_sites_df,
#               amp_long_df =T1_data_amp_long, 
#               wgs_snv_list = t1_snv_list,
#               all_amp_punches =T1_all_amp_punches,
#               punch_coordinates = T1_punch_coordinates,
#               boundary_coordinates = T1_boundary_coordinates,
#               n_sites = 10000, edge_cutoff = 2,
#               n_sequenced_sites = 3E6,
#               state_clocks_template_xml = "xml_files/T1_state_clocks_test_temp.xml", 
#               state_clocks_xml = "xml_files/T1_state_clocks_variable_sites.xml",
#               include_invariable = FALSE)

#To make xmls with all sites + invariable + unsampled

write_li_xml(wgs_all_sites =t1_all_sites_df,
              amp_long_df =T1_data_amp_long, 
              wgs_snv_list = t1_snv_list,
              all_amp_punches =T1_all_amp_punches,
              punch_coordinates = T1_punch_coordinates,
              boundary_coordinates = T1_boundary_coordinates,
              n_sites = nrow(t1_all_sites_df),
              edge_cutoff = 2,
              n_sequenced_sites = 3E9,
              state_clocks_template_xml = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/xml_files/T1_state_clocks_temp.xml", 
              state_clocks_xml = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/xml_files/T1_state_clocks_all_sites.xml",
              include_invariable = TRUE,
              amplicon_only = FALSE,
              sites_file = "data/t1_all_sites_file.csv",
              amplicon_weights_file =  "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/t1_amp_site_weights.csv")


write_li_xml( wgs_all_sites =t2_all_sites_df,
              amp_long_df =T2_data_amp_long, 
              wgs_snv_list = t2_snv_list,
              all_amp_punches =T2_all_amp_punches,
              punch_coordinates = T2_punch_coordinates,
              boundary_coordinates = T2_boundary_coordinates,
              n_sites = nrow(t2_all_sites_df),
              n_sequenced_sites = 3E9,
              edge_cutoff = 1.5,
              state_clocks_template_xml = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/xml_files/T1_state_clocks_temp.xml", 
              state_clocks_xml = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/xml_files/T2_state_clocks_all_sites.xml",
              include_invariable = TRUE,
              amplicon_only = FALSE,
              sites_file = "data/t2_all_sites_file.csv",
              amplicon_weights_file =  "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/t2_amp_site_weights.csv")

set.seed(9182)
write_li_xml( wgs_all_sites =t1_all_sites_df,
              amp_long_df =T1_data_amp_long, 
              wgs_snv_list = t1_snv_list,
              all_amp_punches =T1_all_amp_punches,
              punch_coordinates = T1_punch_coordinates,
              boundary_coordinates = T1_boundary_coordinates,
              n_sites = NULL, edge_cutoff = 2,
              n_sequenced_sites = 3E6,
              state_clocks_template_xml = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/xml_files/T1_state_clocks_amplicon_only_temp.xml", 
              state_clocks_xml = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/xml_files/T1_state_clocks_amplicon_only.xml",
              include_invariable = FALSE,
              amplicon_only = TRUE,
              sites_file = "data/t1_amp_sites_file.csv")

set.seed(2122)
write_li_xml( wgs_all_sites =t2_all_sites_df,
              amp_long_df =T2_data_amp_long, 
              wgs_snv_list = t2_snv_list,
              all_amp_punches =T2_all_amp_punches,
              punch_coordinates = T2_punch_coordinates,
              boundary_coordinates = T2_boundary_coordinates,
              n_sites = NULL, edge_cutoff = 1.5,
              n_sequenced_sites = 3E6,
              state_clocks_template_xml = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/xml_files/T1_state_clocks_amplicon_only_temp.xml", 
              state_clocks_xml = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/xml_files/T2_state_clocks_amplicon_only.xml",
              include_invariable = FALSE,
              amplicon_only = TRUE,
              sites_file = "data/t2_amp_sites_file.csv")

##wgs punches only ####

write_li_xml(wgs_all_sites =t1_all_sites_df,
             amp_long_df =T1_data_amp_long, 
             wgs_snv_list = t1_snv_list,
             all_amp_punches =names(t1_snv_list),
             punch_coordinates = T1_punch_coordinates,
             boundary_coordinates = T1_boundary_coordinates,
             n_sites = nrow(t1_all_sites_df),
             edge_cutoff = 2,
             n_sequenced_sites = 3E9,
             state_clocks_template_xml = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/xml_files/T1_state_clocks_temp.xml", 
             state_clocks_xml = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/xml_files/T1_state_clocks_all_sites_wgs_punches.xml",
             include_invariable = TRUE,
             amplicon_only = FALSE,
             sites_file = "data/t1_all_sites_file_wgs_punches.csv",
             amplicon_weights_file =  "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/t1_amp_site_weights.csv")


write_li_xml( wgs_all_sites =t2_all_sites_df,
              amp_long_df =T2_data_amp_long, 
              wgs_snv_list = t2_snv_list,
              all_amp_punches =names(t2_snv_list),
              punch_coordinates = T2_punch_coordinates,
              boundary_coordinates = T2_boundary_coordinates,
              n_sites = nrow(t2_all_sites_df),
              n_sequenced_sites = 3E9,
              edge_cutoff = 1.5,
              state_clocks_template_xml = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/xml_files/T1_state_clocks_temp.xml", 
              state_clocks_xml = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/xml_files/T2_state_clocks_all_sites_wgs_punches.xml",
              include_invariable = TRUE,
              amplicon_only = FALSE,
              sites_file = "data/t2_all_sites_file_wgs_punches.csv",
              amplicon_weights_file =  "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/t2_amp_site_weights.csv")





# set.seed(31331)   
# state_clocks_xml_file <- paste0("xml_files/T1_state_clocks_wgs_amplicon_300_sites_2.xml", sep = "")
#     
# write_li_xml( wgs_all_sites =t1_all_sites_df,
#                   amp_long_df =T1_data_amp_long, 
#                   wgs_snv_list = t1_snv_list,
#                   all_amp_punches =T1_all_amp_punches,
#                   punch_coordinates = T1_punch_coordinates,
#                   boundary_coordinates = T1_boundary_coordinates,
#                   n_sites = 10000, edge_cutoff = 2,
#                   n_sequenced_sites = 3E6,
#                   state_clocks_template_xml = "xml_files/T1_state_clocks_test_temp.xml", 
#                   state_clocks_xml = state_clocks_xml_file,
#                 n_amplicon_sampled_sites = 300)

# write_alignment_to_state_xml_from_template(state_clocks_template_xml = state_clocks_template_xml,
#                                            alignment_string = T1_alignment_string,
#                                            taxon_traits_string = T1_taxon_traits_string,
#                                            state_clocks_xml = li_T1_state_clocks_xml,
#                                            weights_string = weights_string)



# subclone_labels <- read.csv("data/T1_subclone_labels.csv") %>% 
#     dplyr::mutate(Punch = tolower(Punch)) %>% 
#     dplyr::mutate("WGS" = Punch %in% wgs_punches)
# 
# write_csv(subclone_labels , file = "/Users/mayalewinsohn/Documents/PhD/Bedford_lab/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/subclone_labels.csv")
