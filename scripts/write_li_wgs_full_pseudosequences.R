#write_li_wgs_full_pseudosequences.R

library(tumortree)
library(tidyverse)
library(ggrepel)

#This analysis starts with the raw SNV files provided by Li et al. 
## Note: Not all files are accessible, but pseudosequences can be re-produced from 
## downstream saved SNV files re-labeled to hide contig and position info

##Whole genome estimated to be ~3Gb
###41% GC content

###DATA PROCESSING ####

#tumor 1
## Files not public skip if using github rep
t1_wgs_snvs_files <- list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/data/All_SNV_MBE2022Li", pattern = "t1", full.names = TRUE)

#tumor 2
## Files not public skip if using github rep
t2_wgs_snvs_files <- list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/data/All_SNV_MBE2022Li", pattern = "t2", full.names = TRUE)

## Files not public skip if using github rep
header_info <- read.table(file = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/data/All_SNV_MBE2022Li/head.csv")

#note: files are named w/ .csv extension but are actually tab-delimited

#subclone labels


# t1_subclone_colors <- c("alpha" = "#F4CAC7" , "alpha1" = "#F4ADAE", "alpha2" = "#EFA6B9", "beta" = "#F6E3B1", "beta1" = "#F7C0A2", "gamma"= "#B3D9B0", "other" = "grey")
# t1_subclone_labels <- data.frame("sample" = c("t1z5", "t1l13", "t1z1", "t1f24", "t1z3", "t1f23", "t1f11", "t1l8", "t1f14", "t1f9", "t1l1", "t1l10", "t1l3", "t1f2", "t1l6", "t1f4"),
#                                  "subclone" = c("gamma", "beta", "beta", "beta", "beta1", "beta1", "alpha2", "alpha2", "alpha2", "alpha2", "alpha1", "alpha1", "alpha1", "alpha1", "alpha1", "alpha1"))
# 
# t2_subclone_colors <- c("theta" = "#68AEDE","delta2"="#A3CDAA", "delta1" = "#8BC394", "other" = "#CBCACB")
# t2_subclone_labels <- data.frame("sample" = c("t2z11", "t2z13", "t2z1", "t2f9", "t2f5", "t2f2", "t2f13", "t2z9", "t2z6"),
#                                  "subclone" = c("other", "theta", "theta", "delta2", "delta2", "delta1", "delta1", "delta1", "delta1"))



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
    distinct

#add variant labels
t1_all_sites_df$var_id <- paste0("var", 1:nrow(t1_all_sites_df))

t2_all_sites_df <- t2_all_sites_list %>%
    bind_rows %>%
    distinct

t2_all_sites_df$var_id <- paste0("var", 1:nrow(t2_all_sites_df))

#punch_snvs <- t2_snv_list[[1]]
#all_sites <- t2_all_sites_df

make_vafs <- function(punch_snvs, all_sites, vaf_cutoff = 0.05) {

    
    punch_snvs_vafs <- punch_snvs %>%
        dplyr::select(contig, position, ref_allele, alt_allele, t_ref_count, t_alt_count) %>%
        left_join(all_sites, .) %>%
        dplyr::mutate(t_alt_vaf = t_alt_count /(t_ref_count + t_alt_count)) %>%
        dplyr::mutate(basecall = ifelse(is.na(t_alt_vaf), ref_allele, ifelse(t_alt_vaf <= vaf_cutoff, ref_allele, alt_allele))) %>%
        dplyr::mutate(variant = ifelse(basecall == alt_allele, 1, 0)) %>%
        tibble::add_column(sample = punch_snvs$sample[1] )


    return(punch_snvs_vafs)

}

t1_vaf_list <- purrr::map(t1_snv_list, function(x) make_vafs(punch_snvs = x, all_sites = t1_all_sites_df, vaf_cutoff = 0.05))
t2_vaf_list <- purrr::map(t2_snv_list, function(x) make_vafs(punch_snvs = x, all_sites = t2_all_sites_df, vaf_cutoff = 0.05))


###### MAKE SEQUENCES #####
get_sequence_from_vaf <- function(punch_snvs_vafs)  {
    sequence <- paste(punch_snvs_vafs$basecall, collapse = "")

    return(sequence)
}

t1_vaf_df <- t1_vaf_list %>%
    bind_rows %>%
    dplyr::select(var_id,t_alt_vaf, sample)



t2_vaf_df <- t2_vaf_list %>%
    bind_rows %>%
    dplyr::select(var_id,t_alt_vaf, sample)

write.csv(t1_vaf_df, "../li-application/li_hcc_data/t1_vaf_df.csv")
write.csv(t2_vaf_df, "../li-application/li_hcc_data/t2_vaf_df.csv")

t1_all_sites_df %>% 
    dplyr::select(-contig, -position) %>% 
    write.csv(., "../li-application/li_hcc_data/t1_all_sites_df.csv")

t2_all_sites_df %>% 
    dplyr::select(-contig, -position) %>% 
    write.csv(., "../li-application/li_hcc_data/t2_all_sites_df.csv")

###################### START HERE WITH PUBLIC DATA ######################

t1_vaf_df <- read.csv("../li-application/li_hcc_data/t1_vaf_df.csv")
t2_vaf_df <- read.csv("../li-application/li_hcc_data/t2_vaf_df.csv")

t1_sequences <- purrr::map(t1_vaf_list, function(x) get_sequence_from_vaf(x))
t2_sequences <- purrr::map(t2_vaf_list, get_sequence_from_vaf)


#write fasta file
fasta_file_t1 <- "../li-application/li_hcc_data/li_t1_wgs.fa"
fasta_file_t2 <- "../li-application/li_hcc_data/li_t2_wgs.fa"


write_fasta_file <- function(sequences_list, fasta_file, include_blood = FALSE, all_sites_df = NULL) {
    
    first_line = TRUE
    for (i in 1:length(sequences_list)) {
    
        punch <- names(sequences_list)[i]
    
        write(paste0(">", punch, collapse = ""), file=fasta_file, append=!first_line)
    
        first_line = FALSE
    
        write(sequences_list[[i]], file=fasta_file, append=!first_line)
    }
    
    if (include_blood) {
        
        write(">blood", file=fasta_file, append=!first_line)
        write(paste0(all_sites_df$ref_allele, collapse=""), file=fasta_file, append=!first_line)
    }
    
    
}

write_fasta_file(sequences_list = t1_sequences, fasta_file = fasta_file_t1, include_blood = TRUE, all_sites_df = t1_all_sites_df)
write_fasta_file(sequences_list = t2_sequences, fasta_file = fasta_file_t2, include_blood = TRUE, all_sites_df = t2_all_sites_df)

# write_fasta_file(sequences_list = t1_sequences, fasta_file = fasta_file_t1, include_blood = FALSE, all_sites_df = t1_all_sites_df)
# write_fasta_file(sequences_list = t2_sequences, fasta_file = fasta_file_t2, include_blood = FALSE, all_sites_df = t2_all_sites_df)
### Explore subclonal variation

t1_subclonal_summary_df <-  t1_vaf_df %>% 
    group_by(sample) %>%
    summarise("n_SNVs" = sum(t_alt_vaf > 0.02, na.rm = TRUE),
              "frac_subclonal" = sum(t_alt_vaf < 0.25 & t_alt_vaf > 0.02,
                                     na.rm = TRUE) / sum(t_alt_vaf >= 0.02, na.rm = TRUE)) %>% 
    left_join(., t1_subclone_labels, by = "sample")

ggplot(t1_subclonal_summary_df, aes(x = frac_subclonal, y = n_SNVs)) +
    geom_point(size = 2, aes(color = subclone)) + theme_classic() +
    geom_label_repel(data=subset(t1_subclonal_summary_df, frac_subclonal > 0.4 & n_SNVs >30000), aes(label = sample)) +
    scale_color_manual(values = t1_subclone_colors)


###### ASSIGN STATES #####
t1_li_edge_labels <- data.frame("punch" = c("t1z5", "t1l13", "t1z1", "t1f24", "t1z3", "t1f23", "t1f11", "t1l8", "t1f14", "t1f9", "t1l1", "t1l10", "t1l3", "t1f2", "t1l6", "t1f4"),
                                 "edge" = c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))


t2_li_edge_labels <- data.frame("punch" = c("t2z11", "t2z13", "t2z1", "t2f9", "t2f5", "t2f2", "t2f13", "t2z9", "t2z6"),
                                 "edge" = c(1, 1, 1, 0, 0, 0, 0, 0, 0))

T1_spatial_data <- read.csv("data/T1_spatial_data.csv", header = TRUE)
T2_spatial_data <- read.csv("data/T2_spatial_data.csv", header = TRUE)
######### GET CLOCK RATE #############

expected_origin_time <- 600 #days

t1_n_mutations <- data.frame(sample = names(t1_snv_list),
                             n_mutations = unlist(purrr::map(t1_snv_list, nrow)))

t2_n_mutations <- data.frame(sample = names(t2_snv_list),
                             n_mutations = unlist(purrr::map(t2_snv_list, nrow)))

t1_clock_rate <- mean(t1_n_mutations$n_mutations / expected_origin_time / nrow(t1_all_sites_df))
t2_clock_rate <- mean(t2_n_mutations$n_mutations / expected_origin_time / nrow(t2_all_sites_df))

T1_clock_rate_string <- as.character(round(t1_clock_rate, 5))
T2_clock_rate_string <- as.character(round(t2_clock_rate, 5))

###### COALESCENT MODEL POP SIZ ##

t1_pop_size <- length(t1_sequences) / 0.0001
t2_pop_size <- length(t2_sequences) / 0.0001
######## WRITE ALIGNMENT AND TAXON STRINGS ############
T1_alignment_string <- ""

for (i in 1:length(t1_sequences)) {
    T1_alignment_string <- paste0(T1_alignment_string, paste0('\t', '<sequence id="seq_',names(t1_sequences)[i], '" spec="Sequence" taxon="',names(t1_sequences)[i], '" totalcount="4" value="',
                                                              t1_sequences[i], '"/>', sep = ""), sep = "\n")
}

T2_alignment_string <- ""

for (i in 1:length(t2_sequences)) {
    T2_alignment_string <- paste0(T2_alignment_string, paste0('\t', '<sequence id="seq_',names(t2_sequences)[i], '" spec="Sequence" taxon="'
                                                              ,names(t2_sequences)[i], '" totalcount="4" value="',
                                                              t2_sequences[i], '"/>', sep = ""), sep = "\n")
}

T1_taxon_traits <- c()

for (i in 1:nrow(t1_li_edge_labels)) {

    T1_taxon_traits <- c(T1_taxon_traits, paste0(t1_li_edge_labels$punch[i], "=loc", t1_li_edge_labels$edge[i], sep = ""))
}

T1_taxon_traits_string <- paste(T1_taxon_traits, collapse = ",")

T2_taxon_traits <- c()

for (i in 1:nrow(t2_li_edge_labels)) {
    
    T2_taxon_traits <- c(T2_taxon_traits, paste0(t2_li_edge_labels$punch[i], "=loc", t2_li_edge_labels$edge[i], sep = ""))
}

T2_taxon_traits_string <- paste(T2_taxon_traits, collapse = ",")

######## WRITE XML FILE FROM TEMPLATE ############
write_alignment_to_state_xml_from_template <- function(state_clocks_template_xml,
                                                       alignment_string,
                                                       taxon_traits_string,
                                                       state_clocks_xml,
                                                       clock_rate_string) {

    con = file(state_clocks_template_xml, "r")

    first_line = TRUE

    while(length(x <- readLines(con, n = 1)) > 0) {



        if (grepl("insert_sequence", x, fixed = TRUE)) { #alignment insertion

            write(alignment_string, file=state_clocks_xml, append=!first_line)

        } else if (grepl("insert_trait", x, fixed = TRUE)) { #taxon insertion

            write(gsub("insert_trait", taxon_traits_string , x), file=state_clocks_xml, append = !first_line)

        } else if (grepl("insert_clock_rate", x, fixed = TRUE)) { #taxon insertion

            write(gsub("insert_clock_rate", clock_rate_string , x), file=state_clocks_xml, append = !first_line)

        } else {

            write(x, file=state_clocks_xml, append=!first_line)
        }

        first_line = FALSE


    }
    close(con)

}
#Files for templates and xml files
state_clocks_template_xml <- "../li-application/templates/li_state_clocks_template.xml"

li_T1_state_clocks_xml <- "../li-application/templates/T1_wgs_state_clocks.xml"
li_T2_state_clocks_xml <- "../li-application/templates/T2_wgs_state_clocks.xml"
write_alignment_to_state_xml_from_template(state_clocks_template_xml = state_clocks_template_xml,
                                           alignment_string = T1_alignment_string,
                                           taxon_traits_string = T1_taxon_traits_string,
                                           state_clocks_xml = li_T1_state_clocks_xml,
                                           clock_rate_string = T1_clock_rate_string)

write_alignment_to_state_xml_from_template(state_clocks_template_xml = state_clocks_template_xml,
                                           alignment_string = T2_alignment_string,
                                           taxon_traits_string = T2_taxon_traits_string,
                                           state_clocks_xml = li_T2_state_clocks_xml,
                                           clock_rate_string = T2_clock_rate_string)



