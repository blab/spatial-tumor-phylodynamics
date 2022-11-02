#write_all_ling_state_clocks_xml.R

library(tumortree)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stats)
library(stringr)
library(raster)


ling_data <- read.csv("../ling-application/data/dataset_S2_cleaned_allele_freqs_with_validation.csv")
ling_validation_snvs <- read.csv("../ling-application/data/dataset_S3_wes_only.csv", header = TRUE)
ling_data_all_alleles <- read.csv("../ling-application/data/alleles.csv")

state_clocks_template_xml <- "../ling-application/xmls/state-clocks-template2.xml"
ling_state_clocks_xml <- "../ling-application/xmls/hcc-wes_unidir_state_rep1.xml"
ling_fasta <- "../ling-application/data/ling-punch-sequences.fa"

ling_data$SNV_fixed_status <- gsub(" ", "", ling_data$SNV_fixed_status)

#sum(ling_data$SNV_fixed_status == "Polymorphic")
##set requirements for calling

min_var_allele_freq <- 0.05 #threshold to call presense of variant allele

ref_cutoff <- 10 #threshold for calling variant as absent

#pivot data, this takes some wrangling because there are multiple columns per punch
ling_data_long <- ling_data %>%
    dplyr::select(-contains("N1")) %>% 
    tidyr::pivot_longer(cols = contains("mut_freq"),
                        names_to = c("Punch"),
                        names_pattern = "^(.+?)_mut_freq",
                        values_to = "mut_allele_freq")  %>%
    tidyr::pivot_longer(cols = ends_with("_ref"),
                        names_to = c("Punch1"),
                        names_pattern = "^(.+?)_ref",
                        values_to = "ref") %>%
    dplyr::filter(Punch == Punch1) %>% # filter to only matching punches between pivots to avoid data overload
    dplyr::select(-Punch1) %>%
    tidyr::pivot_longer(cols = ends_with("_mut"),
                        names_to = c("Punch2"),
                        names_pattern = "^(.+?)_mut", values_to = "mut") %>%
    dplyr::filter(Punch == Punch2) %>%
    dplyr::select(-Punch2) %>%
    tidyr::pivot_longer(cols = ends_with("_mut_allele"),
                        names_to = c("Punch3"),
                        names_pattern = "^(.+?)_mut_allele",
                        values_to = "mut_allele") %>%
    dplyr::filter(Punch == Punch3) %>%
    dplyr::select(-Punch3) %>%
    dplyr::filter(! is.na(start)) %>%
    dplyr::arrange(chr, start) %>%
    dplyr::filter(SNV_fixed_status != "SNV_in_CNA_regions" & Punch != "N1") %>% 
    dplyr::mutate(mut_allele_freq = ifelse(mut_allele_freq < min_var_allele_freq & ref < ref_cutoff,
                                           NA,
                                           mut_allele_freq)) %>% 
    dplyr::mutate(mut_id = paste(chr, start, mut_allele, sep = "_"))



#this data defines all sites from WES data, a subset of these were genotyped

# wes_sites <- ling_data_long %>%
#     dplyr::select(chr, start, ref_allele, mut_allele) %>%
#     distinct %>%
#     arrange(chr,start)

#In data there are sometimes two different alleles detected
## Only take the majority allele

wes_sites <- ling_data_long %>% 
  group_by(chr, start, ref_allele, mut_allele, mut_id) %>% 
  summarise(count = n()) %>% 
  arrange(-count) %>% 
  ungroup() %>% 
  distinct(chr, start, .keep_all = TRUE)

#Only include filtered WES sites
ling_data_long <- ling_data_long %>% 
  dplyr::filter(mut_id %in% wes_sites$mut_id)


ling_data_genotyped_long <- ling_data_all_alleles %>%
    tidyr::pivot_longer(cols = -Sites, names_to = c("location"), values_to = "mut_allele_freq") %>%
    separate(col=location, into=c("chr", "start") ,sep = "_") %>%
    rename(Punch=Sites) %>%
    mutate(mut_allele_freq = na_if(mut_allele_freq, "-")) %>%
    mutate(mut_allele_freq = as.numeric(na_if(mut_allele_freq, "F"))) %>%
    mutate(start = as.integer(start)) %>%
    left_join(., wes_sites, by = c("chr", "start"))
    
  
    #dplyr::mutate("base_call" = ifelse(is.na(mut_allele_freq), "N", ifelse((mut_allele_freq >= min_var_allele_freq), mut_allele, ref_allele)))

wes_punches <- unique(ling_data_long$Punch)

comb_data <- ling_data_long %>%
    dplyr::select(Punch, chr, start, mut_allele_freq, ref_allele, mut_allele) %>%
    bind_rows(., ling_data_genotyped_long)%>%
    mutate(Punch = str_replace(Punch, fixed("*") ,"")) %>%
    group_by(Punch, chr, start, ref_allele, mut_allele) %>%
    dplyr::summarise(mut_allele_freq = mean(mut_allele_freq, na.rm= TRUE)) %>%
    filter(Punch != "Z2" & Punch %in% wes_punches)


ling_mat <- ling_data %>% 
  dplyr::filter(SNV_fixed_status != "SNV_in_CNA_regions") %>% 
  dplyr::filter(! is.na(start)) %>% 
  dplyr::select(contains("mut_freq")) %>% 
  as.matrix()

mut_ids <- ling_data %>% 
  dplyr::filter(SNV_fixed_status != "SNV_in_CNA_regions") %>% 
  dplyr::filter(! is.na(start)) %>% 
  dplyr::mutate(mut_id = paste(chr, start, mut_allele, sep = "_"))

rownames(ling_mat) <- mut_ids$mut_id

ling_mat[is.na(ling_mat)] <- 0
ling_clustering_mut <- hclust(dist(ling_mat))
ling_clustering_punch <- hclust(dist(t(ling_mat)))

punch_order <- gsub("_mut_freq", "", ling_clustering_punch$labels[ling_clustering_punch$order])
mut_order <- ling_clustering_mut$labels[ling_clustering_mut$order]

ling_data_long_ordered <- ling_data_long %>% 
  dplyr::mutate(mut_id = paste(chr, start, sep = "_")) %>% 
  mutate(mut_id = factor(mut_id, levels = mut_order), 
         Punch = factor(Punch, levels = punch_order))

# ling_data_long %>% 
#   dplyr::filter(SNV_fixed_status == "Fixed") %>% 
#   ggplot(., aes(x=mut/(ref+mut))) + geom_histogram(aes(fill = Punch)) +
#   facet_wrap(~Punch) + theme_classic() + xlim(c(0, 50))

#Visualize heatmap of variants
# ling_data_long_ordered %>% 
#   mutate(mut_allele_freq = ifelse(mut_allele_freq < 0.05 & ref < 10, NA, mut_allele_freq)) %>% 
#   ggplot(., aes(x=Punch, y = mut_id)) +
#   geom_tile(aes(fill = ifelse(mut_allele_freq > 0.05, "present", "absent")))


  #scale_fill_brewer()

# comb_data %>% 
#   group_by(Punch) %>% 
#   summarise(cov = mean(ref+mut, na.rm = TRUE)) %>% 
#   ggplot(., aes(x=cov, y = Punch)) + geom_point(aes(color = Punch), alpha = 0.5) + theme_classic()


get_punch_sequence <- function(punch, comb_data, wes_sites, min_var_allele_freq = 0.05,
                               collapse = TRUE, weighted = FALSE, wes_punches = NULL) {


    punch_data <- comb_data %>%
        ungroup %>%
        dplyr::filter(Punch == punch) %>%
        dplyr::select(-Punch) %>%
        left_join(wes_sites,., by = c("chr", "start", "ref_allele", "mut_allele")) %>%
        #dplyr::mutate(mut_allele_freq = ifelse(mut_allele_freq <  min_var_allele_freq & ref < ref_count_cutoff, NA, mut_allele_freq)) %>%
        dplyr::mutate("base_call" = ifelse(is.na(mut_allele_freq),
                                           "N",
                                           ifelse((mut_allele_freq >= min_var_allele_freq), mut_allele, ref_allele)))


    if (collapse) {

        punch_sequence <- paste(punch_data$base_call, collapse = "")

    } else {

        punch_sequence <- matrix(punch_data$base_call, nrow = 1, byrow = TRUE)
        dimnames(punch_sequence) <- list(c(punch), as.character(c(1:length(punch_data$base_call))))
    }

    if(weighted) {

      if(is.null(wes_punches)) {

        error("Must provide WES punch labels for weighted option")
      }

      if(punch %in% wes_punches) {

        punch_sequence <- paste(punch_sequence, "ATGC", sep = "")

      } else {

        punch_sequence <- paste(punch_sequence, "NNNN", sep = "")
      }
    }

    return(punch_sequence)
}

punches <- str_sort(unique(comb_data$Punch), numeric = TRUE)
punches <- punches[punches != "Z1"]
punch_sequences <- unlist(purrr::map(punches,
                                     function(p) get_punch_sequence(punch = p, comb_data = comb_data,
                                                                    wes_sites = wes_sites,
                                                                    weighted = FALSE,
                                                                    wes_punches = wes_punches)))

names(punch_sequences) <- punches

#santity check
all(map_dbl(punch_sequences, nchar) == nrow(wes_sites))

#perimeter:12.8cm

## from webDigitizer2 (tar file contains original project)
ling_boundary_points <- read.csv("../ling-application/spatial_measurements/boundary_points.csv")
ling_boundary_points <- bind_rows(ling_boundary_points, ling_boundary_points[1,])

ling_punches_points <- read.csv("../ling-application/spatial_measurements/punches_points.csv")

find_closest_dist_to_boundary <- function(x,y, boundary_points) {
    return(min(raster::pointDistance(p1 = c(x,y), p2 = boundary_points, lonlat = FALSE)))
}

ling_punches_points$dist_to_boundary <- map2_dbl(ling_punches_points$X, ling_punches_points$Y, function(x,y) find_closest_dist_to_boundary(x,y, boundary_points = ling_boundary_points))

edge_cutoff <- 0.35 #cm


ling_punches_points <- ling_punches_points %>%
    mutate("edge" = as.integer(dist_to_boundary < edge_cutoff)) %>%
    filter(Punch %in% punches)

write_csv(ling_punches_points, file="../ling-application/data/edge-center-calls.csv")


################ estimating sampling proportion ######################
total_area <- (17.5^2)*pi / 100 #cm

edge_area <- 128 * 3.5 / 100
edge_fraction_area <- edge_area/ total_area 

punch_diameter <-  0.5 #mm
punch_area <- 0.25^2 *pi
n_edge_punches <- sum(ling_punches_points$edge)
edge_area_sampled <- n_edge_punches*punch_area
fraction_edge_sampled <- edge_area_sampled / (edge_area *100)

n_center_punches <- sum(1 - ling_punches_points$edge)
center_area_sampled <- n_center_punches*punch_area
fraction_center_sampled <-center_area_sampled / ((total_area - edge_area) *100)
############################ PLOTTING #################
scale_bar_length <- 0.5 #cm
scale_bar_position_x <- 3.6
scale_bar_position_y <- 0.5
scale_bar <- data.frame("X" = c(scale_bar_position_x, scale_bar_position_x + scale_bar_length),
                        "Y" = c(scale_bar_position_y, scale_bar_position_y))

ggplot(data = ling_boundary_points, aes(x = X, y = -Y)) +
    geom_path(size= 1.5, color =  "lightgrey") +
    geom_point(data = ling_punches_points, size = 4, aes(color = as.factor(edge))) +
    geom_text(data = ling_punches_points, aes(label = Punch), size = 2) +
    theme_void() +
    scale_color_manual(values = c( "brown3", "cadetblue")) +
    labs(color = "") + ggtitle("Ling et al. HCC Tumor") +
    coord_fixed(ratio = 1) +
    geom_line(data = scale_bar) +
    theme(legend.position = "none") +
    geom_text(data = data.frame("X" = c(scale_bar_position_x + scale_bar_length/2), "Y" = c(scale_bar_position_y - 0.1)), label= c(paste0(scale_bar_length, "cm", sep = " ")), size = 4)


###write xml file

# #Add N's 
# 
# punch_sequences <- paste(punch_sequences, "NNNNNNNNNNNNNNNN", sep = "")
alignment_string <- ""

for (i in 1:length(punches)) {
    alignment_string <- paste0(alignment_string, paste0('\t', '<sequence id="seq_',punches[i], '" spec="Sequence" taxon="',punches[i], '" totalcount="4" value="',
                                                        punch_sequences[i], '"/>', sep = ""), sep = "\n")
}

taxon_traits <- c()

for (i in 1:nrow(ling_punches_points)) {

    taxon_traits <- c(taxon_traits, paste0(ling_punches_points$Punch[i], "=loc", ling_punches_points$edge[i], sep = ""))
}

taxon_traits_string <- paste(taxon_traits, collapse = ",")

con = file(state_clocks_template_xml, "r")

first_line = TRUE

while(length(x <- readLines(con, n = 1)) > 0) {



    if (grepl("insert_sequence", x, fixed = TRUE)) { #alignment insertion

        write(alignment_string, file=ling_state_clocks_xml, append=!first_line)

    } else if (grepl("insert_trait", x, fixed = TRUE)) { #taxon insertion

        write(gsub("insert_trait", taxon_traits_string , x), file=ling_state_clocks_xml, append = !first_line)

    } else if (grepl("insert_weights", x, fixed = TRUE)) { #taxon insertion

      write(gsub("insert_weights", nucleotide_weights_string , x), file=ling_state_clocks_xml, append = !first_line)

    } else {

        write(x, file=ling_state_clocks_xml, append=!first_line)
    }

    first_line = FALSE


}
close(con)



###fasta file
punches <- c(punches, "Normal")


ref_seq <- paste(wes_sites$ref_allele, collapse = "")
punch_sequences <- c(punch_sequences, ref_seq)
con = file(ling_fasta, "w")

first_line = TRUE

for (p in 1:length(punches)) {
  
  if (! (punches[p] == "Z1")) {
    write(paste0(">", punches[p]), file=ling_fasta, append = !first_line)
    first_line = FALSE
    write(punch_sequences[p], file=ling_fasta, append = !first_line)
  }
}


close(con)

#Make nexstrain metadata
metadata <- ling_punches_points %>% 
  rename("strain" = Punch) %>% 
  add_column("date" = "2015-01-01") %>% 
  add_row(strain = "Normal")

write_tsv(metadata, file = "../ling-application/data/punch_metadata.tsv")

fixed_snvs <- ling_data %>% 
  dplyr::filter(SNV_fixed_status == "Fixed") %>% 
  dplyr::mutate(mut_id = paste(chr, start, sep = "_")) %>% 
  distinct(mut_id)

#for figure annotation
length(fixed_snvs$mut_id)

comb_snv_data <- comb_data %>% 
  left_join(.,ling_punches_points, by = "Punch") %>% 
  dplyr::mutate(mut_id = paste(chr, start, sep = "_")) %>% 
  mutate("fixed" = (mut_id %in% fixed_snvs$mut_id))
  

#Combine information into one table to save for histogram plots
write_tsv(comb_snv_data,
          file = "../ling-application/data/combined_variant_data.tsv")
