#sdevo_li_estimates_summary.R
library(tidyverse)
library(ggtree)
library(treeio)
library(beastio)
library(tumortree)
library(smoothr)
library(coda)
library(HDInterval)
library(ggrepel)
#Record of local directory
# log_files <- list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/curr_runs/out",
#                         pattern=".log", 
#                         full.names = TRUE)

log_files <- list.files(path = "../li-application/out",
                        pattern="[0-9].log",
                        full.names = TRUE)

log_files <- log_files[! grepl("chain1", log_files)]
log_files <- log_files[! grepl("T1red", log_files)]
log_files <- log_files[! grepl("bidir", log_files)]

#logs <- purrr::map(log_files, readLog)

process_logs <- function(log_file) {
    print(log_file)
    migration_model <- ifelse(grepl("unidir", basename(log_file)), "unidirectional", "bidirectional")
    tumor_extract <- regmatches(basename(log_file), gregexpr("T[12]", basename(log_file), perl = TRUE))[[1]]
    states_extract <- ifelse(grepl("newstates", basename(log_file)), "newstates", "oldstates")
    clock_extract <- ifelse(grepl("strict", basename(log_file)), "strict", "state-dependent")
    subset_extract <- regmatches(basename(log_file),
                                 gregexpr("(?<=unidir_)[0-9]", basename(log_file), perl = TRUE))[[1]]
    if(length(subset_extract) == 0) {
        subset_extract <- 1
    }
    log <- readLog(log_file, burnin = 0.2)
    log_df <- as.data.frame(log) %>%
        add_column("migration_model" = migration_model,
                   "tumor" = tumor_extract,
                   "states" = states_extract,
                   "subset" = subset_extract,
                   "clock_model" = clock_extract) %>%
        dplyr::mutate(birthRateDiff = birthRateSVCanonical.loc1 - birthRateSVCanonical.loc0,
                      birthRateRatio = birthRateSVCanonical.loc1 / birthRateSVCanonical.loc0)
    
    return(log_df)
}

all_logs_df <- purrr::map(log_files, process_logs) %>% 
    bind_rows

all_logs_birthRate_df <- all_logs_df %>% 
    pivot_longer(., cols = c("birthRateSVCanonical.loc0", "birthRateSVCanonical.loc1"),
                 names_to = "state", 
                 names_prefix = "birthRateSVCanonical.", 
                 values_to = "birthRate")

write_tsv(all_logs_birthRate_df, "../li-application/stats/sdevo_estimates_summary.tsv")
