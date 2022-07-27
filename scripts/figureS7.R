#figureS7.R

##### PLOT POSTERIORS ####

#Record of local directory
# log_files <- list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/curr_runs/out",
#                         pattern=".log", 
#                         full.names = TRUE)

log_files <- list.files(path = "../li-application/out",
                        pattern=".log",
                        full.names = TRUE)


#temporarily get rid of last corrupted log file
log_files <- log_files[! grepl("T1_wgs_newstates_bidir_state_rep2.log", log_files)]

#logs <- purrr::map(log_files, readLog)

process_logs <- function(log_file) {
    print(log_file)
    rep_extract <- regmatches(basename(log_file),
                              gregexpr("(?<=rep)[0-9]", basename(log_file), perl = TRUE))[[1]]
    migration_model <- ifelse(grepl("unidir", basename(log_file)), "unidirectional", "bidirectional")
    tumor_extract <- regmatches(basename(log_file), gregexpr("T[12]", basename(log_file), perl = TRUE))[[1]]
    states_extract <- ifelse(grepl("newstates", basename(log_file)), "newstates", "oldstates")
    clock_extract <- ifelse(grepl("strict", basename(log_file)), "strict", "state-dependent")
    log <- readLog(log_file)
    log_df <- as.data.frame(log) %>%
        add_column("migration_model" = migration_model,
                   "rep" = rep_extract,
                   "tumor" = tumor_extract,
                   "states" = states_extract,
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
#for supplement compare clock models (both unidirectional + oldstates)

## Tumor 1
t1_wgs_posteriors_plot_oldstates_clock_comparison <- all_logs_birthRate_df %>% 
    filter(migration_model == "unidirectional",
           states == "oldstates",
           tumor == "T1") %>% 
    ggplot(., aes(x=birthRate), color = "black") +
    geom_density(aes(fill = state), alpha=0.8) +
    facet_grid(rows = vars(clock_model)) +
    theme_classic() + scale_fill_manual(values=colors_loc) +
    theme(text=element_text(size=20))+
    xlab("Estimated birth rate") +
    theme(legend.position = "none") +ylab("") 
t1_wgs_posteriors_plot_oldstates_clock_comparison
ggsave(plot=t1_wgs_posteriors_plot_oldstates_clock_comparison,
       file ="../figures/t1_li_wgs_posteriors_oldstates_clock_comparison.png", height = 5, width = 5)

## Tumor 20
t2_wgs_posteriors_plot_oldstates_clock_comparison <- all_logs_birthRate_df %>% 
    filter(migration_model == "unidirectional",
           states == "oldstates",
           tumor == "T2") %>% 
    ggplot(., aes(x=birthRate), color = "black") +
    geom_density(aes(fill = state), alpha=0.8) +
    facet_grid(rows = vars(clock_model)) +
    theme_classic() + scale_fill_manual(values=colors_loc) +
    theme(text=element_text(size=20))+
    xlab("Estimated birth rate") +
    theme(legend.position = "none") +ylab("") 
t2_wgs_posteriors_plot_oldstates_clock_comparison
ggsave(plot=t2_wgs_posteriors_plot_oldstates_clock_comparison,
       file ="../figures/t2_li_wgs_posteriors_oldstates_clock_comparison.png", height = 5, width = 5)

