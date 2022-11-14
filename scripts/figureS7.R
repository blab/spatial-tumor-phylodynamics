#figureS7.R

library(tumortree)
library(tidyverse)
##### PLOT POSTERIORS ####

#Record of local directory
# log_files <- list.files(path = "/Volumes/BALAENA/projects/spatial_tumor_growth_simulation/outputs/beast_analysis/state_dependent_clock_model/primary_tumor_analysis/li/curr_runs/out",
#                         pattern=".log", 
#                         full.names = TRUE)

li_log_files <- list.files(path = "../li-application/combined",
                           pattern=".log",
                           full.names = TRUE)

#get reid of inidivudal chains
li_log_files <- li_log_files[! grepl("chain1", li_log_files)]
li_log_files <- li_log_files[! grepl("T1red", li_log_files)]
li_log_files <- li_log_files[! grepl("bidir", li_log_files)]


# ling_log_files <- list.files(path = "../ling-application/out",
#                              pattern=".log",
#                              full.names = TRUE)
# 
# ling_log_files <- ling_log_files[! grepl("chain1", ling_log_files)]
# ling_log_files <- ling_log_files[! grepl("bidir", ling_log_files)]
# 
# log_files <- c(li_log_files, ling_log_files)

colors_loc <- tumortree::get_color_palette(names = c("edge", "center"))
names(colors_loc) <- c("loc1", "loc0")

process_logs <- function(log_file) {
    print(log_file)
    migration_model <- ifelse(grepl("unidir", basename(log_file)), "unidirectional", "bidirectional")
    tumor_extract <- regmatches(basename(log_file), gregexpr("T[1-3]", basename(log_file), perl = TRUE))[[1]]
    if(length(tumor_extract) == 0) {
        tumor_extract = "T3"
    }
    subset_extract <- regmatches(basename(log_file),
                                 gregexpr("(?<=unidir_)[0-9]", basename(log_file), perl = TRUE))[[1]]
    if (length(subset_extract) == 0) {
        subset_extract <- "1"
    }
    states_extract <- ifelse(grepl("newstates", basename(log_file)), "newstates",
                             ifelse(grepl("oristates", basename(log_file)), "oldstates", NA))
    clock_extract <- ifelse(grepl("strict", basename(log_file)), "strict", "state-dependent")
    log <- readLog(log_file)

    log_df <- as.data.frame(log) %>%
        add_column("migration_model" = migration_model,
                   "tumor" = tumor_extract,
                   "states" = states_extract,
                   "subset" = subset_extract,
                   "clock_model" = clock_extract)
    
    if ("birthRateSVCanonical.loc2" %in% colnames(log_df)) {
        
        log_df <- log_df %>% 
            rename("birthRateSVCanonical.loc0" = birthRateSVCanonical.loc1) %>% 
            rename("birthRateSVCanonical.loc1" = birthRateSVCanonical.loc2)
    }
    log_df <- log_df %>% 
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
    mutate("category" = paste(subset, state)) %>% 
    ggplot(., aes(x=birthRate, group = category), color = "black") +
    geom_density(aes(fill = state), alpha=0.7) +
    facet_wrap(~clock_model, nrow = 1) +
    theme_classic() + scale_fill_manual(values=colors_loc) +
    theme(text=element_text(size=20))+
    xlab("Estimated birth rate") +
    theme(legend.position = "none") +ylab("") 

t1_wgs_posteriors_plot_oldstates_clock_comparison
ggsave(plot=t1_wgs_posteriors_plot_oldstates_clock_comparison,
       file ="../figures/t1_li_wgs_posteriors_oldstates_clock_comparison.png", height = 4, width = 6)

## Tumor 2
t2_wgs_posteriors_plot_oldstates_clock_comparison <- all_logs_birthRate_df %>% 
    filter(migration_model == "unidirectional",
           states == "oldstates",
           tumor == "T2") %>% 
    mutate("category" = paste(subset, state)) %>% 
    ggplot(., aes(x=birthRate, group = category), color = "black") +
    geom_density(aes(fill = state), alpha=0.8) +
    facet_wrap(~clock_model, nrow = 1) +
    theme_classic() + scale_fill_manual(values=colors_loc) +
    theme(text=element_text(size=20))+
    xlab("Estimated birth rate") +
    theme(legend.position = "none") +ylab("") 

t2_wgs_posteriors_plot_oldstates_clock_comparison

ggsave(plot=t2_wgs_posteriors_plot_oldstates_clock_comparison,
       file ="../figures/t2_li_wgs_posteriors_oldstates_clock_comparison.png", height = 4, width = 6)


# Get ratios 

#Tumor 1
t1_wgs_clock_comparison_ratios <- all_logs_birthRate_df %>% 
    filter(migration_model == "unidirectional",
           states == "oldstates",
           tumor == "T1") %>% 
    mutate(clock_model = factor(clock_model, levels = c("state-dependent", "strict"))) %>% 
    #mutate("category" = paste0(clock_model, subset)) %>% 
    group_by(clock_model, subset) %>% 
    summarise(mean = mean(birthRateRatio),
              hpd_lower = hdi(birthRateRatio, credMass = 0.95)[1],
              hpd_upper = hdi(birthRateRatio, credMass = 0.95)[2]) %>% 
    ggplot(., aes(x=subset, y=mean), color = "black") +
    geom_point(size =2) +
    facet_wrap(~clock_model, ncol = 2) +
    geom_errorbar(aes(ymin = hpd_lower, ymax= hpd_upper), width = 0, size = 1) +
    theme_classic() +
    
    theme(text=element_text(size=15))+
    ylab("") +
    xlab("") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    #coord_flip() +
    theme(strip.text = element_blank()) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank())


t1_wgs_clock_comparison_ratios
ggsave(plot=t1_wgs_clock_comparison_ratios,
       file ="../figures/t1_li_wgs_posteriors_ratios_clock_comparison.png", height = 3, width = 3)


#Tumor 2
t2_wgs_clock_comparison_ratios <- all_logs_birthRate_df %>% 
    filter(migration_model == "unidirectional",
           states == "oldstates",
           tumor == "T2") %>% 
    mutate(clock_model = factor(clock_model, levels = c("state-dependent", "strict"))) %>% 
    group_by(clock_model, subset) %>% 
    summarise(mean = mean(birthRateRatio),
              hpd_lower = hdi(birthRateRatio, credMass = 0.95)[1],
              hpd_upper = hdi(birthRateRatio, credMass = 0.95)[2]) %>% 
    ggplot(., aes(x=subset, y=mean), color = "black") +
    geom_point(size =2) +
    geom_errorbar(aes(ymin = hpd_lower, ymax= hpd_upper), width = 0, size = 1) +
    theme_classic() +
    
    theme(text=element_text(size=15))+
    ylab("") +
    xlab("") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    #coord_flip() +
    theme(strip.text = element_blank()) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank()) +
    facet_wrap(~clock_model, ncol = 2)


t2_wgs_clock_comparison_ratios
ggsave(plot=t2_wgs_clock_comparison_ratios,
       file ="../figures/t2_li_wgs_posteriors_ratios_clock_comparison.png", height = 3, width = 3)

# # Add Ling et al tumor
# ling_log_files <- list.files(path = "../ling-application/logs",
#                         pattern="_comb.log",
#                         full.names = TRUE)
# ling_all_logs_df <- purrr::map(ling_log_files, process_logs) %>% 
#     bind_rows
# 
# ling_all_logs_birthRate_df <- ling_all_logs_df %>% 
#     pivot_longer(., cols = c("birthRateSVCanonical.loc0", "birthRateSVCanonical.loc1"),
#                  names_to = "state", 
#                  names_prefix = "birthRateSVCanonical.", 
#                  values_to = "birthRate")

# t3_wgs_posteriors_clock_comparison <- all_logs_birthRate_df %>% 
#     filter(migration_model == "unidirectional",
#            #states == "oldstates",
#            tumor == "T3") %>% 
#     ggplot(., aes(x=birthRate), color = "black") +
#     geom_density(aes(fill = state), alpha=0.8) +
#     facet_wrap(~clock_model, nrow = 2) +
#     theme_classic() + scale_fill_manual(values=colors_loc) +
#     theme(text=element_text(size=20))+
#     xlab("Estimated birth rate") +
#     theme(legend.position = "none") +ylab("")  
# 
# t3_wgs_posteriors_clock_comparison
# ggsave(plot=t3_wgs_posteriors_clock_comparison,
#        file ="../figures/t3_ling_wes_posteriors_clock_comparison.png", height = 5, width = 4)
# 
# ##### Plot estimated edge/center ratios
# 
# t3_wgs_clock_comparison_ratios <- all_logs_birthRate_df %>% 
#     filter(migration_model == "unidirectional",
#            #states == "oldstates",
#            tumor == "T3") %>% 
#     mutate(clock_model = factor(clock_model, levels = c("state-dependent", "strict"))) %>% 
#     group_by(clock_model, subset) %>% 
#     summarise(mean = mean(birthRateRatio),
#               hpd_lower = hdi(birthRateRatio, credMass = 0.95)[1],
#               hpd_upper = hdi(birthRateRatio, credMass = 0.95)[2]) %>% 
#     ggplot(., aes(x=subset, y=mean), color = "black") +
#     geom_point(size =2) +
#     geom_errorbar(aes(ymin = hpd_lower, ymax= hpd_upper), width = 0, size = 1) +
#     theme_classic() +
#     
#     theme(text=element_text(size=15))+
#     ylab("") +
#     xlab("") +
#     geom_hline(yintercept = 1, linetype = "dashed") +
#     coord_flip() +
#     facet_wrap(~clock_model, ncol = 1) +
#     theme(strip.text = element_blank()) +
#     theme(axis.text.y = element_blank(), 
#           axis.ticks.y = element_blank())
# 
# 
# t3_wgs_clock_comparison_ratios
# ggsave(plot=t3_wgs_clock_comparison_ratios,
#        file ="../figures/t3_ling_wes_posteriors_ratios_clock_comparison.png", height = 2, width = 3)
#     