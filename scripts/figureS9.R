library(tumortree)
library(tidyverse)

#VAFs calculated / formatted in format_snv_vafs_for_histograms.R

t1_vaf_freqs <- read_tsv(file = "../li-application/li_hcc_data/t1_vaf_freqs.tsv")

t2_vaf_freqs <- read_tsv(file = "../li-application/li_hcc_data/t2_vaf_freqs.tsv")

t1_vaf_df_marked <- t1_vaf_freqs %>%
    dplyr::filter(mark == "truncal")

t2_vaf_df_marked <- t2_vaf_freqs %>%
    dplyr::filter(mark == "truncal")

t1_vaf_df_masked <- t1_vaf_freqs %>%
    dplyr::filter(masked)

t2_vaf_df_masked <- t2_vaf_freqs %>%
    dplyr::filter(masked)

edge_center_colors <- get_color_palette(names = c("edge", "center"))
t1_vaf_histograms <- ggplot(t1_vaf_freqs, (aes(x=t_alt_vaf, fill = ifelse(edgeP == 1, "edge", "center")))) +
    geom_histogram(color = "black", bins = 30) + theme_bw() + facet_wrap(~toupper(sample), ncol = 6) +
    geom_histogram(color = "black", bins = 30, fill = "grey", data = t1_vaf_df_marked, alpha = 0.5) + 
    #geom_vline(xintercept = 0.25, linetype = "dashed") +
    scale_fill_manual(values = edge_center_colors) +
    theme(legend.position = "none") +
    theme(text=element_text(size = 10)) +
    xlab("Variant Allele Frequency")
    #geom_point(data = t1_vaf_df_marked, color = "black", y=0)
t1_vaf_histograms

t1_vaf_histograms_masked <- ggplot(t1_vaf_df_masked, (aes(x=t_alt_vaf, fill = ifelse(edgeP == 1, "edge", "center")))) +
    geom_histogram(color = "black", bins = 30) + theme_bw() + facet_wrap(~toupper(sample), ncol = 6) +
    #geom_histogram(color = "black", bins = 30, fill = "grey", data = t1_vaf_df_masked, alpha = 0.5) + 
    #geom_vline(xintercept = 0.25, linetype = "dashed") +
    scale_fill_manual(values = edge_center_colors) +
    theme(legend.position = "none") +
    theme(text=element_text(size = 10)) +
    xlab("Variant Allele Frequency")
#geom_point(data = t1_vaf_df_marked, color = "black", y=0)
t1_vaf_histograms_masked

ggsave(t1_vaf_histograms, file = "../figures/t1_vaf_histograms.png", height = 5, width = 12)

t2_vaf_histograms <- ggplot(t2_vaf_df, (aes(x=t_alt_vaf, fill = ifelse(edgeP == 1, "edge", "center")))) +
    geom_histogram(color = "black", bins = 30) + theme_bw() +
    geom_histogram(color = "black", fill = "grey", data = t2_vaf_df_marked, alpha = 0.5, bins = 30) +
    facet_wrap(~toupper(sample), ncol = 6) +
    #geom_vline(xintercept = 0.25, linetype = "dashed") +
    scale_fill_manual(values = edge_center_colors) +
    theme(legend.position = "none")  +
    theme(text=element_text(size = 10)) +
    xlab("Variant Allele Frequency")
    #geom_point(data = t2_vaf_df_marked, color = "black", y=0)
t2_vaf_histograms 

t2_vaf_histograms_masked <- ggplot(t2_vaf_df_masked, (aes(x=t_alt_vaf, fill = ifelse(edgeP == 1, "edge", "center")))) +
    geom_histogram(color = "black", bins = 30) + theme_bw() + facet_wrap(~toupper(sample), ncol = 6) +
   # geom_histogram(color = "black", bins = 30, fill = "grey", data = t2_vaf_df_masked, alpha = 0.5) + 
    #geom_vline(xintercept = 0.25, linetype = "dashed") +
    scale_fill_manual(values = edge_center_colors) +
    theme(legend.position = "none") +
    theme(text=element_text(size = 10)) +
    xlab("Variant Allele Frequency")
t2_vaf_histograms_masked 
ggsave(t2_vaf_histograms, file = "../figures/t2_vaf_histograms.png", height = 5*0.666, width = 12)

#TUMOR 3 (Ling et al)

# ling_vaf_data <- read_tsv(file = "../ling-application/data/combined_variant_data.tsv")
# 
# ling_vaf_data_fixed <- ling_vaf_data %>% 
#     filter(fixed, Punch != "Z1", mut_allele_freq > 0.05)  %>% 
#     mutate(Punch = paste0("T3", Punch, sep =""))
# edge_center_colors = get_color_palette(c("edge", "center"))
# 
# t3_vaf_histograms <- ling_vaf_data %>%
#     filter(Punch != "Z1", mut_allele_freq > 0.05) %>% 
#     mutate(Punch = paste0("T3", Punch, sep ="")) %>% 
#     ggplot(., aes(x=mut_allele_freq)) +
#     geom_histogram(aes(fill = ifelse(edge == 1, "edge", "center")), color = "black",bins = 30) +
#     geom_histogram(color = "black", fill = "grey", data = ling_vaf_data_fixed, alpha = 0.5, bins = 30) +
#     facet_wrap(~Punch,  ncol = 6) +
#     theme_bw() +
#     scale_fill_manual(values = edge_center_colors) +
#     theme(legend.position = "none") +
#     #geom_vline(xintercept = 0.25, linetype = "dashed") +
#     theme(text=element_text(size = 25)) +
#     xlab("Variant Allele Frequency")
# t3_vaf_histograms 
# ggsave(t3_vaf_histograms, file = "../figures/t3_vaf_histograms.png", height = 5*1.33, width = 12)
# # t1_vaf_df %>% 
# #     group_by(sample) %>% 
# #     summarize(n_variants = sum(variant, na.rm = TRUE),
# #               max_vaf = max(t_alt_vaf,
# #                             na.rm = TRUE),
# #               mean_vaf = mean(t_alt_vaf, na.rm = TRUE), 
# #               edgeP = mean(edgeP, na.rm = TRUE)) %>% 
# #     ggplot(aes(x=max_vaf, y = n_variants, color = ifelse(edgeP == 1, "edge", "center"))) + geom_point() +
# #     scale_color_manual(values =edge_center_colors) + theme_classic()
# # 
# # t1_vaf_df %>% 
# #     group_by(sample) %>% 
# #     summarize(n_variants = sum(variant, na.rm = TRUE),
# #               max_vaf = max(t_alt_vaf,
# #                             na.rm = TRUE),
# #               mean_vaf = mean(t_alt_vaf, na.rm = TRUE), 
# #               edgeP = mean(edgeP, na.rm = TRUE)) %>% 
# #     ggplot(aes(x=mean_vaf, y = n_variants, color = ifelse(edgeP == 1, "edge", "center"))) + geom_point() +
# #     scale_color_manual(values =edge_center_colors) + theme_classic()
# # # 
# #
# # T1_all_amp_punches <- unique(T1_data_amp_long$Punch)
# # T2_all_amp_punches <- unique(T2_data_amp_long$Punch)
# 
# 
# 
# 
