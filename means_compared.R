library(tidyverse)
library(ggpubr)
library(gtools)
library(hrbrthemes)
library(forcats)
library(scales)

# load in data for Uganda by therapy
ebo_therapy <- readxl::read_xlsx('EBO_clean_therapy.xlsx')

# long format
ebo_therapy_long <- ebo_therapy %>% gather('marker', 'value', -c(Group, Drug,  Assay)) %>% na.omit() %>% 
  dplyr::group_by(marker) %>% dplyr::filter(n_distinct(value) >= 113) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(
    subj = str_extract(Assay, "(^.*)(?=R|C)"),
    country = if_else(grepl("CJ", subj), "SA", "U")
  )

# load in data for SA
cj <- readxl::read_xlsx('CJ_cleaned.xlsx')

# long format
cj_long <- cj %>% gather('marker', 'value', -c(Group, Assay,)) %>% na.omit() %>% 
  dplyr::group_by(marker) %>% dplyr::filter(n_distinct(value) >= 101) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(value = as.numeric(value))


# bind data together
combined_therapy_long <- rbind(ebo_therapy_long, cj_therapy_long)

# calculate means
r_means <- combined_therapy_long %>% dplyr::select(-Assay) %>% 
  dplyr::filter(Group %in% 'Randomization') %>% 
  group_by(Group, marker, Drug) %>% 
  dplyr::summarise(group1_mean = mean(value)) %>% 
  dplyr::rename(group1 = Group)

c_means <- combined_therapy_long %>% dplyr::select(-Assay) %>% 
  dplyr::filter(Group %in% 'Circumcision') %>% 
  group_by(Group, marker, Drug) %>% 
  dplyr::summarise(group2_mean = mean(value)) %>% 
  dplyr::rename(group2 = Group)

paired_samples <- readxl::read_xlsx("olink_pairs.xlsx")

# missing pair samples
missing_pair <- combined_therapy_long %>% dplyr::filter(subj %in% paired_samples$subj) %>% tab(subj, marker) %>% dplyr::filter(N < 2)

combined_therapy_stats <- compare_means(value ~ Group, data = combined_therapy_long %>% dplyr::filter(subj %in% paired_samples$subj) %>% 
                                          dplyr::filter(!marker %in% missing_pair$marker), 
                                        group.by = c('marker', 'Drug'), method = 'wilcox.test', p.adjust.method = "fdr", paired = T) %>%  # FALSE in the paper!!!
  left_join(.,r_means, by = c('marker', 'group1', 'Drug')) %>% 
  left_join(.,c_means, by = c('marker', 'group2', 'Drug')) %>% 
  dplyr::mutate(FC = foldchange(group2_mean, group1_mean))

combined_sig_markers <- combined_therapy_stats %>% dplyr::filter(p.adj < 0.01) %>% pull(marker) %>% unique()

combined_therapy_stats %>% 
  dplyr::filter(marker %in% combined_sig_markers) %>%
  dplyr::mutate(value = ifelse(p.adj<0.01, 1, 0)) %>% 
  dplyr::mutate(value = factor(value),
  Drug = gsub("FTAF", "FTC-TAF", Drug),
  Drug = gsub("TRUVADA", "FTC-TDF", Drug)) %>%
  ggplot(aes(marker, factor(Drug, levels = rev(c('CONTROL', 'FTC-TAF', 'FTC-TDF'))), color = value)) +
  geom_point(aes(size = p.adj)) +
  scale_size("p-values", trans="log10", range=c(10, 0)) +
  scale_color_manual(values = c('#fafafa', 'orange')) +
  scale_fill_manual(values = c('blue', 'red')) +
  # coord_flip() +
  theme_bw() +
  labs(x ='', y = '')+
  theme(
    # legend.position = 'none',
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    panel.grid.major = element_line(linetype = "dashed", colour = "grey", size = 0.1)
  )

# save figure
ggsave('plot_combined_paired.pdf', plot = last_plot(), device = cairo_pdf, height = 3, width = 15)
