library(tidyverse)
library(ggpubr)
library(tabulator)
library(gtools)

# load in data
unblinding_groups <- readxl::read_xlsx("unblinding_groups.xlsx")
olink_all <- readxl::read_xlsx("UNPAIRED_ALL.xlsx")

plasma_annotated <- full_join(unblinding_groups, olink_all) %>% na.omit()

### paired samples
byDrug_comparison <- compare_means(olink_value ~ treatment, group.by = c("protein"), data = plasma_annotated %>% 
                                     filter(
                                       # subj %in% paired_samples$subj, 
                                            group == "C")) 
byDrug_proteins <- byDrug_comparison %>% dplyr::filter(p<0.01) %>% pull(protein) %>% unique()

# boxplot by drug
plasma_annotated %>% dplyr::filter(protein %in% byDrug_proteins, group == "C"
                                     # subj %in% paired_samples$subj
                                     ) %>% 
  ggplot(aes(x = factor(treatment, levels = c("FTC-TAF", "FTC-TDF", "CT")), 
             y = olink_value, 
             fill = factor(treatment, levels = c("FTC-TAF", "FTC-TDF", "CT")))) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.2) +
  expand_limits(y = c(0, 100)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(x = '', y = 'NPx'
       # title = "Significantly different proteins in plasma according to drug (combined, circumcision)"
       ) +
  theme_ipsum() +
  facet_wrap(~ protein, scales = "fixed", nrow = 2) +
  ggthemes::scale_fill_tableau(name=NULL) +
  theme(
    legend.position = 'none',
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  geom_signif(comparisons = list(c("FTC-TAF", "CT"),
                                 c("FTC-TDF", "CT"),
                                 c("FTC-TAF", "FTC-TDF")),
              map_signif_level=TRUE, step_increase = 0.15, textsize= 4, size = 0.5, y_position = 1.4)

annotated_olink %>% dplyr::filter(marker %in% byDrug_proteins, Group == "Circumcision", subj %in% paired_samples$subj) %>% tab(Drug, marker)

# save figures
ggsave("combined_byDrug__.pdf", plot = last_plot(), device = cairo_pdf, width = 4, height = 3, scale = 2)
ggsave("combined_byDrug_alldata.png", plot = last_plot(), device = "png", width = 6, height = 6, scale = 2)
