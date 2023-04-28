library(tidyverse)
library(ggpubr)
library(tabulator)
library(gtools)

# load in protein abundace data
ebo_therapy <- readxl::read_xlsx('~/Documents/olink_chaps/EBO_clean_therapy.xlsx')
cj_therapy <- readxl::read_xlsx('~/Documents/olink_chaps/CJ_clean_therapy.xlsx')

cj_therapy <- cj_therapy %>% 
  dplyr::filter(!Assay == "CJ-051R" | Drug == "TRUVADA") %>% # CJ-051R is duplicated in control and truvada. should oly be truvada
  dplyr::filter(!Assay == "CJ-055R" | Drug == "FTAF") %>% # CJ-055R is duplicated in control and ftaf. should oly be ftaf
  gather('marker', 'value', -c(Group, Drug,  Assay)) %>% na.omit() %>%
  dplyr::group_by(marker) %>% dplyr::filter(n_distinct(value) >= 53) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(
    subj = str_extract(Assay, "(^.*)(?=R|C)"),
    country = if_else(grepl("CJ", subj), "SA", "U")
  )

cj_complete_pairs <- cj_therapy %>% tab(subj, marker) %>% dplyr::filter(N == 2) %>% pull(subj) %>% unique()

ebo_therapy <- ebo_therapy %>% 
  gather('marker', 'value', -c(Group, Drug,  Assay)) %>% na.omit() %>%
  dplyr::group_by(marker) %>% dplyr::filter(n_distinct(value) >= 113) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(
    subj = str_extract(Assay, "(^.*)(?=R|C)"),
    country = if_else(grepl("CJ", subj), "SA", "U")
  )

ebo_complete_pairs <- ebo_therapy %>% tab(subj, marker) %>% dplyr::filter(N == 2) %>% pull(subj) %>% unique()

complete_pairs <- c(cj_complete_pairs, ebo_complete_pairs)

combined_therapy <- rbind(ebo_therapy, cj_therapy) %>% dplyr::filter(subj %in% complete_pairs)

data_long <- combined_therapy

data <- data_long %>% 
  janitor::clean_names() %>%
  # dplyr::filter(drug  == 'FTAF') %>%
  dplyr::filter(drug  == 'TRUVADA') %>%
  dplyr::filter(drug  != 'CONTROL') %>%
  select(-drug, -subj, -country) %>% 
  tidyr::spread(marker, value) %>% 
  select(-assay) %>%
  select_if(~ !any(is.na(.))) # ramove markers with missing values but preserve pairs



data$group <- factor(data$group, levels = c('Randomization', 'Circumcision'))

#rfe does not run unless data is a data.frame
data <- as.data.frame(data)

# find top predictors using randomForest
library(mlbench)
library(caret)
set.seed(12)
control <- rfeControl(functions=rfFuncs, method="cv", number=10)
results <- rfe(data[,2:58], data[,1], sizes=c(2:58), rfeControl=control)
print(results)
# list the chosen features
predictors(results)

##### plot the results #######
# pdf("combined_accuracy_taf.pdf", width = 5, height = 4)
# plot(results, type=c("g", "o"))
# dev.off()

# rank by importance
# prepare training scheme
control <- trainControl(method="repeatedcv", number=10, repeats=5)
# train the model
model <- train(group~., data=data, method="rf", preProcess="scale", trControl=control)
# estimate variable importance
importance <- varImp(model, scale=FALSE)
# summarize importance
print(importance)
# plot importance
d<- importance$importance
d <- rownames_to_column(d, 'protein')
d$protein <- gsub('`', '', d$protein)
levels <- d %>% arrange(desc(Overall)) %>% pull(protein)
d$protein <- factor(d$protein, levels = rev(levels))
d <- d %>% arrange(desc(Overall))

## create separate data fpr each drug
d_taf <- d
d_tdf <- d

theme_set(theme_bw(base_size = 12))
p2 <- ggplot(head(d, 10), aes(y=Overall, x= `protein`)) + #use d_top30 for top 30 proteins; size is 4x6
  geom_segment(aes(y = 0, x = `protein`, yend = Overall, xend = `protein`)) + 
  geom_point(stat='identity', color="orange", size=2) +
  labs(x = '', y = 'Importance', title = 'Top 10 predictors (FTC-TDF)') +
  ylim(0,25) +
  theme(panel.grid = element_blank()) +
  coord_flip() +
  hrbrthemes::theme_ipsum()


ggsave("combined_predictors_tdf_pairs.pdf", plot = p2, device = cairo_pdf, width = 5, height = 5)

# merge data for both drugs

d_taf <- d_taf %>% dplyr::rename(imp_taf = Overall)
taf_top10 <- d_taf %>% arrange(desc(imp_taf)) %>% head(10) %>% pull(protein) %>% as.character()
d_tdf <- d_tdf %>% dplyr::rename(imp_tdf = Overall)
tdf_top10 <- d_tdf %>% arrange(desc(imp_tdf)) %>% head(10) %>% pull(protein) %>% as.character()

top10_shared <- c(taf_top10, tdf_top10) %>% unique()


d_all <- full_join(d_taf, d_tdf)

p3 <- d_all %>% dplyr::filter(protein %in% top10_shared) %>% 
  dplyr::mutate(max_imp = pmax(imp_taf, imp_tdf)) %>% 
  arrange(desc(max_imp)) %>% 
  dplyr::mutate(protein = factor(protein, levels = rev(protein))) %>% 
  head(10) %>% 
  ggplot() + 
  geom_segment(aes(y = 0, x = `protein`, yend = max_imp, xend = `protein`)) + 
  geom_point(aes(x = `protein`, y = imp_taf), stat='identity', size=2, color= "#f94144") + # red
  geom_point(aes(x = `protein`, y = imp_tdf), stat='identity', size=2, color= "#277da1") + #blue
  labs(x = '', y = 'Importance', title = 'Top 10 predictors (O-link)') +
  ylim(0,25) +
  theme(panel.grid = element_blank()) +
  coord_flip() +
  hrbrthemes::theme_ipsum()

#save figure
ggsave("combined_predictors_tdf_taf_pairs.pdf", plot = p3, device = cairo_pdf, width = 5, height = 5)
