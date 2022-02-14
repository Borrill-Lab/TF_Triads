#'TF evolution analysis
#'7/9/21
#'Let's take coding_variants_chisquared_step6_c.R and make the plotting and chi squared into functions
#'which I can call with as many data frame + variables combos as I want
#'This script calls the functions from functions_graph_chisquared.R
#'
#'21/9/21
#'Run with SNPs filtered to remove SNPs with >25% missing values from 2021-09-21.
#'Plot stacked bar plots by TF family, including TF families containing >30 genes.
#'
#'7/10/21
#'Sort by deleterious and stop SNPs and add non-TF genes
#'
#'#'18/10/21
#'Save matching summary count and proportion tables
#'
#'11/01/22
#'Use new input file 2022-01-11. Plot expressed genes only.
#'
#'19/01/22
#'Use new input file 2022-01-19. Plot SNPs NOT in sweep regions and MAF>0.01
#'
#'21/01/22
#'Use new input file 2022-01-21. Plot SNPs represented in the middle 95% of plants, not only in the plants with extreme mutation load.
#'
#'27/01/22
#'Use new input file 2022-01-27. Plot synonymous SNPs which are synonymous ONLY.
#'
#'30/01/22
#'Use new input file 2022-01-30. Reduced graph size to 80*115mm
#'
#'09/02/22
#'Plot TF families which have >=5 SNPs in the graph
#'
#'Built with R 4.0.5 and tidyverse

#PREP
require(tidyverse)
require(lubridate)
setwd("/rds/projects/b/borrillp-phd-students/Catherine/TF_homoeolog_SNPs/scripts")
source("functions_graph_chisquared.R")

#Read in tables
setwd("/rds/projects/b/borrillp-phd-students/Catherine/TF_homoeolog_SNPs/1_raw")
big_TF_families <- read_csv("TF_families_with_over_30_genes.csv")

setwd("/rds/projects/b/borrillp-phd-students/Catherine/TF_homoeolog_SNPs/4_analysis")
combined_table_F25 <- read_csv("coding_variant_tf_allele_frequency_table_MAF0.01_2022-02-10.csv")

#FILTER
#Keep SNPs in TFs in families with over 30 genes
summary(as.factor(combined_table_F25$TF_family.x))
summary(combined_table_F25$Expressed)
summary(combined_table_F25$Sweep)
summary(combined_table_F25$Extreme)
summary(combined_table_F25$Synonly)
# combined_table_big_TF_families <- combined_table_F25 %>%
#   filter(TF_family.x %in% unlist(big_TF_families)) #want to keep non-TFs as well as TFs

#RUN
table <- combined_table_F25 #skip filter step
sift_labs <- c(missense_variant.deleterious = "deleterious missense", missense_variant.tolerated = "tolerated missense", stop_gained.na = "stop gained", synonymous_variant.na = "synonymous")
sift_order <- c("stop_gained.na", "missense_variant.deleterious", "missense_variant.tolerated","synonymous_variant.na")
#Some last minute table tweaking
 table <- table %>%
   unique %>%
#   mutate(SIFT_prediction = replace_na(SIFT_prediction, replace = "na")) %>% #NAs were propagating and messing up the plot
#   mutate(SIFT_allconf = as.factor(str_split_fixed(SIFT_prediction, "_", 2)[, 1])) %>% #split the first word off the front
   mutate(interaction = fct_relevel(interaction(Consequence, SIFT_allconf), sift_order)) #make a new variable for my four categories, in the right order


table <- table %>%
  filter(Expressed == TRUE & Sweep == FALSE & Extreme == FALSE & Synonly == TRUE)
summary(table) #Just check how many we have

summary(table$interaction)
table_30 <- table %>%
  filter(isTF == FALSE | TF_family.x %in% unlist(big_TF_families)) %>% #want to keep non-TFs as well as TFs in big families
  mutate(TF_family.reorder = replace_na(fct_reorder(TF_family.x, interaction, .fun='prop.deleterious'), replace = "non-TF")) #order TF families by the proportion of deleterious & stop mutations. custom function.

summary(table_30) #Check things are working effectively

#Replot without families containing very few SNPs
count_SNP <- table %>%
  group_by(TF_family.x) %>%
  count()
print("Counts of SNPs")
print(count_SNP)
SNP5 <- count_SNP %>%
  filter(n>=5) %>% #At least 5 SNPs per TF family
  select(TF_family.x) %>%
  ungroup() %>%
  unlist()
SNP5

table_SNP5 <- table %>%
  filter(isTF == FALSE | (TF_family.x %in% unlist(big_TF_families) & TF_family.x %in% SNP5)) %>% #want to keep non-TFs as well as TFs in big families
  mutate(TF_family.reorder = replace_na(fct_reorder(TF_family.x, interaction, .fun='prop.deleterious'), replace = "non-TF")) #order TF families by the proportion of deleterious & stop mutations. custom function.
  
# table_SNP5 <- table_30 %>%
#   filter(isTF == FALSE | (TF_family.x %in% SNP5))

summary(table_SNP5) #Check things are working effectively

#First plot
description <- "TF_families_gene30_MAF0.01"
stacked_plot <- return.stacked.plot(table_30,
                  x = "TF_family.reorder",
                  fill = "interaction",
                  description = description,
                  invert = TRUE) + #Put TF families on the vertical axis to match graphs in paper
  scale_fill_brewer(limits = sift_order, labels = sift_labs, palette = "Reds", direction = -1) + #legend labels
  labs(x = "Transcription Factor Family", y = "Proportion of SNPs in group", fill = "")
setwd("/rds/projects/b/borrillp-phd-students/Catherine/TF_homoeolog_SNPs/figures")
  # theme(legend.margin = margin(10,10,10,10, unit = "mm"))

ggsave(file = paste("stacked_plot_of_", description, "_", lubridate::today(),".pdf", sep = ""), plot = stacked_plot, device = "pdf", width = 80, height = 115, units = "mm")

description <- "TF_families_SNP5_MAF0.01"
stacked_plot <- return.stacked.plot(table_SNP5,
                                    x = "TF_family.reorder",
                                    fill = "interaction",
                                    description = description,
                                    invert = TRUE) + #Put TF families on the vertical axis to match graphs in paper
  scale_fill_brewer(limits = sift_order, labels = sift_labs, palette = "Reds", direction = -1) + #legend labels
  labs(x = "Transcription Factor Family", y = "Proportion of SNPs in group", fill = "")
setwd("/rds/projects/b/borrillp-phd-students/Catherine/TF_homoeolog_SNPs/figures")
  # theme(legend.margin = margin(10,10,10,10, unit = "mm"))

ggsave(file = paste("stacked_plot_of_", description, "_", lubridate::today(),".pdf", sep = ""), plot = stacked_plot, device = "pdf", width = 80, height = 115, units = "mm")

description <- "TF_families_SNP5_MAF0.01_legend"
stacked_plot <- return.stacked.plot(table_SNP5,
                                    x = "TF_family.reorder",
                                    fill = "interaction",
                                    description = description,
                                    invert = TRUE, #Put TF families on the vertical axis to match graphs in paper
                                    legend.position = "right") + #add legend here
  scale_fill_brewer(limits = sift_order, labels = sift_labs, palette = "Reds", direction = -1) + #legend labels
  labs(x = "Transcription Factor Family", y = "Proportion of SNPs in group", fill = "")
setwd("/rds/projects/b/borrillp-phd-students/Catherine/TF_homoeolog_SNPs/figures")
  # theme(legend.margin = margin(10,10,10,10, unit = "mm"))
  
ggsave(file = paste("stacked_plot_of_", description, "_", lubridate::today(),".pdf", sep = ""), plot = stacked_plot, device = "pdf", width = 80, height = 115, units = "mm")


#COUNTS TABLE
count_tf_family_table <- table_30 %>%
  group_by(Consequence, SIFT_allconf, TF_family.reorder) %>% #same order of TF family as figures
  count() %>%
  pivot_wider(names_from = TF_family.reorder, values_from = n)
proportion_table <- cbind(count_tf_family_table[1:2], proportions(count_tf_family_table[-(1:2)])) #proportions calculates the proportion of the column




write_csv(count_tf_family_table, file = paste("count_table_tf_family_", lubridate::today(), ".csv", sep = ""), col_names = TRUE)
write_csv(proportion_table, file = paste("proportion_table_tf_family_", lubridate::today(), ".csv", sep = ""), col_names = TRUE)

