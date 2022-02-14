#'TF evolution analysis
#'
#'06/09/21
#'Filter coding variants for different MAF levels and remove variants with too many missing values
#'
#'11/01/22
#'Use new file from 2022-01-11; different allele frequencies and added Expression column
#'
#'19/01/22
#'Use new file from 2022-01-19; added Sweep column
#'
#'21/01/22
#'Use new file from 2022-01-21; added Extreme column
#'
#'26/01/22
#'Use new file from 2022-01-26; added Synonly column
#'
#'30/01/22
#'Use new file from 2022-01-30; Insist on lines being unique.
#'
#'Built with R 4.0.5 and tidyverse

#PREP section
require(tidyverse)
require(readr)
require(lubridate)

setwd("/rds/projects/b/borrillp-phd-students/Catherine/TF_homoeolog_SNPs/4_analysis")
combined_table <- read_csv("coding_variant_tf_allele_frequency_table_MAFall_expression_2022-01-30.csv")

sift_labs <- c(missense_variant.deleterious = "deleterious missense", missense_variant.tolerated = "tolerated missense", stop_gained.na = "stop gained", synonymous_variant.na = "synonymous")
sift_order <- c("stop_gained.na", "missense_variant.deleterious", "missense_variant.tolerated","synonymous_variant.na")
#Some last minute table tweaking
combined_table <- combined_table %>%
#  unique %>%
  mutate(SIFT_prediction = replace_na(SIFT_prediction, replace = "na")) %>% #NAs were propagating and messing up the plot
  mutate(SIFT_allconf = as.factor(str_split_fixed(SIFT_prediction, "_", 2)[, 1])) %>% #split the first word off the front
  mutate(interaction = fct_relevel(interaction(Consequence, SIFT_allconf), sift_order)) #make a new variable for my four categories, in the right order


total_variants <- sum(combined_table$`./.`[1], combined_table$`0/0`[1], combined_table$`0/1`[1], combined_table$`1/1`[1])
threshold <- 0.25

combined_table_F25 <- combined_table %>%
  filter((`./.`/total_variants) <= threshold)

combined_table_MAF0.05 <- combined_table_F25 %>%
  filter(MAF>=0.05)

combined_table_MAF0.01 <- combined_table_F25 %>%
  filter(MAF>=0.01)

print(paste("Variants in triad genes:", nrow(combined_table)))
print(paste("Variants in triad genes with <=25% missing values:", nrow(combined_table_F25)))
print(paste("Variants in triad genes with MAF >=0.05", nrow(combined_table_MAF0.05)))
print(paste("Variants in triad genes with MAF >=0.01", nrow(combined_table_MAF0.01)))

write_csv(combined_table_F25, file = paste("coding_variant_tf_allele_frequency_table_F25_", lubridate::today(),".csv", sep = ""), col_names = TRUE)
write_csv(combined_table_MAF0.05, file = paste("coding_variant_tf_allele_frequency_table_MAF0.05_", lubridate::today(),".csv", sep = ""), col_names = TRUE)
write_csv(combined_table_MAF0.01, file = paste("coding_variant_tf_allele_frequency_table_MAF0.01_", lubridate::today(),".csv", sep = ""), col_names = TRUE)
