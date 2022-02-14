#Filter for expressed genes, step 6f
#'
#'11/01/2022
#'only use those in choulet_URGI_tpm$gene for transcription factors, these are genes with some detectable expression in at least one tissue
#'Comes after combining the other tables (step 6b), before MAF filters (step 6d) and constructing tf family plot (step 7a)
#'
#'19/01/2022
#'Also use this script to add data about selective sweep sites
#'
#'26/1/2022
#'Also use this script to add data about mutation load and "synonymous only" mutations
#'
#'30/01/2022
#'Use new extreme load file from 2022-01-30
#'
#'Built with R 4.0.5 and tidyverse

#PREP section
require(tidyverse)
require(lubridate)

setwd("/rds/projects/b/borrillp-phd-students/Catherine/TF_homoeolog_SNPs/4_analysis")
combined_table_MAFall <- read_csv("coding_variant_tf_allele_frequency_table_MAFall_2022-01-21.csv")
non_sweep_regions <- read_delim("He19_coding_variant_not_in_sweep_regions.vcf", delim = "\t", col_names = FALSE)
allele_frequency_with_load <- read_csv(file = "allele_frequency_coding_variant_extreme_load_2022-01-30.csv")

setwd("/rds/projects/b/borrillp-phd-students/Catherine/expression_and_load_data")
choulet_URGI_tpm <- read.table(file="choulet_URGI_tpm.tsv",header=T)
sites_to_exclude <- read.table(file="sites_to_exclude_2.txt", header = T)

#Expression across all tissues (?)
choulet_URGI_tpm$all <- rowSums(choulet_URGI_tpm[2:ncol(choulet_URGI_tpm)])
#Filter for expressed genes
choulet_URGI_tpm <- choulet_URGI_tpm[choulet_URGI_tpm$all > 0.5,]
print("Number of expressed genes:")
length(choulet_URGI_tpm$gene)
print("number of non-sweep regions")
length(non_sweep_regions$X3)
print("number of synonymous sites with other annotations")
length(sites_to_exclude$x)

#Mutation load
#This script has run the counts with the 2.5% tails of varieties with extreme load at each end removed
#Which SNPs have homozygous sites in the middle 95% dataset?
alleles_in_extreme_load <- allele_frequency_with_load %>%
  mutate(Extreme = (`1/1`==0 | `0/0`==0)) %>% #Only minor alleles were in extreme load plants
  select(ID, Extreme)

#Add information to table, this allows me to choose later whether to filter on this or not!
combined_table_MAFall_expression <- combined_table_MAFall %>%
  mutate(Expressed = Gene %in% choulet_URGI_tpm$gene) %>%
  mutate(Sweep = !(Uploaded_variation %in% non_sweep_regions$X3)) %>%
  inner_join(alleles_in_extreme_load, by = c("Uploaded_variation"="ID")) %>% #Add the 'Extreme' column
  mutate(Synonly =!(Location %in% sites_to_exclude$x))

print("Counts of SNPs in expressed genes")
combined_table_MAFall_expression %>%
  group_by(Expressed) %>%
  count()
combined_table_MAFall_expression %>%
  group_by(isTF, Expressed) %>%
  count()
combined_table_MAFall_expression %>%
  group_by(TF_family.x, Expressed) %>%
  count()

print("Counts of SNPs with load")
combined_table_MAFall_expression %>%
  group_by(isTF, Extreme) %>%
  count()

setwd("/rds/projects/b/borrillp-phd-students/Catherine/TF_homoeolog_SNPs/4_analysis")
write_csv(combined_table_MAFall_expression, file = paste("coding_variant_tf_allele_frequency_table_MAFall_expression_", lubridate::today(),".csv", sep = ""), col_names = TRUE)


