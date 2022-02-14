#'10/08/21
#'SIFT scores
#'Extract and use SIFT scores from the VEP data
#'Let's start trying with a small subset of the data
#'
#'11/08/21
#'This is the data preparation script.
#'Moved in another section to combine with TF table.
#'
#'24/08/21
#'Repeat for coding variant data (but without SIFT score part)
#'
#'03/09/21
#'Actually I do want the SIFT score part.
#'Also specify col types, remove multi-type variants and variants without SIFT scores, and print summary statistics to stdout
#'
#'06/09/21
#'ALSO combine with allele frequency table.
#'
#'21/01/22
#'Use new input file from 2022-01-21.
#'
#'Built with R 4.0.5 and tidyverse


#PREP
#Load packages
require(tidyverse)
require(readr)
require(lubridate)

#Read in file
getwd() #should be in scripts folder
TF_table <- read_delim(file = "../2_triads/triads.txt", delim = "\t", col_names = c("v1.1_ID","TF_family.x","group_num","type","isTF"))
vep_colnames = c("Uploaded_variation","Location","Allele","Gene","Feature","Feature_type", "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "Extra")
vep_coltypes = "ccccfffiiicccc" #c = character, f = factor, i = integer
variant_effect_coding_variant <- read_delim(file = "../3_filtered/variant_effect_coding_variant_all.txt", delim = "\t", skip = 32, col_names = vep_colnames, col_types = vep_coltypes)
allele_frequency_table <- read_csv(file = "../4_analysis/allele_frequency_coding_variant_2022-01-21.csv")

#TIDY
#The SIFT data is in col variant_effect_missense_variant$Extra
#Extract using regular expressions - see https://regexone.com
#Also see https://r4ds.had.co.nz/strings.html#extract-matches 

#Extract SIFT scores and other columns from Extra
sift_table_coding_variant <- variant_effect_coding_variant %>%
  tidyr::extract(Extra, c("Impact", "Canonical"), "IMPACT=([A-Z]+).*CANONICAL=([A-Z]+)", remove = FALSE) %>%
  tidyr::extract(Extra, c("SIFT_prediction","SIFT_score"), "SIFT=([a-z_]+).([\\d|\\.]+).", remove = FALSE)

#Combine with TF table
combined_table <- sift_table_coding_variant %>%
  inner_join(TF_table, by = c("Gene"="v1.1_ID")) #inner join gets rid of the ones that aren't in a triad gene

combined_table_allele_freq <- combined_table %>%
  inner_join(allele_frequency_table, by = c("Uploaded_variation"="ID")) #all the SNPs SHOULD be in this table

#Select pertinent columns
combined_table_narrow <- combined_table_allele_freq %>%
  select(Uploaded_variation, Location, Allele, Gene, Consequence, SIFT_prediction, SIFT_score, isTF, TF_family.x, `./.`,`0/0`,`0/1`,`1/1`,`MAF`)

#Filter
combined_table_F1 <- combined_table_narrow %>%
  filter(Consequence %in% c("synonymous_variant", "missense_variant","stop_gained")) #No multi Consequence types
multi_type_SNPs <- combined_table_narrow %>%
  filter(!(Consequence %in% c("synonymous_variant", "missense_variant","stop_gained")))

combined_table_F2 <- combined_table_F1 %>%
  filter(!(Consequence == "missense_variant" & is.na(SIFT_prediction))) #No missense variants without SIFT scores
no_sift_score <- combined_table_F1 %>%
  filter(Consequence == "missense_variant" & is.na(SIFT_prediction))


#CHECK
print("Summary tables: full and filtered")
summary(combined_table_allele_freq)
summary(combined_table_F2)
print(paste("Total variant effects", nrow(variant_effect_coding_variant)))
print(paste("Variant effects in triad genes: ", nrow(combined_table)))
print(paste("Variant effects in triad genes with allele frequency: ", nrow(combined_table_allele_freq)))
print(paste("Variants with a single effect in triad genes: ", nrow(combined_table_F1)))
print(paste("Variants with a single effect in triad genes, excluding missense variants without SIFT scores: ", nrow(combined_table_F2)))
print("Multi-type SNP variants")
head(multi_type_SNPs, n=5)
summary(multi_type_SNPs)
print("Missense variants with missing SIFT scores")
head(no_sift_score, n=5)
summary(no_sift_score)

#OUTPUT
write_csv(combined_table_allele_freq, file = paste("../4_analysis/coding_variant_tf_allele_frequency_table_", lubridate::today(),".csv", sep = ""), col_names = TRUE)
write_csv(combined_table_F2, file = paste("../4_analysis/coding_variant_tf_allele_frequency_table_MAFall_", lubridate::today(),".csv", sep = ""), col_names = TRUE)
