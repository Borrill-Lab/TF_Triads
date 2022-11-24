# Aim is to find out whether TFs in gene families which are retained less as 1:1:1 triads are lower expressed than genes which are retained as 1:1:1 triads

# Philippa Borrill
# 02-12-2020

library("tidyverse")
library("fields")
library("ggrepel")


## do this analysis with Choulet Chinese Spring data 

# input tpms:
tpm_input_file <- "choulet_URGI_tpm.tsv" # tpm file available from https://opendata.earlham.ac.uk/wheat/under_license/toronto/Ramirez-Gonzalez_etal_2018-06025-Transcriptome-Landscape/expvip/RefSeq_1.1/ByGene/

# load tpms
tpms <- read.table(file=tpm_input_file,sep = "\t", header=T)
head(tpms)
dim(tpms)


# load homoeologs
homologies <- read.csv(file="v1.1_genes_TF_homoeolog_info.csv")
head(homologies)
dim(homologies)

# calculate average per condition (i.e. mean of 2 reps):
head(tpms)

tpms_av <- tpms %>%
  mutate(grain_Z71 = (grain_Z71_rep1+grain_Z71_rep2)/2, 
         grain_Z75 = (grain_Z75_rep1+grain_Z75_rep2)/2,
         grain_Z85 = (grain_Z85_rep1+grain_Z85_rep2)/2,
         leaf_Z10 = (leaf_Z10_rep1+leaf_Z10_rep2)/2,
         leaf_Z23 = (leaf_Z23_rep1+leaf_Z23_rep2)/2,
         leaf_Z71 = (leaf_Z71_rep1+leaf_Z71_rep2)/2,
         root_Z10 = (root_Z10_rep1+root_Z10_rep2)/2,
         root_Z13 = (root_Z13_rep1+root_Z13_rep2)/2,
         root_Z39 = (root_Z39_rep1+root_Z39_rep2)/2,
         spike_Z32 = (spike_Z32_rep1+spike_Z32_rep2)/2,
         spike_Z39 = (spike_Z39_rep1+spike_Z39_rep2)/2,
         spike_Z65 = (spike_Z65_rep1+spike_Z65_rep2)/2,
         stem_Z30 = (stem_Z30_rep1+stem_Z30_rep2)/2, 
         stem_Z32 = (stem_Z32_rep1+stem_Z32_rep2)/2, 
         stem_Z65 = (stem_Z65_rep1+stem_Z65_rep2)/2,
         all_tissues = (grain_Z71 + grain_Z75 + grain_Z85 +
                          leaf_Z10 + leaf_Z23 + leaf_Z71 +
                          root_Z10 + root_Z13 + root_Z39 + 
                          spike_Z32 + spike_Z39 + spike_Z65 +
                          stem_Z30 + stem_Z32 + stem_Z65) / 15)%>%
  select(gene, grain_Z71:all_tissues)
head(tpms_av)
dim(tpms_av)

tpms_av <- tpms_av %>%
  rowwise() %>%
  mutate( median_all_tissues = median(grain_Z71, grain_Z75 , grain_Z85 ,
                            leaf_Z10 , leaf_Z23 , leaf_Z71,
                            root_Z10 , root_Z13 , root_Z39 , 
                            spike_Z32 , spike_Z39 , spike_Z65 ,
                            stem_Z30 , stem_Z32 , stem_Z65) ) %>%
  ungroup()
head(tpms_av)
dim(tpms_av)


# tpms_av$all_tissues is the mean across all tissues NB includes genes expressed <0.5 tpm
tpms_av_HC <- tpms_av[!grepl("LC",tpms_av$gene),] # remove LC genes
head(tpms_av_HC)
dim(tpms_av_HC)

# add in TF family info:
tpms_av_HC_TF <- merge(tpms_av_HC,homologies, by.x="gene",by.y="v1.1_ID")
head(tpms_av_HC_TF)
dim(tpms_av_HC_TF)


# load in % of each TF family in triads
TF_fam_in_triads <- read.csv(file=paste0(Y,"polyploidy\\TFs\\2_TFs_vs_gene_homoeologs_iTAK_annotation\\summary_table_of_triad_groups_TFs_vs_non_TF_per_TF_family.csv"), header=T)
head(TF_fam_in_triads)

TF_triads_only <- TF_fam_in_triads %>%
  filter(type == "1:1:1" & !is.na(TF_family.x))
head(TF_triads_only)
dim(TF_triads_only)
tail(TF_triads_only)


### try with median (could solve if some TFs had extreme expression level which skews the results):

#calc average expression level per TF family (using mean values per tissue):
median_expression_per_TF_family <- tpms_av_HC_TF %>%
  select(gene, all_tissues, TF_family.x) %>% 
  filter(!is.na(TF_family.x )) %>%
  group_by(TF_family.x) %>%
  summarise(expr_level = median(all_tissues),
            n_TF_in_family = length(all_tissues) )
head(median_expression_per_TF_family)


## now merge together the % in triads table with the expression level table

merged_table_median <- merge(TF_triads_only, median_expression_per_TF_family, by="TF_family.x")
head(merged_table_median)

p_all_families_log <- ggplot(data=merged_table_median, aes(x=percent, y=log(expr_level))) +geom_point() +
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95) +
  theme_bw() + labs(x="percentage of TF family in 1:1:1 triads", y="median expression level (log(tpm))")

R2_log <- summary(lm(log(expr_level) ~ percent, data=merged_table_median))$r.squared 
R2_log

p_value <- summary(lm(log(expr_level) ~ percent, data=merged_table_median))$coefficients[2,4] 
p_value

p_all_families_log_annot <- p_all_families_log + annotate("text", x = min(merged_table_median$percent)+30, y = max(log(merged_table_median$expr_level)), 
                                                          label =paste0("R^2 = ",round(R2_log,3), ",  p-value =", round(p_value,4)))

ggsave(p_all_families_log_annot, file="median_expression_level_per_family_mean_per_tissue_log_vs_perc_family_in_triads_ggsave.pdf", 
       width=4, height =4, units=c("in"), device = cairo_pdf)


# try with other families with >30 genes
p_30 <- ggplot(data=merged_table_median[merged_table_median$n >30,], aes(x=percent, y=log(expr_level))) +geom_point() +
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95) +
  theme_bw() + labs(x="percentage of TF family in 1:1:1 triads", y="median expression level (log(tpm))") 

R2_log <- summary(lm(log(expr_level) ~ percent, data=merged_table_median[merged_table_median$n >30,]))$r.squared 

p_value <- summary(lm(log(expr_level) ~ percent, data=merged_table_median[merged_table_median$n >30,]))$coefficients[2,4] 

p_30_log_annot <- p_30 + annotate("text", x = min(merged_table_median$percent)+30, y = max(log(merged_table_median$expr_level)), 
                                  label =paste0("R^2 = ",round(R2_log,3), ",  p-value =", round(p_value,4)))


ggsave(p_30_log_annot, file="median_expression_level_per_family_mean_per_tissue_log_vs_perc_family_in_triads_TF_familes_min_30_genes_ggsave.pdf",
       width=4, height=4, units=c("in"), device = cairo_pdf)

## try adding labels

p_30_log_annot_with_labels <- p_30_log_annot+  geom_label_repel(aes(label=TF_family.x),
                                                                min.segment.length = 0,
                                                                box.padding   = 0.1, 
                                                                point.padding = 0.5,
                                                                segment.color = 'grey50', size=2)

ggsave(p_30_log_annot_with_labels, file="median_expression_level_per_family_mean_per_tissue_log_vs_perc_family_in_triads_TF_familes_min_30_genes_ggsave_with_labels.pdf",
       width=4, height=4, units=c("in"), device = cairo_pdf)


