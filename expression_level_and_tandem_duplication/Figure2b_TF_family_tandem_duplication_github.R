R# Philippa Borrill
# 02-12-2020

# Aim is to find out whether gene families which are retained less as 1:1:1 are more prone to tandem duplication

library("tidyverse")
library("fields")
#install.packages("ggrepel")
library("ggrepel")


# load homoeologs
homologies <- read.csv(file="v1.1_genes_TF_homoeolog_info.csv"))
head(homologies)
dim(homologies)

# find out % of each TF family which is tandem duplicates
TFs <- homologies %>%
  filter(!is.na(TF_family.x))
head(TFs)
dim(TFs)

TFs_tandem_dup_info <- TFs %>%
  mutate(chr = substr(v1.1_ID,1,12),
         geneID = substr(v1.1_ID,13,18), # these rows extract the gene ID numeric part
         geneID_numeric = (as.numeric(geneID))) %>%
  select(TF_family.x,chr, geneID_numeric) %>%
  group_by(TF_family.x, chr) %>% # group by TF family and chromosome
  arrange(geneID_numeric, .by_group=TRUE) %>%
  mutate(diff = geneID_numeric - lag(geneID_numeric, default = first(geneID_numeric)), # this line calculates the difference between subsequent rows which are in the same TF family and have the same chr
         tandem_dup = (diff == "100" | diff == "200"| diff == "300" )) # adjust this line to be diff=="100" for only adjacent; # adjust this line to be diff=="100" | diff=="200" for adjacent plus with 1 in between gene # leaving the line as it is gives you tandem duplicates which can have up to 2 genes between

head(TFs_tandem_dup_info,10)
tail(TFs_tandem_dup_info,10)  

only_tandem_dup_TFs <- TFs_tandem_dup_info %>%
  filter(tandem_dup == TRUE) %>%
  group_by(TF_family.x) %>%
  summarise(num = sum(tandem_dup)) # counts how many genes are tandem dup per TF family NB this is just the duplicated copy i.e. there is another copy of each of these genes

only_tandem_dup_TFs

# now want to normalise for size of TF family

head(TFs)
total_TF <- TFs %>%
  group_by(TF_family.x) %>%
  summarise(total_family_size = sum(isTF))

total_TF


tandem_dup_TF_norm <- merge(total_TF, only_tandem_dup_TFs, by = "TF_family.x", all.x=T)
tandem_dup_TF_norm

final_tandem_dup <- tandem_dup_TF_norm %>%
  rename(num_tandem = num) %>%
  replace_na(list(num_tandem = 0)) %>%
  mutate(perc_tandem_dup = num_tandem/total_family_size*100)

final_tandem_dup
dim(final_tandem_dup)

# calc % of each TF family in triads

per_TF_family_hom_category <- homologies %>%
  group_by(TF_family.x, type, isTF) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  group_by(TF_family.x) %>% 
  mutate(percent = n/sum(n)*100) # add column with percentage of genes in that category per family

per_TF_family_hom_category     


TF_triads_only <- per_TF_family_hom_category %>%
  filter(type == "1:1:1" & !is.na(TF_family.x))
head(TF_triads_only)
dim(TF_triads_only)
tail(TF_triads_only)

# merge together % in triads and % tandem dup

merged_data <- merge(TF_triads_only, final_tandem_dup)
head(merged_data)
dim(merged_data)


## plot some graphs

p_all_families <- ggplot(data=merged_data, aes(x=percent, y=perc_tandem_dup)) +geom_point() +
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95) +
  theme_bw() + labs(x="percentage of TF family in 1:1:1 triads", y="percentage of genes which are tandem duplicates")

R2 <- summary(lm((perc_tandem_dup) ~ percent, data=merged_data))$r.squared 
R2

p_value <- summary(lm((perc_tandem_dup) ~ percent, data=merged_data))$coefficients[2,4] 
p_value

p_all_families_annot <- p_all_families + annotate("text", x = min(merged_data$percent)+30, y = max((merged_data$perc_tandem_dup)-3), 
                                                  label =paste0("R^2 = ",round(R2,3), ",  p-value =", round(p_value,4)))

p_all_families_annot


ggsave(p_all_families_annot, file="perc_tandem_dup_vs_perc_family_in_triads_ggsave.pdf", width=4, height =4,
       units=c("in"), device = cairo_pdf)


# trying families with >30 genes 

p_30 <- ggplot(data=merged_data[merged_data$n>30,], aes(x=percent, y=perc_tandem_dup)) +geom_point() +
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95) +
  theme_bw() + labs(x="percentage of TF family in 1:1:1 triads", y="percentage of genes which are tandem duplicates")

R2 <- summary(lm((perc_tandem_dup) ~ percent, data=merged_data[merged_data$n>30,]))$r.squared 
R2

p_value <- summary(lm((perc_tandem_dup) ~ percent, data=merged_data[merged_data$n>30,]))$coefficients[2,4] 
p_value

p_30_annot <- p_30 + annotate("text", x = min(merged_data$percent)+30, y = max((merged_data$perc_tandem_dup)-3), 
                              label =paste0("R^2 = ",round(R2,3), ",  p-value =", round(p_value,4)))

p_30_annot_NAC_B3_MADSM_labelled <- p_30_annot + geom_label_repel(data = merged_data[merged_data$n>30 & merged_data$TF_family.x %in% c("NAC","MADS-M-type","B3"),], 
                            aes(x=percent, y=perc_tandem_dup, label=TF_family.x),
                            box.padding   = 0.1, 
                            point.padding = 0.5,
                            segment.color = 'grey50')

pdf(file="perc_tandem_dup_vs_perc_family_in_triads_min_30_genes_selected_TF_labelled.pdf", width=4, height =4)
p_30_annot_NAC_B3_MADSM_labelled
dev.off()

ggsave(p_30_annot, file="perc_tandem_dup_vs_perc_family_in_triads_min_30_genes_ggsave.pdf", width=4, height =4,
       units=c("in"), device = cairo_pdf)

