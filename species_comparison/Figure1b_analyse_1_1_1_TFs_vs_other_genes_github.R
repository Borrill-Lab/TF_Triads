# Aim is to analyse TFs that are 1:1:1 compared to other genes in v1.1 annotation
# 
# Philippa Borrill
# 01-06-2020 #updated 10-07-20 for iTAK annotations
# updated 25-10-2021 with statistical test

# load required packages and setwd

library(ggplot2)
library(tidyverse)


##### need to load in required files #####

tidy_all_info <- read.csv("v1.1_genes_TF_homoeolog_info.csv")
head(tidy_all_info)

# make some summary tables to found out of TFs have more 1:1:1 than other genes

summary_table <- tidy_all_info %>%
  group_by(type, isTF) %>%
  summarise(n=n()) 

summary_table

write.csv(file="summary_table_of_triad_groups_TFs_vs_non_TF.csv",as.data.frame(summary_table), row.names = F)

#calc stats
overall_TF_vs_nonTF_contingency <- tidy_all_info %>%
  mutate("is_triad" = type=="1:1:1") %>%
  replace_na(list(is_triad ="FALSE")) %>%
  group_by(isTF,is_triad) %>%
  summarise(n=n())  %>%
  pivot_wider(names_from = isTF, values_from=n, names_prefix = "isTF_")

overall_TF_vs_nonTF_contingency

overall_fisher_res <- fisher.test(overall_TF_vs_nonTF_contingency[,2:3])
overall_fisher_res

# want to know how many TFs in each family are in 1:1:1  triads, and what % of this family is in that category

per_TF_family_hom_category <- tidy_all_info %>%
  group_by(TF_family.x, type, isTF) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  group_by(TF_family.x) %>% 
  mutate(percent = n/sum(n)*100) # add column with percentage of genes in that category per family

per_TF_family_hom_category            

write.csv(file="summary_table_of_triad_groups_TFs_vs_non_TF_per_TF_family.csv",as.data.frame(per_TF_family_hom_category), row.names = F, quote=T)

# rename type to cardinality_abs to match previous dataframe

per_TF_family_hom_category <-  per_TF_family_hom_category %>%
  rename(cardinality_abs = type)

per_TF_family_hom_category

# need to replace implicit NA with "other" in cardinality_abs (i.e. other homoeolog category)
per_TF_family_hom_category$cardinality_abs <- replace_na(per_TF_family_hom_category$cardinality_abs, "other")
per_TF_family_hom_category

per_TF_family_hom_category$cardinality_abs <- as.factor(per_TF_family_hom_category$cardinality_abs)
per_TF_family_hom_category


# try plotting as lollipop plot to get gene numbers showing as colour just 1:1:1

triad_lollipop_num_colour <- per_TF_family_hom_category %>%
  filter(cardinality_abs == "1:1:1" &
           !is.na(TF_family.x)) %>%
  ggplot(aes(x=reorder(TF_family.x, percent), y=percent)) +
  geom_segment( aes(xend=TF_family.x, yend=0)) +
  geom_point(aes(fill=n), colour="black", alpha =1, pch=21, size=3) +
  scale_fill_gradient(low = "white", high = "red", limits=c(0,500)) +
  coord_flip() + geom_hline(yintercept = 55.9, linetype="dashed") +
  
  theme_bw() +
  ylab("Percent") +
  xlab("Transcription Factor Family") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 102))

triad_lollipop_num_colour

pdf(file = "triad_lollipop_num_colour.pdf", width=8, height=8)
triad_lollipop_num_colour
dev.off()

ggsave(triad_lollipop_num_colour, filename = "triad_lollipop_num_colour_ggsave.pdf", device = cairo_pdf, 
       width = 8, height = 8, units = "in")

# only keep gene families >30 genes - for ggsave
head(per_TF_family_hom_category)

triad_lollipop_num_colour_min_30_genes <- per_TF_family_hom_category %>%
  filter(cardinality_abs == "1:1:1" &
           !is.na(TF_family.x) &
           n >30) %>%
  ggplot(aes(x=reorder(TF_family.x, percent), y=percent)) +
  geom_segment( aes(xend=TF_family.x, yend=0)) +
  geom_point(aes(fill=n), colour="black", alpha =1, pch=21, size=3) +
  scale_fill_gradient(low = "white", high = "red", limits=c(0,400)) +
  coord_flip() + geom_hline(yintercept = 55.9, linetype="dashed") +
  
  theme_bw() +
  ylab("Percent") +
  xlab("Transcription Factor Family") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 105))

triad_lollipop_num_colour_min_30_genes

ggsave(triad_lollipop_num_colour_min_30_genes, filename = "triads_lollipop_perc_TF_family_in_category_min_30_genes_ggsave_max400.pdf", device = cairo_pdf, 
       width = 4, height = 5, units = "in")

### now do stats to see if the number in each family in triads is different from non-TFs 

per_TF_family_hom_category %>%
  filter(cardinality_abs == "1:1:1" &
           is.na(TF_family.x))

per_TF_family_hom_category %>%
  filter(cardinality_abs == "1:1:1" &
           !is.na(TF_family.x))

# first need to make a table which has the number per TF family which is in triads, and number not in triads

per_TF_family_hom_category %>% 
  select(TF_family.x, cardinality_abs, isTF, n) %>%
  mutate("is_1_1_1" = cardinality_abs == "1:1:1") %>%
  group_by(TF_family.x, is_1_1_1) %>% 
  mutate(n_1_1_1 = sum(n)) %>%
  head(20)
  
# now summarise to only have 1 row per non 1:1:1 per TF family
data_for_fisher_test <- per_TF_family_hom_category %>% 
  select(TF_family.x, cardinality_abs, isTF, n) %>%
  mutate("is_1_1_1" = cardinality_abs == "1:1:1") %>%
  group_by(TF_family.x, is_1_1_1) %>% 
  summarise(n_1_1_1 = sum(n)) %>%
  pivot_wider(names_from = is_1_1_1, values_from = n_1_1_1,
              names_prefix = "is_triad_") %>% # make into wide format 
  replace_na(list(is_triad_FALSE=0,is_triad_TRUE=0 )) %>%
  mutate(total_n = is_triad_FALSE + is_triad_TRUE)
data_for_fisher_test
  
#try for 1 family:
df_1_family <- data_for_fisher_test %>%
  ungroup() %>%
  filter(TF_family.x == "Alfin-like" | is.na(TF_family.x)) %>%
  select(is_triad_FALSE,is_triad_TRUE)

fisher.result <- fisher.test(df_1_family)
print(fisher.result$p.value)

# now do for each family vs non-TFs and then correct for multiple testing

fisher_res <- data.frame(TF_family = as.character(), pvalue = as.numeric() ) #setup empty vector
fisher_res

for (i in data_for_fisher_test$TF_family.x[1:66]) {
  fisher.result <- fisher.test(data_for_fisher_test %>%
    ungroup() %>%
    filter(TF_family.x == i | is.na(TF_family.x)) %>%
    select(is_triad_FALSE,is_triad_TRUE)
  )
  fisher_res[nrow(fisher_res)+1,] <- c(TF_family = i, pvalue = fisher.result$p.value )
}
fisher_res

fisher_res$padj <- p.adjust(fisher_res$pvalue, method = "fdr", n = length(fisher_res$pvalue))
fisher_res
fisher_res$sig_padj <- fisher_res$padj <0.05
fisher_res

write.csv(fisher_res, file="Fisher_exact_test_TF_vs_nonTF_triads_per_TF_family.csv",row.names = F)
