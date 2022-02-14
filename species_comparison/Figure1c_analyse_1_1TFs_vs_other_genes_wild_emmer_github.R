# Aim is to analyse TFs that are 1:1 compared to other genes in 
#
# Philippa Borrill
# 05-06-2020

# load required packages and setwd

library(ggplot2)
library(tidyverse)

##### need to load in required files #####

tidy_all_info <- read.csv(file="WEW_genes_TF_homoeolog_info.csv")
head(tidy_all_info)

# make some summary tables to found out of TFs have more 1:1:1 than other genes

summary_table <- tidy_all_info %>%
  group_by(isHomoeolog, isTF) %>%
  summarise(n=n()) 

summary_table

write.csv(file="WEW_summary_table_of_homoeolog_groups_TFs_vs_non_TF.csv",as.data.frame(summary_table), row.names = F)

#calc stats
overall_TF_vs_nonTF_contingency <-summary_table %>%
  pivot_wider(names_from = isTF, values_from=n, names_prefix = "isTF_")

overall_TF_vs_nonTF_contingency

overall_fisher_res <- fisher.test(overall_TF_vs_nonTF_contingency[,2:3])
overall_fisher_res

# want to know how many TFs in each family are in 1:1:1  triads, and what % of this family is in that category

per_TF_family_hom_category <- tidy_all_info %>%
  group_by(TF_family, isHomoeolog, isTF) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  group_by(TF_family) %>% 
  mutate(percent = n/sum(n)*100) # add column with percentage of genes in that category per family

per_TF_family_hom_category            

# rename type to cardinality_abs to match previous dataframe

per_TF_family_hom_category <-  per_TF_family_hom_category %>%
  mutate(cardinality_abs = case_when (isHomoeolog =="TRUE" ~ "1:1",
                            isHomoeolog=="FALSE"~ "other"))

per_TF_family_hom_category

# need to replace implicit NA with "other" in cardinality_abs (i.e. other homoeolog category)
per_TF_family_hom_category$cardinality_abs <- as.factor(per_TF_family_hom_category$cardinality_abs)
per_TF_family_hom_category


# try plotting as lollipop plot to get gene numbers showing as colour just 1:1

diad_lollipop_num_colour <- per_TF_family_hom_category %>%
  filter(cardinality_abs == "1:1" &
           !is.na(TF_family)) %>%
  ggplot(aes(x=reorder(TF_family, percent), y=percent)) +
  geom_segment( aes(xend=TF_family, yend=0)) +
  geom_point(aes(fill=n), colour="black", alpha =1, pch=21, size=3) +
  scale_fill_gradient(low = "white", high = "red", limits=c(0,250)) +
  coord_flip() + geom_hline(yintercept = 69.2, linetype="dashed") +
  
  theme_bw() +
  ylab("Percent") +
  xlab("Transcription Factor Family") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 102))

diad_lollipop_num_colour

pdf(file = "WEW_diads_lollipop_perc_TF_family_in_category.pdf", width=8, height=8)
diad_lollipop_num_colour
dev.off()

ggsave(diad_lollipop_num_colour, filename = "WEW_diads_lollipop_perc_TF_family_in_category_ggsave.pdf", device = cairo_pdf, 
       width = 8, height = 8, units = "in")


#adjust for ggsave - extend x axis and make max n = 400 to match triads
over20_genes_diad_lollipop_num_colour <- per_TF_family_hom_category %>%
  filter(cardinality_abs == "1:1" &
           !is.na(TF_family) &
           n >20) %>%
  ggplot(aes(x=reorder(TF_family, percent), y=percent)) +
  geom_segment( aes(xend=TF_family, yend=0)) +
  geom_point(aes(fill=n), colour="black", alpha =1, pch=21, size=3) +
  scale_fill_gradient(low = "white", high = "red", limits=c(0,400)) +
  coord_flip() + geom_hline(yintercept = 69.2, linetype="dashed") +
  
  theme_bw() +
  ylab("Percent") +
  xlab("Transcription Factor Family") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 105))

over20_genes_diad_lollipop_num_colour

ggsave(over20_genes_diad_lollipop_num_colour, filename = "WEW_diads_lollipop_perc_TF_family_in_category_over20_genes_ggsave_max400.pdf", device = cairo_pdf, 
       width = 4, height = 5, units = "in")



### now do stats to see if the number in each family in triads is different from non-TFs 

# first need to make a table which has the number per TF family which is in triads, and number not in triads

per_TF_family_hom_category %>% 
  select(TF_family, cardinality_abs, isTF, n) %>%
  mutate("is_1_1" = cardinality_abs == "1:1") %>%
  group_by(TF_family, is_1_1) %>% 
  mutate(n_1_1 = sum(n)) %>%
  head(20)

# now summarise to only have 1 row per non 1:1:1 per TF family
data_for_fisher_test <- per_TF_family_hom_category %>% 
  select(TF_family, cardinality_abs, isTF, n) %>%
  mutate("is_1_1" = cardinality_abs == "1:1") %>%
  group_by(TF_family, is_1_1) %>% 
  summarise(n_1_1 = sum(n)) %>%
  pivot_wider(names_from = is_1_1, values_from = n_1_1,
              names_prefix = "is_diad_") %>% # make into wide format 
  replace_na(list(is_diad_FALSE=0,is_diad_TRUE=0 )) %>%
  mutate(total_n = is_diad_FALSE + is_diad_TRUE)
data_for_fisher_test

#try for 1 family:
df_1_family <- data_for_fisher_test %>%
  ungroup() %>%
  filter(TF_family == "Alfin-like" | is.na(TF_family)) %>%
  select(is_diad_FALSE,is_diad_TRUE)

fisher.result <- fisher.test(df_1_family)
print(fisher.result$p.value)

# now do for each family vs non-TFs and then correct for multiple testing 

fisher_res <- data.frame(TF_family = as.character(), pvalue = as.numeric() ) #setup empty vector
fisher_res

for (i in data_for_fisher_test$TF_family[1:63]) {
  fisher.result <- fisher.test(data_for_fisher_test %>%
                                 ungroup() %>%
                                 filter(TF_family == i | is.na(TF_family)) %>%
                                 select(is_diad_FALSE,is_diad_TRUE)
  )
  fisher_res[nrow(fisher_res)+1,] <- c(TF_family = i, pvalue = fisher.result$p.value )
}
fisher_res

fisher_res$padj <- p.adjust(fisher_res$pvalue, method = "fdr", n = length(fisher_res$pvalue))
fisher_res
fisher_res$sig_padj <- fisher_res$padj <0.05
fisher_res

write.csv(fisher_res, file="WEW_Fisher_exact_test_TF_vs_nonTF_diads_per_TF_family.csv",row.names = F)
