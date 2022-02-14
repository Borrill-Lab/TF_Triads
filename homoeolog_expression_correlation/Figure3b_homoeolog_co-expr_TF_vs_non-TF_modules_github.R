# Philippa Borrill
# 02-12-2020
# edited 11-05-2021
# add stats test 4-11-2021

# Aim is to find out whether TFs have higher or lower co-expr within their triad than non-TFs

library("tidyverse")
library("fields")
#install.packages("ggrepel")
library("ggrepel")



##### First look at 850 WGCNA network #####
modules <- read.csv(file="genes_with_modules_mergeCutHeight0.15_no_expr_values.csv")
head(modules)

# problem - this is in v1.0 annotation but the rest of my analysis is v1.1

# convert from v1.0 to v1.1 
head(gsub("01G", "02G", modules$X))

modules$X <- (gsub("01G", "02G", modules$X))
head(modules)
dim(modules)

length(unique(modules$X)) # number of HC genes with go terms before removing ones which don't match v1.0 to v1.1


# only keep genes which were >99 % ID > 90% coverage from v1.0 to v1.1 # avoids erroneous transfer of GO terms from v1.0 to v1.1
genes_to_transfer <- read.csv(file=paste0(Y,"\\polyploidy\\TFs\\expression_browser\\genes_to_transfer_qcov90_pident99_same_ID.csv"))
head(genes_to_transfer)

modules <- modules[modules$X %in% genes_to_transfer$gene_v1.1,]
head(modules)
dim(modules)


# now load in triads 

# load homoeologs
homologies <- read.csv(file="v1.1_genes_TF_homoeolog_info.csv")
head(homologies)
dim(homologies)

triad_homologies <- homologies %>%
  filter(type == "1:1:1")
head(triad_homologies)

dim(triad_homologies)

# add module info to triads
head(triad_homologies)
head(modules)

triads_modules <- merge(triad_homologies, modules, by.x= "v1.1_ID", by.y="X")
head(triads_modules)
dim(triads_modules)

# want to select only triads where all 3 hom are assigned in WGCNA
triads_with_3_homoeologs_in_network <-  triads_modules %>%
  group_by(group_num) %>%
  summarise(count = length(group_num)) %>%
  filter(count == 3)
head(triads_with_3_homoeologs_in_network)
dim(triads_with_3_homoeologs_in_network)


triads_modules_3_hom <- triads_modules[triads_modules$group_num %in% triads_with_3_homoeologs_in_network$group_num,]
head(triads_modules_3_hom)
dim(triads_modules_3_hom)
triads_modules_3_hom %>%
  group_by(group_num) %>%
  summarise(count = length(group_num)) # check they have 3 per group

# now calculate per TF family (or non-TF)

#want to add total num of genes in that family (that were in the network as triads)

head(triads_modules_3_hom)
num_triads_per_family <- triads_modules_3_hom %>%
  group_by(TF_family.x) %>%
  summarise(num_triads_in_family = length(v1.1_ID)/3) 
 

### remove module 0 (invariant genes) #####

# how many triads have all 3 homoeologs in the same module

num_triads_same_module_per_TF_family_no_module0 <- triads_modules_3_hom %>%
  group_by(TF_family.x, group_num, bwnetModuleLabels) %>%
  summarise(gene_count = n()) %>%
  filter(gene_count == 3 & bwnetModuleLabels != 0) %>% # only keep rows where the 3 homoeologs are in the same row i.e. 3 hom in same module and exclude module 0
  group_by(TF_family.x) %>% # summarise per TF family
  summarise(num_triads_same_module = sum(gene_count)/3)

head(num_triads_same_module_per_TF_family_no_module0)
print(num_triads_same_module_per_TF_family_no_module0,n=60)
# now merge the num triads same module with the num triads in total per family

head(num_triads_same_module_per_TF_family_no_module0)
head(num_triads_per_family)

merged_data_modules_triads_no_module0 <- merge(num_triads_same_module_per_TF_family_no_module0, num_triads_per_family, all.y=T)
head(merged_data_modules_triads_no_module0)
tail(merged_data_modules_triads_no_module0)

merged_data_modules_triads_no_module0$num_triads_same_module[is.na(merged_data_modules_triads_no_module0$num_triads_same_module)] <- 0  ## replace NA with 0
head(merged_data_modules_triads_no_module0)
tail(merged_data_modules_triads_no_module0)

triad_module_to_plot_no_module0 <- merged_data_modules_triads_no_module0 %>%
  mutate(perc_same_module = num_triads_same_module/num_triads_in_family*100,
         num_triads_diff_module = round(num_triads_in_family - num_triads_same_module),
         isTF = !is.na(TF_family.x))
triad_module_to_plot_no_module0

# do stats
overall_TF_vs_nonTF_contingency <- triad_module_to_plot_no_module0 %>%
  group_by(isTF) %>%
  summarise(count_triad_in_same_module = sum(num_triads_same_module ),
            count_triad_in_diff_module = sum(num_triads_diff_module  )) 
overall_TF_vs_nonTF_contingency

overall_fisher_res <- fisher.test(overall_TF_vs_nonTF_contingency[,2:3])
overall_fisher_res 

# calc average percentage for non-TFs:
perc_same_module_non_TFs_no0 <- triad_module_to_plot_no_module0 %>% # get % for non TFs
  filter(is.na(TF_family.x)) %>%
  select(perc_same_module) %>%
  as.numeric()
perc_same_module_non_TFs_no0

# calc average percentage for TFs:
perc_same_module_TFs_no0 <- triad_module_to_plot_no_module0 %>%
  filter(!is.na(TF_family.x)) %>%
  select(num_triads_same_module, num_triads_in_family) %>%
  summarise("total_TF_in_same_module" = sum(na.omit(num_triads_same_module)),
            "total_TF_in_network" = sum(num_triads_in_family)) %>%
  mutate(perc_same_module <- total_TF_in_same_module/total_TF_in_network)
perc_same_module_TFs_no0

## now to do stats per family
# I have a table which has the number per TF family which is insame module, and number not in same module
head(triad_module_to_plot_no_module0)
tail(triad_module_to_plot_no_module0)


#try for 1 family:
df_1_family <- triad_module_to_plot_no_module0 %>%
  ungroup() %>%
  filter(TF_family.x == "Alfin-like" | is.na(TF_family.x)) %>%
  select(num_triads_same_module, num_triads_diff_module)

fisher.result1 <- fisher.test(df_1_family)
print(fisher.result1$p.value)

# now do for each family vs non-TFs and then correct for multiple testing 

fisher_res <- data.frame(TF_family = as.character(), pvalue = as.numeric() ) #setup empty vector
fisher_res
#i <- "Alfin-like"

for (i in triad_module_to_plot_no_module0$TF_family.x[1:64]) {
  print(i)
  fisher.result <- fisher.test(triad_module_to_plot_no_module0 %>%
                                 ungroup() %>%
                                 filter(TF_family.x == i | is.na(TF_family.x))  %>% 
                                 select(num_triads_same_module, num_triads_diff_module)
  )
  fisher_res[nrow(fisher_res)+1,] <- c(TF_family = i, pvalue = fisher.result$p.value )
}
fisher_res

fisher_res$padj <- p.adjust(fisher_res$pvalue, method = "fdr", n = length(fisher_res$pvalue))
fisher_res
fisher_res$sig_padj <- fisher_res$padj <0.05
fisher_res

write.csv(fisher_res, file="Fisher_exact_test_TF_triads_in_same_module_vs_non_TF_triads_in_same_module.csv",row.names = F)




# plot with fill for family size
plot_TF_same_module_filled_no_module0 <- ggplot(triad_module_to_plot_no_module0[!is.na(triad_module_to_plot_no_module0$TF_family.x),], aes(x=reorder(TF_family.x, perc_same_module), y=perc_same_module)) + 
  geom_segment(aes(xend=TF_family.x, yend=0)) +
  geom_point(aes(fill=num_triads_in_family ), colour="black", alpha =1, pch=21, size =2)+
  scale_fill_gradient(low = "white", high = "red", limits=c(0,120)) +
  coord_flip() + theme_bw() + geom_hline(yintercept = perc_same_module_non_TFs_no0, color="black", linetype="dashed") +
  ylab("Percent of triads with homooelogs in same module") +
  xlab("Transcription Factor Family") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 102))

plot_TF_same_module_filled_no_module0

ggsave(file="Percentage_of_triads_in_same_module_per_TF_family_coloured_num_triads_no_module0_ggsave_all_inc_0perc.pdf", width=8, height=8, units=c("in"), device = cairo_pdf)
plot_TF_same_module_filled_no_module0
dev.off()

# try plotting families with >10 triads (i.e. 30 genes)
triad_module_to_plot_10_triad_families_no_module0 <- triad_module_to_plot_no_module0 %>%
  filter(num_triads_in_family >10 | TF_family.x == "HB-WOX" | TF_family.x == "HB-BELL" | TF_family.x == "MADS-M-type"
  ) # why are HB-WOX and MADS-M-type not in the list? # because they have 0 % in the same module
nrow(triad_module_to_plot_10_triad_families_no_module0)

plot_TF_same_module_filled_10_triads_no_module0 <- ggplot(triad_module_to_plot_10_triad_families_no_module0[!is.na(triad_module_to_plot_10_triad_families_no_module0$TF_family.x),], aes(x=reorder(TF_family.x, perc_same_module), y=perc_same_module)) + 
  geom_segment(aes(xend=TF_family.x, yend=0)) +
  geom_point(aes(fill=num_triads_in_family ), colour="black", alpha =1, pch=21, size =5)+
  scale_fill_gradient(low = "white", high = "red", limits=c(0,133.3), breaks=c(0,33.3, 66.6, 100, 133.3)) +
  coord_flip() + theme_bw() + geom_hline(yintercept = perc_same_module_non_TFs_no0, color="black", linetype="dashed") +
  ylab("Percent of triads with homooelogs in same module") +
  xlab("Transcription Factor Family") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 85))

plot_TF_same_module_filled_10_triads_no_module0

ggsave(file="Percentage_of_triads_in_same_module_per_TF_familes_over_10_triads_coloured_num_triads_no_module0_ggsave_all_inc_0perc.pdf", width=7, height=8, units=c("in"), device = cairo_pdf)
plot_TF_same_module_filled_10_triads_no_module0
dev.off()
