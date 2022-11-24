# Philippa Borrill
# 02-12-2020

# Aim is to find out whether TFs have higher or lower co-expr within their triad than non-TFs

library("tidyverse")
library("fields")
#install.packages("ggrepel")
library("ggrepel")
#install.packages("DescTools")
library("DescTools")

### try using Pearson's correlation ###

# load homoeologs
homologies <- read.csv(file="v1.1_genes_TF_homoeolog_info.csv")
head(homologies)
dim(homologies)


homologies_1_1_1 <- homologies %>%
  filter(type == "1:1:1")
head(homologies_1_1_1)

dim(homologies_1_1_1)


# load expression data

# input tpms:
tpm_input_file <- "choulet_URGI_tpm.tsv" # tpm file available from https://opendata.earlham.ac.uk/wheat/under_license/toronto/Ramirez-Gonzalez_etal_2018-06025-Transcriptome-Landscape/expvip/RefSeq_1.1/ByGene/

# load tpms
tpms <- read.table(file=tpm_input_file,sep = "\t", header=T)
head(tpms)
dim(tpms)


# now I want to only keep genes which are expressed >0.5 tpm in at least 1 condition

# first calculate average per condition (i.e. mean of 2 reps):
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
                          stem_Z30 + stem_Z32 + stem_Z65) / 15)  %>%
  select(gene, grain_Z71:all_tissues)
head(tpms_av)
dim(tpms_av)


# tpms_av$all_tissues is the average across all tissues NB includes genes expressed <0.5 tpm
tpms_av_HC <- tpms_av[!grepl("LC",tpms_av$gene),] #remove LC genes
head(tpms_av_HC)
dim(tpms_av_HC)

# add in TF family info:
tpms_av_HC_TF <- merge(tpms_av_HC,homologies, by.x="gene",by.y="v1.1_ID")
head(tpms_av_HC_TF)
dim(tpms_av_HC_TF)

# select only 1:1:1 triads
triad_tpms <- tpms_av_HC_TF %>%
  filter(type == "1:1:1")
head(triad_tpms)
dim(triad_tpms)

# for one pair
cor(as.numeric(triad_tpms[triad_tpms$gene =="TraesCS1A02G001800",2:16]),as.numeric(triad_tpms[triad_tpms$gene =="TraesCS1A02G001900",2:16]), method="pearson" )

# realise that if one of the homooelogs has 0 expression across all tissues then it messes up the calculation and you get NAs as a result
tpms_id <- triad_tpms %>% mutate(id=1:nrow(.))
over0tpm <- tpms_id %>%
  gather('attribute', 'value', grain_Z71:stem_Z65) %>%
  group_by(id) %>%
  summarize(max_tpm=max(value)) %>%
  right_join(tpms_id, by='id') %>%
  filter(max_tpm >0)
head(over0tpm)
dim(over0tpm)

# keep groups with have 3 homoeolog >0tpm
groups_with_three_hom_over0tpm<- over0tpm %>%
  group_by(group_num)%>%
  summarise(count= length(group_num)) %>%
  filter(count == 3)

groups_with_three_hom_over0tpm 

head(triad_tpms)

# just keep triads which have 3 hom over 0 tpm

triad_tpms_3hom_expr <- triad_tpms[triad_tpms$group_num %in% groups_with_three_hom_over0tpm$group_num,]
head(triad_tpms_3hom_expr)
dim(triad_tpms_3hom_expr)

## now I want to keep triads with >0.5 tpm total expression level

tpms_id <- tpms_av %>% mutate(id=1:nrow(.))
over0.5tpm <- tpms_id %>%
  gather('attribute', 'value', grain_Z71:stem_Z65) %>%
  group_by(id) %>%
  summarize(max_tpm=max(value)) %>%
  right_join(tpms_id, by='id') %>%
  filter(max_tpm >0.5)

head(over0.5tpm)
dim(over0.5tpm)

# now only keep HC genes
HC_expr <- over0.5tpm[!grepl("LC",over0.5tpm$gene),]
head(HC_expr)
dim(HC_expr)

hc_genes_to_use <- data.frame(gene=HC_expr$gene)
head(hc_genes_to_use)
nrow(hc_genes_to_use)


# which triads have at least 1 homoeolog expressed >0.5 tpm?
head(homologies_1_1_1)

# make into long format and add column whether each gene is "hc_genes_to_use"
long_homoeologs <- homologies_1_1_1 %>%
  mutate("gene_to_use" =v1.1_ID %in% hc_genes_to_use$gene ) 
head(long_homoeologs)
dim(long_homoeologs)
nrow(long_homoeologs[long_homoeologs$gene_to_use == "TRUE",])

# does each triad have an expressed gene?
group_expr <- long_homoeologs %>%
  group_by(group_num) %>%
  summarise (n_true = sum(gene_to_use == "TRUE"),
             n_false = sum(gene_to_use == "FALSE"))
head(group_expr)
dim(group_expr)

# select only groups which have at least 1 homoeolog expressed over 0.5 tpm

expressed_groups <- group_expr[group_expr$n_true >0 ,]
head(expressed_groups)
dim(expressed_groups)

# now also select triads where all three homoeologs are >0tpm
triad_tpms_3hom_expr_0.5tpm <- triad_tpms_3hom_expr[triad_tpms_3hom_expr$group_num %in% expressed_groups$group_num,]

dim(triad_tpms_3hom_expr)
dim(triad_tpms_3hom_expr_0.5tpm)

# for each group_num calculate the 3 pairwise correlations and calculate an average. Put into a new table. 

#setup empty dataframe
correlation_per_triad <- data.frame("group_num" = numeric(), "corr1" = numeric(), 
                                   "corr2" = numeric(), "corr3" = numeric())
correlation_per_triad


for (i in (unique(triad_tpms_3hom_expr_0.5tpm$group_num))) { 
  
group_info <- triad_tpms_3hom_expr_0.5tpm %>%
  filter(group_num == i) %>%
  select(grain_Z71:stem_Z65, group_num) 

cor1 <- cor(as.numeric(group_info[1,-16]),as.numeric(group_info[2,-16]), method="pearson" )
cor2 <- cor(as.numeric(group_info[1,-16]),as.numeric(group_info[3,-16]), method="pearson" )
cor3 <- cor(as.numeric(group_info[2,-16]),as.numeric(group_info[3,-16]), method="pearson" )


correlation_per_triad[nrow(correlation_per_triad)+1,] <-  c("group_num" = i, "corr1" = cor1, 
                                 "corr2" = cor2, "corr3" = cor3)
correlation_per_triad


}

head(correlation_per_triad,30)
tail(correlation_per_triad)
dim(correlation_per_triad)



## add in TF family info

head(correlation_per_triad)
head(triad_tpms_3hom_expr_0.5tpm)

group_num_TF_fam <- triad_tpms_3hom_expr_0.5tpm %>%
  select(TF_family.x, group_num, isTF) %>%
  unique()
head(group_num_TF_fam)
dim(group_num_TF_fam)


correl_per_triad_TF_fam <- merge(correlation_per_triad, group_num_TF_fam, by= "group_num")
head(correl_per_triad_TF_fam)
tail(correl_per_triad_TF_fam)
unique(correl_per_triad_TF_fam$TF_family.x)

#make into long data using corr1, corr2 and corr3
correl_per_triad_TF_fam_non_av <- correl_per_triad_TF_fam %>%
  select(group_num, corr1, corr2, corr3, TF_family.x, isTF) %>%
  pivot_longer(cols = corr1:corr3, names_to="corr_source", values_to= "indiv_corr")

head(correl_per_triad_TF_fam)
head(correl_per_triad_TF_fam_non_av)
dim(correl_per_triad_TF_fam_non_av)


# calc median for TF and non-TF
correl_per_triad_TF_fam_non_av %>%
  mutate(FisherZ_correl = FisherZ(indiv_corr)) %>%
  group_by(isTF) %>%
  summarise(median_corr = FisherZInv(median(FisherZ_correl)))


# calc sig diff Man Whitney
head(correl_per_triad_TF_fam_non_av)
wilcox.test(correl_per_triad_TF_fam_non_av$indiv_corr ~ correl_per_triad_TF_fam_non_av$isTF)

correl_per_triad_TF_fam_repl_NA <- correl_per_triad_TF_fam_non_av %>%
  mutate(TF_family = replace_na(as.character(TF_family.x), "nonTF"))
 
head(correl_per_triad_TF_fam_repl_NA)

correl_per_triad_TF_fam_repl_NA$TF_family <- as.factor(correl_per_triad_TF_fam_repl_NA$TF_family)
head(correl_per_triad_TF_fam_repl_NA)
tail(correl_per_triad_TF_fam_repl_NA)

# now write loop for all samples 

i <- "C2H2"

output_man_whitney_df <- data.frame( TF_family=character(), p.value= numeric(), median_corr_TF= numeric(),median_corr_nonTF= numeric())
output_man_whitney_df

unique(correl_per_triad_TF_fam_repl_NA$TF_family)[2:66]

# calc median and p-values per TF family
for (i in unique(correl_per_triad_TF_fam_repl_NA$TF_family)[2:66]) { # only keep 2 to 66 (i.e. excl nonTF)


    print(i)

    man_whitney_correl_loop <- correl_per_triad_TF_fam_repl_NA %>%
      ungroup() %>%
      filter(TF_family== i | TF_family=="nonTF") %>%
      select(indiv_corr, isTF, TF_family)
    
   # head(man_whitney_correl_loop[man_whitney_correl_loop$isTF==TRUE,])
  # tail(man_whitney_correl_loop)
    
    man_whitney <- wilcox.test(man_whitney_correl_loop$indiv_corr ~ man_whitney_correl_loop$isTF)
    #man_whitney$p.value
    
    medians_table <- man_whitney_correl_loop %>%
      group_by(TF_family)  %>%
      summarise(median = FisherZInv(median(FisherZ(indiv_corr))))
    
    median_TF <- medians_table %>%
      filter(TF_family ==i) %>%
      select(median)
    
    median_nonTF <- medians_table %>%
      filter(TF_family =="nonTF") %>%
      select(median)
    
   # as.numeric(median_TF)
  #  as.numeric(median_nonTF)
    
    output_mat <- cbind("TF_family" =i, "p.value" =man_whitney$p.value, "median_corr_TF" = as.numeric(median_TF), "median_corr_nonTF" =   as.numeric(median_nonTF))# add together this information about each sample
    output_df <- as.data.frame(output_mat)
    #head(output_df)
    
    output_man_whitney_df <- rbind(output_man_whitney_df, output_df)
    
  }
  
  head(output_man_whitney_df)
  
  ### Adjust for multiple testing
  
output_man_whitney_df$padj <- p.adjust(output_man_whitney_df$p.value, method = "fdr", n = length(output_man_whitney_df$p.value))

  
head(output_man_whitney_df)

output_man_whitney_df_sig <- output_man_whitney_df %>%
  mutate(is_sig = as.numeric(as.character(p.value)) <0.05,
         is_sig_adj = as.numeric(as.character(padj)) <0.05)

output_man_whitney_df_sig
write.csv(file="man_whitney_pearson_correl_TF_vs_non_TF_triads_per_family_FDR_non_av_per_triad.csv", output_man_whitney_df_sig)


# try plotting boxplot with only TF families >30 genes
TF_families_30_genes<- read.csv(file="TF_families_with_min_30_genes_in_triads.csv"))
TF_families_30_genes
TF_families_30_genes <- rbind(TF_families_30_genes, c(NA, "1:1:1", "lots"))
TF_families_30_genes

data_TF_families_30_genes <- correl_per_triad_TF_fam_non_av%>%
  filter(TF_family.x %in% TF_families_30_genes$TF_family.x )
head(data_TF_families_30_genes)
tail(data_TF_families_30_genes)
dim(data_TF_families_30_genes)
dim(correl_per_triad_TF_fam_non_av)



# try to colour code boxplot by significance
head(data_TF_families_30_genes)
head(output_man_whitney_df_sig)

data_TF_families_30_genes_sig <- merge(data_TF_families_30_genes, output_man_whitney_df_sig, by.x="TF_family.x", by.y = "TF_family", all.x=T)

head(data_TF_families_30_genes_sig)
tail(data_TF_families_30_genes_sig)

#calc median to plot line on graph
median_correl_non_TF_calc <- correl_per_triad_TF_fam_non_av %>%
  filter(isTF==FALSE) %>%
  summarise(median = median(indiv_corr))
median_correl_non_TF <- as.numeric(median_correl_non_TF_calc)

p_boxplot_min_30_genes_sig <- ggplot(data=data_TF_families_30_genes_sig, aes(x=reorder(TF_family.x, indiv_corr), y=indiv_corr, fill=is_sig_adj)) + 
  geom_boxplot() +
  coord_flip() + theme_bw()   + geom_hline(yintercept = median_correl_non_TF, color="black", linetype="dashed") +
  ylab("Pearson's correlation") +
  xlab("Transcription Factor Family") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
  scale_fill_manual(values=c("white","grey")) 

p_boxplot_min_30_genes_sig
ggsave(file="boxplot_correlation_between_homoeologs_sum_0.5tpm_min_30_genes_sig_col_FDR_non_av_per_triad.pdf", height=8, width=8,units=c("in"), device = cairo_pdf)


## now make boxplot for all families with colour

# try plotting boxplot of all data

correl_per_triad_TF_fam_sig <- merge(correl_per_triad_TF_fam_non_av, output_man_whitney_df_sig, by.x="TF_family.x", by.y = "TF_family", all.x=T)

head(correl_per_triad_TF_fam_sig)
tail(correl_per_triad_TF_fam_sig)

p_boxplot <- ggplot(data=correl_per_triad_TF_fam_sig, aes(x=reorder(TF_family.x, indiv_corr), y=indiv_corr, fill=is_sig_adj)) + 
  geom_boxplot() +
  coord_flip() + theme_bw()   + geom_hline(yintercept = median_correl_non_TF, color="black", linetype="dashed") +
  ylab("Pearson's correlation") +
  xlab("Transcription Factor Family") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.02))+
  scale_fill_manual(values=c("white","grey"))

p_boxplot
ggsave("boxplot_correlation_between_homoeologs_sum_0.5tpm_sig_coloured_FDR_non_av_per_triad.pdf", height=8, width=8,units=c("in"), device = cairo_pdf)


