# Aim is to find out whether TFs in 1:1:1 triads have similar or different expression patterns. Are these expression patterns more extreme than for other genes in 1:1: triads? I willclassify triads as balanced, dominant or suppressed like Ricardo did

# Philippa Borrill
# 14-03-2019

# Used https://github.com/Uauy-Lab/WheatHomoeologExpression/blob/master/02.%20Calculate%20triad%20category.ipynb#Definition-of-homoeolog-expression-bias-categories as a model

library("tidyverse")
library("fields")
library("NMF")


## do this analysis with Choulet Chinese Spring data 
# input tpms:
tpm_input_file <- "choulet_URGI_tpm.tsv" # tpm file available from https://opendata.earlham.ac.uk/wheat/under_license/toronto/Ramirez-Gonzalez_etal_2018-06025-Transcriptome-Landscape/expvip/RefSeq_1.1/ByTranscript/
tpm_input_file
# load tpms
tpms <- read.table(file=tpm_input_file,sep = "\t", header=T)
head(tpms)
dim(tpms)


# load homoeologs
homologies <- read.csv(file="v1.1_genes_TF_homoeolog_info.csv")
head(homologies)

# remove genes not in a homoeolog group
homologies <- homologies %>%
  filter(group_num != "NA")
head(homologies)
dim(homologies)

# now only keep 1:1:1 homologies
homologies_1_1_1 <- homologies %>%
  filter(type=="1:1:1")
head(homologies_1_1_1)
tail(homologies_1_1_1)
dim(homologies_1_1_1)

# now I want to only keep genes which are expressed >0.5 tpm in at least 1 condition
dim(tpms)

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

# select only groups which have at least 1 homoeolog expressed

expressed_groups <- group_expr[group_expr$n_true >0 ,]
head(expressed_groups)
dim(expressed_groups)

# now I want to combine the expression data with this list of expressed groups
head(expressed_groups)

head(long_homoeologs)
dim(long_homoeologs)

long_homoeologs_to_use <- long_homoeologs[long_homoeologs$group_num %in% expressed_groups$group_num,]
head(long_homoeologs_to_use)
dim(long_homoeologs_to_use)

# do a test to make a matrix for tpm expression for one sample:
head(tpms_av)
# make the geneIDs into rownames and remove 1st column
row.names(tpms_av) <- tpms_av$gene
tpms_av <- tpms_av[,-1]
head(tpms_av)

tpm_long_homoeologs_to_use <- merge(long_homoeologs_to_use, tpms_av, by.x="v1.1_ID", by.y=0) # add the tpm values for the long_homoeologs_to_use
head(tpm_long_homoeologs_to_use)
tail(tpm_long_homoeologs_to_use)
dim(tpm_long_homoeologs_to_use)

# now rename columns and add colum saying which homoeolog is which

tpm_long_homoeologs_to_use <- tpm_long_homoeologs_to_use %>%
  rename(gene = v1.1_ID, 
         group_id = group_num, 
         TF_family = TF_family.x) %>% 
  mutate(homoeolog = case_when(grepl("A",gene) == TRUE ~ "A", 
                               grepl("B",gene) == TRUE ~ "B", 
                               grepl("D",gene) == TRUE ~ "D"))

head(tpm_long_homoeologs_to_use)
tail(tpm_long_homoeologs_to_use)
dim(tpm_long_homoeologs_to_use)
unique(tpm_long_homoeologs_to_use$homoeolog)

#make new dataframe which matches group_id to TF_family
head(tpm_long_homoeologs_to_use)

group_id_TF <- tpm_long_homoeologs_to_use %>%
  select(TF_family,group_id) %>%
  unique()
head(group_id_TF)
dim(group_id_TF)

# NB a couple of group_id are in two different TF families- need to remove these!!
counts_group_id_TF <- group_id_TF %>%
  group_by(group_id) %>%
  summarise(count = n()) %>%
  left_join(group_id_TF, by='group_id') %>%
  filter(count >1)
counts_group_id_TF
# so I now have a table which shows which triads (group_id) have different TF families assigned
# I need to go into the tpm_long_homoeologs_to_use and correct these...
# I'm not sure there is any way to do this efficiently, there are 54 triads affected 
  

# make a subset of tpm_long_homoeologs_to_use which I will edit: 

tpm_long_homoeologs_to_use_subset <- tpm_long_homoeologs_to_use %>%
  filter(group_id %in% counts_group_id_TF$group_id ) %>%
  arrange(group_id)

dim(tpm_long_homoeologs_to_use_subset)
head(tpm_long_homoeologs_to_use_subset,20)
nrow(tpm_long_homoeologs_to_use_subset)/3

# decided that if all three homoeologs were annotated as a TF but the families disagreed (e.g. MYB/MYB-related or HB-other/DTT)
# I will assign the homoeolog which is different to match with the other 2 homoeologs
# If 2 homoeologs are "NA" (i.e. not TFs) but the other one is a TF I will assign it to be NA
head(tpm_long_homoeologs_to_use_subset,20)
subset_TF_family_group_id <- tpm_long_homoeologs_to_use_subset %>%
  select(gene,TF_family, group_id) %>% 
  group_by(group_id,TF_family) %>%
 summarise(count= n()) %>% # count the number of genes in each group/TF_family combination
 filter(count == 2) 

head(tpm_long_homoeologs_to_use_subset,20)

tpm_long_homoeologs_to_use_subset_cleaned <- merge(tpm_long_homoeologs_to_use_subset, subset_TF_family_group_id, by="group_id")
head(tpm_long_homoeologs_to_use_subset_cleaned)
tpm_long_homoeologs_to_use_subset_cleaned <- tpm_long_homoeologs_to_use_subset_cleaned %>%
  select(-TF_family.x, -count) %>% # get rid of original TF family column 
  rename(TF_family = TF_family.y) %>% # rename new TF_family column to original column (in this one all three homoeologs agree)
  mutate(isTF2 = !is.na(TF_family)) %>%
  select(-isTF) %>%
  rename(isTF = isTF2)
head(tpm_long_homoeologs_to_use_subset_cleaned)
tpm_long_homoeologs_to_use_subset_cleaned %>%
  group_by(isTF) %>%
  summarise(count=n()/3)

 
# make final list of triads to use incorporated the "cleaned" triads which didn't agree about the TF family:
tpm_long_homoeologs_to_use_final <- tpm_long_homoeologs_to_use %>%
filter(!group_id %in% counts_group_id_TF$group_id ) %>%
  bind_rows(tpm_long_homoeologs_to_use_subset_cleaned) %>%
  arrange(group_id)
head(tpm_long_homoeologs_to_use_final)
tail(tpm_long_homoeologs_to_use_final)
dim(tpm_long_homoeologs_to_use_final)
nrow(tpm_long_homoeologs_to_use_final)/3

## ok so I have the final datatable to use: tpm_long_homoeologs_to_use_final

### now calculate relative homoeolog expression for each triad, do this for each sample AND keep original values:
colnames(tpm_long_homoeologs_to_use_final)

list_of_samples <- c(colnames(tpm_long_homoeologs_to_use_final[,7:(ncol(tpm_long_homoeologs_to_use_final)-1)]))
list_of_samples

# make output dataframe:
head(output_df)
output_df_all_samples <- data.frame(A_tpm= numeric(), B_tpm= numeric(), D_tpm= numeric(), A= numeric(), B= numeric(), D= numeric(), 
                                    group_id = numeric(), sample=character())     
output_df_all_samples

for(sample in list_of_samples) {
  print(sample)
  # now select just 1 sample and calculate the relative expression of A, B, D for each triad:
  test_sample <-
    tpm_long_homoeologs_to_use_final %>%
    select(group_id, homoeolog, sample) %>%
    spread(homoeolog,  sample)
  head(test_sample)
  
  # now calculate relative ABD
  test_sample$total <- test_sample$A + test_sample$B + test_sample$D
  head(test_sample)
  
  test_sample$A_rel <- test_sample$A/test_sample$total
  test_sample$B_rel <- test_sample$B/test_sample$total
  test_sample$D_rel <- test_sample$D/test_sample$total
  
  head(test_sample)
  dim(test_sample)
  # only keep triads with a sum >0.5 tpm
  test_sample <- test_sample[test_sample$total >0.5,] 
  
  # make matrix to have group_id as rownames
  test_mat <- as.matrix(test_sample[,c("A_rel","B_rel","D_rel")])
  #head(test_mat)
  rownames(test_mat)<- test_sample$group_id
  colnames(test_mat) <- c("A","B","D")
  #head(test_mat)
 
  output_mat <- cbind("A_tpm" =test_sample$A, "B_tpm" =test_sample$B, "D_tpm"=test_sample$D, test_mat)# add together this information for each triad
  output_df <- as.data.frame(output_mat)
  head(output_df)
  output_df$group_id <- rownames(output_df) # make group_id into a column
  output_df$sample <- sample # add which sample this is as a column
  
  head(output_df)
  
  output_df_all_samples <- rbind(output_df_all_samples, output_df) # puts all samples into a big table
  
}

head(output_df_all_samples)
tail(output_df_all_samples)
dim(output_df_all_samples)

## now I want to add in info about which ones are TFs
head(tpm_long_homoeologs_to_use_final)
TF_info <- tpm_long_homoeologs_to_use_final %>%
  select(TF_family, group_id, isTF) %>%
  unique()
head(TF_info)
TF_info$group_id <- as.character(TF_info$group_id) # convert to character so I can do left-join


output_df_all_samples_TF_info <- left_join(output_df_all_samples,TF_info, by="group_id")

head(output_df_all_samples_TF_info)
tail(output_df_all_samples_TF_info)
dim(output_df_all_samples_TF_info)


#  calculating standard deviations:

head(output_df_all_samples_TF_info)

output_df_all_samples_TF_info2 <- output_df_all_samples_TF_info
output_df_all_samples_TF_info2$A <- as.numeric(as.character(output_df_all_samples_TF_info2$A)) # have to convert factor to character and then to numeric
output_df_all_samples_TF_info2$B <- as.numeric(as.character(output_df_all_samples_TF_info2$B))
output_df_all_samples_TF_info2$D <- as.numeric(as.character(output_df_all_samples_TF_info2$D))

output_df_all_samples_TF_info2$A_tpm <- as.numeric(as.character(output_df_all_samples_TF_info2$A_tpm)) # have to convert factor to character and then to numeric
output_df_all_samples_TF_info2$B_tpm <- as.numeric(as.character(output_df_all_samples_TF_info2$B_tpm))
output_df_all_samples_TF_info2$D_tpm <- as.numeric(as.character(output_df_all_samples_TF_info2$D_tpm))


head(output_df_all_samples_TF_info)
head(output_df_all_samples_TF_info2)

output_df_all_samples_TF_info_sd <- output_df_all_samples_TF_info2 %>%
  rowwise() %>%
  mutate(stdev = sd(c(A,B,D)),
         stdev_tpm = sd(c(A_tpm,B_tpm,D_tpm)))
head(output_df_all_samples_TF_info_sd)
tail(output_df_all_samples_TF_info_sd)         

dim(output_df_all_samples_TF_info_sd)

##  compare stdev for each tissue for each triad for TFs vs non-TFs i.e. exclude all_tissues
# for man-whitney (aka wilcox in R)  using sd from normalised tpm per triad
man_whitney_stdev_14tissues <- output_df_all_samples_TF_info_sd %>%
  ungroup() %>%
  filter(sample!="all_tissues") %>% 
  select(stdev, isTF, group_id, sample)
head(man_whitney_stdev_14tissues)
dim(man_whitney_stdev_14tissues)

# do i have the expected numbers of triads and tissues
summary_table <- man_whitney_stdev_14tissues %>%
  group_by(sample, isTF) %>%
  summarise( n = n()) 


# now write loop for all samples 

output_man_whitney_df <- data.frame( sample=character(), p.value= numeric())
output_man_whitney_df


for(i in unique(output_df_all_samples_TF_info_sd$sample)) {
  print(i)

  man_whitney_stdev_loop <- output_df_all_samples_TF_info_sd %>%
    ungroup() %>%
    filter(sample== i) %>%
    select(stdev, isTF)

man_whitney <- wilcox.test(man_whitney_stdev_loop$stdev ~ man_whitney_stdev_loop$isTF)
#man_whitney$p.value

output_mat <- cbind("sample" =i, "p.value" =man_whitney$p.value)# add together this information about each sample
output_df <- as.data.frame(output_mat)
#head(output_df)

output_man_whitney_df <- rbind(output_man_whitney_df, output_df)

}

head(output_man_whitney_df)

write.csv(file="man_whitney_stdev_TF_vs_non_TF_triads.csv", output_man_whitney_df)

# calc mean and median for each tissue
mean_and_median_per_tissue <- output_df_all_samples_TF_info_sd %>%
  group_by(sample, isTF) %>%
  summarise(mean = mean(stdev),
            median=median(stdev))

write.csv(file="mean_and_median_stdev_relative_tpm_per_tissue_TF_vs_non_TF_triads.csv", mean_and_median_per_tissue,row.names = F)

