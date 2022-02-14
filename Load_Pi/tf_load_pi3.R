setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(dgof)
library(ggpubr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(car)
require(readr)
library(reshape2)
library(lme4)
library(cowplot)

### 1. first get list of all TF and non-TF sites in each site class ###

## filtered and annotated SNP set obtained from Catherine Evans, not excluding positively selected sites

combined_table_MAF0.01 <- read_csv("../data_files/coding_variant_tf_allele_frequency_table_MAF0.01_2022-01-28.csv")
combined_table_MAF0.01 <- combined_table_MAF0.01 %>%
  filter(Expressed == TRUE & Synonly == TRUE)
colnames(combined_table_MAF0.01)[2] <- "site"
combined_table_MAF0.01$site <- gsub(":","_",combined_table_MAF0.01$site)
combined_table_MAF0.01$site <- paste("chr",combined_table_MAF0.01$site,sep="")
combined_table_MAF0.01 <- combined_table_MAF0.01[!duplicated(combined_table_MAF0.01),] ## remove any duplicated lines

## contains per-site pi for all polymorphic sites (from vcftools)
out.sites.pi <- read.table(file="../data_files/filtered_sites.sites.pi",row.names = NULL,header=T)
out.sites.pi$site <- paste(out.sites.pi$CHROM,out.sites.pi$POS,sep="_") ## site id

## obtain site pi for missense, stop codons and synonymous sites in transcription factors (tf)
out.sites.pi_sg_tf <- out.sites.pi[out.sites.pi$site %in% combined_table_MAF0.01[combined_table_MAF0.01$interaction %in% "stop_gained.na" & combined_table_MAF0.01$isTF == TRUE,]$site,] ## Stop gained
out.sites.pi_ms_tf <- out.sites.pi[out.sites.pi$site %in% combined_table_MAF0.01[combined_table_MAF0.01$interaction %in% "missense_variant.deleterious" & combined_table_MAF0.01$isTF == TRUE,]$site,] ## only keep missense sites that are highly Deleterious (sift score <= 0.05 for analyses)
out.sites.pi_to_tf <- out.sites.pi[out.sites.pi$site %in% combined_table_MAF0.01[combined_table_MAF0.01$interaction %in% "missense_variant.tolerated" & combined_table_MAF0.01$isTF == TRUE,]$site,] ## only keep missense sites that are tolerated (sift score > 0.05 for analyses)
out.sites.pi_s_tf <- out.sites.pi[out.sites.pi$site %in% combined_table_MAF0.01[combined_table_MAF0.01$interaction %in% "synonymous_variant.na" & combined_table_MAF0.01$isTF == TRUE,]$site,] ## synonymous

## obtain site pi for missense, stop codons and synonymous sites in background - all genes except tf
out.sites.pi_sg <- out.sites.pi[out.sites.pi$site %in% combined_table_MAF0.01[combined_table_MAF0.01$interaction %in% "stop_gained.na" & combined_table_MAF0.01$isTF == FALSE,]$site,] ## Stop gained
out.sites.pi_ms <- out.sites.pi[out.sites.pi$site %in% combined_table_MAF0.01[combined_table_MAF0.01$interaction %in% "missense_variant.deleterious" & combined_table_MAF0.01$isTF == FALSE,]$site,] ## only keep missense sites that are highly Deleterious (sift score <= 0.05 for analyses)
out.sites.pi_to <- out.sites.pi[out.sites.pi$site %in% combined_table_MAF0.01[combined_table_MAF0.01$interaction %in% "missense_variant.tolerated" & combined_table_MAF0.01$isTF == FALSE,]$site,] ## only keep missense sites that are tolerated (sift score > 0.05 for analyses)
out.sites.pi_s <- out.sites.pi[out.sites.pi$site %in% combined_table_MAF0.01[combined_table_MAF0.01$interaction %in% "synonymous_variant.na" & combined_table_MAF0.01$isTF == FALSE,]$site,] ## synonymous

tfsites <- rbind(out.sites.pi_ms_tf[1:2],out.sites.pi_sg_tf[1:2],out.sites.pi_to_tf[1:2],out.sites.pi_s_tf[1:2])
backgroundsites <- rbind(out.sites.pi_ms[1:2],out.sites.pi_sg[1:2],out.sites.pi_to[1:2],out.sites.pi_s[1:2])

write.table(tfsites,file="tf_sites.txt",row.names = F, col.names = F, quote = F,sep="\t")
write.table(backgroundsites,file="background_sites.txt",row.names = F, col.names = F, quote = F,sep="\t")

## get total number of sites
rbind(cbind(nrow(out.sites.pi_sg_tf),nrow(out.sites.pi_ms_tf),nrow(out.sites.pi_to_tf),nrow(out.sites.pi_s_tf)),cbind(nrow(out.sites.pi_sg),nrow(out.sites.pi_ms),nrow(out.sites.pi_to),nrow(out.sites.pi_s)))

### 2. next, find regions under positive selection in He et al (2019). Different between TF and non-TF, so exclude such regions from comparisons. ###

sites_sel <- read.csv(file="../data_files/sites_sel.csv")

chisq.test(sites_sel[sites_sel$Type %in% "Sweep sites",3:4])
chisq.test(sites_sel[sites_sel$Type %in% "Introgression sites",3:4])
chisq.test(sites_sel[sites_sel$Type %in% "Environmental adaptation sites",3:4])
chisq.test(sites_sel[sites_sel$Type %in% "Improvement selection sites",3:4])

sites_sel$Proportion <- sites_sel$Yes/(sites_sel$Yes + sites_sel$No)

ggplot(data=sites_sel,aes(x=Class,y=Proportion,fill=Class))+
  facet_wrap(~Type,scales = "free_y",ncol=4)+
  geom_col(position = "dodge2")

### 3. estimate nucleotide diversity (PI) for TF and non TFs  ###

## filtered and annotated SNP set obtained from Catherine Evans,  excluding positively selected sites

combined_table_MAF0.01 <- read_csv("../data_files/coding_variant_tf_allele_frequency_table_MAF0.01_2022-01-28.csv")
combined_table_MAF0.01 <- combined_table_MAF0.01 %>%
  filter(Expressed == TRUE & Sweep == FALSE & Synonly == TRUE)
colnames(combined_table_MAF0.01)[2] <- "site"
combined_table_MAF0.01$site <- gsub(":","_",combined_table_MAF0.01$site)
combined_table_MAF0.01$site <- paste("chr",combined_table_MAF0.01$site,sep="")
combined_table_MAF0.01 <- combined_table_MAF0.01[!duplicated(combined_table_MAF0.01),] ## remove any duplicated lines

## contains per-site pi for all polymorphic sites (from vcftools)
out.sites.pi <- read.table(file="../data_files/filtered_sites.sites.pi",row.names = NULL,header=T)
out.sites.pi$site <- paste(out.sites.pi$CHROM,out.sites.pi$POS,sep="_") ## site id

## obtain site pi for missense, stop codons and synonymous sites in transcription factors (tf)
out.sites.pi_sg_tf <- out.sites.pi[out.sites.pi$site %in% combined_table_MAF0.01[combined_table_MAF0.01$interaction %in% "stop_gained.na" & combined_table_MAF0.01$isTF == TRUE,]$site,] ## Stop gained
out.sites.pi_ms_tf <- out.sites.pi[out.sites.pi$site %in% combined_table_MAF0.01[combined_table_MAF0.01$interaction %in% "missense_variant.deleterious" & combined_table_MAF0.01$isTF == TRUE,]$site,] ## only keep missense sites that are highly Deleterious (sift score <= 0.05 for analyses)
out.sites.pi_to_tf <- out.sites.pi[out.sites.pi$site %in% combined_table_MAF0.01[combined_table_MAF0.01$interaction %in% "missense_variant.tolerated" & combined_table_MAF0.01$isTF == TRUE,]$site,] ## only keep missense sites that are tolerated (sift score > 0.05 for analyses)
out.sites.pi_s_tf <- out.sites.pi[out.sites.pi$site %in% combined_table_MAF0.01[combined_table_MAF0.01$interaction %in% "synonymous_variant.na" & combined_table_MAF0.01$isTF == TRUE,]$site,] ## synonymous

## obtain site pi for missense, stop codons and synonymous sites in background - all genes except tf
out.sites.pi_sg <- out.sites.pi[out.sites.pi$site %in% combined_table_MAF0.01[combined_table_MAF0.01$interaction %in% "stop_gained.na" & combined_table_MAF0.01$isTF == FALSE,]$site,] ## Stop gained
out.sites.pi_ms <- out.sites.pi[out.sites.pi$site %in% combined_table_MAF0.01[combined_table_MAF0.01$interaction %in% "missense_variant.deleterious" & combined_table_MAF0.01$isTF == FALSE,]$site,] ## only keep missense sites that are highly Deleterious (sift score <= 0.05 for analyses)
out.sites.pi_to <- out.sites.pi[out.sites.pi$site %in% combined_table_MAF0.01[combined_table_MAF0.01$interaction %in% "missense_variant.tolerated" & combined_table_MAF0.01$isTF == FALSE,]$site,] ## only keep missense sites that are tolerated (sift score > 0.05 for analyses)
out.sites.pi_s <- out.sites.pi[out.sites.pi$site %in% combined_table_MAF0.01[combined_table_MAF0.01$interaction %in% "synonymous_variant.na" & combined_table_MAF0.01$isTF == FALSE,]$site,] ## synonymous

## combine  Deleterious and synonymous mutations for TF
pi_tf <- as.data.frame(rbind(cbind(out.sites.pi_sg_tf$PI,"Stop gained"),cbind(out.sites.pi_ms_tf$PI,"Deleterious missense"),cbind(out.sites.pi_to_tf$PI,"Tolerated missense"),cbind(out.sites.pi_s_tf$PI,"Synonymous"))) 
pi_tf$V1 <- as.numeric(pi_tf$V1)
colnames(pi_tf) <- c("pi","effect")

## combine  Deleterious and synonymous mutations for non-TF
pi_all <- as.data.frame(rbind(cbind(out.sites.pi_sg$PI,"Stop gained"),cbind(out.sites.pi_ms$PI,"Deleterious missense"),cbind(out.sites.pi_to$PI,"Tolerated missense"),cbind(out.sites.pi_s$PI,"Synonymous")))  
pi_all$V1 <- as.numeric(pi_all$V1)

## compare and plot pi for tf and background genes
colnames(pi_all) <- c("pi","effect")
pi_all$class <- "Background"
pi_tf$class <- "TF"
pi_merge <- rbind(pi_tf,pi_all) ## combine tf and background genes
pi_merge$class <- as.factor(pi_merge$class)

## stats to check if significant difference in pi between TF and background, no significant differences in either missense or synonymous sites

wilcox.test(pi_tf[pi_tf$effect %in% "Stop gained",]$pi,pi_all[pi_all$effect %in% "Stop gained",]$pi)
wilcox.test(pi_tf[pi_tf$effect %in% "Deleterious missense",]$pi,pi_all[pi_all$effect %in% "Deleterious missense",]$pi)
wilcox.test(pi_tf[pi_tf$effect %in% "Tolerated missense",]$pi,pi_all[pi_all$effect %in% "Tolerated missense",]$pi)
wilcox.test(pi_tf[pi_tf$effect %in% "Synonymous",]$pi,pi_all[pi_all$effect %in% "Synonymous",]$pi)

pistatvalues <- as.data.frame(rbind(c("Stop gained","W = 237.5,\np-value = 0.11"),c("Deleterious missense","W = 142395,\np-value = 0.99"),c("Tolerated missense","W = 1550269,\np-value = 0.32"),c("Synonymous","W = 1174306,\np-value = 0.07")))
colnames(pistatvalues) <- c("effect","value")
pistatvalues$class <- "Background"

## get total number of sites
rbind(cbind(nrow(out.sites.pi_sg_tf),nrow(out.sites.pi_ms_tf),nrow(out.sites.pi_to_tf),nrow(out.sites.pi_s_tf)),cbind(nrow(out.sites.pi_sg),nrow(out.sites.pi_ms),nrow(out.sites.pi_to),nrow(out.sites.pi_s)))

## now plot the histograms

misdel_tf_hist <- hist(pi_tf[pi_tf$effect %in% "Deleterious missense",]$pi, breaks = hist(pi_all[pi_all$effect %in% "Synonymous",]$pi)$breaks)$counts/sum(hist(pi_tf[pi_tf$effect %in% "Deleterious missense",]$pi)$counts)
mistol_tf_hist <- hist(pi_tf[pi_tf$effect %in% "Tolerated missense",]$pi, breaks = hist(pi_all[pi_all$effect %in% "Synonymous",]$pi)$breaks)$counts/sum(hist(pi_tf[pi_tf$effect %in% "Tolerated missense",]$pi)$counts)
sg_tf_hist <- hist(pi_tf[pi_tf$effect %in% "Stop gained",]$pi, breaks = hist(pi_all[pi_all$effect %in% "Synonymous",]$pi)$breaks)$counts/sum(hist(pi_tf[pi_tf$effect %in% "Stop gained",]$pi)$counts)
syn_tf_hist <- hist(pi_tf[pi_tf$effect %in% "Synonymous",]$pi, breaks = hist(pi_all[pi_all$effect %in% "Synonymous",]$pi)$breaks)$counts/sum(hist(pi_tf[pi_tf$effect %in% "Synonymous",]$pi)$counts)

misdel_tf_hist <- as.data.frame(misdel_tf_hist)
colnames(misdel_tf_hist) <- "Proportion"
misdel_tf_hist$effect <- "Deleterious missense"
misdel_tf_hist$class <- "TF"
misdel_tf_hist$Pi <- hist(pi_all[pi_all$effect %in% "Synonymous",]$pi)$breaks[0:(length(hist(pi_all[pi_all$effect %in% "Synonymous",]$pi)$breaks)-1)]

mistol_tf_hist <- as.data.frame(mistol_tf_hist)
colnames(mistol_tf_hist) <- "Proportion"
mistol_tf_hist$effect <- "Tolerated missense"
mistol_tf_hist$class <- "TF"
mistol_tf_hist$Pi <- hist(pi_all[pi_all$effect %in% "Synonymous",]$pi)$breaks[0:(length(hist(pi_all[pi_all$effect %in% "Synonymous",]$pi)$breaks)-1)]

sg_tf_hist <- as.data.frame(sg_tf_hist)
colnames(sg_tf_hist) <- "Proportion"
sg_tf_hist$effect <- "Stop gained"
sg_tf_hist$class <- "TF"
sg_tf_hist$Pi <- hist(pi_all[pi_all$effect %in% "Synonymous",]$pi)$breaks[0:(length(hist(pi_all[pi_all$effect %in% "Synonymous",]$pi)$breaks)-1)]

syn_tf_hist <- as.data.frame(syn_tf_hist)
colnames(syn_tf_hist) <- "Proportion"
syn_tf_hist$effect <- "Synonymous"
syn_tf_hist$class <- "TF"
syn_tf_hist$Pi <- hist(pi_all[pi_all$effect %in% "Synonymous",]$pi)$breaks[0:(length(hist(pi_all[pi_all$effect %in% "Synonymous",]$pi)$breaks)-1)]

misdel_all_hist <- hist(pi_all[pi_all$effect %in% "Deleterious missense",]$pi, breaks = hist(pi_all[pi_all$effect %in% "Synonymous",]$pi)$breaks)$counts/sum(hist(pi_all[pi_all$effect %in% "Deleterious missense",]$pi)$counts)
mistol_all_hist <- hist(pi_all[pi_all$effect %in% "Tolerated missense",]$pi, breaks = hist(pi_all[pi_all$effect %in% "Synonymous",]$pi)$breaks)$counts/sum(hist(pi_all[pi_all$effect %in% "Tolerated missense",]$pi)$counts)
sg_all_hist <- hist(pi_all[pi_all$effect %in% "Stop gained",]$pi, breaks = hist(pi_all[pi_all$effect %in% "Synonymous",]$pi)$breaks)$counts/sum(hist(pi_all[pi_all$effect %in% "Stop gained",]$pi)$counts)
syn_all_hist <- hist(pi_all[pi_all$effect %in% "Synonymous",]$pi, include.lowest = T, breaks = hist(pi_all[pi_all$effect %in% "Synonymous",]$pi)$breaks)$counts/sum(hist(pi_all[pi_all$effect %in% "Synonymous",]$pi)$counts)

misdel_all_hist <- as.data.frame(misdel_all_hist)
colnames(misdel_all_hist) <- "Proportion"
misdel_all_hist$effect <- "Deleterious missense"
misdel_all_hist$class <- "Background"
misdel_all_hist$Pi <- hist(pi_all[pi_all$effect %in% "Synonymous",]$pi)$breaks[0:(length(hist(pi_all[pi_all$effect %in% "Synonymous",]$pi)$breaks)-1)]

mistol_all_hist <- as.data.frame(mistol_all_hist)
colnames(mistol_all_hist) <- "Proportion"
mistol_all_hist$effect <- "Tolerated missense"
mistol_all_hist$class <- "Background"
mistol_all_hist$Pi <- hist(pi_all[pi_all$effect %in% "Synonymous",]$pi)$breaks[0:(length(hist(pi_all[pi_all$effect %in% "Synonymous",]$pi)$breaks)-1)]

sg_all_hist <- as.data.frame(sg_all_hist)
colnames(sg_all_hist) <- "Proportion"
sg_all_hist$effect <- "Stop gained"
sg_all_hist$class <- "Background"
sg_all_hist$Pi <- hist(pi_all[pi_all$effect %in% "Synonymous",]$pi)$breaks[0:(length(hist(pi_all[pi_all$effect %in% "Synonymous",]$pi)$breaks)-1)]

syn_all_hist <- as.data.frame(syn_all_hist)
colnames(syn_all_hist) <- "Proportion"
syn_all_hist$effect <- "Synonymous"
syn_all_hist$class <- "Background"
syn_all_hist$Pi <- hist(pi_all[pi_all$effect %in% "Synonymous",]$pi)$breaks[0:(length(hist(pi_all[pi_all$effect %in% "Synonymous",]$pi)$breaks)-1)]

all_hist <- rbind(misdel_tf_hist,mistol_tf_hist,sg_tf_hist,syn_tf_hist, misdel_all_hist, mistol_all_hist, sg_all_hist, syn_all_hist)
all_hist$class <- gsub("Background","non-TF",all_hist$class)
pistatvalues$class <- gsub("Background","non-TF",pistatvalues$class)

cairo_pdf(file="pi.pdf",width=8.5,height=4, family = "Arial")
ggplot(data=all_hist,aes(x=Pi,y=Proportion,fill=class))+
  geom_col(position = "dodge2")+
  facet_grid(~factor(effect,levels= c("Stop gained","Deleterious missense","Tolerated missense","Synonymous")))+
  geom_text(data=pistatvalues,aes(x=0.3,y=0.7,label=value))+
  theme_bw()+
  xlab("Nucleotide diversity (π)")+
  ylab("Proportion of sites")+
  theme(panel.grid = element_blank(),legend.title = element_blank())
dev.off()

## get list of Deleterious sites for Deleterious load analysis (site positions will be used to subset vcf file)
write.table(out.sites.pi_ms_tf[1:2],file="missense_Deleterious_tf.txt",row.names = F, col.names = F, quote = F,sep="\t")
write.table(out.sites.pi_sg_tf[1:2],file="missense_stop_tf.txt",row.names = F, col.names = F, quote = F,sep="\t")
write.table(out.sites.pi_to_tf[1:2],file="missense_tolerated_tf.txt",row.names = F, col.names = F, quote = F,sep="\t")
write.table(out.sites.pi_s_tf[1:2],file="missense_syn_tf.txt",row.names = F, col.names = F, quote = F,sep="\t")

write.table(out.sites.pi_ms[1:2],file="missense_Deleterious_background.txt",row.names = F, col.names = F, quote = F,sep="\t")
write.table(out.sites.pi_sg[1:2],file="missense_stop_background.txt",row.names = F, col.names = F, quote = F,sep="\t")
write.table(out.sites.pi_to[1:2],file="missense_tolerated_background.txt",row.names = F, col.names = F, quote = F,sep="\t")
write.table(out.sites.pi_s[1:2],file="missense_syn_background.txt",row.names = F, col.names = F, quote = F,sep="\t")


### 4. Association between pi and allele frequency ### 

missense_stop_background_freq <- read.table(file="missense_stop_background_freq.frq",row.names = NULL)
missense_stop_background_freq$X.ALLELE.FREQ. <- as.numeric(gsub(".:","",missense_stop_background_freq$X.ALLELE.FREQ.))

missense_stop_tf_freq <- read.table(file="missense_stop_tf_freq.frq",row.names = NULL)
missense_stop_tf_freq$X.ALLELE.FREQ. <- as.numeric(gsub(".:","",missense_stop_tf_freq$X.ALLELE.FREQ.))

missense_del_background_freq <- read.table(file="missense_del_background_freq.frq",row.names = NULL)
missense_del_background_freq$X.ALLELE.FREQ. <- as.numeric(gsub(".:","",missense_del_background_freq$X.ALLELE.FREQ.))

missense_del_tf_freq <- read.table(file="missense_del_tf_freq.frq",row.names = NULL)
missense_del_tf_freq$X.ALLELE.FREQ. <- as.numeric(gsub(".:","",missense_del_tf_freq$X.ALLELE.FREQ.))

missense_tolerated_background_freq <- read.table(file="missense_tolerated_background_freq.frq",row.names = NULL)
missense_tolerated_background_freq$X.ALLELE.FREQ. <- as.numeric(gsub(".:","",missense_tolerated_background_freq$X.ALLELE.FREQ.))

missense_tolerated_tf_freq <- read.table(file="missense_tolerated_tf_freq.frq",row.names = NULL)
missense_tolerated_tf_freq$X.ALLELE.FREQ. <- as.numeric(gsub(".:","",missense_tolerated_tf_freq$X.ALLELE.FREQ.))

missense_syn_background_freq <- read.table(file="missense_syn_background_freq.frq",row.names = NULL)
missense_syn_background_freq$X.ALLELE.FREQ. <- as.numeric(gsub(".:","",missense_syn_background_freq$X.ALLELE.FREQ.))

missense_syn_tf_freq <- read.table(file="missense_syn_tf_freq.frq",row.names = NULL)
missense_syn_tf_freq$X.ALLELE.FREQ. <- as.numeric(gsub(".:","",missense_syn_tf_freq$X.ALLELE.FREQ.))


cairo_pdf(file="pi_vs_freq.pdf",height=10,width=6, family = "Arial")
par(mfrow = c(4, 2))
plot(x=out.sites.pi_sg$PI,y=missense_stop_background_freq$X.ALLELE.FREQ.,xlab="Nucleotide diversity (π)",ylab="Allele frequency",main="Stop gained, non-TF",ylim=c(0,1),xlim=c(0,0.55))
plot(x=out.sites.pi_sg_tf$PI,y=missense_stop_tf_freq$X.ALLELE.FREQ.,xlab="Nucleotide diversity (π)",ylab="Allele frequency",main="Stop gained, TF",ylim=c(0,1),xlim=c(0,0.55))

plot(x=out.sites.pi_ms$PI,y=missense_del_background_freq$X.ALLELE.FREQ.,xlab="Nucleotide diversity (π)",ylab="Allele frequency",main="Deleterious missense, non-TF",ylim=c(0,1),xlim=c(0,0.55))
plot(x=out.sites.pi_ms_tf$PI,y=missense_del_tf_freq$X.ALLELE.FREQ.,xlab="Nucleotide diversity (π)",ylab="Allele frequency",main="Deleterious missense, TF",ylim=c(0,1),xlim=c(0,0.55))

plot(x=out.sites.pi_to$PI,y=missense_tolerated_background_freq$X.ALLELE.FREQ.,xlab="Nucleotide diversity (π)",ylab="Allele frequency",main="Tolerated missense, non-TF",ylim=c(0,1),xlim=c(0,0.55))
plot(x=out.sites.pi_to_tf$PI,y=missense_tolerated_tf_freq$X.ALLELE.FREQ.,xlab="Nucleotide diversity (π)",ylab="Allele frequency",main="Tolerated missense, TF",ylim=c(0,1),xlim=c(0,0.55))

plot(x=out.sites.pi_s$PI,y=missense_syn_background_freq$X.ALLELE.FREQ.,xlab="Nucleotide diversity (π)",ylab="Allele frequency",main="Synonymous, non-TF",ylim=c(0,1),xlim=c(0,0.55))
plot(x=out.sites.pi_s_tf$PI,y=missense_syn_tf_freq$X.ALLELE.FREQ.,xlab="Nucleotide diversity (π)",ylab="Allele frequency",main="Synonymous, TF",ylim=c(0,1),xlim=c(0,0.55))
dev.off()

### 5. Mutation LOAD analysis for TF and non-TFs ###

## gff file to estimate gene sizes
gff <- read.table(file="../data_files/Triticum_aestivum.IWGSC.51.cds.gff3",sep="\t")
gff <- tidyr::separate(gff,V9,into=c("ID","Parent","protein_id"),sep=";")
gff$Parent <- gsub("Parent=transcript:","",gff$Parent)

## read in vcf file contain genotype calls
missense_tf.recode <- read.table(file="missense_del_tf.recode.vcf",header=T) ## vcf for deleterious missense for transcription factor genes
missense_background.txt.recode <- read.table(file="missense_del_background.recode.vcf",header=T) ## vcf for deleterious missense mutations for background sites
stop_tf.recode <- read.table(file="missense_stop_tf.recode.vcf",header=T) ## vcf for stop gained mutations for transcription factor genes
stop_background.txt.recode <- read.table(file="missense_stop_background.recode.vcf",header=T) ## vcf for stop gained mutations for background sites
tolerated_tf.recode <- read.table(file="missense_tolerated_tf.recode.vcf",header=T) ## vcf for tolerated missense  mutations for transcription factor genes
tolerated_background.txt.recode <- read.table(file="missense_tolerated_background.recode.vcf",header=T) ## vcf for tolerated deleterious mutations for background sites
syn_tf.recode <- read.table(file="missense_syn_tf.recode.vcf",header=T) ## vcf for synonymous mutations for transcription factor genes
syn_background.txt.recode <- read.table(file="missense_syn_background.recode.vcf",header=T) ## vcf for synonymous mutations for background sites

## file to get transcripts associated with sites
variant_effect_output_all <- read.table(file="../data_files/variant_effect_output_all.txt")
colnames(variant_effect_output_all)[2] <- "site"
variant_effect_output_all$site <- gsub(":","_",variant_effect_output_all$site)
variant_effect_output_all$site <- paste("chr",variant_effect_output_all$site,sep="")
variant_effect_output_all <- variant_effect_output_all[grepl("CANONICAL=YES",variant_effect_output_all$V14),] ## only analyse sites in canonical transcripts

## now sum canonical transcript lengths where tfs are found

translist <- variant_effect_output_all[variant_effect_output_all$V1 %in% c(missense_tf.recode$ID,stop_tf.recode$ID,tolerated_tf.recode$ID,syn_tf.recode$ID),]$V5 ## vector of transcript where Deleterious sites are found
totaltflength <- sum(gff[gff$Parent %in% translist,]$V5 - gff[gff$Parent %in% translist,]$V4)

## now sum canonical transcript lengths where background are found

translist <- variant_effect_output_all[variant_effect_output_all$V1 %in% c(missense_background.txt.recode$ID,stop_background.txt.recode$ID,tolerated_background.txt.recode$ID,syn_background.txt.recode$ID),]$V5 ## vector of transcript where Deleterious sites are found
totalbackgroundlength <- sum(gff[gff$Parent %in% translist,]$V5 - gff[gff$Parent %in% translist,]$V4)

## this loop counts the number of homozygous deleterious missense allele per individual
missense_tf.recode <- missense_tf.recode[-c(1:9)]
value_tf <- ""
for (c in 1:ncol(missense_tf.recode)){
  value_tf[c] <- sum(missense_tf.recode[,c] %in% "1/1")
}

## this loop counts the number of homozygous deleterious missense allele per individual
missense_background.txt.recode <- missense_background.txt.recode[-c(1:9)]
value_background <- ""
for (c in 1:ncol(missense_background.txt.recode)){
  value_background[c] <- sum(missense_background.txt.recode[,c] %in% "1/1")
}

## this loop counts the number of homozygous stop gained  allele per individual
stop_tf.recode <- stop_tf.recode[-c(1:9)]
value_tf_stop <- ""
for (c in 1:ncol(stop_tf.recode)){
  value_tf_stop[c] <- sum(stop_tf.recode[,c] %in% "1/1")
}

## this loop counts the number of homozygous stop gained allele per individual
stop_background.txt.recode <- stop_background.txt.recode[-c(1:9)]
value_background_stop <- ""
for (c in 1:ncol(stop_background.txt.recode)){
  value_background_stop[c] <- sum(stop_background.txt.recode[,c] %in% "1/1")
}

## this loop counts the number of homozygous tolerated missense  allele per individual
tolerated_tf.recode <- tolerated_tf.recode[-c(1:9)]
value_tf_tolerated <- ""
for (c in 1:ncol(tolerated_tf.recode)){
  value_tf_tolerated[c] <- sum(tolerated_tf.recode[,c] %in% "1/1")
}

## this loop counts the number of homozygous tolerated deleterious  allele per individual
tolerated_background.txt.recode <- tolerated_background.txt.recode[-c(1:9)]
value_background_tolerated <- ""
for (c in 1:ncol(tolerated_background.txt.recode)){
  value_background_tolerated[c] <- sum(tolerated_background.txt.recode[,c] %in% "1/1")
}

## this loop counts the number of homozygous synonymous allele per individual
syn_tf.recode <- syn_tf.recode[-c(1:9)]
value_tf_syn <- ""
for (c in 1:ncol(syn_tf.recode)){
  value_tf_syn[c] <- sum(syn_tf.recode[,c] %in% "1/1")
}

## this loop counts the number of homozygous synonymous allele per individual
syn_background.txt.recode <- syn_background.txt.recode[-c(1:9)]
value_background_syn <- ""
for (c in 1:ncol(syn_background.txt.recode)){
  value_background_syn[c] <- sum(syn_background.txt.recode[,c] %in% "1/1")
}

## now combine the load estimates for tf and background genes
value_tf <- as.data.frame(as.numeric(value_tf))
colnames(value_tf) <- "number"
value_background <- as.data.frame(as.numeric(value_background))
colnames(value_background) <- "number"
value_tf_stop <- as.data.frame(as.numeric(value_tf_stop))
colnames(value_tf_stop) <- "number"
value_background_stop <- as.data.frame(as.numeric(value_background_stop))
colnames(value_background_stop) <- "number"
value_tf_tolerated <- as.data.frame(as.numeric(value_tf_tolerated))
colnames(value_tf_tolerated) <- "number"
value_background_tolerated <- as.data.frame(as.numeric(value_background_tolerated))
colnames(value_background_tolerated) <- "number"
value_tf_syn <- as.data.frame(as.numeric(value_tf_syn))
colnames(value_tf_syn) <- "number"
value_background_syn <- as.data.frame(as.numeric(value_background_syn))
colnames(value_background_syn) <- "number"

## aside: individuals at 2.5% tails in any of four classes of mutations in TF transcripts
extremeload <- colnames(stop_tf.recode)[unique(c(which(value_tf$number < quantile(value_tf$number,probs=c(0.025))),which(value_tf$number > quantile(value_tf$number,probs=c(0.975))),which(value_tf_stop$number < quantile(value_tf_stop$number,probs=c(0.025))),which(value_tf_stop$number > quantile(value_tf_stop$number,probs=c(0.975))),which(value_tf_tolerated$number < quantile(value_tf_tolerated$number,probs=c(0.025))),which(value_tf_tolerated$number > quantile(value_tf_tolerated$number,probs=c(0.975))),which(value_tf_syn$number < quantile(value_tf_syn$number,probs=c(0.025))),which(value_tf_syn$number > quantile(value_tf_syn$number,probs=c(0.975)))))]
write.table(extremeload,file="extremeload.txt", row.names = F,col.names = F, quote = F)


## now correcting load estimates by the total length of canonical transcripts
delload <- as.data.frame(rbind(as.matrix(cbind((value_background/totalbackgroundlength),"missense_Background")),as.matrix(cbind((value_tf/totaltflength),"missense_Transcription factors")),as.matrix(cbind((value_background_stop/totalbackgroundlength),"stop_Background")),as.matrix(cbind((value_tf_stop/totaltflength),"stop_Transcription factors")))) ## here number of mutations are scaled by total canonical transcript length
synload <- as.data.frame(rbind(as.matrix(cbind((value_background_syn/totalbackgroundlength),"syn_Background")),as.matrix(cbind((value_tf_syn/totaltflength),"syn_Transcription factors")),as.matrix(cbind((value_background_tolerated/totalbackgroundlength),"tolerated_Background")),as.matrix(cbind((value_tf_tolerated/totaltflength),"tolerated_Transcription factors")))) ## here number of mutations are scaled by total canonical transcript length
colnames(delload) <- c("load","Class")
colnames(synload) <- c("load","Class")
delload <- rbind(delload,synload)
delload$load <- as.numeric(delload$load)

## plot load distributions for for tf and background genes

delload <- tidyr::separate(delload,Class,into=c("type","Class"),sep="_")
delload$type <- gsub("stop","Stop gained",delload$type)
delload$type <- gsub("missense","Deleterious missense",delload$type)
delload$type <- gsub("tolerated","Tolerated missense",delload$type)
delload$type <- gsub("syn","Synonymous",delload$type)
delload$Class <- gsub("Background","non-TF",delload$Class)
delload$Class <- gsub("Transcription factors","TF",delload$Class)

## since distributions are mostly normal, general a linear model and use Tukey post hoc to compare treatments
model <- lm(data=delload,load ~ type * Class)
Anova(model)
TukeyHSD(aov(model), conf.level=.95)

## generate a boxplot of the load
# * means P<0.001
dat <- data.frame(y=c(NA,0.0226,0.158,0.105),type=c("Stop gained","Deleterious missense","Tolerated missense","Synonymous"),Class="TF",label=c("***","***","***","***"))
dat2 <- data.frame(y=c(NA,0.0215,0.15,0.1),type=c("Stop gained","Deleterious missense","Tolerated missense","Synonymous"),Class="TF")

ggsave(file="load2.pdf",width=93, height=118,device = cairo_pdf, units="mm")
ggplot(delload,aes(y=load*1000,x=Class))+
  facet_wrap(~factor(type,levels = c("Stop gained","Deleterious missense","Tolerated missense","Synonymous")),scales = "free")+
  geom_violin(fill="lightgray")+
  geom_boxplot(width=0.2,fill="white")+
  geom_text(data=dat,aes(x=1.5,y=y,label=label))+
  geom_segment(data=dat2,aes(x=1,xend=2,y=y,yend=y))+
  theme_bw()+
  guides(color="none")+
  xlab("") +
  ylab("Number of mutations per kilobase")+
  theme(legend.position="none",panel.grid = element_blank(), axis.text = element_text(colour = "black", size = 9, family="Arial"))+
  theme(axis.title.y = element_text(colour = "black", size = 10, family="Arial"))
dev.off()

delload %>%
  dplyr::group_by(type,Class) %>%
  dplyr::summarise(median(load))

