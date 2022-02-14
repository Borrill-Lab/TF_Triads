#26/01/22
#Script to identify "synonymous only" sites
#Script written by Arun
#

setwd("/rds/projects/b/borrillp-phd-students/Catherine/TF_homoeolog_SNPs/3_vep")

variant_effect_output_all <- read.table(file="variant_effect_output_all.txt")
#variant_effect_output_all$site <- gsub(":","_",variant_effect_output_all$V2)
#variant_effect_output_all$site <- paste("chr",variant_effect_output_all$site,sep="")
variant_effect_output_all$site <- variant_effect_output_all$V2 #Catherine prefers original VEP format
synonly <- variant_effect_output_all[variant_effect_output_all$V7 %in% "synonymous_variant",]$V1
synonly <- variant_effect_output_all[variant_effect_output_all$V1 %in% synonly,]
synonly <- synonly[!(synonly$V7 %in% "synonymous_variant"),] ## identify "synonymous" sites with other annotations

setwd("/rds/projects/b/borrillp-phd-students/Catherine/expression_and_load_data")
write.csv(synonly$site, file="sites_to_exclude_2.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
