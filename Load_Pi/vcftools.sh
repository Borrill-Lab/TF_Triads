## only use sites with minor allele frequency > 0.01

vcftools --vcf all.GP08_mm75_het3_publication01142019.vcf --site-pi --out filtered_sites

## identify number of positively selected and introgressed sites in TFs and non TFs

vcftools --vcf all.GP08_mm75_het3_publication01142019.vcf --recode --positions tf_sites.txt  --out tf_sites
vcftools --vcf all.GP08_mm75_het3_publication01142019.vcf --recode --positions background_sites.txt  --out background_sites

vcftools --vcf tf_sites.recode.vcf --bed all_possel.bed  --out all_possel_tf
vcftools --vcf tf_sites.recode.vcf --bed introgress.bed --out introgress_tf
vcftools --vcf tf_sites.recode.vcf --bed env_sel.bed  --out env_sel_tf
vcftools --vcf tf_sites.recode.vcf --bed improve_sel.bed  --out improve_sel_tf

#all_possel_tf.log:After filtering, kept 16 out of a possible 3142 Sites
#env_sel_tf.log:After filtering, kept 1795 out of a possible 3142 Sites
#improve_sel_tf.log:After filtering, kept 574 out of a possible 3142 Sites
#introgress_tf.log:After filtering, kept 176 out of a possible 3142 Sites

vcftools --vcf background_sites.recode.vcf --bed all_possel.bed  --out all_possel_background
vcftools --vcf background_sites.recode.vcf --bed introgress.bed --out introgress_background
vcftools --vcf background_sites.recode.vcf --bed env_sel.bed  --out env_sel_background
vcftools --vcf background_sites.recode.vcf --bed improve_sel.bed  --out improve_sel_background

#all_possel_background.log:After filtering, kept 421 out of a possible 51619 Sites
#env_sel_background.log:After filtering, kept 31860 out of a possible 51619 Sites
#improve_sel_background.log:After filtering, kept 9776 out of a possible 51619 Sites
#introgress_background.log:After filtering, kept 2815 out of a possible 51619 Sites


## get list of positively selected and introgressed sites

cat all_possel.bed env_sel.bed introgress.bed improve_sel.bed >allsel.bed


## make separate vcfs for each type of variant for TFs and non TFs

vcftools --vcf all.GP08_mm75_het3_publication01142019.vcf --positions missense_Deleterious_background.txt  --freq --out missense_del_background_freq
vcftools --vcf all.GP08_mm75_het3_publication01142019.vcf --positions missense_tolerated_background.txt  --freq --out missense_tolerated_background_freq
vcftools --vcf all.GP08_mm75_het3_publication01142019.vcf --positions missense_stop_background.txt  --freq --out missense_stop_background_freq
vcftools --vcf all.GP08_mm75_het3_publication01142019.vcf --positions missense_syn_background.txt  --freq --out missense_syn_background_freq
vcftools --vcf all.GP08_mm75_het3_publication01142019.vcf --positions missense_syn_tf.txt  --freq --out missense_syn_tf_freq
vcftools --vcf all.GP08_mm75_het3_publication01142019.vcf --positions missense_stop_tf.txt  --freq --out missense_stop_tf_freq
vcftools --vcf all.GP08_mm75_het3_publication01142019.vcf --positions missense_tolerated_tf.txt  --freq --out missense_tolerated_tf_freq
vcftools --vcf all.GP08_mm75_het3_publication01142019.vcf --positions missense_Deleterious_tf.txt  --freq --out missense_del_tf_freq

## make separate vcfs for each type of variant for TFs and non TFs

vcftools --vcf all.GP08_mm75_het3_publication01142019.vcf --recode --positions missense_Deleterious_background.txt  --out missense_del_background
vcftools --vcf all.GP08_mm75_het3_publication01142019.vcf --recode --positions missense_tolerated_background.txt  --out missense_tolerated_background
vcftools --vcf all.GP08_mm75_het3_publication01142019.vcf --recode --positions missense_stop_background.txt  --out missense_stop_background
vcftools --vcf all.GP08_mm75_het3_publication01142019.vcf --recode --positions missense_syn_background.txt  --out missense_syn_background
vcftools --vcf all.GP08_mm75_het3_publication01142019.vcf --recode --positions missense_syn_tf.txt  --out missense_syn_tf
vcftools --vcf all.GP08_mm75_het3_publication01142019.vcf --recode --positions missense_stop_tf.txt  --out missense_stop_tf
vcftools --vcf all.GP08_mm75_het3_publication01142019.vcf --recode --positions missense_tolerated_tf.txt  --out missense_tolerated_tf
vcftools --vcf all.GP08_mm75_het3_publication01142019.vcf --recode --positions missense_Deleterious_tf.txt  --out missense_del_tf

for file in  missense_*_*.recode.vcf; do sed -i 's/#CHROM/CHROM/' $file; done





