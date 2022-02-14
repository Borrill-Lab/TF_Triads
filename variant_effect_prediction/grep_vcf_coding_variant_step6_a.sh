#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --time 00-12:00
#SBATCH --qos castles
#SBATCH --account borrillp-phd-students
#SBATCH --mail-type ALL
#SBATCH --mem 5G
#SBATCH -o /rds/projects/b/borrillp-phd-students/Catherine/TF_homoeolog_SNPs/slurm_output/%x.%N.%j.out # STDOUT
#SBATCH -e /rds/projects/b/borrillp-phd-students/Catherine/TF_homoeolog_SNPs/slurm_output/%x.%N.%j.err # STDERR

set -e
module purge; module load bluebear # this line is required

cd /rds/projects/b/borrillp-phd-students/Catherine/TF_homoeolog_SNPs

cut -f 1 ./3_filtered/variant_effect_coding_variant_all.txt > ./3_filtered/variant_effect_coding_variant_fields1.txt

grep -F -e '#CHROM' --file ./3_filtered/variant_effect_coding_variant_fields1.txt /rds/projects/b/borrillp-phd-students/Catherine/He2019/all.GP08_mm75_het3_publication01142019.vcf > ./4_analysis/He19F_coding_variant_in_triad_genes.vcf

