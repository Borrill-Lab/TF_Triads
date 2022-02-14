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
module load BEDTools/2.29.2-GCC-9.3.0

cd /rds/projects/b/borrillp-phd-students/Catherine/TF_homoeolog_SNPs

bedtools intersect -v -a ./4_analysis/He19F_coding_variant_in_triad_genes_fields10_header.vcf -b /rds/projects/b/borrillp-phd-students/Catherine/expression_and_load_data/allsel.bed > ./4_analysis/He19_coding_variant_not_in_sweep_regions.vcf
# -v outputs records in -a which are not in -b

echo finished


