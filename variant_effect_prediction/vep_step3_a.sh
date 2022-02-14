#!/bin/bash
#SBATCH --ntasks 4
#SBATCH --time 00-12:00
#SBATCH --qos castles
#SBATCH --account borrillp-phd-students
#SBATCH --mail-type ALL
#SBATCH -o /rds/projects/b/borrillp-phd-students/Catherine/TF_homoeolog_SNPs/slurm_output/%x.%N.%j.out
#SBATCH -e /rds/projects/b/borrillp-phd-students/Catherine/TF_homoeolog_SNPs/slurm_output/%x.%N.%j.err


set -e
module purge; module load bluebear # this line is required
module load VEP/99.2-foss-2019b-Perl-5.30.0

cd /rds/projects/b/borrillp-phd-students/Catherine/TF_homoeolog_SNPs

vep --cache --offline --dir_cache "/rds/projects/b/borrillp-phd-students/Catherine/.vep" --cache_version 51 \
 --species triticum_aestivum -v --canonical --sift b --fork 4\
 -i "./2_variants/variants_in_genes_in_triads_He_2019.vcf" -o "./3_vep/variant_effect_output_all.txt"

echo finished variant effect prediction

