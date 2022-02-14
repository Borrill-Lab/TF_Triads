#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --time 00-12:00
#SBATCH --qos castles
#SBATCH --account borrillp-phd-students
#SBATCH --mail-type ALL
#SBATCH --mem 10G
#SBATCH -o /rds/projects/b/borrillp-phd-students/Catherine/TF_homoeolog_SNPs/slurm_output/%x.%N.%j.out # STDOUT
#SBATCH -e /rds/projects/b/borrillp-phd-students/Catherine/TF_homoeolog_SNPs/slurm_output/%x.%N.%j.err # STDERR


set -e
module purge; module load bluebear # this line is required
module load BEDTools/2.29.2-GCC-9.3.0

cd /rds/projects/b/borrillp-phd-students/Catherine/TF_homoeolog_SNPs

bedtools intersect -wa -sorted -a /rds/projects/b/borrillp-phd-students/Catherine/He2019/He19_no_chr.vcf -b ./2_triads/Triticum_aestivum.IWGSC.51.genesintriads2.gff3 > ./2_variants/variants_in_genes_in_triads_He_2019.vcf

echo finished


