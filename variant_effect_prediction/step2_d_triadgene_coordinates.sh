#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --time 5
#SBATCH --qos bbshort
#SBATCH --account borrillp-phd-students
#SBATCH --mail-type ALL
#SBATCH -o /rds/projects/b/borrillp-phd-students/Catherine/TF_homoeolog_SNPs/slurm_output/%x.%N.%j.out # STDOUT
#SBATCH -e /rds/projects/b/borrillp-phd-students/Catherine/TF_homoeolog_SNPs/slurm_output/%x.%N.%j.err # STDERR

set -e
module purge; module load bluebear # this line is required

cd /rds/projects/b/borrillp-phd-students/Catherine/TF_homoeolog_SNPs

#grep \tgene\t ./2_triads/Triticum_aestivum.IWGSC.51.triads.gff3 > ./2_triads/Triticum_aestivum.IWGSC.51.genesintriads2.gff3

awk --field-separator '\t' '$3 == "gene" {print}' ./2_triads/Triticum_aestivum.IWGSC.51.triads.gff3  > ./2_triads/Triticum_aestivum.IWGSC.51.genesintriads2.gff3

