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

echo starting
grep  01:01:01  ./1_raw/v1.1_genes_TF_homoeolog_info.txt > ./2_triads/triads.txt 
echo finished
