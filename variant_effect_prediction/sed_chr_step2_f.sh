#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --time 00-05:00
#SBATCH --qos castles
#SBATCH --account borrillp-phd-students
#SBATCH --mail-type ALL
#SBATCH -o /rds/projects/b/borrillp-phd-students/Catherine/TF_homoeolog_SNPs/slurm_output/%x.%N.%j.out # STDOUT
#SBATCH -e /rds/projects/b/borrillp-phd-students/Catherine/TF_homoeolog_SNPs/slurm_output/%x.%N.%j.err # STDERR


set -e
module purge; module load bluebear # this line is required

cd /rds/projects/b/borrillp-phd-students/Catherine/TF_homoeolog_SNPs

sed 's/chr//' /rds/projects/b/borrillp-phd-students/Catherine/He2019/all.GP08_mm75_het3_publication01142019.vcf > /rds/projects/b/borrillp-phd-students/Catherine/He2019/He19_no_chr.vcf

echo finished
