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

cd /rds/projects/b/borrillp-phd-students/Catherine/TF_homoeolog_SNPs/scripts

module load R/4.0.3-foss-2020b

srun Rscript stacked_plot_by_tf_family_step7_a.R

