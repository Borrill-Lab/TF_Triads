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

##unzip command commented out as already done
#gunzip ./gff_files/Triticum_aestivum.IWGSC.51.gff3.gz
#echo gff unzipped
grep -F --file ./2_triads/triad_genenames.txt ./gff_files/Triticum_aestivum.IWGSC.51.gff3 > ./2_triads/Triticum_aestivum.IWGSC.51.triads.gff3

