#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=6G
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=[your e-mail@email.com]
#SBATCH --job-name="FineMap"

# These are needed modules in UT HPC to get singularity and Nextflow running. Replace with appropriate ones for your HPC.
module load java-1.8.0_40
module load singularity/3.5.3
module load squashfs/4.4

# Define paths
nextflow_path=[full path to your Nextflow executable]

gwaslist=[path to gwas list file]
genotypefolder=[path to genotype .vcf.gz used]
imputationfile=[path to imputation quality file]
output_path=[name of the output path]

# Command
${nextflow_path}/nextflow run FineMap.nf \
--gwaslist ${gwaslist} \
--genotypefolder ${genotypefolder} \
--imputationfile ${imputationfile} \
--outdir ${output_path}  \
-profile slurm,singularity \
-resume
