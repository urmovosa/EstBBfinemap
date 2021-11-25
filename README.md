## Pipeline for running fine-mapping

This pipeline uses summary statistics from EstBB GWAS's, corresponding genotype data for LD reference and sample phenotypes (e.g. case-control status) to run fine-mapping with SuSiE for every GWAS locus.

### What does it do?

For every GWAS summary statistics file it:

- Uses GWAS P-value threshold to extract lead SNPs.
- Uses genomic window to construct loci.
- Uses MAF and INFO score thresholds to filter SNPs by QC.
- Includes only loci which have >= 50 variants after filtering.
- Uses genotype data (.bgen format) from the same samples to construct LD matrix for each locus (by LDstore2).
- Runs SuSiE fine-mapping (function: susie_suff_stat) for each locus and outputs SuSiE results.

### Instructions

#### Requirements

- You need to have Java 8 installed in your HPC environment.
- You need to have Singularity installed in your HPC environment.
- You need SLURM as scheduler.

#### Input

- Space-delimited text file with following information and headers:
    - PhenoName: name of the GWAS phenotype.
    - SumStat: path to the corresponding GWAS summary statistics file. Has to be in the SAIGE output format and can be (b)gzipped.
    - SampleFile: file with phenotype information for samples included to GWAS. Has to contain headers: VKOOD (sample ID) and corresponding phenotype name. This file can contain columns for multiple different GWASs/phenotypes.
- Genotype folder for samples used in GWASs (.bgen format). This is automatically subsetted based on phenotype information file and LD matrices for GWAS loci are constructed, based on exactly the same samples as were in the GWAS.
- Imputation information file. File which contains imputation INFO score for each variant in the genotype data. Has to contain columns: ID, CHR, POS, REF, ALT and INFO and can be gzipped. **NB!** in the future this info should be available in the summary statistics files not in the separate file.

#### Settings

Mandatory arguments:

`--gwaslist`  Space separated file with header (column names: PhenoName, SumStat, SampleFile) and three columns. First column: phenotype name, second column: path to gwas summary statistics file, third column: path to the file which contains measurements for given phenotype (binary: 0, 1, NA; continuous: continuous numbers). Has to contain phenotype as a column name and column named "VKOOD" for sample IDs.

`--genotypefolder`    Folder containing bgen files on which all those GWAS's were ran. File names have to contain the string "chr[1-23]".

`--imputationfile`    Separate file containing imputation INFO score for each SNP in the genotype data.

`--outdir`            Folder where output is written.

Optional arguments:

Data management settings:

Defaults for those are specified in the nextflow.conf, corresponding to SAIGE file format and EstBB sample ID name.

`--ChrCol`  Name of the chromosome column in the summary statistics files.

`--ChrPosCol`   Name of the genomic position column in the summary statistics files.

`--RefAllCol`   Name of the reference allele column in the summary statistics files.

`--EffAllCol`   Name of the effect allele column in the summary statistics files.

`--MafCol`  Name of the MAF column in the summary statistics files.

`--BetaCol` Name of the beta column in the summary statistics files.

`--SeCol`   Name of the se(beta) column in the summary statistics files.

`--PvalCol` Name of the P-value column in the summary statistics files.

`--SampleId`    Name of the sample ID in the sample Sample file.

Filtering:

`--PvalThresh`    GWAS P-value threshold for defining significant loci. Defaults to 5e-8.

`--Win`   Genomic window to extract loci for finemapping. Defaults 1000000bp to either side of lead SNP.

`--MafThresh`   MAF threshold to filter the input GWAS data. Defaults to 0.01.

`--InfoThresh`    INFO score threshold to filter the input GWAS data. Defaults to 0.4.

SuSiE settings:

`--MaxCausalSnps`   Maximal number of causal SNPs tested in SuSiE analysis. Defaults to 10.

#### Command

SLURM template script:

```
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
```

#### Output

Outputs R .rds files for each SuSiE object tested and .txt.gz files with SNP information (including PIPs and credible set membership). It also outputs .html report with per-locus summary and plots.

### TODO

- implement .Rmd report with overview and plots.
- Check the susie_rss in addition to susie_suff_stat.

## Acknowledgements

SuSiE is the work of Wang and colleagues:

[Wang, G., Sarkar, A., Carbonetto, P., & Stephens, M. (2020). A simple new approach to variable selection in regression, with application to genetic fine mapping. Journal of the Royal Statistical Society, Series B 82, 1273–1300. https://doi.org/10.1111/rssb.12388](https://doi.org/10.1111/rssb.12388)

[Website](https://stephenslab.github.io/susieR/index.html)

LDstore2 is the work of Benner and colleagues:

[Benner, C. et al. Prospects of fine-papping trait-associated genomic regions by using summary statistics from genome-wide association studies. Am. J. Hum. Genet. (2017).](https://www.sciencedirect.com/science/article/pii/S0002929717303348?via%3Dihub)

[Website](http://www.christianbenner.com/#)

Some parts of this pipeline were adapted from [FinnGen fine-mapping pipeline](https://github.com/FINNGEN/finemapping-pipeline) by Masahiro Kanai and colleagues.
