library(data.table)
library(tidyverse)
library(susieR)
library(Cairo)

# Parse arguments ----
# Argument parser
option_list <- list( 
    make_option(c("-s", "--SummaryStatisticsFile"), type = "character",
    help = "Summary statistics file. Needs to adhere to specific format and have extension .txt."),
    make_option(c("-l", "--LdFile"), type = "character", 
    help = "SNP-SNP correlation file in .rds format. SNP names have to be in the format chr:pos_A1_A2, where A1 and A2 are alphabetically ordered."),
	make_option(c("-i", "--MaxIter"), default = 100, 
    help = "Number of maximum iterations for SuSiE to run. Defaults to 100."),
    make_option(c("-c", "--MaxCausalSnps"), default = 10, 
    help = "Maximum number of causal SNPs SuSie considers. Defaults to 10")
    make_option(c("-o", "--output"), type = "character", 
    help = "Output rds file, containing the SuSie output object with finemapping results. Has to have an extension .rds.")
    )
parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)

# Report settings ----
message(paste(""))
message(paste("Input summary statistics file:", args$SummaryStatisticsFile))
message(paste("Input LD matrix:", args$LdFile))
message(paste("Output file:", args$output))
message(paste("Number of iterations:", args$MaxIter))
message(paste("Number of causal SNPs:", args$MaxCausalSnps))
message(paste("-----------------------------------"))
message(paste(""))
message(paste(""))
message(paste(""))

# Read in and process summary statistics file ----
message("Read in and process summary statistics file")

sumstats <- fread(args$SummaryStatisticsFile)

# Read in and process LD matrix file ----
message("Read in and process LD matrix file")

LDmat <- readRDS(args$LdFile)
LDmat <- as.data.frame(LDmat)
rownames(LDmat) <- colnames(LDmat)
LDmat <- as.matrix(LDmat)
# Force matrix to be symmetric (probably numeric precision problem)
# NB!!! there is probably a bug in emeraLD
LDmat[lower.tri(LDmat)] <- t(LDmat)[lower.tri(LDmat)]


# Sanity check: are SNPs the same in sumstats and LD matrix
if(!nrow(LDmat) == nrow(sumstats)){stop("Different number of SNPs in inputs!")}
as.data.frame(table(colnames(LDmat) == sumstats$UniqueSnpId))

# Run Susie finemapping ----
message("Run Susie finemapping")
# Convert beta and se(beta) to Z-score
sumstats$Z <- sumstats$BETA/sumstats$SE

fitted_rss <- susie_rss(z = sumstats$Z, 
                        R = LDmat,
                        estimate_residual_variance = TRUE,
                        estimate_prior_variance = TRUE,
                        max_iter = as.numeric(args$MaxIter),
                        L = as.numeric(args$MaxCausalSnps),
                        track_fit = TRUE,
                        check_R = TRUE,
                        check_z = TRUE)
