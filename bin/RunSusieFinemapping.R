library(data.table)
library(dplyr)
library(stringr)
library(susieR)
library(optparse)

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
    help = "Maximum number of causal SNPs SuSie considers. Defaults to 10"),
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

LDmat <- fread(args$LdFile)
LDmat <- as.data.frame(LDmat)
LDmat <- LDmat[, -1]
rownames(LDmat) <- colnames(LDmat)
LDmat <- as.matrix(LDmat)

# Sanity check: are SNPs the same in sumstats and LD matrix
if(!nrow(LDmat) == nrow(sumstats)){stop("Different number of SNPs in inputs!")}
as.data.frame(table(colnames(LDmat) == sumstats$UniqueSnpId))

# Run Susie finemapping ----
message("Run Susie finemapping")
# Convert beta and se(beta) to Z-score
# Test: standardize beta and se(beta)
sumstats$MAF <- sumstats$AF_Allele2
sumstats[MAF > 0.5, ]$MAF <- 1 - sumstats[MAF > 0.5, ]$MAF

sumstats$beta_s <- sumstats$BETA * sqrt(2 * sumstats$MAF * (1 - sumstats$MAF))
sumstats$se_s <- sumstats$SE * sqrt(2 * sumstats$MAF * (1 - sumstats$MAF))

sumstats$Z <- sumstats$beta_s/sumstats$se_s

fitted_rss <- susie_rss(z = sumstats$Z, 
                        R = LDmat,
                        estimate_residual_variance = TRUE,
                        estimate_prior_variance = TRUE,
                        max_iter = as.numeric(args$MaxIter),
                        L = as.numeric(args$MaxCausalSnps),
                        track_fit = TRUE,
                        check_R = TRUE,
                        check_z = TRUE)


fitted_rss <- susie_rss(z = sumstats$Z[1:500], 
                        R = LDmat[1:500, 1:500],
                        estimate_residual_variance = TRUE,
                        estimate_prior_variance = TRUE,
                        max_iter = 100,
                        L = 10,
                        track_fit = TRUE,
                        check_R = FALSE,
                        check_z = TRUE)

saveRDS(fitted_rss, file = args$output)
