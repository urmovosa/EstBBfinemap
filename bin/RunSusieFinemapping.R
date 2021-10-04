library(data.table)
library(tidyverse)
library(susieR)
library(Cairo)

# Parse arguments ----
parser <- ArgumentParser()

parser$add_argument("-s", "--SummaryStatisticsFile", help = "Summary statistics file. Needs to adhere to specific format and have extension .txt.", required = TRUE)
parser$add_argument("-ld", "--LdFile", help = "SNP-SNP correlation file. It has to be calculated with R command cor and have extension .txt. It also has to contain chromosome in the file name, in the format of chr[num].", required = TRUE)
parser$add_argument("-o", "--output", help = "Output rds file, containing the SuSie output object with finemapping results. Has to have extension .rds.", required = TRUE)
parser$add_argument("--MaxIter", default = 100, help = "Number of maximum iterations for SuSie to run.", required = FALSE)
parser$add_argument("--MaxCausalSnps", default = 30, help = "Maximum number of causal SNPs SuSie considers.", required = FALSE)
parser$add_argument("-imp", "--ImpThreshold", default = 0, help = "Imputation INFO score threshold to filter the input for. Defaults to 0 (no filtering).", required = FALSE)
parser$add_argument("-m", "--MafThreshold", default = 0, help = "Minor allele frequency threshold to filter the input for. Defaults to 0 (no filtering).", required = FALSE)

args <- parser$parse_args()

# Report settings ----
message(paste(""))
message(paste("Input summary statistics file:", args$SummaryStatisticsFile))
message(paste("Input LD matrix:", args$LdFile))
message(paste("Output file:", args$output))
message(paste("Number of iterations:", args$MaxIter))
message(paste("Number of causal SNPs:", args$MaxCausalSnps))
message(paste("Imputation info score filter:", args$ImpThreshold))
message(paste("MAF filter:", args$MafThreshold))
message(paste("-----------------------------------"))
message(paste(""))
message(paste(""))
message(paste(""))
