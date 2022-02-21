library(data.table)
library(optparse)
library(dplyr)
library(stringr)

setDTthreads(4)

# Argument parser
option_list <- list( 
    make_option(c("-g", "--gwas_file"), type = "character",
    help = "GWAS summary statistics file."),
    make_option(c("-i", "--info_file"), type = "character", 
    help = "Help file with SNP imputation info."),
    make_option(c("-P", "--P_threshold"), default = 5e-8,
    help = "GWAS P-value threshold. Defaults to 5e-8."),
    make_option(c("-M", "--Maf_threshold"), default = 0.01,
    help = "MAF threshold for filtering SNPs. Defaults to 0.01."),
    make_option(c("-I", "--INFO_threshold"), default = 0.4,
    help = "Info score threshold. Defaults to 0.4."),
    make_option(c("-W", "--win_size"), default = 1000000, 
    help = "Genomic window size for defining loci. Defaults to 1,000,000bp (both sides from lead SNP)."),
    make_option(c("-C", "--chr_col"), type = "character",
    help = "Name of the chromosome column in the summary statistics files."),
   make_option(c("-n", "--pos_col"), type = "character",
    help = "Name of the genomic position column in the summary statistics files."),
    make_option(c("-r", "--ref_allele_col"), type = "character",
    help = "Name of the reference allele column in the summary statistics files."),
    make_option(c("-e", "--eff_allele_col"), type = "character",
    help = "Name of the effect allele column in the summary statistics files."),
    make_option(c("-m", "--maf_col"), type = "character",
    help = "Name of the MAF column in the summary statistics files."),
    make_option(c("-b", "--beta_col"), type = "character",
    help = "Name of the beta column in the summary statistics files."),
    make_option(c("-s", "--se_col"), type = "character",
    help = "Name of the se(beta) column in the summary statistics files."),
    make_option(c("-p", "--p_col"), type = "character",
    help = "Name of the P-value column in the summary statistics files.")
    )
parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)

# Default thresholds
MAF_thresh <- args$Maf_threshold
INFO_thresh <- args$INFO_threshold
WIN <- args$win_size
P_thresh <- args$P_threshold
# Functions
IdentifyLeadSNPs <- function(data, window = 1000000, Pthresh = 5e-8, RemoveHLA = TRUE) {

  data_f <- data[data$P < as.numeric(Pthresh), ]
  data_f$Z <- data_f$beta / data_f$se
  
  # Iteratively identify most significant SNP, and remove all other SNPs in the window
  res <- data_f[-c(1:nrow(data_f)), ]
  
  while (nrow(data_f) > 0) {
    
    lead_snp <- data_f[abs(data_f$Z) == max(data_f$Z), ][1, ] # If multiple SNPs have exactly the same effect, take first
    res <- rbind(res, lead_snp[, -ncol(lead_snp), with = FALSE])
    data_f <- data_f[!(data_f$chr == lead_snp$chr & data_f$pos > lead_snp$pos - window & data_f$pos < lead_snp$pos + window),]
    message(paste("Added:", lead_snp$chr, lead_snp$pos))
  }
  
  # Remove SNPs for which the region overlaps with HLA region (hg19)
  res <- res[!(res$chr == 6 & ((res$pos - window > 28477797 & res$pos - window < 33448354) | (res$pos + window > 28477797 & res$pos + window < 33448354))), ]
  message("SNPs overlapping hg19 MHC region removed!")

  res <- paste0(res$chr, ":", res$pos - window, "-", res$pos + window)
  
  return(res)
}

test_data <- data.table(CHR = c(1, 3, 4, 5, "X"))

AdjustFileFormat <- function(input, chr = "CHR", pos = "POS", ref_allele = "ref", 
effect_allele = "alt", maf = "maf", p_value = "p", beta = "beta", se = "se"){

  input <- as.data.table(input)
  output <- input[, colnames(input) %in% c(chr, pos, ref_allele, effect_allele, maf, beta, se, p_value), with = FALSE]
  output <- output[, match(c(chr, pos, ref_allele, effect_allele, maf, beta, se, p_value), colnames(output)), with = FALSE]
  colnames(output) <- c("chr", "pos", "ref_allele", "effect_allele", "maf", "beta", "se", "P")
  
  output <- as.data.frame(output)
  output[output$maf > 0.5, ]$maf <- 1 - output[output$maf > 0.5, ]$maf

  output$chr <- str_replace(output$chr, "chr", "")
  if (nrow(output[output$chr == "X", ]) > 1){
  output[output$chr == "X", ]$chr <- "23"
  }

  return(output)
}

#########

# Read in summary stats file
sum_stat <- fread(args$gwas_file)

# Convert summary stats into correct format
sum_stat <- AdjustFileFormat(sum_stat, 
                              chr = args$chr_col,
                              pos = args$pos_col,
                              ref_allele = args$ref_allele_col,
                              effect_allele = args$eff_allele_col,
                              maf = args$maf_col,
                              beta = args$beta_col,
                              se = args$se_col,
                              p_value = args$p_col)

fwrite(sum_stat, paste0("standardized_", args$gwas_file), sep = "\t", quote = FALSE, row.names = FALSE)

# Filter:
# MAF
message(paste("Before MAF filter", nrow(sum_stat), "variants."))
sum_stat <- sum_stat[sum_stat$maf > MAF_thresh, ]
message(paste("After MAF filter", nrow(sum_stat), "variants."))

fwrite(sum_stat, paste0("standardized_", args$gwas_file), sep = "\t", quote = FALSE, row.names = FALSE)


if (!nrow(sum_stat[sum_stat$P < P_thresh, ]) > 0){
  message("None of the variants met the specified significance threshold!")

  abi <- data.table(SS = "temp", Region = "temp")
  gwas_file_name <- str_replace(args$gwas_file, "\\..*", "")
  abi$Region <- str_replace(abi$Region, "23:", "X:")
  fwrite(abi[-1, ], paste0(gwas_file_name, "_regions.txt"), sep = "\t")

} else {

  regions <- IdentifyLeadSNPs(sum_stat, window = WIN, Pthresh = P_thresh)
  gwas_file_name <- str_replace(args$gwas_file, "\\..*", "")
  regions <- str_replace(regions, "23:", "X:")
  fwrite(data.table(SS = gwas_file_name, Region = regions), paste0(gwas_file_name, "_regions.txt"), sep = "\t")

}
