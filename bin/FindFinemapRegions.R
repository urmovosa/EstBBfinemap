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
    help = "Genomic window size for defining loci. Defaults to 1,000,000bp (both sides from lead SNP).")
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
  
  # Iteratively identify most significant SNP, and remove all other SNPs inthe window
  res <- data_f[-c(1:nrow(data_f)), ]
  
  while (nrow(data_f) > 0) {
    
    lead_snp <- data_f[data_f$P == min(data_f$P), ]
    res <- rbind(res, lead_snp)
    data_f <- data_f[!(data_f$chr == lead_snp$chr & data_f$pos > lead_snp$pos - window & data_f$pos < lead_snp$pos + window),]
    message(paste("Added:", lead_snp$chr, lead_snp$pos))
  }
  
  # Remove SNPs for which the region overlaps with HLA region (hg19)
  res <- res[!(res$chr == 6 & ((res$pos - window > 28477797 & res$pos - window < 33448354) | (res$pos + window > 28477797 & res$pos + window < 33448354))), ]
  message("SNPs overlapping hg19 MHC region removed!")

  res <- paste0(res$chr, ":", res$pos - window, "-", res$pos + window)
  
  return(res)
}

#########

# Read in summary stats file
sum_stat <- fread(args$gwas_file)

# Read in imputation quality file
imp <- fread(args$info_file)

# Filter:
# MAF
message(paste("Before MAF filter", nrow(sum_stat), "variants."))
sum_stat <- sum_stat[sum_stat$AF_Allele2 > MAF_thresh & 
sum_stat$AF_Allele2 < (1 - MAF_thresh), ]
message(paste("After MAF filter", nrow(sum_stat), "variants."))

sum_stat <- sum_stat[, c(1, 2, 3, 4, 5, 7, 9, 10, 11, 13), with = FALSE]
colnames(sum_stat)[c(1, 2, 10)] <- c("chr", "pos", "P")

if (!nrow(sum_stat[sum_stat$P < P_thresh, ]) > 0){
  message("None of the variants met the specified significance threshold!")

  abi <- data.table(SS = "temp", Region = "temp")
  gwas_file_name <- str_replace(args$gwas_file, "\\..*", "")
  fwrite(abi[-1, ], paste0(gwas_file_name, "_regions.txt"), sep = "\t")

} else {

  regions <- IdentifyLeadSNPs(sum_stat, window = WIN, Pthresh = P_thresh)
  gwas_file_name <- str_replace(args$gwas_file, "\\..*", "")

  fwrite(data.table(SS = gwas_file_name, Region = regions), paste0(gwas_file_name, "_regions.txt"), sep = "\t")
}