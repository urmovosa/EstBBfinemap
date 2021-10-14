library(data.table)
library(optparse)
library(stringr)
library(dplyr)

setDTthreads(4)

# Argument parser
option_list <- list( 
    make_option(c("-g", "--gwas_file"), type = "character",
    help = "GWAS summary statistics file."),
    make_option(c("-i", "--info_file"), type = "character", 
    help = "Help file with SNP imputation info."),
    make_option(c("-r", "--region"), type = "character", 
    help = "Region to filter from summary statistics."),
    make_option(c("-P", "--P_threshold"), default = 5e-8, 
    help = "GWAS P-value threshold. Defaults to 5e-8."),
    make_option(c("-M", "--Maf_threshold"), default = 0.01, 
    help = "MAF threshold for filtering SNPs. Defaults to 0.01."),
    make_option(c("-I", "--INFO_threshold"), default = 0.4, 
    help = "Info score threshold. Defaults to 0.4.")
    )
parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)

# Default thresholds
MAF_thresh <- args$Maf_threshold
INFO_thresh <- args$INFO_threshold
P_thresh <- args$P_threshold
# Functions
ConvUniqSNPName <- function(chr = chr, pos = pos, allele1 = a1, allele2 = a2) {
  .datatable.aware = TRUE

  # Clean "chr" and/or "Chr" from the start of each chromosome name
  message('Removing "chr/Chr" from chromosome name')
  chr <- str_replace(chr, "chr", "")
  chr <- str_replace(chr, "Chr", "")

  # Make chr to character

  chr <- as.character(chr)

  # Check if chromosomes are 1:22, X, Y, MT
  if (!all(chr %in% c(as.character(1:23), "X", "Y", "MT"))) {
    stop("Error: some of the chromosomes are not 1:23, X, Y, MT. Please check!")
  }

  # Combine chr and pos
  chr_pos <- paste(chr, pos, sep = ":")

  # Sort alleles

  allele1 <- as.character(allele1)
  allele2 <- as.character(allele2)

  # Check if deletions are not annotated as empty
  if (any(allele1 == "") | any(allele1 == " ")) {
    stop("Error: some alleles are empty or designated as NA. Those have to be specified as " - " in the data. Please check!")
  }

  # Check if chromosome position is integer
  if (!all(round(pos) == pos)) {
    stop("Error: chromosome position is not integer")
  }

  temp_table <- data.table(chr_pos = chr_pos, all1 = as.character(allele1), all2 = as.character(allele2))

  # If same position happens >1 in the data, add separate indicator
  temp_table$id <- rowidv(temp_table, cols = c("chr_pos"))

  temp_table <- melt(temp_table, measure.vars = c("all1", "all2"))

  temp_table$value <- as.character(temp_table$value)

  # temp_table <- temp_table %>%
  #   group_by(chr_pos, id) %>%
  #   mutate(sort_alleles = paste(sort(unique(value)), collapse = "_"))

  # Attempt to speed up with data.table
  temp_table <- as.data.table(temp_table)
  temp_table[, sort_alleles := paste(sort(unique(as.character(value))), collapse = "_"), by = .(chr_pos, id)]

  temp_table <- temp_table[temp_table$variable != 'all2', ]

  # new SNP ID
  SNPID <- paste(temp_table$chr_pos, temp_table$sort_alleles, sep = "_")
  return(SNPID)
}

chr <- str_replace(args$region, ":.*", "")
reg <- str_replace(args$region, ".*:", "")
start <- str_replace(reg, "-.*", "")
end <- str_replace(reg, ".*-", "")

# start making log file
log_file <- data.table(gwas_file = args$gwas_file, 
region = args$region, 
nr_of_SNPs = 0, 
after_MAF_filter = 0,
after_INFO_filter = 0)

# Read file in
# Read in summary stats file
filter_cmd <- paste0("gunzip -c ", args$gwas_file, " | awk '{ if($1 == ",  as.numeric(chr), " && $2 > ", as.numeric(start), " && $2 < ", as.numeric(end), ") { print }}' ")
sum_stat <- fread(cmd = filter_cmd)
help_head <- fread(cmd = paste0("gunzip -c ", args$gwas_file, " | head -n 2"))
colnames(sum_stat) <- colnames(help_head)

message("Sumstats read!")

log_file$nr_of_SNPs <- nrow(sum_stat)

# Read in imputation quality file
imp <- fread(args$info_file)
message("Imputation file read!")

# Filter:
## Based on region
imp <- imp[as.character(imp$CHR) == chr & as.numeric(imp$POS) > as.numeric(start) & as.numeric(imp$POS) < as.numeric(end), ]
imp <- imp[, c(2, 3, 4, 5, 7), with = FALSE]

## MAF
message(paste("Before MAF filter", nrow(sum_stat), "variants."))
sum_stat <- sum_stat[sum_stat$AF_Allele2 > MAF_thresh & 
sum_stat$AF_Allele2 < (1 - MAF_thresh), ]
message(paste("After MAF filter", nrow(sum_stat), "variants."))

sum_stat <- sum_stat[, c(1, 2, 3, 4, 5, 7, 9, 10, 11, 13), with = FALSE]
colnames(sum_stat)[c(1, 2, 10)] <- c("chr", "pos", "P")

log_file$after_MAF_filter <- nrow(sum_stat)

## Based on INFO score
### Unique SNP IDs
sum_stat$UniqueSnpId <- ConvUniqSNPName(chr = sum_stat$chr, pos = sum_stat$pos, allele1 = sum_stat$Allele1, allele2 = sum_stat$Allele2)
imp$UniqueSnpId <- ConvUniqSNPName(chr = imp$CHR, pos = imp$POS, allele1 = imp$REF, allele2 = imp$ALT)
### Merge and filter
sum_stat <- merge(sum_stat, imp[, c(6, 5), with = FALSE], by = "UniqueSnpId")
sum_stat <- sum_stat[sum_stat$INFO > INFO_thresh, ]
message(paste("After INFO score filter", nrow(sum_stat), "variants."))

log_file$after_INFO_filter <- nrow(sum_stat)

# REF_ALT SNP list
ref_alt_SNP <- paste0(sum_stat$chr, ":", sum_stat$pos, "", sum_stat$Allele1, ",", sum_stat$Allele2)

### Extra check: are there any SNPs remaining in the region which pass the P-value threshold?
# Write out filtered file
gwas_file_name <- str_replace(args$gwas_file, "\\..*", "")
fwrite(sum_stat, paste0(gwas_file_name, "_", args$region, "_region.txt"), sep = "\t", quote = FALSE)
fwrite(as.data.table(ref_alt_SNP), "variants_filter.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
fwrite(log_file, paste0(gwas_file_name, '_', args$region, '.log'), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
