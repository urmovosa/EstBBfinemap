library(data.table)
library(optparse)
library(stringr)
library(dplyr)

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
  if (!all(chr %in% c(as.character(1:22), "X", "Y", "MT"))) {
    stop("Error: some of the chromosomes are not 1:22, X, Y, MT. Please check!")
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

# Argument parser
option_list <- list( 
    make_option(c("-f", "--filtered_vcf_file"), type = "character",
    help = "Filtered vcf file."),
    make_option(c("-r", "--region"), type = "character", 
    help = "Region."),
	make_option(c("-o", "--output"), type = "character", 
    help = "Base name for the output file.")
    )
parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)


# Functions
emeraLD2R <- function(path, bin = "/emeraLD/bin/emeraLD"){
	require(data.table)
	if(!file.exists(bin)) stop(paste0("bin = '", bin, "' file does not exist"))
	function(region, matrix.out = TRUE, info.out = TRUE){
		require(data.table)
		opts <- paste(c("--matrix", "--stdout --extra")[c(matrix.out, TRUE)], collapse = " ")
		pc <- "| tr ':' '\t'"
		chr <- strsplit(region, ":")[[1]][1]
		vcfile <- gsub("\\$chr", chr, path)
		if(!file.exists(vcfile)) stop(paste0(vcfile, " does not exist"))
		out <- suppressMessages(fread(
			input  = paste(bin, "-i", vcfile, "--region", region, opts, pc), 
			header = FALSE, showProgress = FALSE
		))
		info <- NULL
		if(info.out){
			info <- out[,1:5]
			colnames(info) <- c("chr", "pos", "id", "ref", "alt")
		}
		if(matrix.out) out <- as.matrix(out[,-(1:5)]); colnames(out) <- NULL
		list("Sigma" = out, "info" = info)
	}
}

getLD <- emeraLD2R(path = args$filtered_vcf_file)
ld_data <- getLD(region = args$region)

info_fields <- as.data.table(ld_data$info)
info_fields$UniqueSnpId <- ConvUniqSNPName(chr = info_fields$chr, pos = info_fields$pos, allele1 = info_fields$ref, allele2 = info_fields$alt)

LdMat <- as.matrix(ld_data$Sigma)
colnames(LdMat) <- info_fields$UniqueSnpId
rownames(LdMat) <- info_fields$UniqueSnpId

#fwrite(as.data.table(LdMat), paste0(args$region, "_LDmatrix.ld"), sep = "\t", quote = FALSE, row.names = FALSE, colnames = FALSE)
saveRDS(as.data.table(LdMat), file = paste0(args$region, "_LDmatrix.rds"))