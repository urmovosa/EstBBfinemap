# Functions ----

#' Convert SNP information to more unique SNP ID.
#'
#' This function takes SNP chromosome, position, allele information and constructs the new SNP ID in the form chr:pos_FirsAllele_SecondAllele, where alleles are ordered alphabetically. This is useful when comparing data (e.g. effect directions) from SNPs which come from different studies, potentially using different dbSNP builds. Pay attention that function does not take care of differences in genome build (hg18, hg19, hg38). It also does not take care of potential strand issues, therefore care is needed when working with older studies and A/T or C/G SNPs.
#'
#' @param chr Vector with chromosome position for each SNP. Expects 1:22, X, Y, MT. Chr/chr at the beginning is allowed but is removed from output SNP ID.
#' @param pos Vector with genomic positions for each SNP. Function does not check any genomic build (hg18, hg19, hg38) or whether data is in 0-based/1-based format.
#' @param allele1 Vector with first allele for each SNP. Program does not make any checks, so indels are also supported. "Missing" alleles should be denoted as "-" (E.g. alleles: "-/A")
#' @param allele2 Vector with second allele for each SNP.
#'
#' @return Returns vector with new SNP IDs in the format chr:pos_OneAllele_SecondAllele, where alleles are ordered alphabetically. This function might be useful when one needs to compare SNPs from different studies.
#' @export
#'
#' @examples
#' sample_data <- data.frame(
#'   rs = c("rs100", "rs101", "rs102", "rs103", "rs104", "rs105", "rs100"),
#'   chr = c("1", "5", "X", "12", "chr3", "chr3", "1"),
#'   pos = c(10034, 1341351, 13515, 23413153, 1342525, 1342525, 10034),
#'   all1 = c("T", "A", "C", "G", "A", "A", "T"),
#'   all2 = c("C", "T", "G", "A", "-", "AT", "C")
#' )
#'
#' sample_data$SNPID <- ConvUniqSNPName(chr = sample_data$chr, pos = sample_data$pos, allele1 = sample_data$all1, allele2 = sample_data$all2)
ConvUniqSNPName <- function(chr = chr, pos = pos, allele1 = a1, allele2 = a2) {
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(data.table)

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
  