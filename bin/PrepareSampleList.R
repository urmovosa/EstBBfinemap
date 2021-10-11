library(data.table)
library(dplyr)
library(stringr)
library(optparse)

setDTthreads(4)

option_list <- list( 
    make_option(c("-s", "--samples_file"), type = "character",
    help = "File including samples and their phenotypes."),
    make_option(c("-p", "--phenotype"), type = "character", 
    help = "Phenotype name.")
)

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)

pheno <- fread(args$samples_file, na.strings = "NA")
pheno <- pheno[, colnames(pheno) %in% c("VKOOD", args$phenotype), with = FALSE]
colnames(pheno)[2] <- "Phenotype"
pheno <- pheno[!is.na(pheno$Phenotype), ]

fwrite(pheno, paste0(args$phenotype, "_PhenoList.txt"), sep = "\t", quote = FALSE, row.names = FALSE)