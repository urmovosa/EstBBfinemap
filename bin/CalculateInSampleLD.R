library(optparse)

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



# Functions
emeraLD2R <- function(path, bin = "bin/emeraLD"){
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

in_path <- "example/chr$chr.1KG.25K_m.m3vcf.gz"

getLD <- emeraLD2R(path = in_path)

ld_data <- getLD(region = "20:83061-92955")
