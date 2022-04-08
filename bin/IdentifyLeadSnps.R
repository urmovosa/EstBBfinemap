IdentifyLeadSNPs <- function(data, window = 1000000, Pthresh = 5e-8, RemoveHLA = TRUE) {

  data_f <- data[data$P < as.numeric(Pthresh), ]
  
  # Iteratively identify most significant SNP, and remove all other SNPs in the window
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
