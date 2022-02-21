#!/usr/bin/env Rscript

library(data.table)

args = commandArgs(trailingOnly=TRUE)

phen <- fread(args[1])
colnames(phen)[1] <- "Sample"
trait <- "Phenotype"

if (length(colnames(phen)[!colnames(phen) %in% trait]) == 0){break("Phenotype not in pheno file!")} else {

    phen <- phen[, colnames(phen) %in% c("Sample", trait), with = FALSE]
    colnames(phen)[2] <- "Phenotype"
    if (length(phen$Phenotype[phen$Phenotype %in% c(0, 1, NA)]) < length(phen$Phenotype)){
        #message("Phenotype seems to be quantitative, setting var(y)=1")
        var_y <- 1
        #message(paste("var(y) is set to", var_y))
    } else {
            phi <- length(phen$Phenotype[phen$Phenotype == 1])/(length(phen$Phenotype[phen$Phenotype == 1]) + length(phen$Phenotype[phen$Phenotype == 0]))
            phi = vals["1"] / (vals["1"]+vals["0"])
            var_y = phi*(1-phi)
        #message(paste("var(y) is set to", var_y))

    }

}

cat(var_y)
