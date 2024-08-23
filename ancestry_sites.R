#!/bin/env Rscript

library(data.table)
library(dtplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Incorrect arguments", call. = FALSE)

filename <- args[1]

check_genotypes <- function(gt) {
  gt[, unique(unlist(lapply(.SD, function(x) unique(x)))), .SDcols = !c("chrom", "pos", "ref", "alt")]
}

# read the genotype data
gt <- fread(filename, na.strings = "-")
gt[, chrom := paste0("chr", chrom)]

samples <- setdiff(names(gt), c("chrom", "pos", "ref", "alt"))

check_genotypes(gt)

# replace het sites with NA
gt[, (samples) := lapply(.SD, function(x) fifelse(x  %in% c("0|1", "1|0", "0/1"), NA_character_, x)),
.SDcols = samples]

check_genotypes(gt)

# replace diploid genotypes with ref/alt states
gt[, (samples) := lapply(.SD, function(x) { fcase(x %in% c("0|0", "0/0"), 0L,
   					          x %in% c("1|1", "1/1"), 1L,
					          is.na(x), NA_integer_) }), .SDcols = samples]

check_genotypes(gt)

# read MesoNeo metadata file
metadata <- fread("data/neo.impute.1000g.sampleInfo_clusterInfo.txt")

afr <- metadata[popId == "YRI", sampleId]
neand <- "AltaiNeandertal"

afr_freq <- gt[, rowMeans(.SD), .SDcols = afr]
neand_freq <- gt[, rowMeans(.SD), .SDcols = neand]

fixed_sites <- abs(afr_freq - neand_freq) == 1

info_gt <- gt[fixed_sites, ]

fwrite(info_gt, file.path(dirname(filename), paste0("info_", basename(filename), ".gz")), sep = "\t")
