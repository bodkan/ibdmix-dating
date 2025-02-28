#!/bin/env Rscript

library(data.table)
library(dtplyr)
library(magrittr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Incorrect arguments", call. = FALSE)

filename <- args[1]

check_genotypes <- function(gt) {
  gt[, unique(unlist(lapply(.SD, function(x) unique(x)))), .SDcols = !c("chrom", "pos", "ref", "alt")]
}

# read the complete genotype data
gt <- fread(filename, na.strings = "-")
gt[, chrom := paste0("chr", chrom)]

check_genotypes(gt)

# read MesoNeo metadata file
metadata <- fread("data/neo.impute.1000g.sampleInfo_clusterInfo.txt")

afr <- metadata[popId == "YRI", sampleId]
neand <- "AltaiNeandertal"
afr_neand_gt <- gt[, .SD, , .SDcols = c("chrom", "pos", "ref", "alt", "AA", afr, neand)]

check_genotypes(afr_neand_gt)

# replace het sites with NA
afr_neand_gt[, (c(afr, neand)) := lapply(.SD, function(x) fifelse(x  %in% c("0|1", "1|0", "0/1"), NA_character_, x)),
.SDcols = c(afr, neand)]

check_genotypes(afr_neand_gt)

# replace diploid genotypes with ref/alt states
afr_neand_gt[, (c(afr, neand)) := lapply(.SD, function(x) { fcase(x %in% c("0|0", "0/0"), 0L,
   		                         	                  x %in% c("1|1", "1/1"), 1L,
					                          is.na(x), NA_integer_) }), .SDcols = c(afr, neand)]

check_genotypes(afr_neand_gt)

# filter for African-Neanderthal fixed differences
# (Africans = fixed ancestral, Neanderthals = fixed derived)
ancestral <- afr_neand_gt$AA
afr_anc <- rowMeans(afr_neand_gt[, lapply(.SD, function(x) x == ancestral), .SDcols = afr]) == 1
neand_der <- rowMeans(afr_neand_gt[, lapply(.SD, function(x) x != ancestral), .SDcols = neand]) == 1
fixed_sites <- afr_anc & neand_der

# select genotypes at fixed African-Neanderthal differences in Eurasians
info_gt <- gt[fixed_sites, .SD, .SDcols = !c(afr, neand)][, .SD, .SDcols = !patterns("_rnd")]

check_genotypes(info_gt)

# remove non-phased samples
non_phased <- info_gt[, sapply(.SD, function(x) any(grepl("\\/", x)))]
non_phased <- names(non_phased)[non_phased]
info_gt <- info_gt[, .SD, .SDcols = !non_phased]

check_genotypes(info_gt)

# get all samples to use for further processing
ancient_samples <- metadata[dataSource == "thisStudy", sampleId]
modern_samples <- metadata[dataSource == "1000g", sampleId]
samples <- intersect(colnames(info_gt), c(ancient_samples, modern_samples))
#samples <- setdiff(names(info_gt), c("chrom", "pos", "ref", "alt"))
info_gt <- info_gt[, .SD, .SDcols = c("chrom", "pos", "AA", samples)]

# split diploid genotypes into haploid phased chromosome genotypes and define
# genotype states as "an allele here is Neanderthal-like"
neand_gt <- afr_neand_gt$AltaiNeandertal[fixed_sites] %>% { .[!is.na(.)] }
for (i in seq_along(samples)) {
  ind <- samples[i]
  info_gt[, paste0(ind, "_hap1") := as.integer(neand_gt == sub("\\|.*", "", get(ind)))]
  info_gt[, paste0(ind, "_hap2") := as.integer(neand_gt == sub(".*\\|", "", get(ind)))]
  info_gt[, (ind) := NULL]
}
info_gt[, NEA_1 := neand_gt]

fwrite(info_gt, file.path(dirname(filename), paste0("info_", basename(filename), ".gz")), sep = "\t")
