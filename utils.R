library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(rnaturalearth)

library(GenomicRanges)
library(ggbio)
library(plyranges)

read_metadata <- function() {
  raw_info <- read_tsv("neo.impute.1000g.sampleInfo_clusterInfo.txt")

  info <-
    raw_info %>%
    filter(region %in% c("SouthernEurope", "WesternEurope", "WesternAsia", "NorthernEurope", "CentralEasternEurope"),
           country != "Greenland")

  info
}

read_tracts <- function(set, metadata) {
  info <- filter(metadata, groupAge == set)

  raw_tracts <- read_tsv(here::here("Vindija33.19_raw_eurasian_wModern.gz"))

  tracts <- raw_tracts %>%
    mutate(chrom = paste0("chr", chrom)) %>%
    filter(ID %in% unique(info$sampleId)) %>%
    select(ID, chrom, start, end) %>%
    mutate(set = set)

  tracts
}

generate_windows <- function(tracts_gr, gaps_gr, chrom, window_size, step_size) {
  tracts_gr <- tracts_gr %>% filter(seqnames == chrom)
  chrom_length <- seqlengths(tracts_gr)[chrom]

  windows <- slidingWindows(IRanges(start = 1, end = chrom_length), width = window_size, step = step_size)
  windows_gr <- GRanges(seqnames = chrom, ranges = unlist(windows))
  windows_gr$id <- factor(seq_len(length(windows_gr)))
  seqlengths(windows_gr) <- seqlengths(tracts_gr)[chrom]
  genome(windows_gr) <- genome(tracts_gr)

  to_remove <- queryHits(findOverlaps(windows_gr, gaps_gr))
  windows_gr$gap <- FALSE
  windows_gr[to_remove]$gap <- TRUE

  windows_gr
}

compute_ancestry <- function(windows_gr, cov) {
  chr_coverage <- cov[[unique(seqnames(windows_gr))]]

  # count overlaps between windows and tracts
  average_coverage_per_window <- sapply(
    seq_along(windows_gr),
    function(i) {
      start_idx <- start(windows_gr[i])
      end_idx <- end(windows_gr[i])
      mean(chr_coverage[start_idx:end_idx])
  })

  mcols(windows_gr)$coverage <- average_coverage_per_window
  mcols(windows_gr)$midpoint <- (start(windows_gr) + end(windows_gr)) / 2

  windows_gr
}
