library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(parallel)

library(GenomicRanges)
library(plyranges)

read_metadata <- function() {
  raw_info <- read_tsv("data/neo.impute.1000g.sampleInfo_clusterInfo.txt")

  info <-
    raw_info %>%
    filter(region %in% c("SouthernEurope", "WesternEurope", "WesternAsia", "NorthernEurope", "CentralEasternEurope"),
           country != "Greenland")

  info
}

read_tracts <- function(set, metadata) {
  info <- filter(metadata, groupAge == set)

  raw_tracts <- read_tsv(here::here("data/Vindija33.19_raw_eurasian_wModern.gz"))

  tracts <- raw_tracts %>%
    mutate(chrom = paste0("chr", chrom)) %>%
    filter(ID %in% unique(info$sampleId)) %>%
    select(ID, chrom, start, end) %>%
    mutate(set = set)

  tracts
}

generate_windows <- function(gaps_gr, window_size, step_size) {
  autosomes_gr <- GenomeInfoDb::getChromInfoFromUCSC("hg19") %>%
    filter(grepl("chr\\d+$", chrom)) %>%
    mutate(start = 1) %>%
    dplyr::rename(end = size) %>%
    makeGRangesFromDataFrame()
  seqinfo(autosomes_gr) <- seqinfo(tracts_gr)

  windows_grl <- slidingWindows(autosomes_gr, width = window_size, step = step_size)
  for (i in seq_along(windows_grl)) {
    windows_grl[[i]]$midpoint <- (start(windows_grl[[i]]) + end(windows_grl[[i]])) / 2
    windows_grl[[i]]$gap <- FALSE
    to_remove <- queryHits(findOverlaps(windows_grl[[i]], gaps_gr))
    if (length(to_remove) > 0)
      windows_grl[[i]][to_remove]$gap <- TRUE
  }

  unlist(windows_grl)
}

compute_ancestry <- function(tracts_gr, windows_gr) {
  autosomes <- getChromInfoFromUCSC("hg19") %>%
    filter(grepl("chr\\d+$", chrom)) %>%
    mutate(start = 1) %>%
    dplyr::rename(end = size) %>%
    makeGRangesFromDataFrame()

  # first compute coverage...
  cov <- coverage(tracts_gr)
  # ... then convert that to proportions
  cov <- lapply(cov, function(x) x / length(unique(tracts_gr$ID)))

  ancestry_list <- mclapply(as.character(unique(seqnames(windows_gr))),
                            function(chrom) {
    chrom_coverage <- cov[[chrom]]
    chrom_gr <- windows_gr[seqnames(windows_gr) == chrom]

    # count overlaps between windows and tracts
    average_coverage_per_window <- sapply(
      seq_along(chrom_gr),
      function(i) {
        start_idx <- start(chrom_gr[i])
        end_idx <- end(chrom_gr[i])
        mean(chrom_coverage[start_idx:end_idx])
    })

    mcols(chrom_gr)$coverage <- average_coverage_per_window

    chrom_gr
  }, mc.cores = detectCores())

  ancestry_grl <- GRangesList(ancestry_list)

  unlist(ancestry_grl)
}

plot_desert_ancestry <- function(ancestry_gr, deserts_gr, chrom) {
  desert_df <- deserts_gr %>% filter(seqnames == chrom) %>% as_tibble

  ancestry_df <- as_tibble(ancestry_gr) %>%
    filter(seqnames == chrom) %>%
    select(chrom = seqnames, start, end, ancient, modern, gap, midpoint)

  ancestry_df %>%
    filter(start >= (desert_df$start * 0.9) & end <= (desert_df$end * 1.1)) %>%
    {
      ggplot(data = .) +
        geom_rect(data = desert_df, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "red", alpha = 0.1) +
        geom_line(aes(midpoint, ancient), color = "orange") +
        geom_point(data = filter(., ancient > 0), aes(midpoint, ancient, color = "ancient individuals"), size = 0.8) +
        geom_line(aes(midpoint, modern), color = "blue") +
        geom_point(data = filter(., modern > 0), aes(midpoint, modern, color = "present-day individuals"), size = 0.8) +
        geom_vline(data = desert_df, aes(xintercept = start, linetype = "desert boundary"), color = "red") +
        geom_vline(data = desert_df, aes(xintercept = end, linetype = "desert boundary"), color = "red") +
        scale_color_manual(values = c("orange", "blue")) +
        guides(color = guide_legend("", override.aes = list(size = 5)),
               linetype = guide_legend("")) +
        labs(x = "genomic coordinate [bp]", y = "proportion of Neanderthal ancestry") +
        scale_x_continuous(labels = scales::comma) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 0.01)) +
        coord_cartesian(ylim = c(0, 0.1)) +
        scale_linetype_manual(values = "dashed") +
        theme_minimal() +
        theme(legend.position = "bottom") +
        ggtitle(paste("Archaic ancestry desert on chromosome", gsub("chr", "", .$chrom[1])))
    }

}

plot_desert_correlation <- function(ancestry_gr, chrom) {
  ancestry_df <- as_tibble(ancestry_gr) %>% filter(seqnames == chrom)

  ggplot() +
    geom_smooth(data = filter(ancestry_df, within_desert, ancient > 0, modern > 0), aes(ancient, modern, color = within_desert),
                color = "red", fill = "black", method = "lm", linetype = "dashed", linewidth = 0.8, alpha = 0.35) +
    geom_point(data = filter(ancestry_df, !within_desert, ancient > 0, modern > 0), aes(ancient, modern, color = desert, shape = "outside desert"),
               color = "lightgray", alpha = 0.5) +
    geom_point(data = filter(ancestry_df, within_desert, ancient > 0, modern > 0), aes(ancient, modern, color = within_desert, shape = "within desert"), color = "black") +
    geom_abline(slope = 1, linetype = "dashed") +
    geom_hline(aes(color = "modern", yintercept = mean(ancestry_df$modern)), linetype = "dashed", color = "blue") +
    geom_vline(aes(color = "ancient", xintercept = mean(ancestry_df$ancient)), linetype = "dashed", color = "orange") +
    scale_x_log10(breaks = c(0.0001, mean(ancestry_df$ancient), 1), labels = scales::percent_format(accuracy = 0.01), limits = c(0.00001, 1)) +
    scale_y_log10(breaks = c(0.0001, mean(ancestry_df$modern), 1), labels = scales::percent_format(accuracy = 0.01), limits = c(0.00001, 1)) +
    labs(x = "Neanderthal ancestry proportion\nin ancient Eurasians [log scale]",
         y = "Neanderthal ancestry proportion\nin present-day Eurasians [log scale]") +
    coord_fixed() +
    theme_minimal() +
    theme(legend.position = "bottom") +
    guides(shape = guide_legend("window", override.aes = list(alpha = 1, size = 3)),
           linetype = "none") +
    scale_shape_manual(values = c(4, 20))
}
