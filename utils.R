library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(parallel)

library(GenomicRanges)
library(plyranges)


# i/o functions -------------------------------------------------------------------------------

read_metadata <- function() {
  raw_info <- read_tsv("data/neo.impute.1000g.sampleInfo_clusterInfo.txt")

  info <-
    raw_info %>%
    filter(region %in% c("SouthernEurope", "WesternEurope", "WesternAsia", "NorthernEurope", "CentralEasternEurope"),
           country != "Greenland") %>%
    filter(groupAge != "Archaic") %>%
    mutate(ageAverage = ifelse(groupAge == "Modern", 0, ageAverage)) %>%
    mutate(coverage = ifelse(groupAge == "Modern", Inf, coverage))

  info
}

read_tracts <- function(set) {
  metadata <- read_metadata()

  info <- filter(metadata, groupAge == set)

  raw_tracts <- read_tsv(here::here("data/Vindija33.19_raw_eurasian_wModern.gz"))

  tracts <- raw_tracts %>%
    mutate(chrom = paste0("chr", chrom)) %>%
    filter(ID %in% unique(info$sampleId)) %>%
    select(ID, chrom, start, end) %>%
    mutate(length = end - start, set = set)

  tracts
}

# archaic deserts -----------------------------------------------------------------------------

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

compute_ancestry <- function(tracts_gr, windows_gr, keep_gaps = FALSE) {
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
    if (!keep_gaps)
      mcols(chrom_gr)$coverage[mcols(chrom_gr)$gap] <- NA

    chrom_gr
  }, mc.cores = detectCores())

  ancestry_grl <- GRangesList(ancestry_list)


  unlist(ancestry_grl)
}

plot_desert_ancestry <- function(ancestry_gr, deserts_gr, chrom, full = FALSE) {
  desert_df <- deserts_gr %>% filter(seqnames == chrom) %>% as_tibble

  ancestry_df <- as_tibble(ancestry_gr) %>%
    filter(seqnames == chrom) %>%
    select(chrom = seqnames, start, end, ancient, modern, gap, midpoint)

  if (!full) {
    ancestry_df <- ancestry_df %>%
      filter(start >= (desert_df$start * 0.9) & end <= (desert_df$end * 1.1))
  }

  p <- ancestry_df %>%
    {
      ggplot(data = .) +
        geom_rect(data = desert_df, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "red", alpha = 0.05) +

        geom_line(aes(midpoint, ancient), color = "blue") +

        geom_line(aes(midpoint, modern), color = "orange") +

        geom_vline(data = desert_df, aes(xintercept = start, linetype = "desert boundary"), color = "red") +
        geom_vline(data = desert_df, aes(xintercept = end, linetype = "desert boundary"), color = "red") +

        geom_hline(yintercept = mean(.$ancient, na.rm = TRUE), color = "blue", linetype = "dashed", alpha = 0.5) +
        geom_hline(yintercept = mean(.$modern, na.rm = TRUE), color = "orange", linetype = "dashed", alpha = 0.5) +

        scale_color_manual(values = c("blue", "orange")) +
        guides(color = guide_legend("", override.aes = list(size = 5)),
               linetype = guide_legend("")) +
        labs(x = "genomic coordinate [bp]", y = "proportion of Neanderthal ancestry") +
        scale_x_continuous(labels = scales::comma) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 0.01)) +
        coord_cartesian(ylim = c(0, 0.5)) +
        scale_linetype_manual(values = "dashed") +
        theme_minimal() +
        theme(legend.position = "bottom", text = element_text(size = 13)) +
        ggtitle(paste("Archaic ancestry desert on chromosome", gsub("chr", "", .$chrom[1])))
    }
    if (!full) {
      p <- p +
        geom_point(data = filter(ancestry_df, ancient > 0), aes(midpoint, ancient, color = "ancient individuals"), size = 0.8) +
        geom_point(data = filter(ancestry_df, modern > 0), aes(midpoint, modern, color = "present-day individuals"), size = 0.8)
    }
    p
}

plot_desert_correlation <- function(ancestry_gr, chrom) {
  ancestry_df <- as_tibble(ancestry_gr) %>% filter(seqnames == chrom)

  rho <- ancestry_gr %>% filter(within_desert, seqnames == chrom) %>% { cor(.$modern, .$ancient) }

  ggplot() +
    geom_smooth(data = filter(ancestry_df, within_desert, ancient > 0, modern > 0), aes(modern, ancient, color = within_desert),
                formula = y ~ x, color = "red", fill = "black", method = "lm", linetype = "dashed", linewidth = 0.8, alpha = 0.35) +

    geom_point(data = filter(ancestry_df, !within_desert, ancient > 0, modern > 0), aes(modern, ancient, color = desert, shape = "outside desert"),
               color = "lightgray", alpha = 0.5) +
    geom_point(data = filter(ancestry_df, within_desert, ancient > 0, modern > 0), aes(modern, ancient, color = within_desert, shape = "within desert"), color = "black") +

    geom_abline(slope = 1, linetype = "dashed") +

    geom_vline(aes(color = "modern", xintercept = mean(ancestry_df$modern, na.rm = TRUE)), linetype = "dashed", color = "blue") +
    geom_hline(aes(color = "ancient", yintercept = mean(ancestry_df$ancient, na.rm = TRUE)), linetype = "dashed", color = "orange") +

    scale_x_log10(breaks = c(0.0001, mean(ancestry_df$modern, na.rm = TRUE), 1), labels = scales::percent_format(accuracy = 0.01), limits = c(0.00001, 1)) +
    scale_y_log10(breaks = c(0.0001, mean(ancestry_df$ancient, na.rm = TRUE), 1), labels = scales::percent_format(accuracy = 0.01), limits = c(0.00001, 1)) +
    labs(x = "Neanderthal ancestry proportion\nin present-day Eurasians [log scale]",
         y = "Neanderthal ancestry proportion\nin ancient Eurasians [log scale]") +
    coord_fixed() +
    theme_minimal() +
    theme(legend.position = "bottom", text = element_text(size = 13)) +
    guides(shape = guide_legend("window", override.aes = list(alpha = 1, size = 3)),
           linetype = "none") +
    scale_shape_manual(values = c(4, 20)) +
    ggtitle("", subtitle = paste0("Pearson correlation within desert = ",
                                  formatC(rho, format = "f", digits = 3)))
}

plot_desert_ancestry2 <- function(ancestry_gr, deserts_gr, chrom) {
  desert_df <- deserts_gr %>% filter(seqnames == chrom) %>% as_tibble

  ancestry_df <- as_tibble(ancestry_gr) %>%
    filter(seqnames == chrom) %>%
    select(chrom = seqnames, start, end, modern, chen, gap, midpoint)

  ancestry_df %>%
    filter(start >= (desert_df$start * 0.9) & end <= (desert_df$end * 1.1)) %>%
    {
      ggplot(data = .) +
        geom_rect(data = desert_df, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "red", alpha = 0.1) +

        geom_line(aes(midpoint, chen), color = "orange") +
        geom_point(data = filter(., chen > 0), aes(midpoint, chen, color = "Chen et al."), size = 0.8) +

        geom_line(aes(midpoint, modern), color = "blue") +
        geom_point(data = filter(., modern > 0), aes(midpoint, modern, color = "Alba et al."), size = 0.8) +

        geom_vline(data = desert_df, aes(xintercept = start, linetype = "desert boundary"), color = "red") +
        geom_vline(data = desert_df, aes(xintercept = end, linetype = "desert boundary"), color = "red") +

        geom_hline(yintercept = mean(.$chen), color = "orange", linetype = "dashed", alpha = 0.5) +
        geom_hline(yintercept = mean(.$modern), color = "blue", linetype = "dashed", alpha = 0.5) +

        scale_color_manual(values = c("Chen et al." = "orange", "Alba et al." = "blue")) +
        guides(color = guide_legend("", override.aes = list(size = 5)), linetype = guide_legend("")) +

        labs(x = "genomic coordinate [bp]", y = "proportion of Neanderthal ancestry") +
        scale_x_continuous(labels = scales::comma) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 0.01)) +
        # coord_cartesian(ylim = c(0, 0.1)) +
        scale_linetype_manual(values = "dashed") +
        theme_minimal() +
        theme(legend.position = "bottom") +
        ggtitle(paste("Archaic ancestry desert on chromosome", gsub("chr", "", .$chrom[1])))
    }

}

# LD-based admixture dating -------------------------------------------------------------------

# Define positions of archaic ancestry informative sites
generate_info_sites <- function(tracts_gr, interval) {
  sites_grl <- lapply(seqlevels(tracts_gr), function(chrom) {
    positions <- seq(from = 1, to = seqlengths(tracts_gr)[chrom], by = interval)

    gr <- GRanges(seqnames = chrom, ranges = IRanges(start = positions, end = positions))
    mcols(gr)$index <- seq_len(length(gr))

    gr
  }) %>% GRangesList()
  seqlevels(sites_grl) <- seqlevels(tracts_gr)
  seqlengths(sites_grl) <- seqlengths(tracts_gr)

  sites_grl
}

# Define list of pairs of sites at given distances
# (one element of the list for each distance bin)
collect_pairs <- function(sites_grl, distances, ncores = parallel::detectCores()) {

  chroms <- sapply(sites_grl, function(x) as.character(unique(seqnames(x))))

  chr_pairs <- lapply(chroms, function(chrom) {

    sites_gr <- sites_grl[chroms == chrom, ] %>% unlist

    pairs <- parallel::mclapply(distances, function(distance) {

      pair1 <- c()
      pair2 <- c()

      # iterate through each site one by one...
      for (i in sites_gr$index) {
        index1 <- i
        # ... and find the index of the first site that is at a given distance
        index2 <- sites_gr[start(sites_gr) >= start(sites_gr[i]) + distance]$index[1]

        if (is.na(index2)) {
          if (seqlengths(sites_gr)[chrom] < start(sites_gr[i]) + distance  + distance / 10)
            break
          else
            next
        }

        # otherwise record the indices of the pair of sites and proceed with searching
        # for the next pair
        pair1 <- c(pair1, index1)
        pair2 <- c(pair2, index2)
      }

      list(pair1 = pair1, pair2 = pair2)

    }, mc.cores = ncores)

    pairs

  })

  names(chr_pairs) <- chroms
  chr_pairs
}

# Compute covariances of allele states at pairs of sites
compute_tract_covariances <- function(tracts_gr, sites_grl, pairs) {
  lapply(seqlevels(sites_grl), function(chrom) {

    chrom_sites_gr <- sites_grl[seqlevels(sites_grl) == chrom, ] %>% unlist

    parallel::mclapply(unique(tracts_gr$name), function(name) {

      ind_tracts_gr <- tracts_gr %>% filter(name == !!name, seqnames == chrom)
      ind_sites_gr <- chrom_sites_gr

      # mark sites falling within an introgressed tract
      tract_overlaps <- queryHits(findOverlaps(ind_sites_gr, ind_tracts_gr))
      mcols(ind_sites_gr)$neand <- FALSE
      mcols(ind_sites_gr[tract_overlaps])$neand <- TRUE
      mcols(ind_sites_gr)$neand <- as.integer(mcols(ind_sites_gr)$neand)

      covariances <- sapply(seq_along(distances), function(i) {
        sites1 <- ind_sites_gr[pairs[[chrom]][[i]]$pair1]$neand
        sites2 <- ind_sites_gr[pairs[[chrom]][[i]]$pair2]$neand
        cov(sites1, sites2)
      })

      tibble(
        chrom = chrom,
        name = name,
        sample_age = unique(ind_tracts_gr$sample_age),
        distance = distances,
        covariance = covariances
      )
    }, mc.cores = detectCores()) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
}


# Compute covariances of allele states at pairs of sites
compute_match_covariances <- function(info_gt, pairs, metadata) {
  archaic_name <- "NEA_1"

  lapply(unique(info_gt$chrom), function(chrom) {

    chrom_info_gt <- info_gt[, .SD[chrom %in% ..chrom, ], .SDcols = !c("chrom", "pos")]

    parallel::mclapply(colnames(chrom_info_gt), function(name) {

      ind_matches <- chrom_info_gt[, .(match = .SD[, get(name) == .SD[, get(archaic_name)]])]

      covariances <- sapply(seq_along(distances), function(i) {
        sites1 <- ind_matches[pairs[[chrom]][[i]]$pair1]$match
        sites2 <- ind_matches[pairs[[chrom]][[i]]$pair2]$match
        cov(sites1, sites2)
      })

      tibble(
        chrom = chrom,
        name = name,
        sample_age = filter(metadata, name == !!name)$sample_age,
        distance = distances,
        covariance = covariances
      )
    }, mc.cores = detectCores()) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
}

fit_exponential <- function(cov_df) {
  grid_df <- expand_grid(chrom = unique(cov_df$chrom), name = unique(cov_df$name))

  fit_df <- lapply(1:nrow(grid_df), function(i) {
    name <- grid_df[i, ]$name
    chrom <- grid_df[i, ]$chrom

    data_df <- filter(cov_df, name == !!name, chrom == !!chrom)
    lambda <- tryCatch({
      nls_res <- nls(covariance ~ SSasymp(distance, Asym, R0, lrc), data = data_df)
      exp(coef(nls_res)["lrc"])
    }, error = function(e) NA
    )

    if (!is.na(lambda)) {
      df <- tibble(
        distance = data_df$distance,
        covariance = predict(nls_res, newdata = data_df[, "distance"])
      )
    } else
      df <- NULL

    r <- 1e-8

    tibble(
      name = name,
      chrom = chrom,
      sample_age = data_df$sample_age[1],
      lambda = lambda,
      t_gens_before = lambda / r,
      t_admix = t_gens_before * gen_time + sample_age,
      fit = list(df)
    )
  }) %>%
    do.call(rbind, .) %>%
    unnest(fit)
}
