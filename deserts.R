library(readr)
library(dplyr)
library(ggplot2)

library(GenomicRanges)
library(ggbio)

raw_info <- read_tsv("neo.impute.1000g.sampleInfo_clusterInfo.txt")
raw_tracts <- read_tsv("EURASIA_tracts_archaics_raw.gz")

# info <- read_tsv("neo.impute.1000g.sampleInfo_clusterInfo.txt") %>% filter(dataSource == "thisStudy")
info <- raw_info %>%
  filter(dataSource == "1000g") %>%
  filter(region %in% c("SouthernEurope", "WesternEurope", "NorthernEurope"))

tracts <- raw_tracts %>%
  filter(ID %in% unique(info$sampleId))

glimpse(info)
glimpse(tracts)

tracts <- tracts %>% filter(chrom == 7) %>% select(ID, chrom, start, end)

p_tracts <-
  tracts %>%
  ggplot(aes(x = start, xend = end, y = ID, yend = ID)) +
    geom_segment(linewidth = 1) +
    geom_vline(xintercept = c(106300000, 124700000), color = "red") +
    xlim(1, max(tracts$end)) +
    theme_bw() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank()
    ); p_tracts


# testing tracts set --------------------------------------------------------------------------

tracts_gr <- GRanges(
  seqnames = "chr7",
  ranges = IRanges(c(2, 4, 10, 12), c(6, 6, 17, 15)),
  ID = c("ind1", "ind2", "ind1", "ind3")
)
seqlengths(tracts_gr) <- 20

tracts_gr

ggplot(tracts_gr) +
  geom_rect(aes(group = ID, fill = ID)) +
  scale_x_continuous(breaks = seq_len(seqlengths(tracts_gr))) +
  coord_cartesian(xlim = c(1, seqlengths(tracts_gr))) +
  theme(panel.grid = element_blank())

orig_par <- par(no.readonly = TRUE)

par(mfrow = c(2, 1))

cov <- coverage(tracts_gr)
cov <- cov[[1]]
plot(cov, type = "o", xlim = c(1, seqlengths(tracts_gr)))

# prop_cov <- cov / length(unique(tracts_gr$ID))
# prop_cov
# plot(prop_cov[[1]], type = "o", xlim = c(1, seqlengths(tracts_gr)))

runcov <- runmean(cov, 5)
runcov
plot(as.numeric(runcov[[1]]), type="o", xlim = c(0, seqlengths(tracts_gr)))

par(orig_par)

chrom_length <- seqlengths(tracts_gr)
window_size <- 5
step_size <- 3

windows <- slidingWindows(IRanges(start = 1, end = chrom_length), width = window_size, step = step_size)
windows_gr <- GRanges(seqnames = "chr7", ranges = unlist(windows))
windows_gr$id <- factor(seq_len(length(windows_gr)))
seqlengths(windows_gr) <- 20
windows_gr

# sanity check by plotting overlapping windows sequentially as tiles
autoplot(windows_gr, aes(group = id), color = NA) + xlim(1, seqlengths(windows_gr))

gaps_gr <- GRanges(seqnames = "chr7", ranges = IRanges(start = 8, end = 11))

to_remove <- queryHits(findOverlaps(windows_gr, gaps_gr))
windows_gr$gap <- FALSE
windows_gr[to_remove]$gap <- TRUE

# sanity check by plotting overlapping windows sequentially as tiles
autoplot(windows_gr, aes(group = id, fill = gap), color = NA) + xlim(1, seqlengths(windows_gr))

# count overlaps between windows and tracts
cov <- countOverlaps(windows_gr, tracts_gr)

sliding_window_coverage <- tibble(
  start = start(windows_gr),
  end = end(windows_gr),
  coverage = cov,
  gap = windows_gr$gap
) %>%
  mutate(midpoint = (start + end) / 2,
         coverage = ifelse(gap, NA, coverage))

plot(sliding_window_coverage$midpoint, sliding_window_coverage$coverage,
     type = "o", xlim = c(1, seqlengths(tracts_gr)))



# empirical gaps ------------------------------------------------------------------------------

library(rtracklayer)
library(plyranges)

mySession <- browserSession()
genome(mySession) <- "hg19"
query <- ucscTableQuery(mySession, table = "gap")

# gap table columns: https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=map&hgta_track=gap&hgta_table=gap&hgta_doSchema=describe+table+schema
gaps <- getTable(query)
gaps <- gaps %>% dplyr::filter(grepl("chr\\d+$", chrom)) %>% as_tibble()

gaps_gr <- makeGRangesFromDataFrame(gaps, keep.extra.columns = TRUE, ignore.strand = TRUE)

gaps_gr %>%
  filter(seqnames %in% c("chr1", "chr7"), type %in% c("centromere", "heterochromatin")) %>%
  autoplot(aes(fill = type), color = NA) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom")

gaps_gr %>%
  filter(seqnames %in% c("chr1", "chr7"), type %in% c("centromere", "heterochromatin")) %>%
  reduce() %>%
  autoplot(color = NA) +
  theme_bw() +
  theme(panel.grid = element_blank())

# empirical tracts ----------------------------------------------------------------------------

tracts_gr <- GRanges(
  seqnames = tracts$chrom,
  ranges = IRanges(start = tracts$start, end = tracts$end),
  ID = tracts$ID
)

coverage_gr <- coverage(tracts_gr)

cov <- coverage_gr[[1]]
# prop_cov <- cov / length(unique(tracts$ID))

coverage_df <- tibble(chrom = 7, pos = seq_along(cov), coverage = as.numeric(cov))
coverage_df

# sliding windows along the chromosome --------------------------------------------------------

chrom_length <- max(end(tracts_gr))
window_size <- 200e3
step_size <- 50e3

windows <- slidingWindows(IRanges(start = 1, end = chrom_length), width = window_size, step = step_size)
windows_gr <- GRanges(seqnames = "7", ranges = unlist(windows))

# count overlaps between windows and tracts
overlap_counts <- countOverlaps(windows_gr, tracts_gr)

# compute the proportion of tracts in each window
num_individuals <- length(unique(tracts$ID))
proportions <- overlap_counts / num_individuals

sliding_window_coverage <- tibble(
  start = start(windows_gr),
  end = end(windows_gr),
  coverage = proportions
)

ggplot(sliding_window_coverage, aes(x = start, y = coverage)) +
  geom_line()


# windows along the chromosome ----------------------------------------------------------------

window_size <- 200e3

win_coverage <- coverage_df %>%
  mutate(window = (pos - 1) %/% window_size + 1) %>%
  group_by(chrom, window) %>%
  summarise(
    win_start = min(pos),
    win_end = max(pos),
    coverage = mean(coverage)
  ) %>%
  ungroup() %>%
  select(chrom, win_start, win_end, coverage)

p_coverage <-
  ggplot(win_coverage, aes(x = win_start, y = coverage)) +
    geom_line() +
    geom_vline(xintercept = c(106300000, 124700000), color = "red") +
    expand_limits(y = 0) +
    xlim(1, max(coverage_df$pos))

cowplot::plot_grid(p_tracts, p_coverage, nrow = 2)

