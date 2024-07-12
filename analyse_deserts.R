args <- commandArgs(trailingOnly = TRUE)
set <- args[1]
if (!set %in% c("Modern", "Ancient"))
  stop("Invalid 'set' value", call. = FALSE)

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(rnaturalearth)

library(GenomicRanges)
library(ggbio)
library(plyranges)


# IBDmix desert coordinates from Table S8 of Chan et al., 202 ---------------------------------
# (https://ars.els-cdn.com/content/image/1-s2.0-S0092867420300593-mmc1.pdf)
# (https://ars.els-cdn.com/content/image/1-s2.0-S0092867420300593-figs4_lrg.jpg)

deserts_df <- tribble(
  ~chrom, ~start_ss, ~end_ss,   ~start_ibdmix, ~end_ibdmix,
  "chr1", 102200000, 114900000, 105400000,     120600000,
  "chr3", 76500000,  90500000,  74100000,      89300000,
  "chr7", 106300000, 124700000, 106200000,     123200000,
  "chr8", 53900000,  66000000,  49400000,      66500000
) %>% pivot_longer(cols = c(starts_with("start_"), starts_with("end_")),
                   names_to = c(".value", "method"), names_sep = "_")

# combined Neanderthal and Denisovan deserts
# (https://www.science.org/doi/suppl/10.1126/science.aad9416/suppl_file/vernot-sm.pdf)
ss2_deserts_gr <- tribble(
  ~chrom, ~start, ~end,
  "chr1", 104000000, 114900000,
  "chr3", 76500000, 90500000,
  "chr7", 113600000, 124700000,
  "chr8", 54500000, 65400000
) %>% makeGRangesFromDataFrame()

ss_deserts_gr <- filter(deserts_df, method == "ss") %>% makeGRangesFromDataFrame()
ibdmix_deserts_gr <- filter(deserts_df, method == "ibdmix") %>% makeGRangesFromDataFrame()
deserts_gr <- pintersect(ss_deserts_gr, ibdmix_deserts_gr)
deserts2_gr <- pintersect(ss2_deserts_gr, ibdmix_deserts_gr)

ss_deserts_gr$width <- width(ss_deserts_gr) / 1e6
ss2_deserts_gr$width <- width(ss2_deserts_gr) / 1e6
ibdmix_deserts_gr$width <- width(ibdmix_deserts_gr) / 1e6
deserts_gr$width <- width(deserts_gr) / 1e6
deserts2_gr$width <- width(deserts2_gr) / 1e6

ss_deserts_gr
ss2_deserts_gr
ibdmix_deserts_gr
deserts_gr
deserts2_gr

desert_coords <- deserts2_gr %>% filter(seqnames == "chr7") %>% { c(start(.), end(.)) }

# load Alba's tracts data and MesoNeo metadata ------------------------------------------------

if (is.na(set)) set <- "Ancient"

# raw_tracts <- read_tsv("EURASIA_tracts_archaics_raw.gz")
raw_tracts <- read_tsv("Vindija33.19_raw_eurasian_wModern.gz")
raw_info <- read_tsv("neo.impute.1000g.sampleInfo_clusterInfo.txt")

info <- raw_info %>%
  filter(groupAge == set) %>%
  filter(region %in% c("SouthernEurope", "WesternEurope", "WesternAsia", "NorthernEurope", "CentralEasternEurope"),
         country != "Greenland")

world <- ne_countries(scale = "medium", returnclass = "sf")
sf::st_agr(world) <- "constant"
bbox <- st_as_sfc(st_bbox(c(xmin = -25, xmax = 65, ymin = 25, ymax = 70), crs = st_crs(world)))
western_eurasia <- st_crop(st_make_valid(world), bbox)

info %>%
  filter(!is.na(latitude) & !is.na(longitude)) %>%
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs(4326) %>%
  ggplot() +
    geom_sf(data = western_eurasia) +
    geom_sf(aes(color = region)) +
    coord_sf(crs = 3035)

tracts <- raw_tracts %>%
  mutate(chrom = paste0("chr", chrom)) %>%
  filter(ID %in% unique(info$sampleId)) %>%
  filter(chrom == "chr7") %>%
  select(ID, chrom, start, end)

# glimpse(info)
# glimpse(tracts)

nrow(info)
length(unique(tracts$ID))

p_tracts <-
  tracts %>%
  ggplot(aes(x = start, xend = end, y = ID, yend = ID)) +
    geom_segment(linewidth = 1) +
    # geom_vline(xintercept = c(106300000, 124700000), color = "red") +
    geom_vline(xintercept = desert_coords, color = "red") +
    xlim(1, max(tracts$end)) +
    labs(x = "position along a chromosome [bp]", y = "each row = tracts in an individual") +
    theme_bw() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank()
    )
p_tracts

ggplot2::ggsave(paste0("tracts_", set, ".pdf"), p_tracts, width = 10, height = 7)



# testing tracts set --------------------------------------------------------------------------

tracts_gr <- GRanges(
  seqnames = "chr7",
  ranges = IRanges(c(2, 4, 12, 4, 12, 5), c(6, 6, 17, 7, 15, 7)),
  ID = c("ind1", "ind2", "ind1", "ind3", "ind3", "ind4")
)
seqlengths(tracts_gr) <- 20

tracts_gr

# ggplot(tracts_gr) +
#   geom_rect(aes(group = ID, fill = ID)) +
#   scale_x_continuous(breaks = seq_len(seqlengths(tracts_gr))) +
#   coord_cartesian(xlim = c(1, seqlengths(tracts_gr))) +
#   theme_minimal() +
#   theme(panel.grid = element_blank())

orig_par <- par(no.readonly = TRUE)

pdf("windows_test.pdf", 10, 10)

par(mfrow = c(4, 1))

plot(NA, xlim = c(1, seqlengths(tracts_gr)), ylim = c(1, length(unique(tracts_gr$ID))), ylab = "individual")
segments(x0 = start(tracts_gr), x1 = end(tracts_gr), y0 = as.numeric(factor(tracts_gr$ID)), y1 = as.numeric(factor(tracts_gr$ID)),
         col = factor(tracts_gr$ID))

cov <- coverage(tracts_gr)
cov <- cov[[1]]
cov <- cov / length(unique(tracts_gr$ID))
plot(seq_along(cov), cov, type = "o", xlim = c(1, seqlengths(tracts_gr)), ylab = "coverage per site", ylim = c(0, 1))

# prop_cov <- cov / length(unique(tracts_gr$ID))
# prop_cov
# plot(prop_cov[[1]], type = "o", xlim = c(1, seqlengths(tracts_gr)))

# runcov <- runmean(cov, 5)
# runcov
# plot(seq_along(runcov), runcov, type = "o", xlim = c(0, seqlengths(tracts_gr)), ylab = "runmean() coverage")

chrom_length <- seqlengths(tracts_gr)
window_size <- 5
step_size <- 3

windows <- slidingWindows(IRanges(start = 1, end = chrom_length), width = window_size, step = step_size)
windows_gr <- GRanges(seqnames = "chr7", ranges = unlist(windows))
windows_gr$id <- factor(seq_len(length(windows_gr)))
seqlengths(windows_gr) <- seqlengths(tracts_gr)
windows_gr

# sanity check by plotting overlapping windows sequentially as tiles
# autoplot(windows_gr, aes(group = id), color = NA) + xlim(1, seqlengths(windows_gr))

gaps_gr <- GRanges(seqnames = "chr7", ranges = IRanges(start = 8, end = 11))

to_remove <- queryHits(findOverlaps(windows_gr, gaps_gr))
windows_gr$gap <- FALSE
windows_gr[to_remove]$gap <- TRUE

plot(NA, xlim = c(1, seqlengths(windows_gr)), ylim = c(1, length(windows_gr)), ylab = "sliding window number")
segments(x0 = start(windows_gr), x1 = end(windows_gr), y0 = as.numeric(windows_gr$id), y1 = as.numeric(windows_gr$id),
         col = windows_gr$gap + 1)

# sanity check by plotting overlapping windows sequentially as tiles
# autoplot(windows_gr, aes(group = id, fill = gap), color = NA) + xlim(1, seqlengths(windows_gr))

# count overlaps between windows and tracts
average_coverage_per_window <- function(windows_gr, cov) {
  sapply(seq_along(windows_gr), function(i) {
    start_idx <- start(windows_gr[i])
    end_idx <- end(windows_gr[i])
    mean(cov[start_idx:end_idx])
  })
}
mcols(windows_gr)$coverage <- average_coverage_per_window(windows_gr, as.numeric(cov))
mcols(windows_gr)$midpoint <- (start(windows_gr) + end(windows_gr)) / 2

plot(windows_gr$midpoint, windows_gr$coverage,
     ylab = "mean coverage in sliding window", type = "o", xlim = c(1, seqlengths(tracts_gr)), ylim = c(0, 1))

par(orig_par)

dev.off()


# empirical tracts ----------------------------------------------------------------------------

library(BSgenome.Hsapiens.UCSC.hg19)

tracts_gr <- GRanges(
  seqnames = tracts$chrom,
  ranges = IRanges(start = tracts$start, end = tracts$end),
  ID = tracts$ID
)
seqlengths(tracts_gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)[unique(as.character(seqnames(tracts_gr)))]
genome(tracts_gr) <- "hg19"

coverage_gr <- coverage(tracts_gr)

cov <- coverage_gr[[1]]
cov <- cov / length(unique(tracts$ID))



# empirical gaps ------------------------------------------------------------------------------

library(ggbio)
library(rtracklayer)
library(plyranges)

mySession <- browserSession()
genome(mySession) <- "hg19"
query <- ucscTableQuery(mySession, table = "gap")

# gap table columns: https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=map&hgta_track=gap&hgta_table=gap&hgta_doSchema=describe+table+schema
gaps <- getTable(query)
gaps <- gaps %>% dplyr::filter(grepl("chr\\d+$", chrom)) %>% as_tibble() %>% filter(chrom == "chr7")

gaps_gr <- makeGRangesFromDataFrame(gaps, starts.in.df.are.0based = TRUE, keep.extra.columns = TRUE, ignore.strand = TRUE)
seqlengths(gaps_gr) <- seqlengths(tracts_gr)
genome(gaps_gr) <- genome(tracts_gr)

gaps_gr %>%
  filter(seqnames %in% c("chr1", "chr7")) %>% #, type %in% c("centromere", "heterochromatin")) %>%
  autoplot(aes(fill = type), color = NA) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom")

gaps_gr %>%
  filter(seqnames %in% c("chr1", "chr7")) %>%  #, type %in% c("centromere", "heterochromatin")) %>%
  reduce() %>%
  autoplot(color = NA) +
  theme_bw() +
  theme(panel.grid = element_blank())




# sliding windows along the chromosome --------------------------------------------------------

chrom_length <- seqlengths(tracts_gr)
window_size <- 200e3
step_size <- 50e3

windows <- slidingWindows(IRanges(start = 1, end = chrom_length), width = window_size, step = step_size)
windows_gr <- GRanges(seqnames = "chr7", ranges = unlist(windows))
windows_gr$id <- factor(seq_len(length(windows_gr)))
seqlengths(windows_gr) <- seqlengths(tracts_gr)
genome(windows_gr) <- genome(tracts_gr)
windows_gr

# sanity check by plotting overlapping windows sequentially as tiles
# autoplot(windows_gr, aes(group = id), color = NA) + xlim(1, seqlengths(windows_gr))

to_remove <- queryHits(findOverlaps(windows_gr, gaps_gr))
windows_gr$gap <- FALSE
windows_gr[to_remove]$gap <- TRUE

plot(NA, xlim = c(1, seqlengths(windows_gr)), ylim = c(1, length(windows_gr)), ylab = "sliding window number")
segments(x0 = start(windows_gr), x1 = end(windows_gr), y0 = as.numeric(windows_gr$id), y1 = as.numeric(windows_gr$id),
         col = windows_gr$gap + 1)

plotting_cutoff <- 2e6

plot(NA, xlim = c(1, plotting_cutoff), ylim = c(1, length(windows_gr[start(windows_gr) < plotting_cutoff])), ylab = "sliding window number")
segments(x0 = start(windows_gr), x1 = end(windows_gr), y0 = as.numeric(windows_gr$id), y1 = as.numeric(windows_gr$id),
         col = windows_gr$gap + 1)

plot(NA, xlim = c(seqlengths(windows_gr) - plotting_cutoff, seqlengths(windows_gr)), ylim = c(length(windows_gr) - length(windows_gr[start(windows_gr) < plotting_cutoff]), length(windows_gr)), ylab = "sliding window number")
segments(x0 = start(windows_gr), x1 = end(windows_gr), y0 = as.numeric(windows_gr$id), y1 = as.numeric(windows_gr$id),
         col = windows_gr$gap + 1)

# sanity check by plotting overlapping windows sequentially as tiles
# autoplot(windows_gr, aes(group = id, fill = gap), color = NA) + xlim(1, seqlengths(windows_gr))

# count overlaps between windows and tracts
average_coverage_per_window <- function(windows_gr, cov) {
  sapply(seq_along(windows_gr), function(i) {
    start_idx <- start(windows_gr[i])
    end_idx <- end(windows_gr[i])
    mean(cov[start_idx:end_idx])
  })
}
mcols(windows_gr)$coverage <- average_coverage_per_window(windows_gr, as.numeric(cov))
mcols(windows_gr)$midpoint <- (start(windows_gr) + end(windows_gr)) / 2

pdf(paste0("windows_", set, ".pdf"), 10, 7)

plot(windows_gr$midpoint, windows_gr$coverage,
     ylab = "mean coverage in sliding window", type = "l", xlim = c(1, seqlengths(tracts_gr)), ylim = c(0, 1))
# abline(v = c(106300000, 124700000), col = "red")
abline(v = desert_coords, col = "red")

dev.off()

pdf(paste0("desert_windows_", set, ".pdf"), 10, 7)

windows_gr %>%
  # filter(start >= 100000000 & end <= 130000000) %>% {
  filter(start >= desert_coords[1] * 0.9 & end <= desert_coords[2] * 1.1) %>% {
plot(.$midpoint, .$coverage, ylab = "mean coverage in sliding window", type = "o", ylim = c(0, 0.3), pch = 20)
# abline(v = c(106300000, 124700000), col = "red")
abline(v = desert_coords, col = "red")
}

dev.off()



saveRDS(windows_gr, paste0("windows_", set, ".rds"))


unlink("Rplots.pdf")


# modern vs ancient windows -------------------------------------------------------------------

if (file.exists("windows_Ancient.rds") && file.exists("windows_Modern.rds")) {

win_ancient <- readRDS("windows_Ancient.rds")
win_ancient$set <- "ancient"
win_modern <- readRDS("windows_Modern.rds")
win_modern$set <- "modern"

win <- rbind(
  dplyr::as_tibble(win_ancient) %>% select(chrom = seqnames, start, end, midpoint, coverage, set, gap),
  dplyr::as_tibble(win_modern) %>% select(chrom = seqnames, start, end, midpoint, coverage, set, gap)
) %>%
  pivot_wider(names_from = "set", values_from = "coverage") %>%
  # mutate(desert = start >= 106300000 & end <= 124700000)
  mutate(desert = start >= desert_coords[1] & end <= desert_coords[2])

p1 <- win %>%
  filter(start >= desert_coords[1] * 0.9 & end <= desert_coords[2] * 1.1) %>%
  {
    ggplot(data = .) +
    geom_line(aes(midpoint, ancient), color = "orange") +
    geom_point(data = filter(., ancient > 0), aes(midpoint, ancient, color = "ancient"), size = 0.8) +
    geom_line(aes(midpoint, modern), color = "blue") +
    geom_point(data = filter(., modern > 0), aes(midpoint, modern, color = "modern"), size = 0.8) +
    geom_vline(xintercept = desert_coords, color = "red", linetype = "dashed") +
    scale_color_manual(values = c("orange", "blue")) +
    guides(color = guide_legend("data set", override.aes = list(size = 5))) +
    labs(x = "genomic coordinate [bp]", y = "proportion of Neanderthal ancestry") +
    scale_x_continuous(labels = scales::comma) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    ggtitle(paste("Archaic ancestry desert on chromosome", gsub("chr", "", win$chrom[1])))
  }

mean(win$modern)
mean(win$ancient)

win %>% filter(desert) %>% { cor(.$ancient, .$modern) }

# desert_win <- win %>% filter(desert)
# res <- lm(modern ~ ancient, data = desert_win)
# plot(desert_win$ancient, desert_win$modern, xlim = c(0, 0.01), ylim = c(0, 0.01))
# abline(a = 0, b = 1)
# abline(res, col = "red")

data_range <- c(win$ancient, win$modern) %>% .[. > 0] %>% range

p2 <- ggplot() +
  geom_point(data = filter(win, !desert, ancient > 0, modern > 0), aes(ancient, modern, color = desert, shape = "outside desert"),
             color = "lightgray", alpha = 0.5) +
  geom_point(data = filter(win, desert, ancient > 0, modern > 0), aes(ancient, modern, color = desert, shape = "within desert"), color = "black") +
  geom_smooth(data = filter(win, desert, ancient > 0, modern > 0), aes(ancient, modern, color = desert),
              color = "black", fill = "black", method = "lm", linetype = "dashed", linewidth = 0.8, alpha = 0.35) +
  geom_abline(slope = 1, linetype = "dashed") +
  geom_hline(aes(color = "modern", yintercept = mean(win$modern)), linetype = "dashed", color = "blue") +
  geom_vline(aes(color = "ancient", xintercept = mean(win$ancient)), linetype = "dashed", color = "orange") +
  scale_x_log10(labels = scales::percent_format(accuracy = 0.01)) +
  scale_y_log10(labels = scales::percent_format(accuracy = 0.01)) +
  labs(x = "Neanderthal ancestry proportion\nin ancient Eurasians [%, log scale]",
       y = "Neanderthal ancestry proportion \ninpresent-day Eurasians [%, log scale]") +
  coord_fixed(xlim = data_range, ylim = data_range) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  guides(shape = guide_legend("window", override.aes = list(alpha = 1, size = 3)),
         linetype = "none") +
  scale_shape_manual(values = c(4, 20)); p2

cowplot::plot_grid(p1, p2, nrow = 1, rel_widths = c(1, 0.7))

win %>%
  filter(desert) %>%
  summarise(
    mean(ancient == 0 & modern > 0),
    mean(ancient > 0 & modern == 0),
    mean((ancient == 0 & modern == 0) | (ancient > 0 & modern > 0))
  ) %>%
  pivot_longer(cols = everything())

win %>% filter(desert) %>% filter(modern == 0) %>% {.$ancient * 100} %>% summary

}
