library(GenomicAlignments)
library(GenomicRanges)
library(tidyverse)
bamfiles <- list.files("~/Documents/hff_tb40_timeseries/results/alignments/", pattern = "cmv_tb40e.bam$", full.names = T)
bam_data <- purrr::map(bamfiles, readGAlignmentPairs, strandMode = 2, param=ScanBamParam(what="mapq"))
sample_names <- factor(c("04_hpi", "12_hpi", "24_hpi", "48_hpi", "72_hpi", "mock"))

# bam_data as GRanges
gr <- purrr::map(bam_data[1:5], granges, use.names = T)
names(gr) <- sample_names[1:5]
gr <- GRangesList(gr)
cmv_seqinfo <- Seqinfo(seqnames="KF297339.1",
                       seqlengths=237683,
                       isCircular=FALSE,
                       genome="cmv")
seqinfo(gr) <- cmv_seqinfo


# Normalize 5' ends to coverage
gr.normalized_fiveprimes <- purrr::map(gr, normalized_fiveprime_coverage)
df_normalized_fiveprime <- gr.normalized_fiveprimes |>
  map_depth(.depth = 2, as.list) |>
  as_tibble() |>
  pivot_longer(cols = everything(), names_to = "sample", values_to = "data") |>
  mutate(strand = factor(names(data), levels = c("+", "-"))) |>
  unnest_longer(data, indices_to = "chr") |>
  mutate(chr = factor(chr)) |>
  mutate(sample = factor(sample, levels = levels(sample_names)))

# apply peak finding
df_peaks <- df_normalized_fiveprime |>
  rowwise() |>
  mutate(peaks = list(findtsr_peaks(data, bg_size = 51, height_above_bg =  5))) |>
  mutate(peaks = list(GRanges(seqnames = chr, ranges = peaks, strand = strand))) |>
  group_by(sample)

# Trim low map quality regions
get_average_mapq_by_region <- function(bf, peaks) {
  overlaps <- findOverlaps(bf, peaks)

  avg_mapq <- sapply(seq_along(peaks), function(i) {
    # Get alignments that overlap with region i
    region_aligns <- bf[queryHits(overlaps)[subjectHits(overlaps) == i]]

    if (length(region_aligns) > 0) {
      mapq_values <- region_aligns@first@elementMetadata$mapq
      mapq_values <- mapq_values[!is.na(mapq_values)]
      if (length(mapq_values) > 0) {
        return(mean(mapq_values))
      }
    }
    return(NA)
  })

}

avg_mapq <- map2(c(bam_data[1:5], bam_data[1:5]), df_peaks$peaks, get_average_mapq_by_region)
df_peaks$peaks <- map2(df_peaks$peaks, avg_mapq, function(i,j) {
  i$avg_mapq <- j
  return(i)
})

df_peaks <- df_peaks |>
  rowwise() |>
  mutate(peaks = list(peaks[peaks$avg_mapq > 10]))

# Consolidate peaks
peaks_combined <- df_peaks$peaks |>
  GRangesList() |>
  unlist() |>
  GenomicRanges::reduce(min.gapwidth = 10)

# Apply Views
df <- df_normalized_fiveprime |>
  rowwise() |>
  mutate(peak_views = list(
    Views(data, ranges(peaks_combined[strand(peaks_combined) == as.character(strand) &
                                      seqnames(peaks_combined) == as.character(chr)])))) |>
  mutate(max_pos = list(which.max(peak_views))) |>
  mutate(peak_views2 = list(Views(data, IRanges(start = max_pos - 50, width = 101)))) |>
  mutate(peak_scores = list(peak_views2 |> as_tibble() |> group_by(group) |> mutate(pos = -50:50)))

df2 <- df |>
  ungroup() |>
  select(sample, chr, strand, peak_scores) |>
  unnest(peak_scores)

ggplot(df2) + geom_line(aes(pos, value, group = c(group))) + facet_grid(sample ~ strand, scales = "free")

######################
#
#####################

x <- map2(transpose(gr.normalized_fiveprimes), peaks_combined2.split[1:2], function(cov_by_strand, views_by_strand)
  map(cov_by_strand, function(chromosome) Views(chromosome, views_by_strand)))

max_value <- map_depth(x, 2, ~map(.x, max))
max_pos <- map_depth(x, 2, ~map(.x, which.max))

rtracklayer::export(fiveprime_cov.pos, "notebook/test_pos_norm.bw")
rtracklayer::export(fiveprime_cov.neg*-1, "notebook/test_neg_norm.bw")
rtracklayer::export(peaks, "notebook/peaks.bed")
