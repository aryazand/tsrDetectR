library("GenomicAlignments")
bamfile.4hpi <- BamFile("notebook/4_hpi_flavo_pro-seq_cmv_tb40e.bam")
bamfile.12hpi <- BamFile("notebook/12_hpi_flavo_pro-seq_cmv_tb40e.bam")
bamfiles <- list(bamfile.4hpi, bamfile.12hpi)

x <- purrr::map(bamfiles, readGAlignmentPairs, strandMode = 2)
x <- purrr::map(x, keepSeqlevels, c("KF297339.1"))

gr <- purrr::map(x, granges, use.names = T)
names(gr) <- c("4_hpi", "12_hpi")
gr.normalized_fiveprimes <- purrr::map(gr, normalized_fiveprime_coverage)
peaks <- purrr::map_depth(gr.normalized_fiveprimes, 2, ~map(.x, findtsr_peaks, bg_size = 51, height_above_bg =  5))
peaks <- purrr::transpose(peaks)
peaks$`+` <- purrr::map_depth(peaks$`+`, 2, ~GRanges(seqnames = "KF297339.1", ranges = .x, strand = "+"))
peaks$`-` <- purrr::map_depth(peaks$`-`, 2, ~GRanges(seqnames = "KF297339.1", ranges = .x, strand = "-"))
peaks <- transpose(peaks) |>
  purrr::map(unlist) |>
  purrr::map(GRangesList) |>
  purrr::map(unlist) |>
  purrr::map(unname)

peaks_combined <- peaks |>
  GRangesList() |>
  unlist()

peaks_combined2 <- GenomicRanges::reduce(peaks_combined, min.gapwidth = 10, with.revmap = TRUE)
revmap <- mcols(peaks_combined2)$revmap

mcols(peaks_combined2) <- relist(peaks_combined[unlist(revmap)], revmap) |>
  purrr::map(~mcols(.x) |> data.frame() |> dplyr::summarise_all(sum, na.rm=T)) |>
  dplyr::bind_rows()

peaks_combined2.split <- split(peaks_combined2, strand(peaks_combined2))

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



