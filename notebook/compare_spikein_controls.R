library(GenomicAlignments)
library(GenomicRanges)
library(tidyverse)


library("biomartr")
spodoptera_genome_file <- getAssemblyStats(organism = "Spodoptera frugiperda",
                                           db = "refseq", reference = TRUE)
spodoptera_genome <- read_genome(spodoptera_genome_file)
chr_names <- names(spodoptera_genome) |> str_split_i(pattern = " ", i = 1)


bamfiles <- list.files("//wsl.localhost/Ubuntu/home/arya/Projects/hff_tb40e_timeseries/results/alignments/", pattern = "spodoptera.bam$", full.names = T)
bam_data <- purrr::map(bamfiles, readGAlignmentPairs, strandMode = 2)
sample_names <- bamfiles |>
  str_split_i(pattern = "/", i = 11)|>
  str_split_i(pattern = "_", i = 1) |>
  factor(levels = c("4", "12", "24", "48", "72", "mock"))

gr <- map(bam_data[1:5], granges, use.names = T, on.discordant.seqnames="drop")
names(gr) <- sample_names[1:5]
gr <- map(gr, keepSeqlevels, value = chr_names)


binned_genome <- tileGenome(seqlengths = seqlengths(gr[[1]]), tilewidth = 1000)
overlaps <- map(gr, ~countOverlaps(binned_genome, .x)) |>
  bind_cols()

overlaps |> sample_n(1000) |> GGally::ggpairs() + scale_y_log10() + scale_x_log10()
