## code to prepare `cmv_proseq_sample` dataset goes here
cmv_proseq_sample <- rtracklayer::import.bed(
  "data-raw/SRR16301852_cmv.bed",
  genome = "cytomegalovirus",
  seqinfo = Seqinfo(seqnames = "KF297339.1", seqlengths = 237683))

usethis::use_data(cmv_proseq_sample, overwrite = TRUE)
