## code to prepare `cmv_proseq_hff_12hpi` dataset goes here
cmv_proseq_hff_12hpi <- rtracklayer::import.bed(
  "data-raw/SRR16301856_cmv.bed",
  genome = "cytomegalovirus",
  seqinfo = GenomeInfoDb::Seqinfo(seqnames = "KF297339.1", seqlengths = 237683))

usethis::use_data(cmv_proseq_hff_12hpi, overwrite = TRUE)
