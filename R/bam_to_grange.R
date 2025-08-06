#' Title
#'
#' @param x GAlignmentPairs object
#' @param keep which metadata to keep
#'
#' @returns GenomicRanges object
#' @export
#'
#' @examples
#' x <- system.file("inst", "extdata", "04_hpi_flavo_pro-seq_cmv_tb40e_sampled.bam")
#' bam_to_grange(x, keep = "mapq")
bam_to_grange <- function(x, keep = NULL) {

  gr <- granges(x, use.names = TRUE)

  if("mapq" %in% keep) {
    mcols(gr)$mapq <- x@first@elementMetadata$mapq
  }

  return(gr)

}
