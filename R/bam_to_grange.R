#' Convert a GAlignments object to GRanges object
#'
#' @param x GAlignments or GAlignmentPairs object
#' @param keep which metadata to keep
#'
#' @returns GenomicRanges object
#' @import GenomicAlignments
#' @export
#'
#' @examples
#' reads <- GenomicAlignments::GAlignments(
#'   names = c("a","b","c","d","e","f","g"),
#'   seqnames = Rle(c(rep(c("chr1", "chr2"), 3), "chr1")),
#'   pos = c(1400, 2700, 3400, 7100, 4000, 3100, 5200),
#'   cigar = c("500M", "100M", "300M", "500M", "300M",
#'             "50M200N50M", "50M150N50M"),
#'   strand = strand(rep("+", 7)),
#'   mapq = sample(seq(42), 7))
#'
#' bam_to_grange(reads, keep = "mapq")
bam_to_grange <- function(x, keep = NULL) {

  checkmate::assert(
    checkmate::checkClass(x, "GAlignments"),
    checkmate::checkClass(x, "GAlignmentPairs"),
    combine = "or"
  )

  gr <- GenomicRanges::granges(x, use.names = TRUE)

  if(!is.null(keep))
    if(inherits(x, "GAlignments")) {
      mcols(gr)[keep] <- x@elementMetadata[keep]
    } else if(inherits(x, "GAlignmentPairs")) {
      mcols(gr)[keep] <- x@first@elementMetadata[keep]
  }

  return(gr)

}
