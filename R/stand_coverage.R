
#' Calculate strand specific coverage
#'
#' @description
#' `strand_coverage` takes a GRanges object and returns a list of two Rle objects, the coverage on the positive strand and the coverage on the negative strand.
#'
#' @param x Input GRanges object. Each interval must have either `+` or `-` assigned to its strand property
#'
#' @returns A list of two Rle objects, the coverage on the positive strand and the coverage on the negative strand.
#'
#' @import checkmate
#' @import BiocGenerics
#' @import GenomicRanges
#'
#' @export
#' @examples
#' # load a sample GRanges objects
#' cmv_proseq_sample
#'
#' strand_coverage(cmv_proseq_sample)
strand_coverage <- function(x) {

  # Check input parameters -----------
  checkmate::assert(
    checkmate::checkClass(x, "GRanges"),
    checkmate::checkNumeric(GenomicRanges::width(x), lower = 1),
    checkmate::checkNumeric(BiocGenerics::score(x), lower = 0, any.missing = F),
    checkmate::checkNumber(x@seqinfo@seqlengths)
  )

  if (any(BiocGenerics::strand(x) == "*")) {
    stop("all intervals in x must be strand-specific")
  }

  ###

  # create grange covering each position in genome ----------------------------------

  x.pos <- x[GenomicRanges::strand(x) == "+"] |> GenomicRanges::coverage()
  x.neg <- x[GenomicRanges::strand(x) == "-"] |> GenomicRanges::coverage()

  x <- list(x.pos, x.neg)
  return(x)
}
