
# TODO:
# - add description
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
