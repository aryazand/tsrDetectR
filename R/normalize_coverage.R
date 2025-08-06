#' Get normalized 5' end coverage for fragments mapped to genome
#'
#' @param pos.cov Rle object for positive strand
#' @param neg.cov Rle object for negative strand
#' @param method method by which to normalize. One of "total" (default) or
#'   "numeric"
#' @param normalization_factor if method is "numeric", then this is the factor
#'   by which coverage values are normalized
#'
#' @returns a list of Rle objects. First item in list represents normalized
#'   coverage of positive strand. Second item of list represents normalized
#'   coverage of negative strand.
#' @export
#'
#' @examples
#' set.seed(100)
#'
#' gr <- GenomicRanges::GRanges(
#'   seqnames = "chr1",
#'   ranges = IRanges::IRanges(
#'     start = sample(1:900, 100),
#'     width = sample(1:100, 100)),
#'   strand = sample(c("+", "-"), 100, replace = TRUE))
#'
#' seqlengths(gr) <- 1000
#' x <- strand_coverage(gr)
#' normalize_coverage(x[[1]], x[[2]])
normalize_coverage <- function(pos.cov, neg.cov,
                               method = "total",
                               normalization_factor = NULL) {

  # check input
  checkmate::assert(
    checkmate::checkClass(pos.cov, "Rle"),
    checkmate::checkClass(neg.cov, "Rle"),
    checkmate::checkChoice(method, choices = c("total", "numeric")),
    checkmate::checkNumeric(normalization_factor, lower = 0, len = 1, null.ok = TRUE),
    checkmate::checkTRUE(normalization_factor != 0),
    checkmate::checkTRUE(method == "numeric" & !is.null(normalization_factor))
  )

  if(method == "total") {
    nf <- sum(pos.cov) + sum(neg.cov)
  } else if (method == "numeric") {
    nf <- normalization_factor
  }

  pos.cov <- pos.cov/nf
  neg.cov <- neg.cov/nf

  return(list(pos.cov, neg.cov))

}
