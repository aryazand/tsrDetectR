#' Get normalized 5' end coverage for fragments mapped to genome
#'
#' @param gr GRanges object with seqlength property included
#'
#' @returns a list of Rle objects. First item in list represents normalized
#'   coverage of positive strand. Second item of list represents normalized
#'   coverage of negative strand.
#' @export
#'
#' @description at each position the normalized 5' coverage is calculated as
#'   follows:
#'   \deqn{x^2/y}
#'   where \eqn{x} is the 5' coverage and \eqn{y} is the total fragment coverage
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
#' normalized_fiveprime_coverage(gr)
normalized_fiveprime_coverage <- function(gr) {

  # gr must have seqlengths to properly calculate coverage
  stopifnot(!is.null(GenomeInfoDb::seqlengths(gr)))
  gr.seqinfo <-  GenomeInfoDb::seqinfo(gr)

  # Get coverage from fragments in a strand specific manner
  cov.pos <- GenomicRanges::coverage(gr[strand(gr) == "+"])
  cov.neg <- GenomicRanges::coverage(gr[strand(gr) == "-"])

  # Create GRange of 5' ends
  gr.fiveprime <- GenomicRanges::GRanges(
    seqnames = GenomeInfoDb::seqnames(gr),
    ranges = IRanges::IRanges(
      start = ifelse(strand(gr) == "+", start(gr), end(gr)),
      width = 1
    ),
    strand = BiocGenerics::strand(gr)
  )

  GenomeInfoDb::seqinfo(gr.fiveprime) <- gr.seqinfo

  # Get
  fiveprime_cov.pos <- GenomicRanges::coverage(gr.fiveprime[BiocGenerics::strand(gr.fiveprime) == "+"])
  fiveprime_cov.neg <- GenomicRanges::coverage(gr.fiveprime[BiocGenerics::strand(gr.fiveprime) == "-"])

  normalized_cov.pos <- purrr::map2(fiveprime_cov.pos, cov.pos, ~(.x^2)/(.y+0.1))
  normalized_cov.neg <- purrr::map2(fiveprime_cov.neg, cov.neg, ~(.x^2)/(.y+0.1))

  GenomeInfoDb::seqinfo(fiveprime_cov.pos) <- gr.seqinfo
  GenomeInfoDb::seqinfo(fiveprime_cov.neg) <- gr.seqinfo

  normalized_coverage <- list(`+` = IRanges::RleList(normalized_cov.pos), `-` = IRanges::RleList(normalized_cov.neg))

  return(normalized_coverage)

}
