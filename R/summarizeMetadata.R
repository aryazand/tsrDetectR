#' Summarize GRanges metadata based on overlap with another GRange
#'
#' @param gr a granges object
#' @param grouping_gr a granges object
#' @param metadata charater indicating which metadata column in gr to summarize
#' @param func the function to use when summarizing metadata column. function must be able to
#' @param ... additional arguments to pass to func
#'
#' @returns grange with summarized metadata
#' @export
#'
#' @examples
#'
#' # Create dummy data
#' gr1 <- GRanges(
#'   seqnames = c(rep("chr1", 11)), strand = "+",
#'   ranges = IRanges(c(1000, 3000, 3600, 4000, 4000, 5000, 5400,
#'                      2000, 3000, 7000, 7500),
#'                    width = c(500, 500, 300, 500, 900, 500, 500,
#'                              900, 500, 600, 300))
#' )
#
#' gr2a <- resize(gr1, width = 50, fix = "start")
#' mcols(gr2a)$feature = sample(seq(100), length(gr2a))
#' gr2b <- resize(gr1, width = 50, fix = "center")
#' mcols(gr2b)$feature = sample(seq(100), length(gr2b))
#' gr2 <- c(gr2a, gr2b)
#'
#' # Apply function
#' summarizeMetadata(gr2, gr1, metadata = "feature", func = mean)

summarizeMetadata <- function(gr, grouping_gr, metadata, func = mean, ...) {

  overlaps <- GenomicRanges::findOverlaps(grouping_gr, gr)

  # For each range in grouping_gr, apply a function to the metadata of all
  # ranges in gr overlap with the range in the grouping_gr
  x <- seq_along(grouping_gr) |>
    lapply(function(.x) overlaps@from == .x) |>
    lapply(function(.x) overlaps@to[.x]) |>
    lapply(function(.x) mcols(gr)[[metadata]][.x]) |>
    lapply(function(.x) func(.x))

  # Ensure all element lengths of x
  stopifnot(all(lengths(x) == 1))

  # convert x into a simple vector
  mcols(grouping_gr)[[metadata]] <- unlist(x)

  return(grouping_gr)
}
