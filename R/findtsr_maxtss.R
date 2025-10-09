#' Identify fixed windows that are transcription-start regions (TSRs) based on maximum value within each windows
#'
#' @param x an Rle object
#' @param w width of window
#' @param background size of window to calculate background value. This should be significantly larger than w. To call a window as a TSR, the maximum value within a window must be greater than median of background window.
#' @param threshold the minimum value for maximum value within a window to call a TSR
#'
#' @return an IRanges object of all TSR windows
#' @export
#'
#' @examples
#' x <- strand_coverage(cmv_proseq_sample)
#' x <- x[[1]][[1]]
#' findtsr_maxtss(x = x, w = 21, background = 501, threshold = 10)
findtsr_maxtss <- function(x, w, background = 501, threshold = 0, rel_threshold = 0) {

  # Calculate max for each window
  window.max <- S4Vectors::runq(x, k=w, i = w, endrule = "constant")

  # Calculate threshold value
  threshold.rle = S4Vectors::Rle(values = threshold, lengths = length(x))

  # Calculate relative threshold
  rel_x = x/sum(x)
  rel_threshold.rle = S4Vectors::Rle(values = rel_threshold, lengths = length(x))

  # Calculate background value
  background.rle = S4Vectors::Rle(stats::runmed(x, k = background, endrule = "constant"))

  # Calculate max including overlapping window
  w2 = (w-1)*2 + 1
  overlap.max = S4Vectors::runq(x, k=w2, i = w2, endrule = "constant")

  # Create a IRanges of the maxtss TSRs
  maxtss.pos <- which(x == overlap.max &
                      x >= threshold.rle &
                      x > background.rle &
                      rel_x >= rel_threshold.x)

  maxtss <- IRanges::IRanges(
    start = maxtss.pos,
    width = 1,
    score = as.vector(x)[maxtss.pos],
    background.score = as.vector(background.rle)[maxtss.pos]
  )
  maxtss <- IRanges::resize(maxtss, width = w, fix = "center")

  return(maxtss)

}
