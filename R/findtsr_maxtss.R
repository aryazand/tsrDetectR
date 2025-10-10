#' Identify fixed windows that are transcription-start regions (TSRs) based on maximum value within each windows
#'
#' @param x an Rle object
#' @param w width of window
#' @param background size of window to calculate background value. This should be significantly larger than w. To call a window as a TSR, the maximum value within a window must be greater than median of background window.
#' @param threshold windows are scored by the maximum value within the window. The score must be greater than threshold value for a window to be considered a TSR
#' @param rel_threshold windows are scored by the maximum value within the window. The percentile rank of the relative (relative to all values within the Rle object) must be greater than rel_threshold for the window to be considered a TSR.
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

  # Calculate relative threshold
  rel_x = S4Vectors::Rle((dplyr::cume_dist(as.vector(x))))

  # Calculate background value
  background.rle = S4Vectors::Rle(stats::runmed(as.vector(x), k = background, endrule = "constant"))

  # Calculate max including overlapping window
  w2 = (w-1)*2 + 1
  overlap.max = S4Vectors::runq(x, k=w2, i = w2, endrule = "constant")

  # Create a IRanges of the maxtss TSRs
  maxtss.pos <- which(x == overlap.max &
                      x > background.rle &
                      x >= threshold &
                      rel_x >= rel_threshold)

  maxtss <- IRanges::IRanges(
    start = maxtss.pos,
    width = 1,
    score = as.vector(x)[maxtss.pos],
    background.score = as.vector(background.rle)[maxtss.pos]
  )
  maxtss <- IRanges::resize(maxtss, width = w, fix = "center")

  return(maxtss)

}
