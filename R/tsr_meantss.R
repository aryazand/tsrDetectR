#' Identify fixed windows that are transcription-start regions (TSRs) based on sum value within each windows
#'
#' @param x an Rle object
#' @param w width of window
#' @param background size of window to calculate background value. This should be significantly larger than w. To call a window as a TSR, the mean value within a window must be greater than median value of background window.
#' @param threshold the minimum value for sum value of a window to call a TSR
#'
#' @return an IRanges object of all TSR windows
#' @export
#'
#' @examples
#' x <- strand_coverage(cmv_proseq_sample)
#' x <- x[[1]][[1]]
#' tsr_meantss(x = x, w = 21, background = 501, threshold = 10)
tsr_meantss <- function(x, w, background = 501, threshold = 0) {

  # Calculate max for each window
  window.mean <- S4Vectors::runmean(x, k=w, endrule = "constant")

  # Calculate threshold value
  threshold.rle = S4Vectors::Rle(values = threshold, lengths = length(x))

  # Calculate background value
  background.rle = S4Vectors::Rle(stats::runmed(x, k = background, endrule = "constant"))

  # Calculate max including overlapping window
  w2 = (w-1)*2 + 1
  overlap.mean = S4Vectors::runmean(x, k=w2, endrule = "constant")

  # Create a IRanges of the maxtss TSRs
  meantss.pos <- which(window.mean >= overlap.mean &
                       window.mean > threshold.rle &
                       window.mean > background.rle)

  meantss <- IRanges::IRanges(
    start = meantss.pos,
    width = 20,
    score = as.vector(window.mean)[meantss.pos]
  )
  meantss <- IRanges::resize(meantss, width = w, fix = "center") |>
    IRanges::reduce()

  return(meantss)

}
