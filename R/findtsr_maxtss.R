#' Identify fixed windows that are transcription-start regions (TSRs) based on maximum value within each windows
#'
#' @param x an Rle object
#' @param w width of window
#' @param background size of window to calculate background value. This should be significantly larger than w. To call a window as a TSR, the maximum value within a window must be greater than median of background window.
#' @param threshold windows are scored by the maximum value within the window. The score must be greater than threshold value for a window to be considered a TSR
#' @param rel_threshold windows are scored by the maximum value within the window. The percentile rank of the relative (relative to all values within the Rle object) must be greater than rel_threshold for the window to be considered a TSR.
#'
#' @return an IRanges object of all TSR windows, with metadata columns \code{score} (the maximum value within the window) and \code{background.score} (the background value for the window)
#' @export
#'
#' @examples
#'
#' ##################################
#' # CREATE SAMPLE GENOME COVERAGE
#' ##################################
#' chr_length = 1000
#' n_tss = 10
#' tss_positions = sample(chr_length, n_tss)
#' non_tss_positions = seq(chr_length)[-tss_positions]
#'
#' # Assume transcription starts are a poisson process for simplification
#' avg_background_coverage = 3
#' avg_tss_coverage = 21
#' fiveprime_coverage <- rpois(chr_length, lambda = avg_background_coverage)
#' fiveprime_coverage[tss_positions] <- rpois(n_tss, lambda = avg_tss_coverage)
#' fiveprime_coverage <- Rle(fiveprime_coverage)
#'
#' #####################
#' # FIND PEAKS
#' #####################
#' min_peak = 10
#' peaks <- findtsr_maxtss(fiveprime_coverage, w = 21, background = 101, threshold = min_peak)
#' predicted_tss_positions <- start(resize(peaks, width = 1, fix = "center"))
#'
#' plot(fiveprime_coverage, type = "l")
#' abline(h = min_peak, lty = 2, col = "red")
#' points(tss_positions, as.vector(fiveprime_coverage)[tss_positions], col = "red")
#' points(predicted_tss_positions, as.vector(fiveprime_coverage)[predicted_tss_positions], pch = 20, col = "green")
#' print(paste("Sensitivity:", sum(tss_positions %in% predicted_tss_positions)/n_tss))
#' print(paste("Specificity:", sum(!(non_tss_positions %in% predicted_tss_positions))/length(non_tss_positions)))
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

  # Determine which position meet criteria for maxtss
  is_max_overlap <- x == overlap.max
  is_above_background <- x > background.rle
  is_above_threshold <- x >= threshold
  is_above_rel_threshold <- rel_x >= rel_threshold

  maxtss.pos <- which(is_max_overlap &
                      is_above_background &
                      is_above_threshold &
                      is_above_rel_threshold)

  # Create an IRanges of maxtss
  maxtss <- IRanges::IRanges(
    start = maxtss.pos,
    width = 1,
    score = as.vector(x)[maxtss.pos],
    background.score = as.vector(background.rle)[maxtss.pos]
  )
  maxtss <- IRanges::resize(maxtss, width = w, fix = "center")

  return(maxtss)

}
