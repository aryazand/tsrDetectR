#' Find TSRs using a peaks based approach
#'
#' @param x An S4Vectors::Rle object
#' @param bg_size The window around each position to assess background
#' @param peak_pat The regex pattern used to identify peaks. This is passed to pracma::findpeaks() function
#' @param height_above_bg The minimum (absolute) height a peak has to have to be recognized as such. This is passed to pracma::findpeaks() function
#' @param min_peak_distance the minimum distance (in indices) peaks have to have to be counted. This is passed to pracma::findpeaks() function
#' @param min_gapwidth the distance between two peaks to merge them into a single peak
#'
#' @returns GRanges object
#' @export
#'
#' @examples
#' x <- Rle(c(rep(0,10),2,0,0,0,0,3,3,6,3,3,0,0,0,0,2,rep(0,10)))
#' findtsr_peaks(x, bg_size = 15, height_above_bg = 2)
findtsr_peaks <- function(x, bg_size = 101,
                          peak_pat = "[+]{1,}([0]{1,}[+]{1,}){0,}[-]{1,}([0]{1,}[-]{1,}){0,}",
                          height_above_bg = 1, min_peak_distance = 1,
                          min_gapwidth = 2) {

  # Calculate background value at each position of the genome
  # background value is subtracted
  bg.rle = S4Vectors::Rle(stats::runmed(as.vector(x), k = bg_size, endrule = "median"))
  x <- x - bg.rle
  x[x < 0] <- 0

  # Find peaks using regular expression
  peaks <- pracma::findpeaks(x = as.vector(x),
                              peakpat = peak_pat,
                              minpeakheight = height_above_bg,
                              minpeakdistance = min_peak_distance,
                              sortstr = F)

  # Format at GRanges
  peaks <- IRanges::IRanges(start = peaks[,3], end = peaks[,4])

  # with the above algorithm the start of the peak is the baseline before the
  # actual peak, and the end of the peak is when it returns to the baseline.
  # Thus the minimum lenght of a peak is 3 nucleotides.Here we will trim those
  # first and last nucleotides that are actually at the baseline
  peaks <- IRanges::narrow(peaks, start = 2, end = -2)

  # Combine peaks within a gap of 2 of each other
  peaks <- IRanges::reduce(peaks, min.gapwidth = min_gapwidth)

  return(peaks)

}
