#' Apply a calculation through a sliding window across an rle object
#'
#' @description
#' This function converts an Rle object to a vector, applies a sliding window function to it using zoo::rollapply.
#'
#' @param x Input Rle-object
#' @param width width of the sliding window
#' @param FUN the function to apply to the window
#' @param ... any additional arguments to pass to zoo::rollapply
#'
#' @returns Rle object same length as x
#' @export
#'
#' @examples
#' cmv_proseq_sample
#'
#' x <- strand_coverage(cmv_proseq_sample)
#' runfun(x[[1]][[1]], width = 10, FUN = sum, fill = NA)
runfun <- function(x, width, FUN, ...) {

  # TODO: More description on documentation
  # TODO: Improve efficiency by working on Rle object instead of converting to vector

  x.vector <- as.vector(x)
  x.metric <- zoo::rollapply(x.vector, width, FUN, ...)
  S4Vectors::Rle(x.metric)

}
