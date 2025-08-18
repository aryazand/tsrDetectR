#' Get normalization factor for spikein controls
#'
#' @param grlist a named GRangesList object representing the mapped reads for spike-in 
#'
#' @returns data.frame
#'
#' @export
#' @examples
#' 
#' n = 1000
#' chrLength = 1000
#' chrNames = c("chr1", "chr2")
#' mygenome = Seqinfo(
#'   seqnames = chrNames, 
#'   seqlengths = rep(chrLength, length(chrNames))) 
#' w = 10
#'
#' gr1 <- GRanges(sample(chrNames,n, replace = TRUE), 
#' IRanges(sample(seq(chrLength-w), n, replace = TRUE), width = w), 
#' strand=sample(c("+", "-"), n, replace = TRUE),
#' seqinfo = mygenome)
#'
#' gr2 <- GRanges(sample(chrNames,n, replace = TRUE), 
#' IRanges(sample(seq(chrLength-w), n, replace = TRUE), width = w), 
#' strand=sample(c("+", "-"), n, replace = TRUE),
#' seqinfo = mygenome)
#'
#' grlist <- GRangesList(gr1, gr2)
#' names(grlist) <- c("sample_1", "sample_2")
#'
#' spikein_normalization_factor(grlist)
spikein_normalization_factor <- function(grlist) {

  checkmate::assert(
    checkmate::checkClass(grlist, "GRangesList"),
    checkmate::checkNames(names(grlist))
  )

  # Get Average Read Depth
  read_depth <- sapply(grlist, length)
  df <- data.frame(
    sample = names(read_depth),
    read_depth = read_depth,
    NormalizationFactor = mean(read_depth)/read_depth
  ) 

  return(df)

}

#' Count coverage within genomic tiles 
#'
#' @param grlist a named GRangesList object representing the mapped reads for spike-in 
#' @param tilewidth a positive integer giving the size of genomic tiles
#' @param average.method how to assess the average coverage within each genomic tile. Choice of "median" or "mean".  Or "none" if not averaging should be done.  
#'
#' @returns data.frame
#'
#' @export
#' @examples
#' 
#' n = 1000
#' chrLength = 1000
#' chrNames = c("chr1", "chr2")
#' mygenome = Seqinfo(
#'   seqnames = chrNames, 
#'   seqlengths = rep(chrLength, length(chrNames))) 
#' w = 10
#'
#' gr1 <- GRanges(sample(chrNames,n, replace = TRUE), 
#' IRanges(sample(seq(chrLength-w), n, replace = TRUE), width = w), 
#' strand=sample(c("+", "-"), n, replace = TRUE),
#' seqinfo = mygenome)
#'
#' gr2 <- GRanges(sample(chrNames,n, replace = TRUE), 
#' IRanges(sample(seq(chrLength-w), n, replace = TRUE), width = w), 
#' strand=sample(c("+", "-"), n, replace = TRUE),
#' seqinfo = mygenome)
#'
#' grlist <- GRangesList(gr1, gr2)
#' names(grlist) <- c("sample_1", "sample_2")
#'
#' counts_in_genomicTiles(grlist)
counts_in_genomicTiles <- function(grlist, tilewidth = 50, average.method = "median") {
  
  checkmate::assert(
    checkmate::checkClass(grlist, "GRangesList"),
    checkmate::checkNames(names(grlist)),
    checkmate::checkInt(tilewidth),
    checkmate::checkChoice(average.method, c("none", "median", "mean")),
    combine = "and"
  )

  ####################################
  # Prepare GRangesList for analysis
  ####################################\

  # Convert to tibble  
  gr.df <- grlist |> 
    purrr::map(list) |> 
    tibble::as_tibble() |> 
    tidyr::pivot_longer(tidyselect::everything(), names_to = "sample", values_to = "GRanges")

  ####################################
  # Strand specific tiles
  ####################################

  tiles <- GenomicRanges::tileGenome(
    seqlengths = GenomeInfoDb::seqlengths(grlist[[1]]), 
    tilewidth = tilewidth
  )
  tiles.pos <- tiles
  strand(tiles.pos) <- "+"
  tiles.neg <- tiles
  strand(tiles.neg) <- "-"
  tiles <- c(tiles.pos, tiles.neg)
  
  ####################################
  # Count number of 5' ends within genomic tiles
  ####################################

  tileCounts <- gr.df$GRanges |> 
    purrr::map(~GenomicRanges::countOverlaps(tiles, .x, maxgap = 0, minoverlap = 1)) |> 
    dplyr::bind_cols() |> 
    stats::setNames(gr.df$sample)

  if(average.method != "none") {
    tileCounts <- avg_counts_in_genomicTiles(tileCounts, average.method)
  }

  return(tileCounts)
}

#' Use spike-in control to calculate a regression function to normalize pro-seq data
#'
#' @param grlist a named GRangesList object representing the mapped reads for spike-in. Must 
#' @param tilewidth a positive integer giving the size of genomic tiles
#' @param average.method how to assess the average coverage within each genomic tile. Choice of "median" or "mean".  Or "none" if not averaging should be done.  
#' @param average.min when applying the regression function, this is used to filter out counts that are unlikely to be meaningful (i.e. count of 0)
#'
#' @returns data.frame
#'
#' @export
#' @examples
#' n = 1000
#' chrLength = 1000
#' chrNames = c("chr1", "chr2")
#' mygenome = Seqinfo(
#'   seqnames = chrNames, 
#'   seqlengths = rep(chrLength, length(chrNames))) 
#' w = 10
#'
#' gr1 <- GRanges(sample(chrNames,n, replace = TRUE), 
#'  IRanges(sample(seq(chrLength-w), n, replace = TRUE), width = w), 
#'  strand=sample(c("+", "-"), n, replace = TRUE),
#'  seqinfo = mygenome)
#'
#' gr2 <- GRanges(sample(chrNames,n, replace = TRUE), 
#'  IRanges(sample(seq(chrLength-w), n, replace = TRUE), width = w), 
#'  strand=sample(c("+", "-"), n, replace = TRUE),
#'  seqinfo = mygenome)
#'
#' grlist <- GRangesList(gr1, gr2)
#' names(grlist) <- c("sample_1", "sample_2")
#' grlist <- resize(grlist, width = 1, fix = "start")
#'
#' spikein_regression(grlist)
spikein_regression <- function(grlist, tilewidth = 50, average.method = "median", average.min = 1) {

  checkmate::assert(
    checkmate::checkClass(grlist, "GRangesList"),
    checkmate::checkTRUE(unique(all(width(grlist) == 1))), 
    checkmate::checkTRUE(unique(all(strand(grlist) != "*"))), 
    checkmate::checkNames(names(grlist)),
    checkmate::checkInt(tilewidth),
    checkmate::checkChoice(average.method, c("median", "mean")),
    checkmate::checkNumeric(average.min),
    combine = "and"
  )

  tileCounts <- counts_in_genomicTiles(grlist, average.method = average.method)
  reg_model <- regression_counts_in_genomicTiles(tileCounts, "average")

  return(reg_model)
}

avg_counts_in_genomicTiles <- function(tileCounts, average.method = "median") {
  
  checkmate::assert(
    checkmate::checkClass(tileCounts, "data.frame"),
    checkmate::checkChoice(average.method, c("median", "mean")),
    combine = "and"
  )
  
  # Get central tendency of counts within each tile
  if(average.method == "median") {
    tileCounts <- tileCounts |> 
      dplyr::mutate(average = matrixStats::rowMedians(as.matrix(tileCounts)))
  } else if(average.method == "mean") {
    tileCounts <- tileCounts |> 
      dplyr::mutate(average = matrixStats::rowMeans2(as.matrix(tileCounts)))    
  }

  return(tileCounts)
}

regression_counts_in_genomicTiles <- function(tileCounts, average_col, average.min = 0) {
  
  checkmate::assert(
    checkmate::checkClass(tileCounts, "data.frame"),
    checkmate::checkString(average_col),
    checkmate::checkTRUE(average_col %in% colnames(tileCounts)),
    checkmate::checkNumeric(average.min),
    combine = "and"
  )

  # Apply Regression Function
  reg_models <- tileCounts[tileCounts[[average_col]] > average.min, ] |> 
    tidyr::pivot_longer(cols = !dplyr::matches(average_col), names_to = "sample_name", values_to = "value")
  
  reg_models <- reg_models |> 
    split(reg_models$sample_name) |> 
    lapply(function(.x) stats::lm(average ~ value, data = .x)) |> 
    lapply(list) |> 
    tibble::as_tibble() |> 
    tidyr::pivot_longer(cols = tidyselect::everything(), names_to = "sample_name", values_to = "model")

  return(reg_models)
}

