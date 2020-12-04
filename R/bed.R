
##----------------------------------------------------------------------------##
## parse BED files


#' read_bed
#'
#' @param x string path to the peak file
#'
#' @import rtracklayer
#' @export
#'
read_bed <- function(x) {
  rtracklayer::import(x, format = "BED")
}



#' peak.xls
#'
#' @param x string path to the peak xls, output fo MACS2
#' @export
#'
read_peak_xls <- function(x) {
  xlines <- readLines(x, n = 50)
  xlines <- xlines[grep("^#", xlines)] # head lines
  xlines <- xlines[grep(":", xlines)] # records
  xlines <- gsub("^#", "", xlines) # remove #
  as.data.frame(stringr::str_split(xlines, ": ", simplify = TRUE))
}


#' read_narrowpeak
#'
#' @param x string path to the peak file
#'
#' @import rtracklayer
#' @export
#'
read_narrowpeak <- function(x) {
  ext <- c(signalValue = "numeric", pValue = "numeric",
           qValue = "numeric", peak = "integer")
  rtracklayer::import(x, format = "BED", extraCols = ext)
}


#' intersect_bed2
#'
#' @param blist list, bed files
#' @import GenomicRanges
#' @export
intersect_bed2 <- function(blist){
  stopifnot(length(blist) >= 2)
  bed1 <- blist[[1]]
  bed2 <- blist[[2]]
  # bed 2 gr
  gr1 <- read_narrowpeak(bed1)
  gr2 <- read_narrowpeak(bed2)
  gr12 <- GenomicRanges::findOverlaps(gr1, gr2)

  # intersect
  n1 <- paste("a", seq_len(length(gr1) - length(gr12)), sep = "")
  n2 <- paste("b", seq_len(length(gr12)), sep = "")
  n3 <- paste("c", seq_len(length(gr2) - length(gr12)), sep = "")

  ## list
  x <- list(rep1 = c(n1, n2), rep2 = c(n2, n3))
  # out
  return(x)
}


#' intersect_bed3
#'
#' @param blist list, intersect 3 bed files
#' @import GenomicRanges
#' @export
intersect_bed3 <- function(blist){
  stopifnot(length(blist) >= 3)
  bed1  <- blist[[1]]
  bed2  <- blist[[2]]
  bed3  <- blist[[3]]
  # bed to gr
  gr1   <- read_narrowpeak(bed1)
  gr2   <- read_narrowpeak(bed2)
  gr3   <- read_narrowpeak(bed3)
  # intersect
  gr12  <- findOverlaps(gr1, gr2, ignore.strand = TRUE)
  gr13  <- findOverlaps(gr1, gr3, ignore.strand = TRUE)
  gr23  <- findOverlaps(gr2, gr3, ignore.strand = TRUE)
  gr123 <- GenomicRanges::findOverlapsfindOverlaps(gr1[gr12@from],
                                                   gr3,
                                                   ignore.strand = TRUE)
  # numbers
  n12   <- length(gr12)
  n13   <- length(gr13)
  n23   <- length(gr23)
  #
  n123  <- length(gr123)
  n12   <- n12 - n123
  n13   <- n13 - n123
  n23   <- n23 - n123
  n1    <- length(gr1) - n12 - n13 - n123
  n2    <- length(gr2) - n12 - n23 - n123
  n3    <- length(gr3) - n13 - n23 - n123
  # overlap
  out <- c(n1, n2, n3, n12, n13, n23, n123)
  names(out) <- c("n1", "n2", "n3", "n12", "n13", "n23", "n123")
  # out list
  p <- lapply(seq_len(7), function(i){
    paste(letters[i], seq_len(out[i]), sep = "")
  })
  names(p) <- names(out)
  # combine
  x <- list(
    rep1 = c(p$n1, p$n12, p$n13, p$n123),
    rep2 = c(p$n2, p$n12, p$n23, p$n123),
    rep3 = c(p$n3, p$n13, p$n23, p$n123))
  return(x)
}


#' intersect_bed4
#'
#' @param blist list, intersect 4 bed files
#' @import GenomicRanges
#' @export
intersect_bed4 <- function(blist){
  stopifnot(length(blist) >= 4)
  bed1 = blist[[1]]
  bed2 = blist[[2]]
  bed3 = blist[[3]]
  bed4 = blist[[4]]
  # bed to gr
}
