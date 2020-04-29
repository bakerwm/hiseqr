

## Functions for parsing atac directories

# 1. Summary
# 2. Map reads, reads + peaks
# 3. Mito pct, table + figure
# 4. Fragment length
# 5. TSS enrich
# 6. Correlation between replicates
# 7. Peaks overlap
# 8. IDR


#' is_atac_dir
#'
#' @param x path, the directory of a ATACseq output
#' @param log boolen, whether print the log status, default: FALSE
#'
#' @export
#'
is_atac_dir <- function(x, log = FALSE){
  # main dir
  flag0   <- dir.exists(x)

  # directory exists
  required_dirs <- c("raw_data", "clean_data", "align",
                     "bam_files", "peak", "qc", "report")
  flag1   <- dir.exists(file.path(x, required_dirs))

  # check files
  frip    <- file.path(x, "qc", "FRiP.txt")
  lendist <- file.path(x, "qc", "length_distribution.txt")
  align   <- list.files(file.path(x, "align"), "*.align.txt$", full.names = TRUE)
  peak    <- list.files(file.path(x, "peak"), "*.narrowPeak", full.names = TRUE)

  align   <- ifelse(length(align) == 0, NA, align)
  peak    <- ifelse(length(peak) == 0, NA, peak)
  flag2   <- file.exists(c(frip, lendist, align, peak))

  # report
  df      <- data.frame(id = c("main", required_dirs, "FRiP.txt",
                              "length_distribution.txt", "align.txt",
                              "narrowPeak"),
                       status = c(flag0, flag1, flag2),
                       stringsAsFactors = FALSE)
  if(isTRUE(log)) {
    print(paste(format(df$id, width = 30), ":", df$status, sep = " "))
  }

  # return
  return(all(df$status))
}



#' atac_dirs
#'
#' @param x path to the input directories
#' @param log boolen, whether print the log status, default: FALSE
#'
#' @export
#'
atac_dirs <- function(x, log = FALSE) {
  # all directories
  d1 <- list.dirs(x, full.names = TRUE, recursive = FALSE)
  d2 <- c(x, d1) # include itself

  # atac dirs
  flag <- sapply(d2, is_atac_dir)
  d2 <- d2[flag]

  # warn
  if (length(d2) == 0){
    warning(paste0("ATACseq directories not detected: ", x))
  }

  # return
  return(d2)

}



#' atac_align_txt
#'
#' @param x path to the file, align.txt
#'
#' @export
#'
atac_files <- function(x, group = "align") {
  is_atac_dir(x)

  # align.txt
  align <- list.files(file.path(x, "align"), "*align.txt", full.names = TRUE)

  # peak
  peak  <- list.files(file.path(x, "peak"), "*narrowPeak", full.names = TRUE)
  peak_xls <- gsub(".narrowPeak$", ".xls", peak)

  # length distribution
  lendist <- file.path(x, "qc", "length_distribution.txt")

  # frip
  frip <- file.path(x, "qc", "FRiP.txt")

  # bam file, proper paired
  bam <- list.files(file.path(x, "bam_files"), "*.bam$", full.names = TRUE)
  bam_rmdup  <- bam[grep("rmdup.bam", bam)]
  bam_proper <- bam[grep("proper_pair.bam", bam)]
  if(length(align) == 0) align = NA
  if(length(peak) == 0) peak = NA
  if(length(lendist) == 0) lendist = NA
  if(length(frip) == 0) frip = NA
  if(length(bam_rmdup) == 0) bam_rmdup = NA
  if(length(bam_proper) == 0) bam_proper = NA

  # output
  out <- c(align, peak, peak_xls, frip, lendist, bam_rmdup, bam_proper)
  names(out) <- c("align", "peak", "peak_xls", "frip", "length_distribution", "bam_rmdup", "bam_proper")

  if(group == "all") {
    group <- names(out)
  }

  return(out[group])
}




##----------------------------------------------------------------------------##
## merge replicates

#' is_atac_merge_dir
#'
#' @param x path, the directory of a ATACseq output
#' @param log boolen, whether print the log status, default: FALSE
#'
#' @export
#'
is_atac_merge_dir <- function(x, log = FALSE){
  # main dir
  flag0   <- dir.exists(x)

  # directory exists
  required_dirs <- c("config", "align", "bam_files", "bw_files", "peak",
                     "qc", "report")
  flag1   <- dir.exists(file.path(x, required_dirs))

  # check files
  smp_list  <- file.path(x, "config", "sample_list.txt")
  count_tab <- file.path(x, "qc", "cor.bam.counts.tab")
  frip      <- file.path(x, "qc", "FRiP.txt")
  idr       <- file.path(x, "qc", "idr.txt")
  peak      <- list.files(file.path(x, "peak"), "*.narrowPeak",
                          full.names = TRUE)
  peak <- ifelse(length(peak) == 0, NA, peak[1])
  flag2   <- file.exists(c(smp_list, count_tab, frip, idr, peak))

  # report
  df      <- data.frame(id = c("main", required_dirs, "sample_list",
                               "count_tab", "frip", "idr", "narrowPeak"),
                        status = c(flag0, flag1, flag2),
                        stringsAsFactors = FALSE)
  if(isTRUE(log)) {
    print(paste(format(df$id, width = 30), ":", df$status, sep = " "))
  }

  # return
  # return(all(df$status))
  return(sum(df$status) > 11) # allow 2 missing
}




#' atac_merge_dirs
#'
#' @param x path to the input directories
#' @param log boolen, whether print the log status, default: FALSE
#'
#' @export
#'
atac_merge_dirs <- function(x, log = FALSE) {
  # all directories
  d1 <- list.dirs(x, full.names = TRUE, recursive = FALSE)
  d2 <- c(x, d1) # include itself

  # atac dirs
  flag <- sapply(d2, is_atac_merge_dir)
  d2 <- d2[flag]

  # warn
  if (length(d2) == 0){
    warning(paste0("ATACseq directories not detected: ", x))
  }

  # return
  return(d2)
}



#' atac_merge_files
#'
#' @param x path to the file, align.txt
#'
#' @export
#'
atac_merge_files <- function(x, group = "align") {
  is_atac_merge_dir(x)

  # align.txt
  smp_list_file <- file.path(x, "config", "sample_list.txt")
  smp_list <- readLines(smp_list_file)

  # align.txt
  align_list   <- sapply(smp_list, function(f) atac_files(f, "align"))
  lendist_list <- sapply(smp_list, function(f) atac_files(f, "length_distribution"))
  frip_list    <- sapply(smp_list, function(f) atac_files(f, "frip"))
  peak_list    <- sapply(smp_list, function(f) atac_files(f, "peak"))
  count_tab    <- file.path(x, "qc", "cor.bam.counts.tab")
  idr_png      <- file.path(x, "qc", "idr.txt.png")
  idr_log      <- file.path(x, "qc", "idr.log")

  out <- list("align"     = align_list,
              "length_distribution"  = lendist_list,
              "frip"      = frip_list,
              "peak"      = peak_list,
              "count_tab" = count_tab,
              "idr_png"   = idr_png,
              "idr_log"   = idr_log)

  if(group == "all") {
    group <- names(out)
  }

  return(out[group])
}









#' to-do
#'
#' TSS
#' IDR
#' cor: counts.tab
#'


