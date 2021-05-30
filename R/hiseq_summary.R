#' Functions for HiSeq summary
#'
#' Summary data in project_dir.
#' Processing data
#' Prepare for plot
#' Plotting
#'
#' @name hiseq_summary


#' @describeIn summary_hiseq_r1
#'
#' organize data for project, saved in json format
#' saved in: project_dir/report/data/
#'
#' files/png
#' plot.rds
#' config.json
#' ...
#'
#' @param x path to the ATACseq project
#'
#' @export
summary_hiseq <- function(prj_dir, hiseq_type = "r1") {
  # load data
  pd <- read_hiseq(prj_dir)
  if(! is(pd, "list")) {
    warning(paste0("Not hiseq dir: ", prj_dir))
    return(NULL)
  }
  x_dirs <- list_hiseq_dir(prj_dir, hiseq_type)
  #--02.stat+rds
  # trim, align, chrM, peak
  # qc:
  # frip, frag_size, bam_cor, peak_overlap, PCA
  #--02.1 trim stat
  px <- read_hiseq_stat(x_dirs, keys = TRUE)
  #--03 save to file
  f_rds <- file.path(pd$args$report_dir, "00.project_stat.rds")
  saveRDS(px, f_rds)
  #--04 output
  px
}


#' @describeIn check x is hiseq or not
#' @param x string, path to dir
#'
#' @export
is_hiseq_dir <- function(x, hiseq_type = TRUE) {
  sapply(x, function(i) {
    pd     <- read_hiseq(i)
    i_type <- ifelse(is(pd, "list"), pd$hiseq_type, "")
    i_is_hiseq <- ifelse(is(pd, "list"), pd$is_hiseq, FALSE)
    i_is_hiseq <- i_is_hiseq && (! endsWith(i, "config"))
    if(isTRUE(hiseq_type)) {
      i_is_hiseq
    } else if(is(hiseq_type, "character")) {
      endsWith(i_type, hiseq_type) | startsWith(i_type, hiseq_type)
    }
  })
}


#' list_hiseq_file
#'
#' @param x path to the input directories
#' @param keys character, attributes/names in hiseq_dir, default: ["bam"]
#' @export
list_hiseq_file <- function(x, keys = "bam", hiseq_type = "r1") {
  # get project dirs
  p_dirs <- list_hiseq_dir(x, hiseq_type)
  lapply(p_dirs, function(f) {
    pd <- read_hiseq(f)
    # check
    if(is(pd, "list")) {
      if(keys %in% names(pd$args)) {
        pd$args[[keys]]
      }
    }
  }) %>%
    unlist() # !!!! whether or not
}


#' list_hiseq_dir:
#'
#' @param x path to the input directories
#' @param log boolen, whether print the log status, default: FALSE
#'
#' @export
list_hiseq_dir <- function(x, hiseq_type = "r1"){
  #--01.load data
  pd <- read_hiseq(x)
  if(! is(pd, "list")) {
    warning("Not a hiseq dir: x")
    return(NULL)
  }
  #--02.get dirs
  # getattr, hasattr ?
  if(is_hiseq_dir(x, "r1")) {
    x_dirs <- x
    purrr::keep(x_dirs, is_hiseq_dir(x_dirs, hiseq_type))
  } else if(is_hiseq_dir(x, "rn")) {
    x_dirs <- c(x, pd$args$rep_list) # r1 + rn
    purrr::keep(x_dirs, is_hiseq_dir(x_dirs, hiseq_type))
  } else if(is_hiseq_dir(x, "rx")) {
    # report: r1+rn+rx
    if(is_hiseq_dir(x, "atac")) { # ATACseq
      x_dirs <- list.dirs(x, recursive = FALSE)
      x_dirs <- c(x, x_dirs)
      x_dirs <- purrr::keep(x_dirs, is_hiseq_dir(x_dirs, hiseq_type = TRUE))
    } else {
      if(grepl("^chip|^cnt|^cnr", pd$hiseq_type, ignore.case = TRUE)) {
        rn_dirs <- c(pd$args$ip_dir, pd$args$input_dir) # chip/cnr
      } else if(is_hiseq_dir(x, "rnaseq")) {
        rn_dirs <- c(pd$args$wildtype_dir, pd$args$mutant_dir) # RNAseq
      } else {
        rn_dirs <- c()
      }
      # for r1
      r1 <- lapply(rn_dirs, list_hiseq_dir, hiseq_type="r1") %>%
        unlist()
      # for rn
      rn <- rn_dirs
      # for rx
      x_dirs <- c(r1, rn, x)
    }
    # output
    purrr::keep(x_dirs, is_hiseq_dir(x_dirs, hiseq_type))
  } else if(pd$is_hiseq_rp) {
    x
  } else {
    c()
  }
}




