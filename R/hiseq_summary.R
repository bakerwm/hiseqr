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


