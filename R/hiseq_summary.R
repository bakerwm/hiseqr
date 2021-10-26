#' Functions for HiSeq summary
#'
#' Summary data in project_dir.
#' Processing data
#' Prepare for plot
#' Plotting
#'
#' @name hiseq_summary


#' @describeIn hiseq_summary
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
hiseq_summary <- function(x, hiseq_type = "auto") {
  # parameters
  pd <- read_hiseq(x)
  if(hiseq_type == "auto") {
    hiseq_type <- pd$hiseq_type # auto
  }
  keys <- list(
    k1 = c("peak", "enrich", "report"), # r1, r2, rx
    k2 = c("trim", "align"), # r1 only
    k3 = c("align", "frip", "lendist"), # r1, r2
    k4 = c("cor") # r2, rx
  )
  if(is_hiseq_dir(x)) {
    x_dirs <- list_hiseq_dir(x, hiseq_type)
    # r1: trim, align, peak, frip, lendist, enrich, report
    # rn: align,peak, frip, lendist, enrich, cor, report
    # rx: peak, enrich, cor, report
    ##
    x1 <- purrr::keep(x_dirs, function(i) is_hiseq_dir(i, "_r1"))
    xn <- purrr::keep(x_dirs, function(i) is_hiseq_dir(i, "_rn"))
    xx <- purrr::keep(x_dirs, function(i) is_hiseq_dir(i, "_rx"))
    # data
    p1 <- read_hiseq_stat(x_dirs, keys[["k1"]], add_tag = TRUE)
    p2 <- read_hiseq_stat(x1, keys[["k2"]], add_tag = TRUE)
    p3 <- read_hiseq_stat(c(x1, xn), keys[["k3"]], add_tag = TRUE)
    p4 <- read_hiseq_stat(c(xn, xx), keys[["k4"]], add_tag = TRUE)
    px <- c(p1, p2, p3, p4)
    # px <- read_hiseq_stat(x_dirs, keys = TRUE, add_tag = TRUE)
    # save to rds
    report_dir <- list_hiseq_file(x, "report_dir")
    x_rds <- file.path(report_dir, "00.project_stat.rds")
    saveRDS(px, x_rds)
    # output
    px
  }
}


