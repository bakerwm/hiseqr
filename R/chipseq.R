# Functions for ChIPseq pipeline


#' chipseq_report_single
#'
#' @param input directory to the sample
#' @param output directory to the html file
#' @param template the template, default from hiseqr
#'
#' @export
chipseq_report <- function(input, output) {
  ## input
  input <- normalizePath(input)
  if(! dir.exists(output)) dir.create(output, recursive = TRUE)
  output <- normalizePath(output)
  ## output
  outhtml <- file.path(output, "ChIPseq_report.html")

  if(is_chipseq_single_dir(input)){
    template <- system.file('chipseq', 'chipseq_report_single.Rmd', package = "hiseqr")
  } else if(is_chipseq_merge_dir(input)) {
    template <- system.file('chipseq', 'chipseq_report_merge.Rmd', package = "hiseqr")
  } else if(is_chipseq_multiple_dir(input)) {
    template <- system.file("chipseq", "chipseq_report_multiple.Rmd", package = "hiseqr")
  } else {
    warning("unknown input:")
    stop("Stop report")
  }

  ## copy template to output
  template_to <- file.path(output, basename(template))
  file.copy(template, template_to)

  rmarkdown::render(input       = template_to,
                    output_file = outhtml,
                    params      = list(input_dir = input))
}



#' read ChIPseq directory
#'
#' input dir
#' output files, config
#'
#' @param x string path
#'
#' @export
read_chipseq <- function(x) {
  x_type <- is_hiseq_dir(x)
  if(is.null(x_type)) {
    warning(paste0("Not a ChIPseq dir: ", x))
  } else {
    message(paste0("Reading: ", x_type))
  }

  # read args
  pk1 <- list_arguments_file(x)
  pk1 <- unlist(pk1)
  if(! is.null(pk1)) {
    args <- load_pickle(pk1[1])
  } else {
    args <- list()
  }

  list(chipseq_type = x_type,
       file_list = file_to_list(x, recursive = TRUE),
       args = args)
}



#' @param x string, path to dir
#'
#' @export
is_chipseq_dir <- function(x) {
  x_type <- is_hiseq_dir(x)

  if(! is.null(x_type)) {
    grepl("^chipseq", x_type)
  }
}


#' @param x string, path to dir
#'
#' @export
is_chipseq_single_dir <- function(x) {
  x_type <- is_hiseq_dir(x)
  x_type == "chipseq_r1" # sample:1, replicate:n
}


#' @param x string, path to dir
#'
#' @export
is_chipseq_merge_dir <- function(x) {
  x_type <- is_hiseq_dir(x)
  x_type %in% c("chipseq_rn") # sample:1, replicate:n
}


#' @param x string, path to dir
#'
#' @export
is_chipseq_multiple_dir <- function(x) {
  x_type <- is_hiseq_dir(x)
  x_type %in% c("chipseq_from_design", "chipseq_rx")
}



#' list_chipseq_single_dirs
#'
#' @param x path to the input directories
#' @param log boolen, whether print the log status, default: FALSE
#'
#' @export
list_chipseq_single_dirs <- function(x){
  if(is_chipseq_dir(x)) {
    if(is_chipseq_single_dir(x)) {
      x
    } else if(is_chipseq_merge_dir(x)) {
      rep_list <- read_chipseq(x)$args$rep_list
      sapply(rep_list, list_chipseq_single_dirs, simplify = TRUE)
    } else if(is_chipseq_multiple_dir(x)) {
      px <- read_chipseq(x)
      sapply(c(px$args$ip_dir, px$args$input_dir), list_chipseq_single_dirs, simplify = TRUE)
    }
  }
}


#' chipseq_merge_dirs
#'
#' @param x path to the input directories
#'
#' @export
list_chipseq_merge_dirs <- function(x){
  if(is_chipseq_dir(x)) {
    if(is_chipseq_single_dir(x)) {
      NULL
    } else if(is_chipseq_merge_dir(x)) {
      x
    } else if(is_chipseq_multiple_dir(x)) {
      px <- read_chipseq(x)
      sapply(c(px$args$ip_dir, px$args$input_dir), list_chipseq_merge_dirs, simplify = TRUE)
    }
  }
}


#' chipseq_merge_dirs
#'
#' @param x path to the input directories
#'
#' @export
list_chipseq_multiple_dirs <- function(x){
  if(is_chipseq_dir(x)) {
    if(is_chipseq_multiple_dir(x)) {
      x
    }
  }
}




## read align from chipseq seq
#' @export
get_chipseq_align_stat <- function(x) {
  # search for single dir
  rep_list <- list_chipseq_single_dirs(x)
  lapply(rep_list, function(i){
    px <- read_chipseq(i)
    read_align1(px$args$align_stat)}) %>%
    dplyr::bind_rows()
}



#' @export
get_chipseq_frip_stat <- function(x) {
  # search for single dir
  rep_list <- list_chipseq_single_dirs(x)
  lapply(rep_list, function(i){
    px <- read_chipseq(i)
    read.delim(px$args$frip_txt, sep = "\t") %>%
      dplyr::mutate(id = px$args$smp_name,
                    FRiP = paste0(round(FRiP * 100, 2), "%")) %>%
      dplyr::select(id, total_reads, peak_reads, FRiP)}) %>%
    dplyr::bind_rows()
}


#' @export
get_chipseq_lendist_stat <- function(x) {
  # search for single dir
  rep_list <- list_chipseq_single_dirs(x)
  lapply(rep_list, function(i){
    px <- read_chipseq(i)
    read_frag(px$args$lendist_txt)
  }) %>%
    dplyr::bind_rows()
}


#' @export
get_chipseq_peak_stat <- function(x) {
  # search for single dir
  rep_list   <- list_chipseq_single_dirs(x)
  merge_list <- list_chipseq_merge_dirs(x)
  input_list <- sort(c(rep_list, merge_list))
  lapply(input_list, function(i){
    px <- read_chipseq(i)

    data.frame(
      id = px$args$smp_name,
      count = length(readLines(px$args$peak))
    )
  }) %>%
    dplyr::bind_rows()
}


#' @export
get_chipseq_report <- function(x) {
  # search for single dir
  rep_list      <- list_chipseq_single_dirs(x)
  merge_list    <- list_chipseq_merge_dirs(x)
  multiple_list <- list_chipseq_multiple_dirs(x)
  input_list    <- sort(c(rep_list, merge_list, multiple_list))

  report_list <- sapply(input_list, function(i){
    px <- read_chipseq(i)
    report_dir <- px$args$report_dir
    f_html  <- list.files(report_dir, "*.html", full.names = TRUE)
    if(length(f_html) > 0) {
      f_html[1]
    }
  }, simplify = TRUE, USE.NAMES = FALSE) %>%
    unlist()

  report_list[!is.null(report_list)]
}

