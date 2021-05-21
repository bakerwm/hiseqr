
## Functions for parsing atac directories

# 1. Summary
# 2. Map reads, reads + peaks
# 3. Mito pct, table + figure
# 4. Fragment length
# 5. TSS enrich
# 6. Correlation between replicates
# 7. Peaks overlap
# 8. IDR
# x <- '~/work/devel_pipeline/hiseq/atac/results/ATACseq/ctl_rep1/'


#' atac_report_single
#'
#' @param input directory to the sample
#' @param output directory to the html file
#' @param template the template, default from hiseqr
#'
#' @export
atac_report <- function(input, output) {
  ## input
  input <- normalizePath(input)
  if(! dir.exists(output)) dir.create(output, recursive = TRUE)
  output <- normalizePath(output)
  ## output
  outhtml <- file.path(output, "ATACseq_report.html")

  if(is_atac_single_dir(input)){
    template <- system.file('atacseq', 'atac_report_single.Rmd',
                            package = "hiseqr")
  } else if(is_atac_merge_dir(input)) {
    template <- system.file('atacseq', 'atac_report_merge.Rmd',
                            package = "hiseqr")
  } else if(is_atac_multiple_dir(input)) {
    template <- system.file("atacseq", "atac_report_multiple.Rmd",
                            package = "hiseqr")
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




#' read atac directory
#'
#' input dir
#' output files, config
#'
#' @param x string path
#'
#' @export
read_atac <- function(x) {
  x_type <- is_hiseq_dir(x)
  if(isTRUE(is_atac_dir(x))) {
    message(paste0("Reading: ", x_type))
  } else {
    warning(paste0("Not a ATACseq dir: ", x))
  }

  # read args
  pk1 <- list_arguments_file(x)
  pk1 <- unlist(pk1)
  if(! is.null(pk1)) {
    args <- load_config(pk1[1])
  } else {
    args <- list()
  }

  list(atacseq_type = x_type,
       file_list = file_to_list(x, recursive = TRUE),
       args = args)
}


#' @param x string, path to dir
#'
#' @export
is_atac_dir <- function(x) {
  x_type <- is_hiseq_dir(x)

  if(! is.null(x_type)) {
    grepl("^atac", x_type)
  } else {
    FALSE
  }
}


#' @param x string, path to dir
#'
#' @export
is_atac_single_dir <- function(x) {
  x_type <- is_hiseq_dir(x)
  x_type %in% c("atacseq_r1") # sample:1
}


#' @param x string, path to dir
#'
#' @export
is_atac_merge_dir <- function(x) {
  x_type <- is_hiseq_dir(x)
  x_type %in% c("atacseq_rn") # replicate:n
}


#' @param x string, path to dir
#'
#' @export
is_atac_multiple_dir <- function(x) {
  x_type <- is_hiseq_dir(x)
  x_type %in% c("atacseq_rx")
}



#' list_atac_single_dirs
#'
#' @param x path to the input directories
#' @param log boolen, whether print the log status, default: FALSE
#'
#' @export
list_atac_single_dirs <- function(x){
  if(is_atac_dir(x)) {
    if(is_atac_single_dir(x)) {
      x
    } else if(is_atac_merge_dir(x)) {
      rep_list <- read_atac(x)$args$rep_list
      rep_list <- unique(rep_list) # in case, dupilcates
      sapply(rep_list, list_atac_single_dirs, simplify = TRUE)
    } else if(is_atac_multiple_dir(x)) {
      m_dirs <- list_atac_merge_dirs(x)
      s_dirs <- lapply(m_dirs, list_atac_single_dirs)
      unlist(s_dirs)
    }
  }
}


#' atacseq_merge_dirs
#'
#' @param x path to the input directories
#'
#' @export
list_atac_merge_dirs <- function(x){
  if(is_atac_dir(x)) {
    if(is_atac_single_dir(x)) {
      NULL
    } else if(is_atac_merge_dir(x)) {
      x
    } else if(is_atac_multiple_dir(x)) {
      px      <- read_atac(x)
      px_dirs <- list.dirs(x, full.names = TRUE, recursive = FALSE)
      m_dirs  <- lapply(px_dirs, list_atac_merge_dirs)
      m_dirs  <- unlist(m_dirs)
      m_dirs[file.exists(m_dirs)]
    }
  }
}


#' atacseq_merge_dirs
#'
#' @param x path to the input directories
#'
#' @export
list_atac_multiple_dirs <- function(x){
  if(is_atac_dir(x)) {
    if(is_atac_multiple_dir(x)) {
      x
    }
  }
}


## read align from atacseq seq
#' @export
get_atac_align_stat <- function(x) {
  # search for single dir
  rep_list <- list_atac_single_dirs(x)
  lapply(rep_list, function(i){
    px <- read_atac(i)
    read_align1(px$args$align_stat)}) %>%
    dplyr::bind_rows()
}


#' @export
get_atac_frip_stat <- function(x) {
  # search for single dir
  rep_list <- list_atac_single_dirs(x)
  lapply(rep_list, function(i){
    px <- read_atac(i)
    read.delim(px$args$frip_txt, sep = "\t") %>%
      dplyr::mutate(id = px$args$smp_name,
                    FRiP = paste0(round(FRiP * 100, 2), "%")) %>%
      dplyr::select(id, total_reads, peak_reads, FRiP)}) %>%
    dplyr::bind_rows()
}


#' @export
get_atac_lendist_stat <- function(x) {
  # search for single dir
  rep_list <- list_atac_single_dirs(x)
  lapply(rep_list, function(i){
    px <- read_atac(i)
    read_frag(px$args$lendist_txt)
  }) %>%
    dplyr::bind_rows()
}


#' @export
get_atac_peak_stat <- function(x) {
  # search for single dir
  rep_list   <- list_atac_single_dirs(x)
  merge_list <- list_atac_merge_dirs(x)
  input_list <- sort(c(rep_list, merge_list))
  lapply(input_list, function(i){
    px <- read_atac(i)
    peak_file <- px$args$peak
    if(file.exists(px$args$peak)) {
      npeak <- length(readLines(px$args$peak))
    } else {
      npeak <- 0
    }

    data.frame(
      id = px$args$smp_name,
      count = npeak
    )
  }) %>%
    dplyr::bind_rows()
}


#' @export
get_atac_report <- function(x) {
  # search for single dir
  rep_list      <- list_atac_single_dirs(x)
  merge_list    <- list_atac_merge_dirs(x)
  # multiple_list <- list_atac_multiple_dirs(x)
  input_list    <- sort(c(rep_list, merge_list)) #, multiple_list))

  report_list <- lapply(input_list, function(i){
    px <- read_atac(i)
    report_dir <- px$args$report_dir
    f_html  <- list.files(report_dir, "*.html", full.names = TRUE)
    if(length(f_html) > 0) {
      f_html[1]
    }
  })

  report_list <- unlist(report_list)
}


#' to-do
#'
#' TSS
#' IDR
#' cor: counts.tab
#'


