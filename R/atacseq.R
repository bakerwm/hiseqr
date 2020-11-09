
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
    args <- load_pickle(pk1[1])
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


#'
#' #' atac_align_txt
#' #'
#' #' @param x path to the file, align.txt
#' #'
#' #' @export
#' #'
#' atac_files <- function(x, group = "align") {
#'   is_atac_dir(x)
#'
#'   # align.txt
#'   align <- list.files(file.path(x, "align"), "*align.txt", full.names = TRUE)
#'
#'   # peak
#'   peak  <- list.files(file.path(x, "peak"), "*narrowPeak", full.names = TRUE)
#'   peak_xls <- gsub(".narrowPeak$", ".xls", peak)
#'
#'   # length distribution
#'   lendist <- file.path(x, "qc", "length_distribution.txt")
#'
#'   # frip
#'   frip <- file.path(x, "qc", "FRiP.txt")
#'
#'   # bam file, proper paired
#'   bam <- list.files(file.path(x, "bam_files"), "*.bam$", full.names = TRUE)
#'   bam_rmdup  <- bam[grep("rmdup.bam", bam)]
#'   bam_proper <- bam[grep("proper_pair.bam", bam)]
#'   if(length(align) == 0) align = NA
#'   if(length(peak) == 0) peak = NA
#'   if(length(lendist) == 0) lendist = NA
#'   if(length(frip) == 0) frip = NA
#'   if(length(bam_rmdup) == 0) bam_rmdup = NA
#'   if(length(bam_proper) == 0) bam_proper = NA
#'
#'   # output
#'   out <- c(align, peak, peak_xls, frip, lendist, bam_rmdup, bam_proper)
#'   names(out) <- c("align", "peak", "peak_xls", "frip", "length_distribution", "bam_rmdup", "bam_proper")
#'
#'   if(group == "all") {
#'     group <- names(out)
#'   }
#'
#'   return(out[group])
#' }






#' ### is_hiseq_dir
#' #' is_atac_dir
#' #'
#' #' @param x path, the directory of a ATACseq output
#' #' @param log boolen, whether print the log status, default: FALSE
#' #'
#' #' @export
#' #'
#' is_atac_dir <- function(x, log = FALSE){
#'   # main dir
#'   flag0   <- dir.exists(x)
#'
#'   # directory exists
#'   required_dirs <- c("raw_data", "clean_data", "align",
#'                      "bam_files", "peak", "qc", "report")
#'   flag1   <- dir.exists(file.path(x, required_dirs))
#'
#'   # check files
#'   frip    <- file.path(x, "qc", "FRiP.txt")
#'   lendist <- file.path(x, "qc", "length_distribution.txt")
#'   align   <- list.files(file.path(x, "align"), "*.align.txt$", full.names = TRUE)
#'   peak    <- list.files(file.path(x, "peak"), "*.narrowPeak", full.names = TRUE)
#'
#'   align   <- ifelse(length(align) == 0, NA, align)
#'   peak    <- ifelse(length(peak) == 0, NA, peak)
#'   flag2   <- file.exists(c(frip, lendist, align, peak))
#'
#'   # report
#'   df      <- data.frame(id = c("main", required_dirs, "FRiP.txt",
#'                                "length_distribution.txt", "align.txt",
#'                                "narrowPeak"),
#'                         status = c(flag0, flag1, flag2),
#'                         stringsAsFactors = FALSE)
#'   if(isTRUE(log)) {
#'     print(paste(format(df$id, width = 30), ":", df$status, sep = " "))
#'   }
#'
#'   # return
#'   return(all(df$status))
#' }
#'


#' #' atac_dirs
#' #'
#' #' @param x path to the input directories
#' #' @param log boolen, whether print the log status, default: FALSE
#' #'
#' #' @export
#' #'
#' atac_dirs <- function(x, log = FALSE) {
#'   # all directories
#'   d1 <- list.dirs(x, full.names = TRUE, recursive = FALSE)
#'   d2 <- c(x, d1) # include itself
#'
#'   # atac dirs
#'   flag <- sapply(d2, is_atac_dir)
#'   d2 <- d2[flag]
#'
#'   # warn
#'   if (length(d2) == 0){
#'     warning(paste0("ATACseq directories not detected: ", x))
#'   }
#'
#'   # return
#'   return(d2)
#' }







#' is_atac_merge_dir
#'
#' @param x path, the directory of a ATACseq output
#' @param log boolen, whether print the log status, default: FALSE
#'
#' @export
#'
# is_atac_merge_dir <- function(x, log = FALSE){
#   # main dir
#   flag0   <- dir.exists(x)
#
#   # directory exists
#   required_dirs <- c("config", "align", "bam_files", "bw_files", "peak",
#                      "qc", "report")
#   flag1   <- dir.exists(file.path(x, required_dirs))
#
#   # check files
#   smp_list  <- file.path(x, "config", "sample_list.txt")
#   count_tab <- file.path(x, "qc", "cor.bam.counts.tab")
#   frip      <- file.path(x, "qc", "FRiP.txt")
#   idr       <- file.path(x, "qc", "idr.txt")
#   peak      <- list.files(file.path(x, "peak"), "*.narrowPeak",
#                           full.names = TRUE)
#   peak <- ifelse(length(peak) == 0, NA, peak[1])
#   flag2   <- file.exists(c(smp_list, count_tab, frip, idr, peak))
#
#   # report
#   df      <- data.frame(id = c("main", required_dirs, "sample_list",
#                                "count_tab", "frip", "idr", "narrowPeak"),
#                         status = c(flag0, flag1, flag2),
#                         stringsAsFactors = FALSE)
#   if(isTRUE(log)) {
#     print(paste(format(df$id, width = 30), ":", df$status, sep = " "))
#   }
#
#   # return
#   # return(all(df$status))
#   return(sum(df$status) > 11) # allow 2 missing
# }
#




#
# atac_merge_dirs <- function(x, log = FALSE) {
#   # all directories
#   d1 <- list.dirs(x, full.names = TRUE, recursive = FALSE)
#   d2 <- c(x, d1) # include itself
#
#   # atac dirs
#   flag <- sapply(d2, is_atac_merge_dir)
#   d2 <- d2[flag]
#
#   # warn
#   if (length(d2) == 0){
#     warning(paste0("ATACseq directories not detected: ", x))
#   }
#
#   # return
#   return(d2)
# }



#' #' atac_merge_files
#' #'
#' #' @param x path to the file, align.txt
#' #'
#' #' @export
#' #'
#' atac_merge_files <- function(x, group = "align") {
#'   is_atac_merge_dir(x)
#'
#'   # align.txt
#'   smp_list_file <- file.path(x, "config", "sample_list.txt")
#'   smp_list <- readLines(smp_list_file)
#'
#'   # align.txt
#'   align_list   <- sapply(smp_list, function(f) atac_files(f, "align"))
#'   lendist_list <- sapply(smp_list, function(f) atac_files(f, "length_distribution"))
#'   frip_list    <- sapply(smp_list, function(f) atac_files(f, "frip"))
#'   peak_list    <- sapply(smp_list, function(f) atac_files(f, "peak"))
#'   count_tab    <- file.path(x, "qc", "cor.bam.counts.tab")
#'   idr_png      <- file.path(x, "qc", "idr.txt.png")
#'   idr_log      <- file.path(x, "qc", "idr.log")
#'
#'   out <- list("align"     = align_list,
#'               "length_distribution"  = lendist_list,
#'               "frip"      = frip_list,
#'               "peak"      = peak_list,
#'               "count_tab" = count_tab,
#'               "idr_png"   = idr_png,
#'               "idr_log"   = idr_log)
#'
#'   if(group == "all") {
#'     group <- names(out)
#'   }
#'
#'   return(out[group])
#' }






#'
#' ## read align from atac seq
#' #' @export
#' get_atac_align_stat <- function(x) {
#'   x_type <- is_hiseq_dir(x)
#'
#'   # load data
#'   if(isTRUE(is_atac_dir(x))) {
#'     px <- read_atac(x)
#'   } else {
#'     return(NULL)
#'   }
#'
#'   # check, align
#'   if(x_type == "atacseq_s1r1") {
#'     hiseqr::read_align1(px$args$align_stat)
#'   } else if(x_type == "atacseq_s1rn") {
#'     rep_list <- file.path(px$args$outdir, fq_name(px$args$fq1))
#'     dplyr::bind_rows(lapply(rep_list, get_atac_align_stat))
#'   } else if(x_type %in% c("atacseq_from_design", "atacseq_snrn")) {
#'     rep_list <- list_atac_single_dirs(x)
#'     dplyr::bind_rows(lapply(rep_list, get_atac_align_stat))
#'   } else {
#'     warning(paste0("Unknown dir, atac_* expected, ", x_tyep, " got"))
#'   }
#' }
#'
#'
#'
#' #' @export
#' get_atac_frip_stat <- function(x) {
#'   x_type <- is_hiseq_dir(x)
#'
#'   # load data
#'   if(isTRUE(is_atac_dir(x))) {
#'     px <- read_atac(x)
#'   } else {
#'     return(NULL)
#'   }
#'
#'   # check, align
#'   if(x_type == "atacseq_s1r1") {
#'     f <- px$args$frip_txt
#'     read.delim(f, sep = "\t") %>%
#'       dplyr::mutate(id = px$args$smp_name,
#'                     FRiP = as.character(FRiP)) %>%
#'       dplyr::select(id, total_reads, peak_reads, FRiP)
#'   } else if(x_type == "atacseq_s1rn") {
#'     rep_list <- file.path(px$args$outdir, fq_name(px$args$fq1))
#'     dplyr::bind_rows(lapply(rep_list, get_atac_frip_stat))
#'   } else if(x_type %in% c("atacseq_from_design", "atacseq_snrn")) {
#'     rep_list <- list_atac_single_dirs(x)
#'     dplyr::bind_rows(lapply(rep_list, get_atac_frip_stat))
#'   } else {
#'     warning(paste0("Unknown dir, atac_* expected, ", x_tyep, " got"))
#'   }
#' }
#'
#'
#'
#' #' @export
#' get_atac_lendist_stat <- function(x) {
#'   x_type <- is_hiseq_dir(x)
#'
#'   # load data
#'   if(isTRUE(is_atac_dir(x))) {
#'     px <- read_atac(x)
#'   } else {
#'     return(NULL)
#'   }
#'
#'   # check, align
#'   if(x_type == "atacseq_s1r1") {
#'     f <- px$args$lendist_txt
#'     read.delim(f, sep = "\t", col.names = c("length", "count")) %>%
#'       dplyr::mutate(id = px$args$smp_name)
#'   } else if(x_type == "atacseq_s1rn") {
#'     rep_list <- file.path(px$args$outdir, fq_name(px$args$fq1))
#'     dplyr::bind_rows(lapply(rep_list, get_atac_lendist_stat))
#'   } else if(x_type %in% c("atacseq_from_design", "atacseq_snrn")) {
#'     rep_list <- list_atac_single_dirs(x)
#'     dplyr::bind_rows(lapply(rep_list, get_atac_lendist_stat))
#'   } else {
#'     warning(paste0("Unknown dir, atac_* expected, ", x_tyep, " got"))
#'   }
#' }
#'
#'
#' #' @export
#' get_atac_report <- function(x) {
#'   x_type <- is_hiseq_dir(x)
#'
#'   # load data
#'   if(isTRUE(is_atac_dir(x))) {
#'     px <- read_atac(x)
#'   } else {
#'     return(NULL)
#'   }
#'
#'   # check, align
#'   if(x_type == "atacseq_s1r1") {
#'     rpt_dir <- px$args$report_dir
#'     f_html  <- list.files(rpt_dir, "*.html", full.names = TRUE)
#'     if(length(f_html) > 0) {
#'       f_html[1]
#'     }
#'   } else if(x_type == "atacseq_s1rn") {
#'     rpt_dir <- px$args$report_dir
#'     f_html  <- list.files(rpt_dir, "*.html", full.names = TRUE)
#'     if(length(f_html) > 0) {
#'       f_html <- f_html[1]
#'     }
#'     # sub-dir
#'     rep_list <- file.path(px$args$outdir, fq_name(px$args$fq1))
#'     f_html2 <- sapply(rep_list, get_atac_report, simplify = TRUE)
#'     if(length(f_html) > 0) {
#'       f_html <- f_html[1]
#'     }
#'     f_html2 <- unlist(f_html2)
#'     c(f_html, f_html2)
#'   } else if(x_type %in% c("atacseq_from_design", "atacseq_snrn")) {
#'     rpt_dir <- px$args$report_dir
#'     f_html  <- list.files(rpt_dir, "*.html", full.names = TRUE)
#'     if(length(f_html) > 0) {
#'       f_html <- f_html[1]
#'     }
#'     # sub-dir
#'     merge_list <- list_atac_merge_dirs(x)
#'     f_html2 <- sapply(merge_list, get_atac_report, simplify = TRUE)
#'     f_html2 <- unlist(f_html2)
#'     c(f_html, f_html2)
#'   } else {
#'     warning(paste0("Unknown dir, atac_* expected, ", x_tyep, " got"))
#'   }
#' }



#' #' @export
#' get_atac_peak_stat <- function(x) {
#'   x_type <- is_hiseq_dir(x)
#'
#'   # load data
#'   if(isTRUE(is_atac_dir(x))) {
#'     px <- read_atac(x)
#'   } else {
#'     return(NULL)
#'   }
#'
#'   # check, align
#'   if(x_type == "atacseq_s1r1") {
#'     data.frame(
#'       id = px$args$smp_name,
#'       count = length(readLines(px$args$peak))
#'     )
#'   } else if(x_type == "atacseq_s1rn") {
#'     rep_list <- file.path(px$args$outdir, fq_name(px$args$fq1))
#'     dplyr::bind_rows(lapply(rep_list, get_atac_peak_stat))
#'   } else if(x_type %in% c("atacseq_from_design", "atacseq_snrn")) {
#'     rep_list <- list_atac_single_dirs(x)
#'     dplyr::bind_rows(lapply(rep_list, get_atac_peak_stat))
#'   } else {
#'     warning(paste0("Unknown dir, atac_* expected, ", x_tyep, " got"))
#'   }
#' }



#' to-do
#'
#' TSS
#' IDR
#' cor: counts.tab
#'


