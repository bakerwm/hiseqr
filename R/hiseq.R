#' Functions for HiSeq pipeline
#'
#' Reading files, directories, pipeline, ...
#' Processing data
#' Prepare for plot
#' Plotting
#'


#' hiseq_report
#'
#' @param input directory to the sample
#' @param output directory to the html file
#' @param template the template, default from hiseqr
#'
#' @export
hiseq_report <- function(input, output) {
  ## input
  input <- normalizePath(input)
  if(! dir.exists(output)) dir.create(output, recursive = TRUE)
  output <- normalizePath(output)
  ## output
  outhtml <- file.path(output, "HiSeq_report.html")

  ## determine the template
  hiseq_type = get_hiseq_type(input)

  # hiseq (hiseq, atacseq, chipseq, rnaseq, ...)
  if(is_hiseq_dir(input)) {
    if(startsWith(hiseq_type, "hiseq_")) {
      hiseq_dir = "hiseq"
    } else if(startsWith(hiseq_type, "atacseq_")) {
      hiseq_dir = "atacseq"
    } else if(startsWith(hiseq_type, "chipseq_")) {
      hiseq_dir = "chipseq"
    } else if(startsWith(hiseq_type, "rnaseq_")) {
      hiseq_dir = "rnaseq"
    } else {
      hiseq_dir = "hiseq"
    }
  }

  # subtype
  if(is_hiseq_single_dir(input)) {
    hiseq_subtype = "hiseq_report_single.Rmd"
  } else if(is_hiseq_merge_dir(input)) {
    hiseq_subtype = "hiseq_report_merge.Rmd"
  } else if(is_hiseq_multiple_dir(input)) {
    hiseq_subtype = "hiseq_report_multiple.Rmd"
  } else {
    hiseq_subtype = "tmp"
  }

  # template
  template <- system.file(hiseq_dir, hiseq_subtype, package = "hiseqr")
  stopifnot(file.exists(template))

  ## copy template to output
  template_to <- file.path(output, basename(template))
  file.copy(template, template_to, overwrite = TRUE)

  rmarkdown::render(input       = template_to,
                    output_file = outhtml,
                    params      = list(input_dir = input))
}


#' read_hiseq_trim
#'
#' parsing the trimming status
#'
#' @param x path to hiseq, single
#' @import dplyr
#' @import readr
#'
#' @export
read_hiseq_trim <- function(x) {
  if(is_hiseq_single_dir(x)) {
    px <- read_hiseq(x)

    # version-1
    # from: qc/00.trim_summary.json
    j <- px$args$trim_summary_json

    if(file.exists(j)) {
      jsonlite::read_json(j) %>%
        as.data.frame.list() %>%
        dplyr::select(id, input, output, out_pct)
    } else {
      data.frame(
        id = px$args$smp_name,
        input = 0,
        output = 0,
        out_pct = 100
      )
    }
  } else if(is_hiseq_merge_dir(x)) {
    # for replicates
    px     <- read_hiseq(x)
    p_list <- px$args$rep_list # single dirs
    lapply(p_list, read_hiseq_trim) %>%
      dplyr::bind_rows()
  } else if(is_hiseq_multiple_dir(x)) {
    # for replicates
    px      <- read_hiseq(x)
    p1_list <- px$args$ip_dir
    p2_list <- px$args$input_dir
    lapply(c(p1_list, p2_list), read_hiseq_trim) %>%
      dplyr::bind_rows()
  }
}


#' read_hiseq_align
#'
#' @param x path to hiseq, single
#' @import readr
#' @import dplyr
#'
#' @export
read_hiseq_align <- function(x){
  # hiseq_type <- get_hiseq_type(x)

  # if(endsWith(hiseq_type, "_r1")) {
  if(is_hiseq_single_dir(x)) {
    px <- read_hiseq(x)

    # version-1
    # from qc/01.alignment_summary.json
    # version-2
    # from alignment/*.json
    j1 <- px$args$align_summary_json
    j2 <- px$args$align_json

    if(file.exists(j1)) {
      m1 <- jsonlite::read_json(j1)
    } else if(file.exists(j2)) {
      m1 <- jsonlite::read_json(j2)
    } else {
      m1 <- NULL
    }
    df <- as.data.frame.list(m1)

    ## arrange columns
    cols <- c("fqname", "index_name", "total", "map", "unique", "multiple",
              "unmap", "nodup", "chrM", "spikein")

    ## pick columns
    cols <- cols[cols %in% colnames(df)]

    ## version-1
    if(all(c("nodup", "chrM", "spikein") %in% cols)) {
      df %>%
        dplyr::select(all_of(cols)) %>%
        dplyr::mutate(map_pct  = round(map / total * 100, 2),
                      dup_pct  = round((map - nodup) / total * 100, 2),
                      chrm_pct = round(chrM / total * 100, 2))
    } else {
      ## version-2
      df %>%
        dplyr::select(all_of(cols)) %>%
        dplyr::mutate(nodup    = map,
                      chrM     = 0,
                      spikein  = 0,
                      map_pct  = round(map / total * 100, 2),
                      dup_pct  = round((map - nodup) / total * 100, 2),
                      chrm_pct = round(chrM / total * 100, 2))
    }

    # } else if(endsWith(hiseq_type, "_rn")) {
  } else if(is_hiseq_merge_dir(x)) {
    # for replicates
    px     <- read_hiseq(x)
    p_list <- px$args$rep_list # single dirs
    tmp    <- lapply(p_list, read_hiseq_align)
    dplyr::bind_rows(tmp)

    # } else if(endsWith(hiseq_type, "_rx")) {
  } else if(is_hiseq_multiple_dir(x)) {
    # for ip, input
    px <- read_hiseq(x)

    # get dirs for replicates
    p1     <- px$args$ip_dir
    p2     <- px$args$input_dir
    p_list <- unlist(c(p1, p2))
    tmp    <- lapply(p_list, read_hiseq_align)
    dplyr::bind_rows(tmp)
  }
}


#' read_hiseq_peak_stat
#'
#' number of peaks
#' peak reads
#'
#' @param x path to the directory
#'
#' @export
read_hiseq_peak_stat <- function(x) {
  # search for single dir
  rep_list   <- list_hiseq_single_dirs(x)
  merge_list <- list_hiseq_merge_dirs(x)
  input_list <- sort(c(rep_list, merge_list))

  lapply(input_list, function(i){
    px <- read_hiseq(i)
    data.frame(
      id    = px$args$smp_name,
      count = length(readLines(px$args$peak))
    )
  }) %>%
    dplyr::bind_rows()
}


#' read_hiseq_lendist_stat
#'
#' length distribution
#'
#' @param x path to the directory
#'
#' @export
read_hiseq_lendist_stat <- function(x) {
  # search for single dir
  rep_list   <- list_hiseq_single_dirs(x)
  merge_list <- list_hiseq_merge_dirs(x)
  input_list <- sort(c(rep_list, merge_list))

  lapply(input_list, function(i){
    px <- read_hiseq(i)
    f  <- px$args$lendist_txt
    if(file.exists(f)) {
      read_frag(f)
    } else {
      NULL
    }
  }) %>%
    dplyr::bind_rows()
}



#'
#' Fraction in peak
#'
#' @param x path to the directory
#'
#' @export
read_hiseq_frip_stat <- function(x) {
  # search for single dir
  rep_list   <- list_hiseq_single_dirs(x)
  merge_list <- list_hiseq_merge_dirs(x)
  input_list <- sort(c(rep_list, merge_list))

  lapply(input_list, function(i){
    px <- read_hiseq(i)
    f  <- px$args$frip_txt

    if(file.exists(f)) {
      readr::read_delim(f, delim = "\t") %>%
        dplyr::mutate(id   = px$args$smp_name,
                      FRiP = paste0(round(FRiP * 100, 2), "%")) %>%
        dplyr::select(id, total, peak_reads, FRiP)
    } else {
      NULL
    }
  }) %>%
    dplyr::bind_rows()
}


#' read HiSeq directory
#'
#' input dir
#' output files, config
#'
#' @param x string path
#'
#' @export
read_hiseq <- function(x) {
  if(is_hiseq_dir(x)) {
    message(paste0("Reading HiSeq directory: ", x))
  } else {
    warning(paste0("Not a HiSeq directory: ", x))
    return(NULL)
  }

  # read args
  pk   <- list_arguments_file(x)
  args <- load_pickle(pk[1])

  list(hiseq_type = get_hiseq_type(x),
       files      = file_to_list(x, recursive = TRUE),
       args       = args)
}





#' @param x string, path to dir
#'
#' @export
is_hiseq_single_dir <- function(x) {
  hiseq_type <- get_hiseq_type(x)
  ifelse(is.null(hiseq_type), FALSE, endsWith(hiseq_type, "_r1"))
}


#' @param x string, path to dir
#'
#' @export
is_hiseq_merge_dir <- function(x) {
  hiseq_type <- get_hiseq_type(x)
  ifelse(is.null(hiseq_type), FALSE, endsWith(hiseq_type, "_rn"))
}


#' @param x string, path to dir
#'
#' @export
is_hiseq_multiple_dir <- function(x) {
  hiseq_type <- get_hiseq_type(x)
  ifelse(is.null(hiseq_type), FALSE, endsWith(hiseq_type, "_rx"))
}




#' is_hiseq_dir
#'
#' check config.pickle
#'
#' hiseq_type, ...
#'
is_hiseq_dir <- function(x) {
  st <- get_hiseq_type(x)
  length(st) == 1 #
}



#' hiseq_type
#'
#' @param x path, the directory of RNAseq output
#'
#' @export
get_hiseq_type <- function(x) {
  # check pickle file
  pk_files <- list_arguments_file(x)

  if(is.null(pk_files)) {
    return(NULL)
  }

  # choose the first one
  pk <- load_pickle(pk_files[1])

  # candidate seqtype
  st <- c("hiseq_type", "rnaseq_type", "atacseq_type", "align_type",
          "chipseq_type")

  ht <- lapply(st, function(i){
    if(i %in% names(pk)) {
      pk[[i]]
    }
  })

  ht <- unlist(ht)

  ht
}




#' list_arguments_file
#'
#' @param x path to the directory.
#' expect
#' x/config/arguments.pickle
#' x/{feature}/config/arguments.pickle
#' x/config/{feature}/config/arguments.pickle,
#'
#' filename could be config.pickle
#'
#' @export
list_arguments_file <- function(x) {
  # arguments file
  fname <- c("arguments.pickle", "config.pickle")

  # version-1: x/config/...
  f1 <- file.path(x, "config", fname)

  # version-2: x/{feature}/config/...
  f2 <- file.path(x, "*", "config", fname)

  # version-3: x/config/{feature}/config/arguments.txt
  f3 <- file.path(x, "config",  "*", "config", fname)

  # check exists
  f_list <- c(f1, f2, f3)
  Sys.glob(f_list)
}



#' load pickle file
#'
#' @param x path to pickle,
#'
#' @import reticulate
#' @export
#'
load_pickle <- function(x) {
  # python version issue !!
  if (file.exists(x) & endsWith(x, ".pickle")) {
    pd <- reticulate::import("pandas")
    tag <- pd$read_pickle(x)
  } else {
    tag <- NULL
  }
}





#' list_hiseq_single_dirs: hiseq_r1
#'
#' @param x path to the input directories
#' @param log boolen, whether print the log status, default: FALSE
#'
#' @export
list_hiseq_single_dirs <- function(x){
  if(is_hiseq_single_dir(x)) {
    x
  } else if(is_hiseq_merge_dir(x)) {
    rep_list <- read_hiseq(x)$args$rep_list
    sapply(rep_list, list_hiseq_single_dirs, simplify = TRUE)
  } else if(is_hiseq_multiple_dir(x)) {
    hiseq_type <- get_hiseq_type(x)
    ## for CnR, ChIPseq
    if(grepl("^hiseq_|^chipseq_|^rnaseq_|^cnt_", hiseq_type)) {
      px     <- read_hiseq(x)
      p_list <- unlist(c(px$args$ip_dir, px$args$input_dir))
      sapply(p_list, list_hiseq_single_dirs, simplify = TRUE)
    } else if(grepl("^atacseq_", hiseq_type)) {
      d1 <- list.dirs(x, recursive = FALSE)
      d2 <- unlist(sapply(d1, is_hiseq_single_dir))
      names(d2[d2])
    }
  }
}


#' hiseq_merge_dirs
#'
#' @param x path to the input directories
#'
#' @export
list_hiseq_merge_dirs <- function(x){
  if(is_hiseq_single_dir(x)) {
    NULL
  } else if(is_hiseq_merge_dir(x)) {
    x
  } else if(is_hiseq_multiple_dir(x)) {
    # px <- read_hiseq(x)
    # sapply(c(px$args$ip_dir, px$args$input_dir), list_hiseq_merge_dirs,
    #        simplify = TRUE)
    hiseq_type <- get_hiseq_type(x)
    ## for CnR, ChIPseq
    if(grepl("^hiseq_|^chipseq_|^rnaseq_|^cnt_", hiseq_type)) {
      px     <- read_hiseq(x)
      p_list <- unlist(c(px$args$ip_dir, px$args$input_dir))
      sapply(p_list, list_hiseq_merge_dirs, simplify = TRUE)
    } else if(grepl("^atacseq_", hiseq_type)) {
      d1 <- list.dirs(x, recursive = FALSE)
      d2 <- unlist(sapply(d1, is_hiseq_merge_dir))
      names(d2[d2])
    }
  }
}


#' hiseq_merge_dirs
#'
#' @param x path to the input directories
#'
#' @export
list_hiseq_multiple_dirs <- function(x){
  if(is_hiseq_multiple_dir(x)) {
    x
  }
}




#' list all hiseq report files
#'
#' @export
list_hiseq_report <- function(x) {
  # search for single dir
  rep_list      <- list_hiseq_single_dirs(x)
  merge_list    <- list_hiseq_merge_dirs(x)
  multiple_list <- list_hiseq_multiple_dirs(x)
  input_list    <- sort(c(rep_list, merge_list, multiple_list))

  report_list <- sapply(input_list, function(i){
    px <- read_hiseq(i)
    report_dir <- px$args$report_dir
    f_html  <- list.files(report_dir, "*.html", full.names = TRUE)
    if(length(f_html) > 0) {
      f_html[1]
    }
  }, simplify = TRUE, USE.NAMES = FALSE) %>%
    unlist()

  report_list[!is.null(report_list)]
}












