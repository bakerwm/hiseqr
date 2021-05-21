#' Functions for HiSeq pipeline
#'
#' Reading files, directories, pipeline, ...
#' Processing data
#' Prepare for plot
#' Plotting
#'
#' @name hiseq


#' hiseq_report
#'
#' @param input directory to the sample
#' @param output directory to the html file
#' @param template the template, default from hiseqr
#'
#' @export
hiseq_report <- function(input, output, template_rmd = NULL) {
  ## input
  input  <- normalizePath(input)
  if(! dir.exists(output)) {
    dir.create(output, recursive = TRUE)
  }
  output <- normalizePath(output)
  ## output
  outhtml <- file.path(output, "HiSeq_report.html")

  ## determine the template
  hiseq_type = get_hiseq_type(input)

  # hiseq (hiseq, atacseq, chipseq, rnaseq, ...)
  if(is_hiseq_dir(input)) {
    if(startsWith(hiseq_type, "hiseq_")) {
      hiseq_dir = "hiseq"
    } else if(startsWith(hiseq_type, "cnr_")) {
      hiseq_dir = "cnr"
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
  if(is.null(template_rmd)) {
    template <- system.file(hiseq_dir, hiseq_subtype, package = "hiseqr")
  } else {
    template <- template_rmd
  }
  stopifnot(file.exists(template))

  ## copy template to output
  template_to <- file.path(output, basename(template))
  file.copy(template, template_to, overwrite = TRUE)
  rmarkdown::render(input       = template_to,
                    output_file = outhtml,
                    params      = list(input_dir = input))
}


#' @describeIn  read_hiseq_trim_stat
#'
#' parsing the trimming status
#'
#' @param x path to hiseq, single
#'
#' @import dplyr
#' @import readr
#'
#' @export
read_hiseq_trim_stat <- function(x) {
  dirs <- list_hiseq_single_dirs(x)
  dirs <- purrr::discard(dirs, is.null)
  if(length(dirs) > 0) {
    lapply(dirs, function(i) {
      px <- read_hiseq(i)
      px$args$trim_summary_json %>%
        jsonlite::read_json() %>%
        as.data.frame.list() %>%
        dplyr::select(id, input, output, out_pct, rm_pct)
    }) %>%
      dplyr::bind_rows() %>%
      dplyr::mutate(across(!id, as.numeric)) %>%
      dplyr::rename(clean_pct = out_pct,
                    short_pct = rm_pct) %>%
      dplyr::mutate(clean_pct = round(clean_pct, 2),
                    short_pct = round(short_pct, 2))
  }
}





#' read_hiseq_align
#'
#' @param x path to hiseq, single
#'
#' @description output columns
#' "fqname", "index_name", "total", "map", "unique", "multiple",
#' "unmap", "nodup", "chrM", "spikein"
#'
#' @import readr
#' @import dplyr
#'
#' @export
read_hiseq_align_stat <- function(x){
  dirs <- list_hiseq_single_dirs(x)
  dirs <- purrr::discard(dirs, is.null)
  if(length(dirs) > 0) {
    lapply(dirs, function(i) {
      px <- read_hiseq(i)
      f1 <- px$args$align_summary_json # ver-1: qc/01.alignment_summary.json
      f2 <- px$args$align_json # ver-2: alignment/*.json
      if(file.exists(f1)) {
        cols <- c("fqname", "index_name", "total", "map", "unique", "multiple",
                  "unmap", "nodup", "chrM", "spikein")
        jsonlite::read_json(f1) %>%
          as.data.frame.list %>%
          dplyr::select(all_of(cols)) %>%
          dplyr::mutate(map_pct  = round(map / total * 100, 2),
                        dup_pct  = round((map - nodup) / total * 100, 2),
                        chrm_pct = round(chrM / total * 100, 2),
                        unique_pct = round(unique / total * 100, 2),
                        unmap_pct  = round(unmap / total * 100, 2))
      } else if(file.exists(f2)) {
        cols <- c("fqname", "index_name", "total", "map", "unique", "multiple",
                  "unmap")
        jsonlite::read_json(f2) %>%
          as.data.frame.list %>%
          dplyr::select(all_of(cols)) %>%
          dplyr::mutate(nodup    = map,
                        chrM     = 0,
                        spikein  = 0,
                        map_pct  = round(map / total * 100, 2),
                        dup_pct  = round((map - nodup) / total * 100, 2),
                        chrm_pct = round(chrM / total * 100, 2),
                        unique_pct = round(unique / total * 100, 2),
                        unmap_pct  = round(unmap / total * 100, 2))
      }
    }) %>%
      dplyr::bind_rows()
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
    px     <- read_hiseq(i)
    n_peak <- tryCatch(
      error = function(cnd) 0,
      length(readLines(px$args$peak))
    )
    data.frame(
      id    = px$args$smp_name,
      count = n_peak
    )
  }) %>%
    dplyr::bind_rows()
}





#' @describeIn read_hiseq_lendist_stat
#'
#' length distribution
#'
#' @param x path to the directory
#'
#' @export
read_hiseq_lendist_stat <- function(x) {
  # search for single dir
  r1_list <- list_hiseq_single_dirs(x)
  rn_list <- list_hiseq_merge_dirs(x)
  dirs    <- sort(c(r1_list, rn_list))
  dirs    <- purrr::discard(dirs, is.null)
  if(length(dirs) > 0) {
    lapply(dirs, function(i){
      px <- read_hiseq(i)
      f  <- px$args$lendist_txt
      if(file.exists(f)) {
        read_text(f)
      }
    }) %>%
      dplyr::bind_rows()
  }
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
      df <- read_text(f) %>%
        dplyr::mutate(id   = px$args$smp_name,
                      FRiP = paste0(round(FRiP * 100, 2), "%"))
      if("total_reads" %in% names(df)) {
        df$total <- df$total_reads
      }
      df %>% dplyr::select(id, total, peak_reads, FRiP)
    }
  }) %>%
    dplyr::bind_rows()
}










#' @param x string, path to dir
#'
#' @export
is_hiseq_single_dir <- function(x) {
  ht <- get_hiseq_type(x)
  ifelse(is.null(ht), FALSE, endsWith(ht, "_r1"))
}


#' @param x string, path to dir
#'
#' @export
is_hiseq_merge_dir <- function(x) {
  ht <- get_hiseq_type(x)
  ifelse(is.null(ht), FALSE, endsWith(ht, "_rn"))
}


#' @param x string, path to dir
#'
#' @export
is_hiseq_multiple_dir <- function(x) {
  ht <- get_hiseq_type(x)
  ifelse(is.null(ht), FALSE, endsWith(ht, "_rx"))
}



#' is_hiseq_dir
#'
#' check config.pickle
#'
#' hiseq_type, ...
#'
is_hiseq_dir <- function(x) {
  ht <- get_hiseq_type(x)
  is(ht, "character")
}


#' hiseq_type
#'
#' @param x path, the directory of RNAseq output
#'
#' @export
get_hiseq_type <- function(x) {
  px <- read_hiseq(x)
  if(is(px, "list")) {
    px$hiseq_type
  }
}







#' list_hiseq_single_dirs: hiseq_r1
#'
#' @param x path to the input directories
#' @param log boolen, whether print the log status, default: FALSE
#'
#' @export
list_hiseq_single_dirs <- function(x){
  if(is_hiseq_dir(x)) {
    if(is_hiseq_single_dir(x)) {
      x
    } else {
      if(is_hiseq_merge_dir(x)) {
        px   <- read_hiseq(x)
        px$args$rep_list
      } else if(is_hiseq_multiple_dir(x)) {
        px <- read_hiseq(x)
        hiseq_type <- get_hiseq_type(x)
        if(grepl("^chipseq_|^cnt_|^cnr_|^hiseq_", hiseq_type)) {
          rn <- c(px$args$ip_dir, px$args$input_dir)
        } else if(grepl("^rnaseq_", hiseq_type)) {
          rn <- c(px$args$wildtype_dir, px$args$mutant_dir)
        } else if(grepl("^atacseq_", hiseq_type)) {
          rn <- list.dirs(x, recursive = FALSE)
          rn <- purrr::keep(rn, is_hiseq_merge_dir)
        } else {
          rn <- NULL
        }
        lapply(rn, list_hiseq_single_dirs) %>% unlist
      } else {
        NULL
      }
    }
  }
}




#' hiseq_merge_dirs
#'
#' @param x path to the input directories
#'
#' @export
list_hiseq_merge_dirs <- function(x){
  if(is_hiseq_dir(x)) {
    if(is_hiseq_merge_dir(x)) {
      x
    } else {
      if(is_hiseq_multiple_dir(x)) {
        px <- read_hiseq(x)
        hiseq_type <- get_hiseq_type(x)
        if(grepl("^chipseq_|^cnt_|^cnr_", hiseq_type)) {
          dirs <- c(px$args$ip_dir, px$args$input_dir)
        } else if(grepl("^rnaseq_", hiseq_type)) {
          dirs <- c(px$args$wildtype_dir, px$args$mutant_dir)
        } else if(grepl("^atacseq_", hiseq_type)) {
          dirs <- list.dirs(x, recursive = FALSE)
        } else {
          dirs <- NULL
        }
        purrr::keep(dirs, is_hiseq_merge_dir)
      }
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










#' read HiSeq directory
#'
#' input dir
#' output files: pickle, toml, json
#'
#' @param x string path
#'
#' @export
read_hiseq <- function(x) {
  if(is(x, "character")) {
    x  <- x[1]
    pk <- list_arguments_file(x)
    if(length(pk)) {
      args <- load_config(pk[1])
      ht_list <- c("hiseq_type", "rnaseq_type", "atacseq_type",
                   "align_type", "chipseq_type")
      hiseq_type <- sapply(ht_list, function(ht) {
        if(ht %in% names(args)) {
          args[[ht]]
        }
      }, simplify = TRUE) %>%
        unlist %>%
        head
      list(
        hiseq_type = hiseq_type,
        files      = path_to_list(x[1], recursive = TRUE),
        args       = args
      )
    }
  }
}


#' @describeIn load_config Parsing config files
#'
#' Python version
#'
#' load pickle|toml file
#'
#' @param x path to pickle,
#'
#' @import reticulate
#' @export
#'
load_config <- function(x) {
  # python version issue
  if(is(x, "character")) {
    x <- x[1] # only first item
    if (file.exists(x) & endsWith(x, ".pickle")) {
      reticulate::use_condaenv("hiseq", required = TRUE) # !!! switch to toml ?
      pd <- reticulate::import("pandas")
      pd$read_pickle(x)
    } else if(endsWith(x, ".toml")) {
      configr::read.config(x)
    }
  }
}



#' @describeIn list_arguments_file Search config files, .pickle, .toml
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
  # fname <- c("arguments.pickle", "config.pickle", "config.toml")
  fname <- c("config.toml", "config.json", "config.pickle", "arguments.pickle")
  # version-1: x/config/...
  f1 <- file.path(x, "config", fname)
  # version-2: x/...
  f2 <- file.path(x, fname)
  # # version-2: x/{feature}/config/..., deprecated
  # f2 <- file.path(x, "*", "config", fname)
  # # version-3: x/config/{feature}/config/arguments.txt, deprecated
  # f3 <- file.path(x, "config",  "*", "config", fname)
  # check exists
  f_list <- c(f1, f2)
  Sys.glob(f_list)
}




##-- deprecated functions ------------------------------------------------------
#'
#' #' read_hiseq_align
#' #'
#' #' @param x path to hiseq, single
#' #' @import readr
#' #' @import dplyr
#' #'
#' #' @export
#' read_hiseq_align_stat <- function(x){
#'   # hiseq_type <- get_hiseq_type(x)
#'
#'   # if(endsWith(hiseq_type, "_r1")) {
#'   if(is_hiseq_single_dir(x)) {
#'     px <- read_hiseq(x)
#'
#'     # version-1
#'     # from qc/01.alignment_summary.json
#'     # version-2
#'     # from alignment/*.json
#'     j1 <- px$args$align_summary_json
#'     j2 <- px$args$align_json
#'
#'     if(file.exists(j1)) {
#'       m1 <- jsonlite::read_json(j1)
#'     } else if(file.exists(j2)) {
#'       m1 <- jsonlite::read_json(j2)
#'     } else {
#'       m1 <- NULL
#'     }
#'     df <- as.data.frame.list(m1)
#'
#'     ## arrange columns
#'     cols <- c("fqname", "index_name", "total", "map", "unique", "multiple",
#'               "unmap", "nodup", "chrM", "spikein")
#'
#'     ## pick columns
#'     cols <- cols[cols %in% colnames(df)]
#'
#'     ## version-1
#'     if(all(c("nodup", "chrM", "spikein") %in% cols)) {
#'       df %>%
#'         dplyr::select(all_of(cols)) %>%
#'         dplyr::mutate(map_pct  = round(map / total * 100, 2),
#'                       dup_pct  = round((map - nodup) / total * 100, 2),
#'                       chrm_pct = round(chrM / total * 100, 2))
#'     } else {
#'       ## version-2
#'       df %>%
#'         dplyr::select(all_of(cols)) %>%
#'         dplyr::mutate(nodup    = map,
#'                       chrM     = 0,
#'                       spikein  = 0,
#'                       map_pct  = round(map / total * 100, 2),
#'                       dup_pct  = round((map - nodup) / total * 100, 2),
#'                       chrm_pct = round(chrM / total * 100, 2))
#'     }
#'
#'     # } else if(endsWith(hiseq_type, "_rn")) {
#'   } else if(is_hiseq_merge_dir(x)) {
#'     # for replicates
#'     px     <- read_hiseq(x)
#'     p_list <- px$args$rep_list # single dirs
#'     tmp    <- lapply(p_list, read_hiseq_align)
#'     dplyr::bind_rows(tmp)
#'
#'     # } else if(endsWith(hiseq_type, "_rx")) {
#'   } else if(is_hiseq_multiple_dir(x)) {
#'     # for ip, input
#'     px <- read_hiseq(x)
#'
#'     # get dirs for replicates
#'     p1     <- px$args$ip_dir
#'     p2     <- px$args$input_dir
#'     p_list <- unlist(c(p1, p2))
#'     tmp    <- lapply(p_list, read_hiseq_align)
#'     dplyr::bind_rows(tmp)
#'   }
#' }
#'
