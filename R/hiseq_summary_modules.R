#' Functions for HiSeq summary
#'
#' Summary data in project_dir.
#' Processing data
#' Prepare for plot
#' Plotting
#'


#'  read_hiseq_stat
#' @param x path to hiseq, single
#' @param keys character, quality control groups
#'  options: c("trim", "align", "peak", "lendist", "frip", "report",
#'  "enrich", "cor")
#' @param add_tag add SHA-256 value of dirnname(x), first-7 characters
#'
#'
#'
#' @export
read_hiseq_stat <- function(x, keys = "align", add_tag = FALSE) {
  if(!inherits(x, "character")) {
    message(glue::glue("read_hiseq_stat() faild, x is {class(x)}, expect character"))
    return(NULL)
  }
  # filter
  x <- purrr::keep(x, is_hiseq_dir)
  if(length(x) == 0) {
    message("read_hiseq_stat() failed, not enough x")
    return(NULL)
  }
  # for keys
  k_list <- c(
    "trim", "align", "peak", "lendist", "frip", "report",
    "enrich", "cor")
  if(isTRUE(keys)) {
    keys <- k_list
  }
  keys <- purrr::keep(keys, function(i) i %in% k_list)
  k_rm <- purrr::discard(keys, function(i) i %in% k_list)
  if(length(k_rm) > 0) {
    k_rm_str <- paste(k_rm, collapse = ",")
    message(glue::glue("unknown keys: {k_rm_str}"))
  }
  if(length(keys) == 0) {
    k_str <- paste(k_list, collapse = ",")
    message(glue::glue("no keys found, expect: {k_str}"))
    return(NULL)
  }

  out <- sapply(keys, function(k) {
    f  <- paste0("read_hiseq_", k, "_stat")
    fn <- tryCatch(
      error = function(cnd) {
        warning(paste0("unknown keys: ", k))
        NULL
      },
      match.fun(f)
    )
    # do the things
    tmp2 <- lapply(x, function(i) {
      tmp <- fn(i) # result
      if(isTRUE(add_tag) & inherits(tmp, "data.frame")) {
        tag <- substr(hash_string(dirname(i)), 1, 7)
        tmp$label <- paste0(basename(i), ".", tag)
      }
      tmp
    })
    names(tmp2) <- x
    tmp2 <- purrr::discard(tmp2, is.null)
    # merge data.frame
    if(all(sapply(tmp2, is.data.frame))) {
      tmp2 <- dplyr::bind_rows(tmp2)
    }
    tmp2
  }, USE.NAMES = TRUE, simplify = FALSE)
  # # assign names
  # names(out) <- keys
  out
}


#'  read_hiseq_trim_stat
#'
#' @param x path to hiseq, single
#'
#'
#' @export
read_hiseq_trim_stat <- function(x) {
  p_out <- lapply(x, function(i) {
    if(is_hiseq_dir(i)) {
      j <- list_hiseq_file(i, "trim_summary_json", hiseq_type = TRUE)
      tryCatch(
        {
          # as.data.frame(jsonlite::read_json(j))
          read_hiseq_trim_json(j)
        },
        error = function(cond) {
          return(NULL)
        }
      )
    }
  })
  # single
  if(length(x) == 1 & inherits(p_out, "list")) {
    p_out <- p_out[[1]]
  } else if(length(x) > 1) {
    names(p_out) <- x
  }
  # return
  p_out
}


#'  read_hiseq_align_stat
#'
#' @param x path to hiseq, single
#'
#'
#' @export
read_hiseq_align_stat <- function(x) {
  p_out <- lapply(x, function(i) {
    if(is_hiseq_dir(i)) {
      j <- list_hiseq_file(i, "align_summary_json", hiseq_type = TRUE)
      tryCatch(
        {
          read_hiseq_align_json(j)
        },
        error = function(cond) {
          return(NULL)
        }
      )
    }
  })
  # single
  if(length(x) == 1 & inherits(p_out, "list")) {
    p_out <- p_out[[1]]
  } else if(length(x) > 1) {
    names(p_out) <- x
  }
  # return
  p_out
}



#'  read_hiseq_trim_json
#'
#' @param x json file
#'
#'
#' @export
read_hiseq_trim_json <- function(x) {
  # filtering, only .json file
  x <- purrr::keep(x, function(i) {
    file.exists(i) & endsWith(i, ".json")
  })
  if(length(x) < 1) {
    message("failed reading align_json")
    return(NULL)
  }
  lapply(x, function(f) {
    df <- lapply(f, function(i) {
      jsonlite::read_json(i) %>%
        as.data.frame()
    }) %>%
      dplyr::bind_rows()
  }) %>%
    dplyr::bind_rows()
}



#'  read_hiseq_align_json
#'
#' @param x json file
#'
#'
#' @export
read_hiseq_align_json <- function(x) {
  # filtering, only .json file
  x <- purrr::keep(x, function(i) {
    file.exists(i) & endsWith(i, ".json")
  })
  if(length(x) < 1) {
    message("failed reading align_json")
    return(NULL)
  }
  lapply(x, function(f) {
    df <- lapply(f, function(i) {
      jsonlite::read_json(i) %>%
        as.data.frame()
    }) %>%
      dplyr::bind_rows()
    # column - common
    t_cols <- c("name", "total", "map", "unique", "multi", "unmap")
    t1 <- c("chrM", "spikein") # atac, cnr
    t2 <- c("rRNA", "spikein") # rnaseq
    t3 <- c("dup", "nodup")
    if(all(t1 %in% names(df))) {
      t_cols <- c(t_cols, t1)
    } else if(all(t2 %in% names(df))) {
      t_cols <- c(t_cols, t2)
    } else {
      warning(paste0("unknown hiseq: ", x))
    }
    if(all(t3 %in% names(df))) {
      t_cols <- c(t_cols, t3)
    }
    dplyr::select(df, all_of(t_cols))
  }) %>%
    dplyr::bind_rows()
}


#' read_hiseq_peak_stat
#' @param x path to the directory
#'
#' @export
read_hiseq_peak_stat <- function(x) {
  p_out <- lapply(x, function(i) {
    out <- NULL
    if(is_hiseq_dir(i)) {
      # load json
      j_list <- sapply(c("peak", "macs2_peak", "seacr_peak"), function(j) {
        list_hiseq_file(i, j)
      })
      j_list <- purrr::discard(j_list, is.null)
      if(length(j_list) > 0) {
        n_peak <- tryCatch(
          {
            length(readLines(unlist(j_list[[1]])[1]))
          },
          error = function(cnd) {
            return(0)
          })
        # data.frame
        out <- data.frame(
          id    = list_hiseq_file(i, "smp_name"),
          count = n_peak
        )
      }
    }
    # return
    out
  })
  # single
  if(length(x) == 1 & inherits(p_out, "list")) {
    p_out <- p_out[[1]]
  } else if(length(x) > 1) {
    names(p_out) <- x
  }
  # return
  p_out
}


#' read_hiseq_lendist_stat
#'
#' @param x path to the directory
#'
#' @export
read_hiseq_lendist_stat <- function(x) {
  p_out <- lapply(x, function(i) {
    if(is_hiseq_dir(i)) {
      # load json
      j_list <- sapply(c("lendist_csv", "lendist_json"), function(j) {
        list_hiseq_file(i, j)
      })
      j_list <- purrr::discard(j_list, is.null)
      if(length(j_list) > 0) {
        tryCatch(
          {
            read_hiseq_lendist_json(unlist(j_list[[1]]))
          },
          error = function(cond) {
            return(NULL)
          }
        )
      }
    }
  })
  # single
  if(length(x) == 1 & inherits(p_out, "list")) {
    p_out <- p_out[[1]]
  } else if(length(x) > 1) {
    names(p_out) <- x
  }
  # return
  p_out
}


#' read_hiseq_lendist_json
#'
#' @param x path to the csv/json file
#'
#' @export
read_hiseq_lendist_json <- function(x) {
  # filtering, only .csv file
  x <- purrr::keep(x, function(i) {
    file.exists(i) & (endsWith(i, ".csv") | endsWith(i, ".json"))
  })
  lapply(x, function(f) {
    tryCatch(error = function(cnd) NULL, read_text(f))
  }) %>%
    dplyr::bind_rows()
}


#' read_hiseq_frip_stat
#'
#' @param x path to the directory
#'
#' @export
read_hiseq_frip_stat <- function(x) {
  p_out <- lapply(x, function(i) {
    out <- NULL
    if(is_hiseq_dir(i)) {
      # load json
      j_list <- sapply(c("frip_json", "frip_toml", "frip_txt"), function(j) {
        list_hiseq_file(i, j)
      })
      j_list <- purrr::discard(j_list, is.null)
      if(length(j_list) > 0) {
        out <- read_hiseq_frip_json(unlist(j_list[[1]]))
      }
    } else {
      out <- list(tss = NULL, genebody = NULL)
    }
    # return
    out
  })
  # single
  if(length(x) == 1 & inherits(p_out, "list")) {
    p_out <- p_out[[1]]
  } else if(length(x) > 1) {
    names(p_out) <- x
  }
  # return
  p_out
}


#' read_hiseq_frip_json
#'
#' @param x json file
#'
#' @export
read_hiseq_frip_json <- function(x) {
  lapply(x, function(j) {
    y <- tryCatch(error = function(cnd) NULL, jsonlite::read_json(j))
    if(is(y, "list")) {
      as.data.frame(y)
    }
  }) %>%
    dplyr::bind_rows()
}


#' read_hiseq_enrich_stat
#' @param x path to the directory
#'
#' @export
read_hiseq_enrich_stat <- function(x) {
  p_out <- lapply(x, function(i) {
    if(is_hiseq_dir(i)) {
      out <- list(
        tss = list_hiseq_file(i, "tss_enrich_png"),
        genebody = list_hiseq_file(i, "genebody_enrich_png")
      )
    } else {
      out <- list(tss = NULL, genebody = NULL)
    }
  })
  # single
  if(length(x) == 1 & inherits(p_out, "list")) {
    p_out <- p_out[[1]]
  } else if(length(x) > 1) {
    names(p_out) <- x
  }
  # return
  p_out
}


#' read_hiseq_cor_stat
#'
#' @param x path to the directory
#'
#' @export
read_hiseq_cor_stat <- function(x) {
  p_out <- lapply(x, function(i) {
    if(is_hiseq_dir(i)) {
      qc_dir <- list_hiseq_file(i, "qc_dir")
      if(inherits(qc_dir, "character")) {
        out <- list(
          cor_heatmap     = list_hiseq_file(i, "bam_cor_heatmap_png"),
          cor_pca         = list_hiseq_file(i, "bam_cor_pca_png"),
          cor_fingerprint = list_hiseq_file(i, "bam_fingerprint_png"),
          peak_overlap    = list_hiseq_file(i, "peak_overlap_png"),
          peak_idr        = list_hiseq_file(i, "peak_idr_png"), # error
          peak_idr_list   = list.files(qc_dir, ".*peak_idr.*png", full.names = TRUE)
        )
      }
    }
  })
  # single
  if(length(x) == 1 & inherits(p_out, "list")) {
    p_out <- p_out[[1]]
  } else if(length(x) > 1) {
    names(p_out) <- x
  }
  # return
  p_out
}


#' read_hiseq_report_stat
#'
#' @param x path to the directory
#'
#' @export
read_hiseq_report_stat <- function(x) {
  p_out <- lapply(x, function(i) {
    if(is_hiseq_dir(i)) {
      list_hiseq_file(i, "report_html")
    }
  })
  # single
  if(length(x) == 1) {
    p_out <- p_out[[1]]
  } else if(length(x) > 1) {
    names(p_out) <- x
  }
  p_out
}


#' read_rnaseq_deseq
#'
#' @param x character path to the rnaseq_rx
#'
#' @export
read_rnaseq_deseq <- function(x, sig_list = TRUE) {
  if(is_hiseq_dir(x, "rnaseq_rx")) {
    tryCatch(
      {
        f  <- get_fix_xls(x) # transcripts_deseq2.fix.xls
        df <- readr::read_delim(f, "\t", col_types = readr::cols())
        s_all <- unique(df$sig)
        if(isTRUE(sig_list)) sig_list <- s_all
        sig_list <- purrr::keep(sig_list, function(i) i %in% s_all)
        # subset
        dplyr::filter(df, sig %in% sig_list)
      },
      error = function(cond) {
        stderr = glue::glue(
          "Failed: reading *.fix.xls file in {x}"
        )
        message(stderr)
        # return(NA)
      },
      warning = function(cond) {
        stdout = glue::glue(
          "x caused a warning: {x}",
          "Here's the message: {cond}",
          .sep = "\n"
        )
        message(stdout)
        # return(NULL)
      },
      finally={
        msg = glue::glue(
          "Reading deseq data from directory: {x}"
        )
        message(msg)
      }
    )
  }
}
