#' Functions for HiSeq summary
#'
#' Summary data in project_dir.
#' Processing data
#' Prepare for plot
#' Plotting
#'
#' @name hiseq_summary_modules



#' @describeIn  read_hiseq_stat
#'
#'
#'
#' input: trim_json/trim_stat_json
#' output: name, total, clean, too_short, too_short2, ...
#'
#'
#' @param x path to hiseq, single
#' @param keys character, quality control groups
#'  options: c("trim", "align", "peak", "lendist", "frip", "report")
#'
#' @import dplyr
#' @import readr
#'
#' @export
read_hiseq_stat <- function(x, keys = "align") {
  # check arguments
  pd <- read_hiseq(x)
  if(! is(pd, "list")) {
    warning(paste0("Not a hiseq directory: ", x))
    return(NULL)
  }
  # tss_enrich
  # genebody_enrich
  # bam_cor: PCA, heatmap
  # bam_fingerprint
  # peak_overlap
  # peak_idr
  k_list <- c(
    "trim", "align", "peak", "lendist", "frip", "report",
    "enrich", "cor")
  if(isTRUE(keys)) {
    keys <- k_list
  }
  keys   <- purrr::keep(keys, function(i) i %in% k_list)
  k_rm   <- purrr::discard(keys, function(i) i %in% k_list)
  if(length(k_rm) > 0) {
    warning(paste(c("unknown keys skipped:", k_rm), collapse = ", "))
  }
  if(length(keys) == 0) {
    warning(paste(c("no keys, expect", k_list), collapse = ", "))
    return(NULL)
  }
  # output
  out <- lapply(keys, function(k) {
    # build function
    f  <- paste0("read_hiseq_", k, "_stat")
    fn <- tryCatch(error = function(cnd) {
      warning(paste0("unknown keys: ", k))
      NULL
    },
    match.fun(f))
    # do the things
    fn(x) %>% unique()
  })
  # assign names
  names(out) <- keys
  out
}


##----------------------------------------------------------------------------##
## sub-modules for summary

#' @describeIn  read_hiseq_trim_stat
#'
#' Only for r1 directory (atac, ... )
#' or trim dir
#' parsing the trimming status
#'
#' input: trim_json/trim_stat_json
#' output: name, total, clean, too_short, too_short2, ...
#'
#'
#' @param x path to hiseq, single
#'
#' @import dplyr
#' @import readr
#'
#' @export
read_hiseq_trim_stat <- function(x) {
  f <- list_hiseq_file(x, "trim_summary_json", "r1")
  if(length(f) > 0) {
    lapply(f, function(i) {
      jsonlite::read_json(i) %>%
        as.data.frame()
    }) %>%
      dplyr::bind_rows()
  }
}




#' @describeIn  read_hiseq_align_stat
#'
#' Only for r1 directory (atac, ... )
#' or trim dir
#' parsing the alignment status
#'
#' input: trim_json/trim_stat_json
#' output: name, total, clean, too_short, too_short2, ...
#'
#'
#' @param x path to hiseq, single
#'
#' @import dplyr
#' @import readr
#'
#' @export
read_hiseq_align_stat <- function(x) {
  # pd <- read_hiseq(x)
  f <- list_hiseq_file(x, "align_summary_json", "r1")
  if(length(f) > 0) {
    df <- lapply(f, function(i) {
      jsonlite::read_json(i) %>%
        as.data.frame()
    }) %>%
      dplyr::bind_rows()
    # column - common
    t_cols <- c("name", "total", "map", "unique", "multi", "unmap")
    t1 <- c("chrM", "spikein") # atac, cnr
    t2 <- c("rRNA", "spikein") # rnaseq
    if(all(t1 %in% names(df))) {
      t_cols <- c(t_cols, t1)
    } else if(all(t2 %in% names(df))) {
      t_cols <- c(t_cols, t2)
    } else {
      warning(paste0("unknown hiseq: ", x))
    }
    dplyr::select(df, all_of(t_cols))
  }
}






#' @describeIn read_hiseq_peak_stat
#'
#' number of peaks
#' peak reads
#'
#' @param x path to the directory
#'
#' @export
read_hiseq_peak_stat <- function(x) {
  lapply(x, function(f) {
    if(is_hiseq_dir(f)) {
      pd <- read_hiseq(f)
      # search for: trim_dir
      if(startsWith(pd$hiseq_type, "callpeak")) {
        t <- c("macs2_peak", 'seacr_peak') # peak_dir
      } else {
        t <- c("peak", 'peak_seacr') # pipeline
      }
      t <- purrr::keep(t, function(i) {i %in% names(pd$args)})
      t <- t[1] # first one
      j <- pd$args[[t]]
      n_peak <- tryCatch(error = function(cnd) 0, length(readLines(j)))
      # output
      if(length(j) > 0) {
        data.frame(
          id    = pd$args$smp_name,
          count = n_peak
        )
      }
    }
  }) %>%
    dplyr::bind_rows()
}



#' @describeIn read_hiseq_lendist_stat
#'
#' Gather length distributino of hiseq: r1, rn
#'
#' @param x path to the directory
#'
#' @export
read_hiseq_lendist_stat <- function(x) {
  lapply(x, function(f) {
    if(is_hiseq_dir(f)) {
      pd <- read_hiseq(f)
      # search for: trim_dir
      if(startsWith(pd$hiseq_type, "qc")) {
        t <- c("lendist_csv", "lendist_json") # qc_dir
      } else {
        t <- c("lendist_csv") # pipeline
      }
      t <- purrr::keep(t, function(i) {i %in% names(pd$args)})
      t <- t[1] # first one
      j <- pd$args[[t]]
      # output\
      tryCatch(error = function(cnd) NULL, read_text(j))
    }
  }) %>%
    dplyr::bind_rows()
}



#' @describeIn read_hiseq_frip_stat
#'
#' Fraction in peak
#'
#' @param x path to the directory
#'
#' @export
read_hiseq_frip_stat <- function(x) {
  lapply(x, function(f) {
    if(is_hiseq_dir(f)) {
      pd <- read_hiseq(f)
      # search for: trim_dir
      if(startsWith(pd$hiseq_type, "qc")) {
        t <- c("frip_json", "frip_toml", "frip_txt") # qc_dir
      } else {
        t <- c("frip_json") # pipeline
      }
      t <- purrr::keep(t, function(i) {i %in% names(pd$args)})
      t <- t[1] # first one
      j <- pd$args[[t]]
      # output
      l <- tryCatch(error = function(cnd) NULL, jsonlite::read_json(j))
      if(is(l, "list")) {
        as.data.frame(l)
      }
    }
  }) %>%
    dplyr::bind_rows()
}




#' @describeIn read_hiseq_report_stat
#'
#' Return the report/HiSeq_report.html
#'
#' @param x path to the directory
#'
#' @export
read_hiseq_report_stat <- function(x) {
  sapply(x, function(f) {
    if(is_hiseq_dir(f)) {
      pd <- read_hiseq(f)
      t <- "report_dir"
      j <- ifelse(t %in% names(pd$args), pd$args[[t]], "")
      # for html file
      if(is(j, "character") & file.exists(j)) {
        h <- list.files(j, "*html$", full.names = TRUE)
        if(length(h) > 1) {
          h <- h[1]
          message("More than 1 html detected")
        }
        h
      }
    }
  })
}




#' @describeIn read_hiseq_enrich_stat
#'
#' Return the report/HiSeq_report.html
#'
#' @param x path to the directory
#'
#' @export
read_hiseq_enrich_stat <- function(x) {
  lapply(x, function(f) {
    if(is_hiseq_dir(f)) {
      pd <- read_hiseq(f)
      t <- "qc_dir"
      j <- ifelse(t %in% names(pd$args), pd$args[[t]], "")
      # for html file
      if(is(j, "character") & file.exists(j)) {
        # TSS enrichment
        e_tss      <- list.files(j, "tss_enrich.*png", full.names = TRUE)
        if(length(e_tss) > 1) {
          e_tss <- e_tss[1]
          message("More than 1 tss_enrich png detected")
        }
        # Genebody enrichment
        e_genebody <- list.files(j, "genebody_enrich.*png", full.names = TRUE)
        if(length(e_genebody) > 1) {
          e_genebody <- e_genebody[1]
          message("More than 1 genebody_enrich png detected")
        }
        # output
        list(
          tss      = e_tss,
          genebody = e_genebody
        )
      }
    }
  })
}



#' @describeIn read_hiseq_cor_stat
#'
#' Return the cor heatmap/pca plots
#'
#' @param x path to the directory
#'
#' @export
read_hiseq_cor_stat <- function(x) {
  lapply(x, function(f) {
    if(is_hiseq_dir(f)) {
      pd <- read_hiseq(f)
      t <- "qc_dir"
      j <- ifelse(t %in% names(pd$args), pd$args[[t]], "")
      # for html file
      if(is(j, "character") & file.exists(j)) {
        # cor
        cor_heatmap     <- list.files(j, "cor_heatmap.*png", full.names = TRUE)
        cor_pca         <- list.files(j, "cor_PCA.*png", full.names = TRUE)
        cor_fingerprint <- list.files(j, "fingerprint.png", full.names = TRUE)
        peak_overlap    <- list.files(j, "peak_overlap.png", full.names = TRUE)
        peak_idr        <- list.files(j, "peak_dir.*png", full.names = TRUE)
        # filter
        if(length(cor_heatmap) > 1) cor_heatmap <- cor_heatmap[1]
        if(length(cor_pca) > 1) cor_pca <- cor_pca[1]
        if(length(cor_fingerprint) > 1) cor_fingerprint <- cor_fingerprint[1]
        if(length(peak_overlap) > 1) peak_overlap <- peak_overlap[1]
        list(
          heatmap = cor_heatmap,
          pca     = cor_pca,
          fingerprint = cor_fingerprint,
          overlap = peak_overlap,
          idr     = peak_idr
        )
      }
    }
  })
}


# tss_enrich
# genebody_enrich
# bam_cor: PCA, heatmap
# bam_fingerprint
# peak_overlap
# peak_idr



#' @describeIn read_rnaseq_deseq
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









