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
    fn(x)
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
  lapply(x, function(f) {
    if(is_hiseq_dir(f, "r1")) {
      pd <- read_hiseq(f)
      # search for: trim_dir
      if(startsWith(pd$hiseq_type, "trim")) {
        j <- pd$args$trim_stat_json # trim_json
      } else {
        # pipeline:
        c_dir <- file.path(pd$args$clean_dir, pd$args$smp_name)
        pc    <- read_hiseq(c_dir)
        j     <- pc$args$trim_json
      }
      if(length(j) > 0) {
        # read data
        jsonlite::read_json(j) %>%
          as.data.frame()
      }
    }
  }) %>%
    dplyr::bind_rows()
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
  lapply(x, function(f) {
    if(is_hiseq_dir(f, "r1")) {
      pd <- read_hiseq(f)
      # search for: align_dir
      if(startsWith(pd$hiseq_type, "align")) {
        j <- pd$args$align_json # align_json
      } else {
        t <- c("align_summary_json", 'align_json')
        t <- purrr::keep(t, function(i) {i %in% names(pd$args)})
        t <- t[1] # first one
        j <- pd$args[[t]]
        ## to-be-continue
        # } else {
        #   # pipeline:
        #   c_dir <- file.path(pd$args$clean_dir, pd$args$smp_name)
        #   pc    <- read_hiseq(c_dir)
        #   j     <- pc$args$trim_json
      }
      j <- purrr::keep(j, file.exists)
      if(length(j) > 0) {
        # read data
        # re-arrange columns
        df <- jsonlite::read_json(j) %>%
          as.data.frame()
        # column order
        t_cols <- c("name", "total", "map", "unique", "multi", "unmap", "nodup", "chrM")
        dplyr::select(df, all_of(t_cols))
      }
    }
  }) %>%
    dplyr::bind_rows()
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
        t <- c("lendist_txt", "lendist_json") # qc_dir
      } else {
        t <- c("lendist_txt") # pipeline
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


#------------------------------------------------------------------------------#
# figures for summary
#
#' @describeIn plot_hiseq_trim Create bar_plot for trim stat
#'
#' output from get_rnaseq_trim_stat(),
#' including columns:
#' id, raw, clean, clean_pct, short_pct
#'
#'
#' @param data data.frame From get_trim_stat
#' @param fish string Name of the fish, use `fishualize::fish_palettes()` list
#' all available fish names
#'
#' @export
plot_hiseq_trim <- function(x,
                            hiseq_type = "r1",
                            position = "fill",
                            label = NULL,
                            fish = "Trimma_lantana") {
  if(is(x, "data.frame")) {
    df <- x
  } else if(is(x, "character")) {
    # check input dirs
    x_dirs <- list_hiseq_dir(x, hiseq_type)
    df <- tryCatch(
      error = function(cnd) NULL,
      read_hiseq_trim_stat(x_dirs))
  } else {
    df <- data.frame()
  }
  # input: total, clean
  # output: pct
  if(is.null(df)) {
    warning("Failed reading data")
    return(NULL)
  } else if(nrow(df) == 0) {
    warning("No trim data detected")
    return(NULL)
  }
  # processing data.frame
  t_cols <- c("name", "total", "clean")
  if(!all(t_cols %in% names(df))) {
    warning("missing columns [name, total, clean]")
    return(NULL)
  }
  df2 <- df %>%
    dplyr::mutate(too_short = total - clean) %>%
    dplyr::select(name, total, clean, too_short) %>%
    tidyr::pivot_longer(
      names_to = "group",
      values_to = "count",
      clean:too_short) %>%
    dplyr::group_by(name) %>%
    dplyr::mutate(pct = round(count / total * 100, 1),
                  count = round(count / 1e6, 1))
  # plot
  if(position == "fill") {
    label <- "pct"
  } else {
    label <- NA
  }
  # position: fill/identity
  if(position == "fill") {
    p <- df2 %>%
      hiseqr::bar_plot(x = "pct", y = "name",
                       fill = "group",
                       label = label,
                       direction = "horizontal") +
      xlab("Percentage of reads (%)")
  } else {
    p <- df2 %>%
      hiseqr::bar_plot(x = "count", y = "name",
                       fill = "group",
                       label = "count",
                       position = "identity",
                       direction = "horizontal") +
      xlab("Number of reads (M)")
  }
  p +
    fishualize::scale_fill_fish(discrete = TRUE, option = fish) +
    ggtitle("Trim reads") +
    theme(legend.position = "top")
}




#' @describeIn plot_hiseq_align Create bar_plot for align_stat
#'
#' Output from get_rnaseq_align_stat()
#' including columns:
#' fqname, map, unique, multiple
#'
#' @param data data.frame From get_rnaseq_align_stat()
#' @param mode integer map=1, unique/multiple=2, default: 1
#' @param fish string Name of the fish, use `fishualize::fish_palettes()` list
#'  all availabel fish names
#'
#' @export
plot_hiseq_align <- function(x,
                             mode = 0,
                             columns = NA,
                             add_label = FALSE,
                             position = "fill",
                             hiseq_type = "r1",
                             title = "Alignment",
                             fish = "Trimma_lantana", ...) {
  if(is(x, "data.frame")) {
    df <- x
  } else if(is(x, "character")) {
    # check input dirs
    x_dirs <- list_hiseq_dir(x, hiseq_type)
    df <- tryCatch(
      error = function(cnd) NULL,
      read_hiseq_align_stat(x_dirs))
  } else {
    df <- data.frame()
  }
  # input: total, clean
  # output: pct
  if(is.null(df)) {
    warning("Failed reading data")
    return(NULL)
  } else if(nrow(df) == 0) {
    warning("No trim data detected")
    return(NULL)
  }
  # choose columns
  g1 <- switch(
    mode,
    c("map", "unmap"),
    c("unique", "multi", "unmap"),
    c("chrM", "map", "unmap"),
    c("total"))
  # choose columns
  if(length(g1) > 0) {
    t_cols <- g1
  } else if(is(columns, "character")) {
    t_cols <- columns
  } else {
    t_cols <- c("unique", "multi", "unmap")
  }
  if(! all(t_cols %in% names(df))) {
    warning("missing columns [name, total, map]")
    return(NULL)
  }
  t_cols <- purrr::discard(t_cols, is.null)
  if(! all(t_cols %in% names(df))) {
    t_line <- paste(names(df), collapse = ", ")
    stop("`columns` missing, expect: ", t_line)
  }
  # choose id
  n_cols <- c("id", "name", "fqname")
  n_cols <- purrr::keep(n_cols, function(i) i %in% names(df))
  if(length(n_cols) == 0) {
    stop("name column not found: [id, name, fqname]")
  }
  name_col <- n_cols[1] # first one
  # subset data.frame
  df1 <- df %>%
    dplyr::rename(id = !!name_col) %>%
    dplyr::select(all_of(c("id", "total", t_cols))) %>%
    tidyr::pivot_longer(
      names_to = "group",
      values_to = "count",
      -c(id, total)) %>%
    dplyr::mutate(group = factor(group, levels = t_cols)) %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(pct = round(count / sum(count) * 100, 2))
  # plot
  label_col <- ifelse(isTRUE(add_label), "count", NA)
  df1 %>%
    hiseqr::bar_plot(
      x = "count", y = "id",
      direction = "horizontal",
      fill = "group",
      label = label_col,
      position = position) +
    fishualize::scale_fill_fish(discrete = TRUE, option = fish) +
    ylab(NULL) +
    xlab(NULL) +
    # xlab("Percentage%") +
    theme(legend.position = "top") +
    ggtitle(title)
}








#' @describeIn plot_hiseq_peak Create bar_plot for align_stat
#'
#' Output from get_rnaseq_align_stat()
#' including columns:
#' fqname, map, unique, multiple
#'
#' @param data data.frame From get_rnaseq_align_stat()
#' @param mode integer map=1, unique/multiple=2, default: 1
#' @param fish string Name of the fish, use `fishualize::fish_palettes()` list
#'  all availabel fish names
#'
#' @export
plot_hiseq_peak <- function(x,
                             hiseq_type = "r1",
                             title = "No. of peaks",
                             fish = "Trimma_lantana") {
  if(is(x, "data.frame")) {
    df <- x
  } else if(is(x, "character")) {
    # check input dirs
    x_dirs <- list_hiseq_dir(x, hiseq_type)
    df <- tryCatch(
      error = function(cnd) NULL,
      read_hiseq_peak_stat(x_dirs))
  } else {
    df <- data.frame()
  }
  # input: total, clean
  # output: pct
  if(nrow(df) == 0) {
    warning("No peak data detected")
    return(NULL)
  }
  # subset data.frame
  df %>%
    dplyr::mutate(id = id) %>%
    bar_plot(x = "count", y = "id",
             direction = "horizontal",
             label = "count") +
    ggtitle(title)
}




#' @describeIn plot_hiseq_peak Create bar_plot for align_stat
#'
#' Output from get_rnaseq_align_stat()
#' including columns:
#' fqname, map, unique, multiple
#'
#' @param data data.frame From get_rnaseq_align_stat()
#' @param mode integer map=1, unique/multiple=2, default: 1
#' @param fish string Name of the fish, use `fishualize::fish_palettes()` list
#'  all availabel fish names
#'
#' @export
plot_hiseq_lendist <- function(x,
                            hiseq_type = TRUE,
                            title = "No. of peaks",
                            fish = "Trimma_lantana") {
  if(is(x, "data.frame")) {
    df <- x
  } else if(is(x, "character")) {
    # check input dirs
    x_dirs <- list_hiseq_dir(x, hiseq_type)
    df <- tryCatch(
      error = function(cnd) NULL,
      read_hiseq_stat(x_dirs, "lendist")) %>%
      dplyr::bind_rows()
  } else {
    df <- data.frame()
  }
  # input: total, clean
  # output: pct
  if(nrow(df) == 0) {
    warning("No peak data detected")
    return(NULL)
  }
  # group by sample
  df %>%
    dplyr::mutate(sample = gsub(".rmdup$|_rep\\d", "", id)) %>%
    hiseqr::fragsize_plot() +
    facet_wrap(.~sample, ncol = 2)
}







#' @describeIn plot_rnaseq_sig_count for significant genes
#'
#' @param x data.frame From deseq_dir, fix.xls
#'
#' @export
plot_rnaseq_sig_count <- function(x, fish = "Trimma_lantana") {
  if(is(x, "data.frame")) {
    df <- x
  } else if(is(x, "character")) {
    pd <- read_hiseq(x)
    if(is(pd, "list")) {
      #!!!! to-do !!!!#
      df <- data.frame()
    } else {
      df <- data.frame()
    }
  } else {
    df <- data.frame()
  }
  # check input
  if(nrow(df) == 0) {
    stop("not rnaseq_rx dir")
  }
  # data %>%
  #   group_by(sig) %>%
  #   summarize(count = n()) %>%
  #   bar_plot(x = "sig", y = "count", fill = "sig", label = "count") +
  #   fishualize::scale_fill_fish(discrete = TRUE, option = fish) +
  #   ggtitle("DE genes")
}





#' @describeIn plot_hiseq_bam_cor
#'
#' generate bam cor plots for hiseq dirs
#' @param x path to the directory
#'
#' @export
plot_hiseq_bam_cor <- function(x) {
  xdirs <- list_hiseq_dir(x, "rn")
  if(length(xdirs) == 0) {
    return(NULL)
  }
  p_list <- lapply(xdirs, function(f) {
    # for each rn dirs
    f_qc_dir <- list_hiseq_file(f, "qc_dir", "rn")
    # title
    smp_name <- list_hiseq_file(f, "smp_name", "rn")
    # search for cor.matrix
    m <- list.files(f_qc_dir, "06.*cor.matrix", full.names = TRUE)
    m <- purrr::keep(m, file.exists)
    if(length(m) > 0) {
      ma <- readr::read_delim(m[1], "\t", comment = "#", trim_ws = TRUE,
                              col_types = readr::cols()) %>%
        tibble::column_to_rownames("X1") %>%
        as.matrix()
      rownames(ma) <- .fix_prefix(rownames(ma), "r") # fix prefix
      colnames(ma) <- .fix_prefix(colnames(ma), "r") # fix suffix
    } else {
      df <- NULL
    }
    # wrap long title
    t <- gsub("_", "-", smp_name) %>%
      stringi::stri_wrap(width = 10, whitespace_only = FALSE) %>%
      paste(collapse = "\n")
    # generate plot
    ggcorrplot::ggcorrplot(ma, type = "lower", lab = TRUE) +
      ggtitle(t)
  })
  patchwork::wrap_plots(p_list) +
    patchwork::plot_layout(guides = "collect", ncol = 2)
}







#' remove the longest common prefix from strings
#'
#' RNAseq_wt_rep1, RNAseq_wt_rep2 => r1, r2
#'
.fix_prefix <- function(x, replacement = "r") {
  # longest common prefix/suffix
  lcp <- Biobase::lcPrefix(x, ignore.case = TRUE)
  lcs <- Biobase::lcSuffix(x, ignore.case = TRUE)
  x <- gsub(lcp, replacement, x, ignore.case = TRUE)
  x <- gsub(lcs, "", x, ignore.case = TRUE)
  as.character(x)
}






