




#' @describeIn plot_hiseq_trim Create bar_plot for trim stat
#'
#' output from get_rnaseq_trim_stat(),
#' including columns:
#' id, raw, clean, clean_pct, short_pct
#'
#'
#' @param x data.frame From get_trim_stat
#' @param fish string Name of the fish, use `fishualize::fish_palettes()` list
#' all available fish names
#'
#' @export
plot_hiseq_trim <- function(x, fish = "Trimma_lantana") {
  if(is(x, "data.frame")) {
    col_required <- c("name", "input", "output")
    if(! all(col_required %in% names(data))) {
      on.exit("`data` required columns missing: ",
              paste(col_required, collapse = ", "))
    }
    df <- x
  } else if(is(x, "character")) {
    df <- read_hiseq_trim_stat(x)
  } else {
    on.exit("`data` expect data.frame, failed")
  }
  df %>%
    dplyr::select(name, out_pct, rm_pct) %>%
    dplyr::rename(clean_pct = out_pct,
                  short_pct = rm_pct) %>%
    tidyr::pivot_longer(names_to  = "group",
                        values_to = "count",
                        c(clean_pct, short_pct)) %>%
    bar_plot(x = "count", y = "name", fill = "group", label = "count",
             direction = "horizontal") +
    fishualize::scale_fill_fish(discrete = TRUE, option = fish) +
    ggtitle("Trim reads") +
    theme(legend.position = "top")
}




#' @describeIn plot_rnaseq_align Create bar_plot for align_stat
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
                             mode = 1,
                             columns = NULL,
                             title = "Alignment",
                             add_label = FALSE,
                             position = "fill",
                             fish = "Trimma_lantana",
                             ...) {
  if(is(x, "data.frame")) {
    df <- x
  } else if(is(x, "character")) {
    df <- read_hiseq_align_stat(x)
  } else {
    on.exit("`data` expect data.frame, failed")
  }
  #--pick columns--#
  # id
  i <- c("fqname", "id", "name") # different version
  i <- purrr::keep(i, function(a) a %in% colnames(df))
  # columns: mode,columns
  g <- switch(mode,
              c("map", "unmap"),
              c("unique", "multi", "unmap"),
              c("chrM", "spikein", "map", "unmap"),
              c("map", "rRNA", "spikein", "unmap"),
              c("total"))
  if(is.null(columns)) {
    columns <- g
  }
  # validate:
  if(! all(columns %in% colnames(df))) {
    on.exit("`columns`, `mode` required, unknown column found.")
  }
  # variables
  axis_title_x <- ifelse(
    position == "fill",
    "Percentage%",
    "Count")
  x_col <- ifelse(position == "fill", "pct", "count")
  label_col <- ifelse(isTRUE(add_label), x_col, NA)
  # subset data.frame
  df2 <- df %>%
    dplyr::rename(name = !!i) %>%
    dplyr::select(all_of(c("name", columns))) %>%
    tidyr::pivot_longer(names_to = "group", values_to = "count", -1) %>%
    dplyr::mutate(group = factor(group, levels = columns)) %>%
    dplyr::group_by(name) %>%
    dplyr::mutate(pct = round(count / sum(count) * 100, 2))
  # plot
  df2 %>%
    bar_plot(x = x_col, y = "name", direction = "horizontal",
             fill = "group", label = label_col, ...) +
    fishualize::scale_fill_fish(discrete = TRUE, option = fish) +
    ggtitle(title) +
    ylab(NULL) +
    xlab(axis_title_x) +
    theme(legend.position = "top")
}




#' @describeIn plot_rnaseq_featureCounts
#'
#' @param x character path to RNAseq r1, rn
#'
#' @export
plot_rnaseq_featureCounts <- function(x, fish = "Trimma_lantana") {
  # sense
  c1 <- list_hiseq_file(x, "count_sens", "r1")
  if(is(c1, "character")) {
    c1s <- paste0(c1, ".summary")
    dc1 <- lapply(c1s, read_fc_summary) %>%
      dplyr::bind_rows() %>%
      dplyr::mutate(strand = "sense")
  } else {
    dc1 <- NULL
  }
  ## anti
  c2 <- list_hiseq_file(x, "count_anti", "r1")
  if(is(c2, "character")) {
    c2s <- paste0(c2, ".summary")
    dc2 <- lapply(c2s, read_fc_summary) %>%
      dplyr::bind_rows() %>%
      dplyr::mutate(strand = "anti")
  } else {
    dc2 <- NULL
  }
  # combine
  dc <- dplyr::bind_rows(dc1, dc2)
  if(nrow(dc) > 0) {
    df <- dc %>%
      dplyr::filter(count > 0) %>%
      dplyr::group_by(sample, strand) %>%
      dplyr::mutate(total = sum(count)) %>%
      dplyr::filter(Status == "Assigned") %>%
      tidyr::pivot_wider(names_from = "strand", values_from = "count") %>%
      dplyr::mutate(na = total - sense - anti) %>%
      tidyr::pivot_longer(names_to = "strand", values_to = "count",
                          sense:na) %>%
      dplyr::mutate(strand = factor(strand,
                                    levels = c("sense", "anti", "na")))

    # plot
    hiseqr::bar_plot(df, x = "count", y = "sample", label = NA,
                     fill = "strand", direction = "horizontal",
                     position = "fill") +
      fishualize::scale_fill_fish(discrete = TRUE, option = fish) +
      theme(legend.position = "top")
  }
}




#' @describeIn plot_rnaseq_sig_count for significant genes
#'
#' @param x character path to the rnaseq_rx
#'
#' @export
plot_rnaseq_sig_count <- function(x, fish = "Trimma_lantana") {
  if(is(x, "character")) {
    if(is_hiseq_dir(x, "rnaseq_rx")) {
      fix_xls <- get_fix_xls(x)
    } else if(endsWith(x, ".fix.xls")) {
      fix_xls <- x
    } else {
      on.exit(glue::glue("illegal: {x}"))
    }
    df <- readr::read_delim(fix_xls, "\t", col_types = readr::cols())
    if(nrow(df) > 0) {
      df %>%
        dplyr::group_by(sig) %>%
        dplyr::summarize(count = n()) %>%
        hiseqr::bar_plot(x = "sig", y = "count", fill = "sig", label = "count") +
        fishualize::scale_fill_fish(discrete = TRUE, option = fish) +
        ggplot2::ggtitle("Number of DE genes") +
        ggplot2::theme_bw()
    }
  } else {
    warning(glue::glue("illegal, expect: .fix.xls or rnaseq_rx dir; {x}"))
  }
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












#' #------------------------------------------------------------------------------#
#' # figures for summary
#' #
#' #' @describeIn plot_hiseq_trim Create bar_plot for trim stat
#' #'
#' #' output from get_rnaseq_trim_stat(),
#' #' including columns:
#' #' id, raw, clean, clean_pct, short_pct
#' #'
#' #'
#' #' @param data data.frame From get_trim_stat
#' #' @param fish string Name of the fish, use `fishualize::fish_palettes()` list
#' #' all available fish names
#' #'
#' #' @export
#' plot_hiseq_trim <- function(x,
#'                             hiseq_type = "r1",
#'                             position = "fill",
#'                             label = NULL,
#'                             fish = "Trimma_lantana") {
#'   if(is(x, "data.frame")) {
#'     df <- x
#'   } else if(is(x, "character")) {
#'     # check input dirs
#'     x_dirs <- list_hiseq_dir(x, hiseq_type)
#'     df <- tryCatch(
#'       error = function(cnd) NULL,
#'       read_hiseq_trim_stat(x_dirs))
#'   } else {
#'     df <- data.frame()
#'   }
#'   # input: total, clean
#'   # output: pct
#'   if(is.null(df)) {
#'     warning("Failed reading data")
#'     return(NULL)
#'   } else if(nrow(df) == 0) {
#'     warning("No trim data detected")
#'     return(NULL)
#'   }
#'   # processing data.frame
#'   t_cols <- c("name", "total", "clean")
#'   if(!all(t_cols %in% names(df))) {
#'     warning("missing columns [name, total, clean]")
#'     return(NULL)
#'   }
#'   df2 <- df %>%
#'     dplyr::mutate(too_short = total - clean) %>%
#'     dplyr::select(name, total, clean, too_short) %>%
#'     tidyr::pivot_longer(
#'       names_to = "group",
#'       values_to = "count",
#'       clean:too_short) %>%
#'     dplyr::group_by(name) %>%
#'     dplyr::mutate(pct = round(count / total * 100, 1),
#'                   count = round(count / 1e6, 1))
#'   # plot
#'   if(position == "fill") {
#'     label <- "pct"
#'   } else {
#'     label <- NA
#'   }
#'   # position: fill/identity
#'   if(position == "fill") {
#'     p <- df2 %>%
#'       hiseqr::bar_plot(x = "pct", y = "name",
#'                        fill = "group",
#'                        label = label,
#'                        direction = "horizontal") +
#'       xlab("Percentage of reads (%)")
#'   } else {
#'     p <- df2 %>%
#'       hiseqr::bar_plot(x = "count", y = "name",
#'                        fill = "group",
#'                        label = "count",
#'                        position = "identity",
#'                        direction = "horizontal") +
#'       xlab("Number of reads (M)")
#'   }
#'   p +
#'     fishualize::scale_fill_fish(discrete = TRUE, option = fish) +
#'     ggtitle("Trim reads") +
#'     theme(legend.position = "top")
#' }




#' #' @describeIn plot_hiseq_align Create bar_plot for align_stat
#' #'
#' #' Output from get_rnaseq_align_stat()
#' #' including columns:
#' #' fqname, map, unique, multiple
#' #'
#' #' @param data data.frame From get_rnaseq_align_stat()
#' #' @param mode integer map=1, unique/multiple=2, default: 1
#' #' @param fish string Name of the fish, use `fishualize::fish_palettes()` list
#' #'  all availabel fish names
#' #'
#' #' @export
#' plot_hiseq_align <- function(x,
#'                              mode = 0,
#'                              columns = NA,
#'                              add_label = FALSE,
#'                              position = "fill",
#'                              hiseq_type = "r1",
#'                              title = "Alignment",
#'                              fish = "Trimma_lantana", ...) {
#'   if(is(x, "data.frame")) {
#'     df <- x
#'   } else if(is(x, "character")) {
#'     # check input dirs
#'     x_dirs <- list_hiseq_dir(x, hiseq_type)
#'     df <- tryCatch(
#'       error = function(cnd) NULL,
#'       read_hiseq_align_stat(x_dirs))
#'   } else {
#'     df <- data.frame()
#'   }
#'   # input: total, clean
#'   # output: pct
#'   if(is.null(df)) {
#'     warning("Failed reading data")
#'     return(NULL)
#'   } else if(nrow(df) == 0) {
#'     warning("No trim data detected")
#'     return(NULL)
#'   }
#'   # choose columns
#'   g1 <- switch(
#'     mode,
#'     c("map", "unmap"),
#'     c("unique", "multi", "unmap"),
#'     c("chrM", "map", "unmap"),
#'     c("total"))
#'   # choose columns
#'   if(length(g1) > 0) {
#'     t_cols <- g1
#'   } else if(is(columns, "character")) {
#'     t_cols <- columns
#'   } else {
#'     t_cols <- c("unique", "multi", "unmap")
#'   }
#'   if(! all(t_cols %in% names(df))) {
#'     warning("missing columns [name, total, map]")
#'     return(NULL)
#'   }
#'   t_cols <- purrr::discard(t_cols, is.null)
#'   if(! all(t_cols %in% names(df))) {
#'     t_line <- paste(names(df), collapse = ", ")
#'     stop("`columns` missing, expect: ", t_line)
#'   }
#'   # choose id
#'   n_cols <- c("id", "name", "fqname")
#'   n_cols <- purrr::keep(n_cols, function(i) i %in% names(df))
#'   if(length(n_cols) == 0) {
#'     stop("name column not found: [id, name, fqname]")
#'   }
#'   name_col <- n_cols[1] # first one
#'   # subset data.frame
#'   df1 <- df %>%
#'     dplyr::rename(id = !!name_col) %>%
#'     dplyr::select(all_of(c("id", "total", t_cols))) %>%
#'     tidyr::pivot_longer(
#'       names_to = "group",
#'       values_to = "count",
#'       -c(id, total)) %>%
#'     dplyr::mutate(group = factor(group, levels = t_cols)) %>%
#'     dplyr::group_by(id) %>%
#'     dplyr::mutate(pct = round(count / sum(count) * 100, 2))
#'   # plot
#'   label_col <- ifelse(isTRUE(add_label), "count", NA)
#'   df1 %>%
#'     hiseqr::bar_plot(
#'       x = "count", y = "id",
#'       direction = "horizontal",
#'       fill = "group",
#'       label = label_col,
#'       position = position) +
#'     fishualize::scale_fill_fish(discrete = TRUE, option = fish) +
#'     ylab(NULL) +
#'     xlab(NULL) +
#'     # xlab("Percentage%") +
#'     theme(legend.position = "top") +
#'     ggtitle(title)
#' }




