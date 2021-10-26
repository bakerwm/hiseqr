


#------------------------------------------------------------------------------#

#' @describeIn read_rnaseq_deseq2
#' for upgraded version deseq, including shrink results
#'
#' @param x character path to the rnaseq_rx
#'
#' @export
read_rnaseq_deseq2 <- function(x, shrink_method = "standard",
                               sig_type = "all") {
  if(is_hiseq_dir(x, "rnaseq_rx")) {
    tryCatch(
      {
        if(is_hiseq_dir(x, "rnaseq_rx")) {
          list_hiseq_file(x, "deseq_dir", "_rx") %>%
            deseq_qc_res(shrink_method = shrink_method) %>%
            filt_sig_gene(type = sig_type, force = FALSE)
        }
      },
      error = function(cond) {
        stderr = glue::glue(
          "Failed: reading *.fix.xls file in {x}"
        )
        message(stderr)
        return(NULL)
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




#' @describeIn plot_rnaseq_sig_count for significant genes
#'
#' @param x character path to the rnaseq_rx
#'
#' @export
plot_rnaseq_sig_count2 <- function(x, shrink_method = "standard",
                                   skip_not = TRUE, return_data = FALSE,
                                   fish = "Trimma_lantana") {
  if(is(x, "character")) {
    if(is_hiseq_dir(x, "rnaseq_rx")) {
      df <- read_rnaseq_deseq2(x, shrink_method = shrink_method)
      if(inherits(df, "data.frame")) {
        if(nrow(df) > 0) {
          df2 <- as.data.frame(table(df[["sig"]])) %>%
            dplyr::mutate(shrink = shrink_method)
          if(skip_not) {
            df2 <- df2 %>%
              dplyr::filter(!Var1 == "not")
          }
          p <- df2 %>%
            hiseqr::bar_plot(x = "Freq",  y = "Var1", group = "shrink",
                             label = "Freq",
                             fill = "Var1", direction = "horizontal") +
            fishualize::scale_fill_fish(discrete = TRUE, option = fish) +
            ggplot2::ggtitle(glue::glue(
              "Number of DE genes, (shrink = {shrink_method})"
            )) +
            ggplot2::theme_bw()
          if(return_data) {
            df2
          } else {
            p
          }
        }
      }
    }
  } else {
    warning(glue::glue("illegal, expect: rnaseq_rx dir; {x}"))
  }
}



#------------------------------------------------------------------------------#

#' @describeIn plot_hiseq_trim Create bar_plot for trim stat
#'
#' output from get_rnaseq_trim_stat(),
#' including columns:
#' id, raw, clean, clean_pct, short_pct
#'
#'
#' @param x data.frame From get_trim_stat
#' @param fish string Name of the fish, use `fishualize::fish_palettes()` list
#' @param direction string ["horizontal", "vertical"]
#' @param position string ["stack", "fill"]
#' all available fish names
#'
#' @export
plot_hiseq_trim <- function(x, fish = "Trimma_lantana",
                            direction = "horizontal", position = "fill") {
  if(is(x, "data.frame")) {
    col_required <- c("name", "input", "output")
    if(! all(rlang::has_name(x, col_required))) {
      warning(glue::glue("failed, missing columns {read_hiseq_trim}"))
      return(NULL)
    }
    df <- x
  } else if(is_hiseq_dir(x, "trim")) {
    df <- read_hiseq_trim(x) # for single trim
  } else if(is_hiseq_dir(x, TRUE)) {
    df <- read_hiseq_trim_stat(x) # for pipeline
  } else {
    warning("`data` expect data.frame, failed")
    return(NULL)
  }
  df %>%
    # dplyr::mutate(clean_pct = round(output / input * 100, 1),
    #               short_pct = round(100 - clean_pct, 1)) %>%
    dplyr::mutate(clean = round(output / 1e6, 2),
                  short = round((input - output) / 1e6, 2)) %>%
    dplyr::select(name, clean, short) %>%
    tidyr::pivot_longer(c(clean, short),
                        names_to  = "group",
                        values_to = "count") %>%
    bar_plot(x = "count", y = "name", fill = "group", label = NA,
             direction = direction, position = position) +
    fishualize::scale_fill_fish(discrete = TRUE, option = fish) +
    ggtitle("Trim reads") +
    theme(legend.position = "top")
}



#' @describeIn plot_hiseq_trim2
#' for all available columns
#' Create bar_plot for trim stat
#'
#' output from get_rnaseq_trim_stat(),
#' including columns:
#' id, raw, clean, clean_pct, short_pct
#'
#' @param x data.frame From get_trim_stat
#' @param fish string Name of the fish, use `fishualize::fish_palettes()` list
#' @param direction string ["horizontal", "vertical"]
#' @param position string ["stack", "fill"]
#' all available fish names
#'
#' @export
plot_hiseq_trim2 <- function(x, fish = "Trimma_lantana",
                             direction = "horizontal", position = "fill") {
  if(is(x, "data.frame")) {
    col_required <- c("name", "input", "output")
    if(! all(rlang::has_name(x, col_required))) {
      warning(glue::glue("failed, missing columns {read_hiseq_trim}"))
      return(NULL)
    }
    df <- x
  } else if(is_hiseq_dir(x, "trim")) {
    df <- read_hiseq_trim(x) # for single trim
  } else if(is_hiseq_dir(x, TRUE)) {
    df <- read_hiseq_trim_stat(x) # for pipeline
  } else {
    warning("`data` expect data.frame, failed")
    return(NULL)
  }
  df %>%
    tidyr::pivot_longer(-c(name, input, percent),
                        names_to  = "group",
                        values_to = "count") %>%
    dplyr::group_by(name) %>%
    # dplyr::mutate(pct = round(count / sum(count) * 100, 1),
    #               group = forcats::fct_relevel(group, "output")) %>%
    dplyr::mutate(count = round(count / 1e6, 2),
                  group = forcats::fct_relevel(group, "output")) %>%
    bar_plot(x = "count", y = "name", fill = "group", label = NA,
             direction = direction, position = position) +
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
    df <- read_hiseq_align_stat(x) # pipeline
    if(is.null(df)) {
      df <- read_hiseq_align(x) # alignment only
    }
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
    facet_wrap(.~sample, ncol = 2) +
    theme(
      legend.position = "none"
    )
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













