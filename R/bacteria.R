#' Functions for bacteria
#'


#' bacteria
#' 
#' @description A tmp function for documentation
bacteria <- function(x) {
  x
}

#' hiseq_kraken2_report
#'
#' @param indir character Path to the of kraken2 results
#' @param outdir character path to the html report
#'
#' @export
hiseq_kraken2_report <- function (indir, outdir, preview = FALSE) {
  indir  <- normalizePath(indir)
  outdir <- normalizePath(outdir)
  report_template <- system.file("qc", "hiseq_kraken2_report.Rmd",
                                 package = "hiseqr")
  if(! dir.exists(outdir)){
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }
  output_html <- file.path(outdir, "kraken2_report.html")
  rmarkdown::render(input = report_template,
                    output_file = output_html,
                    params = list(input_dir = indir))
  if (preview) {
    utils::browseURL(output_html)
  }
}



#' read_kraken2_report
#'
#' @param x character path to the directory of kraken2 results
#' @param tax_level character the level of tax, default; [G]
#' @param topN int the number of taxon to show
#'
#' choose topN taxon across all samples
#' assign the rest as "other"
#'
#' A rank code, indicating (U)nclassified, (R)oot, (D)omain, (K)ingdom,
#' (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. Taxa that are
#' not at any of these 10 ranks have a rank code that is formed by using the
#' rank code of the closest ancestor rank with a number indicating the distance
#' from that rank. E.g., "G2" is a rank code indicating a taxon is between
#' genus and species and the grandparent taxon is at the genus rank
#'
#'
#' @export
read_kraken2_report <- function(x, tax_level = "G", topN = 5) {
  # header: .report
  s  <- c('pct', 'reads_in_clade', 'reads_in_tax', 'code', 'taxid', 'name')
  #-- choose the topN taxid --#
  t <- lapply(x, function(f) {
    df <- readr::read_delim(f, "\t",
                            trim_ws   = TRUE,
                            col_types = readr::cols(),
                            col_names = s)
    # topN taxon, by reads_in_clade
    df %>%
      dplyr::filter(code %in% tax_level) %>%
      dplyr::arrange(desc(reads_in_clade)) %>%
      dplyr::mutate(rank = row_number()) %>%
      dplyr::filter(rank <= topN)
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::arrange(rank) %>%
    dplyr::pull(name) %>%
    unique() %>%
    head(topN)
  #--Choose topN across samples --#
  ts <- c(t, "other", "unclassified")
  # reprocessing files, 0=unclassified
  df2 <- lapply(x, function(f) {
    df <- readr::read_delim(f, "\t", trim_ws = TRUE,
                            col_types = readr::cols(),
                            col_names = s)
    # n_unclassified
    n_un   <- dplyr::filter(df, name == "unclassified") %>%
      dplyr::pull("reads_in_clade")
    # n_root
    n_root <- dplyr::filter(df, name == "root") %>%
      dplyr::pull("reads_in_clade")
    # other
    n_topN <- dplyr::filter(df, name %in% t) %>%
      dplyr::pull("reads_in_clade") %>%
      sum()
    n_other <- n_root - n_topN
    # all reads
    n_total <- n_un + n_root
    df_other <- tibble::tibble(
      pct = n_other / n_root * 100,
      reads_in_clade = n_other,
      reads_in_tax   = n_other,
      code           = "X",
      taxid          = 0,
      name           = "other"
    )
    # combine
    df1 <- dplyr::bind_rows(df, df_other) %>%
      dplyr::filter(name %in% ts) %>%
      dplyr::mutate(sample  = gsub(".kraken2|.report$", "", basename(f)),
                    hit_pct = round(reads_in_clade / (n_un + n_root) * 100, 2),
                    hit     = reads_in_clade,
                    total   = n_un + n_root) %>%
      dplyr::select(sample, name, reads_in_clade, hit_pct, hit, total)
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(name = factor(name, levels = ts))
}



#' plot_kraken2_report
#' Generate heatmap for kraken2 report
#'
#' @param x character, path to the kraken2 output
#'
#' A rank code, indicating (U)nclassified, (R)oot, (D)omain, (K)ingdom,
#' (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. Taxa that are
#' not at any of these 10 ranks have a rank code that is formed by using the
#' rank code of the closest ancestor rank with a number indicating the distance
#' from that rank. E.g., "G2" is a rank code indicating a taxon is between
#'
#'
#' @export
plot_kraken2_report <- function(x, tax_level = "G", topN = 5,
                                switch_axis = FALSE) {
  f_list <- list.files(x, ".kraken2.report", full.names = TRUE)
  df <- read_kraken2_report(f_list, tax_level, topN)
  # heatmap
  if(isTRUE(switch_axis)) {
    p <- df %>%
      ggplot(aes(name, sample, fill = hit_pct))
  } else {
    p <- df %>%
      dplyr::mutate(name = forcats::fct_rev(name)) %>%
      ggplot(aes(sample, name, fill = hit_pct))
  }
  p +
    geom_tile() +
    scale_x_discrete(position = "top") +
    scale_fill_gradient(limits = c(0, 100), low = "grey80", high = "red") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = .5, hjust = 0, size = 10),
      axis.text.y = element_text(color = "black", face = "bold", size = 10),
      axis.title  = element_blank()
    )
}