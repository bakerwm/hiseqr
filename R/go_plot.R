#' Functions for generating plots for GO, KEGG analysis
#'
#'
#'
#' Further plots are compatiable, with (gene_list, ...)
#'
#'
#'
#' force "entrezid" in analysis pipeline
#' force "readable=T" for output
#'
#' @name go_plot



#' create plots
#' @param x list of objects of groupGO, enrichGO
#' @param ..., support, fold_change, text_width,
#' fold_change
#' text_width
#' show_category
#'
#'
#' @export
go_barplot <- function(x, ...) {
  #--Default values: BEGIN
  dots <- rlang::list2(
    show_category = 12,
    text_width    = 40,
  )
  dots <- purrr::list_modify(dots, !!!list(...))
  #--Default values: END
  if(is_go_result(x)) {
  barplot(x, dorp = TRUE, showCategory = dots$show_category, order = TRUE, ...) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = dots$text_width)) +
    xlab("Number of genes") +
    ggtitle(x@ontology) +
    theme(plot.title = element_text(hjust = 0.5))
  } else {
    warning("`x` not groupGOResult, go_barplot() skipped...")
    NULL
  }
}







#' create plots
#' @param x object of enrichGO
#' @param ..., support, fold_change, text_width,
#' fold_change
#' text_width
#' show_category
#'
#' @export
go_dotplot <- function(x, ...) {
  #--Default values: BEGIN
  dots <- rlang::list2(
    show_category = 12,
    text_width    = 40,
  )
  dots <- purrr::list_modify(dots, !!!list(...))
  #--Default values: END
  if(is_go_result(x)) {
    clusterProfiler::dotplot(x, showCategory = dots$show_category, ...) +
      scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = dots$text_width)) +
      xlab("Gene Ratio") +
      ggtitle(x@ontology) +
      theme(axis.text.y = element_text(size  = 10),
            plot.title  = element_text(hjust = 0.5))
  } else {
    warning("`x` not enrichGOResult, go_dotplot() skipped...")
    NULL
  }
}










#' create plots
#' @param x object of enrichGO
#' @param ..., support, fold_change, text_width,
#' fold_change
#' text_width
#' show_category
#'
#' @export
go_cnetplot <- function(x, ...) {
  #--Default values: BEGIN
  dots <- rlang::list2(
    fold_change   = NULL,
    show_category = 12,
    layout        = "nicely"
  )
  dots <- purrr::list_modify(dots, !!!list(...))
  #--Default values: END
  if(is_go_result(x)) {
    clusterProfiler::cnetplot(x,
                              showCategory = dots$show_category,
                              foldChange   = dots$fold_change,
                              layout       = dots$layout,
                              categorySize = "pvalue") +
      ggtitle(x@ontology) +
      theme(text = element_text(size = 8),
            plot.title = element_text(hjust = .5))
  } else {
    warning("`x` not enrichGOResult, go_cnetplot() skipped...")
    NULL
  }
}








#' create plots
#' @param x object of enrichGO
#' @param ..., support, fold_change, text_width,
#'
#' orgdb
#' show_category
#' layout, kk, nicely,
#'
#' updated: 2020-12-15, pairwise_termsim(ego)
#'
#'
#' @export
go_emapplot <- function(x, ...) {
  #--Default values: BEGIN
  dots <- rlang::list2(
    orgdb         = NULL,
    fold_change   = NULL,
    show_category = 12,
    layout        = "nicely" # kk
  )
  dots <- purrr::list_modify(dots, !!!list(...))
  #--Default values: END
  if(is_go_result(x)) {
    if(nrow(x) > 1) {
      ont <- tryCatch(
        x@ontology,
        error = function(cnd) "NULL"
      )
      if(is(dots$orgdb, "OrgDb") && ont %in% c("BP", "CC", "MF")) {
        d  <- GOSemSim::godata(dots$orgdb, ont = ont)
        x2 <- enrichplot::pairwise_termsim(x, method = "Wang", semData = d)
      } else {
        x2 <- enrichplot::pairwise_termsim(x)
      }
      p1 <- enrichplot::emapplot(x2,
                                 showCategory = dots$show_category,
                                 layout       = dots$layout,
                                 cex_category = .8,
                                 cex_line     = .4) +
        ggtitle(ont) +
        theme(text = element_text(size = 10))
      p2 <- enrichplot::emapplot_cluster(x2,
                                         showCategory = dots$show_category,
                                         layout       = dots$layout,
                                         cex_category = .8,
                                         cex_line     = .4) +
        ggtitle(ont) +
        theme(text = element_text(size = 10))
      list(
        plot         = p1,
        plot_cluster = p2
      )
    }
  } else {
    warning("`x` not enrichGOResult, go_emapplot() skipped...")
    NULL
  }
}





#' create plots
#' @param x object of enrichGO
#'
#' updated: 2020-12-15, pairwise_termsim(ego)
#'
#'
#' @export
go_emapplot_cluster <- function(x, ...) {
  dots <- rlang::list2(
    layout = "nicely"
  )
  if(is_go_result(x)) {
    if(nrow(x) > 1) {
      x2 <- enrichplot::pairwise_termsim(x) # updated 2020-12-15
      enrichplot::emapplot_cluster(x2,
                                   cex_category = .8,
                                   cex_line     = .4,
                                   layout       = dots$layout) +
        ggtitle(x@ontology) +
        theme(text = element_text(size = 10))
    }
  } else {
    warning("`x` not enrichGOResult, go_emapplot_cluster() skipped...")
    NULL
  }
}









#' create plots
#' @param x object of gseaResult
#' @param ..., support, fold_change, text_width,
#' fold_change
#' text_width
#' show_category
#'
#'
#' @export
go_gsea_plot <- function(x, ...) {
  if(class(x) == "gseaResult") {
    if(nrow(x) > 0) {
      p_list <- lapply(seq_len(nrow(x)), function(i){
        enrichplot::gseaplot2(x,
                              geneSetID    = i,
                              title        = x$Description[i],
                              pvalue_table = TRUE)
      })
      names(p_list) <- seq_len(nrow(x)) # paste0("gsea.", seq_len(nrow(x)))
      p_list
    }
  } else {
    warning("`x` not gseaGOResult, go_gsea_plot() skipped...")
    NULL
  }
}













#' create plots
#' @param x object of enrichGO
#' @param ..., support, fold_change, text_width,
#' fold_change
#' text_width
#' show_category
#'
#'
#' @export
go_heatplot <- function(x, ...) {
  #--Default values: BEGIN
  dots <- rlang::list2(
    fold_change   = NULL,
    show_category = 12
  )
  dots <- purrr::list_modify(dots, !!!list(...))
  #--Default values: END
  if(is_go_result(x)) {
    enrichplot::heatplot(x,
                         foldChange   = dots$fold_change,
                         showCategory = dots$show_category) +
      ggtitle(x@ontology) +
      theme(text = element_text(size = 10))
  } else {
    warning("`x` not enrichGOResult, go_heatplot() skipped...")
    NULL
  }
}










#' wego plot
#'
#' @param x list of objects of groupGO, enrichGO
#'
#' @return
go_wego_plot <- function(x, ...) {
  if(is(x, "list") & all(purrr::map_lgl(x, is_go_result))) {
    message("Generating wego plot")
    if(is(x, "gseaResult")) {
      message("to-do: wego plot for GSEA result, not available now, skipped")
    }
  } else {
    on.exit("`x` expect, list of GO Result")
    return(NULL)
  }
  # prepare data
  x_list <- lapply(x, function(d){
    if(is_go_result(d)) {
      d@result %>%
        dplyr::mutate(ontology = d@ontology)
    }
  })
  # for group/enrich
  df <- dplyr::bind_rows(x_list) %>%
    dplyr::mutate(label = paste0(Description, "-", ID)) %>%
    tidyr::separate("GeneRatio", c("num", "total"), sep = "/") %>%
    dplyr::mutate(pct = Count / as.numeric(total) * 100)
  # blank data
  if(nrow(df) == 0) {
    warning("No data found for wego plot")
    return(NULL)
  }
  # for group
  if("pvalue" %in% colnames(df)) {
    df <- df %>%
      dplyr::group_by(ontology) %>%
      dplyr::top_n(10, wt = Count) %>%
      dplyr::filter(Count > 0) %>%
      dplyr::arrange(desc(Count), pvalue) %>%
      dplyr::filter(row_number() <= 10)
  } else {
    df <- df %>%
      dplyr::group_by(ontology) %>%
      dplyr::top_n(10, wt = Count) %>%
      dplyr::filter(Count > 0) %>%
      dplyr::arrange(desc(Count), ID) %>%
      dplyr::filter(row_number() <= 10)
  }
  # assign colors, for ontology
  df$color <- plyr::mapvalues(df$ontology,
                              c("BP", "CC", "MF"),
                              scales::hue_pal()(3))
  # change ontology names
  df <- df %>%
    dplyr::ungroup() %>%
    dplyr::mutate(ontology = recode(ontology,
                                    "BP" = "Biological Process",
                                    "CC" = "Cell Component",
                                    "MF" = "Molecular Function"))
  # plot
  # second y axis
  coef <- max(df$pct) / max(df$Count)
  # main plot
  ggplot(df, aes(x = reorder(Description, -Count),
                 y = Count,
                 fill = ontology)) +
    geom_col() +
    xlab(NULL) +
    geom_text(aes(label = Count), vjust = 1.1) +
    facet_wrap(ontology~., scales = "free_x") +
    scale_y_continuous(sec.axis = sec_axis(~.*coef,
                                           name = "Percent of Genes (%)"),
                       name = "Number of Genes",
                       breaks = scales::pretty_breaks()) +
    theme_bw() +
    theme(axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 40, hjust = 1, size = 8),
          axis.text.y = element_text(size = 10),
          legend.position = "none",
          panel.grid = element_blank(),
          strip.background = element_rect(fill = "grey90",
                                          color = "black"))
}
