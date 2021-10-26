#' Functions for generating plots for GO, KEGG analysis
#'
#' Further plots are compatiable, with (gene_list, ...)
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
  #----------------------------------------------------------------------------#
  #-- Check arguments
  dots <- rlang::list2(...)
  args <- rlang::list2(
    show_category = 12,
    text_width    = 30, # label_format,
    font_size     = 12, # font.size
    x_axis        = "Count", # GeneRatio
    color_by      = "p.adjust", # color: pvalue, qvalue
    order         = TRUE,
    drop          = TRUE
  )
  dots <- purrr::list_modify(args, !!!dots)
  for(name in names(dots)) {
    assign(name, dots[[name]])
  }
  #-- convert arguments
  dots <- purrr::list_modify(
    dots,
    height       = x,
    showCategory = show_category,
    x            = x_axis,
    label_format = text_width,
    font.size    = font_size,
    color        = color_by, # pvalue, qvalue
    colorBy      = color_by
  )
  #----------------------------------------------------------------------------#
  #-- update arguments
  if(!inherits(x, "enrichResult")) {
    warning(glue::glue(
      "x is {class(x)}, expect 'enrichResult'"
    ))
    return(NULL)
  }
  #----------------------------------------------------------------------------#
  #-- enrichplot barplot.enrichResult, not support !!!
  rlang::exec(barplot, !!!dots) +
    ggtitle(x@ontology) +
    theme(plot.title = element_text(hjust = 0.5))
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
  #----------------------------------------------------------------------------#
  #-- Check arguments
  dots <- rlang::list2(...)
  args <- rlang::list2(
    show_category = 12,
    x_axis        = "geneRatio", # by: "geneRatio", "Percentage" and "count"
    color_by      = "p.adjust", # color, colorBy, colorBy, pvalue, qvalue
    show_catetory = 10,
    by            = "geneRatio", #
    size_by       = "geneRatio", # size; by: "geneRatio", "Percentage" and "count"
    font_size     = 12, # font.size
    order_by      = "x", # orderBy The order of the x-axis
    text_width    = 30, # label_format
    # split         = NULL, # ONTOLOGY,
    # decreasing    = TRUE
  )
  dots <- purrr::list_modify(args, !!!dots)
  for(name in names(dots)) {
    assign(name, dots[[name]])
  }
  #-- convert arguments
  dots <- purrr::list_modify(
    dots,
    object       = x,
    x            = x_axis,
    label_format = text_width,
    color        = color_by, # pvalue, qvalue
    colorBy      = color_by,
    by           = size_by, # GeneRatio
    orderBy      = order_by,
    font.size    = font_size
  )
  #----------------------------------------------------------------------------#
  #-- update arguments
  #-- enrichResult, gseaResult, compareClusterResult
  if(!is_go_result(x)) {
    warning(glue::glue(
      "x is {class(x)}, expect 'enrichResult'"
    ))
    return(NULL)
  }
  #----------------------------------------------------------------------------#
  # do.call(dotplot, dots) +
  rlang::exec(dotplot, !!!dots) +
    xlab("Gene Ratio") +
    ggtitle(x@ontology) +
    theme(axis.text.y = element_text(size  = 10),
          plot.title  = element_text(hjust = 0.5))
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
  #----------------------------------------------------------------------------#
  #-- Check arguments
  dots <- rlang::list2(...)
  args <- rlang::list2(
    fold_change   = NULL, # foldChange
    layout        = "nicely", #  'star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl' or 'lgl'
    show_category = 12
  )
  dots <- purrr::list_modify(args, !!!dots)
  #-- convert arguments
  dots <- purrr::list_modify(
    dots,
    x            = x,
    foldChange   = dots$fold_change,
    showCategory = dots$show_category,
    font.size    = dots$font_size
  )
  #----------------------------------------------------------------------------#
  #-- update arguments
  #-- enrichResult, gseaResult, compareClusterResult
  if(!is_go_result(x)) {
    warning(glue::glue(
      "x is {class(x)}, expect 'enrichResult'"
    ))
    return(NULL)
  }
  #----------------------------------------------------------------------------#
  rlang::exec(clusterProfiler::cnetplot, !!!dots) +
    xlab("Gene Ratio") +
    ggtitle(x@ontology) +
    theme(plot.title  = element_text(hjust = 0.5))
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
#' @export
go_emapplot <- function(x, ...) {
  #----------------------------------------------------------------------------#
  #-- Check arguments
  dots <- rlang::list2(...)
  args <- rlang::list2(
    fold_change   = NULL, # foldChange
    layout        = "nicely", #  'star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl' or 'lgl'
    show_category = 12,
    orgdb         = NULL,
    sim_method    = "Wang" # "Resnik", "Lin", "Rel", "Jiang" , "Wang" and "JC"
  )
  dots <- purrr::list_modify(args, !!!dots)
  #-- convert arguments
  dots <- purrr::list_modify(
    dots,
    x            = x,
    foldChange   = dots$fold_change,
    showCategory = dots$show_category,
    font.size    = dots$font_size
  )
  #----------------------------------------------------------------------------#
  #-- update arguments
  #-- enrichResult, gseaResult, compareClusterResult
  if(!is_go_result(x)) {
    warning(glue::glue(
      "x is {class(x)}, expect 'enrichResult'"
    ))
    return(NULL)
  }
  #----------------------------------------------------------------------------#
  #-- run
  ont <- tryCatch(
    x@ontology,
    error = function(cnd) return(NULL)
  )
  #-- similarity matrix
  if(inherits(dots$orgdb, "OrgDb") && ont %in% c("BP", "CC", "MF")) {
    d <- GOSemSim::godata(dots$orgdb, ont = ont)
    x2 <- enrichplot::pairwise_termsim(x, method = dots$sim_method, semData = d)
  } else {
    x2 <- enrichplot::pairwise_termsim(x)
  }
  #-- plot1
  p1 <- enrichplot::emapplot(
    x2,
    showCategory = dots$show_category,
    layout       = dots$layout,
    cex_category = .8,
    cex_line     = .4
  ) +
    ggtitle(ont) +
    theme(text = element_text(size = 10))
  ## deprecated: version 4.0.5
  # #-- plot2: cluster
  # p2 <- enrichplot::emapplot_cluster(
  #   x2,
  #   showCategory = dots$show_category,
  #   layout       = dots$layout,
  #   cex_category = .8,
  #   cex_line     = .4
  # ) +
  #   ggtitle(ont) +
  #   theme(text = element_text(size = 10))
  list(
    plot         = p1#,
    # plot_cluster = p2
  )
}


#' create plots
#' @param x object of enrichGO
#'
#' updated: 2020-12-15, pairwise_termsim(ego)
#' updated: 2021-09-01, deprecated, see treeplot()
#'
#' @export
go_emapplot_cluster <- function(x, ...) {
  #----------------------------------------------------------------------------#
  #-- Check arguments
  dots <- rlang::list2(...)
  args <- rlang::list2(
    fold_change   = NULL, # foldChange
    layout        = "nicely", #  'star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl' or 'lgl'
    show_category = 12,
    cex_category  = .8,
    cex_line      = .4,
  )
  dots <- purrr::list_modify(args, !!!dots)
  #-- convert arguments
  dots <- purrr::list_modify(
    dots,
    x            = x,
    foldChange   = dots$fold_change,
    showCategory = dots$show_category,
    font.size    = dots$font_size
  )
  #----------------------------------------------------------------------------#
  #-- update arguments
  #-- enrichResult, gseaResult, compareClusterResult
  if(!is_go_result(x)) {
    warning(glue::glue(
      "x is {class(x)}, expect 'enrichResult'"
    ))
    return(NULL)
  }
  #----------------------------------------------------------------------------#
  if(nrow(x) > 1) {
    x2 <- enrichplot::pairwise_termsim(x) # updated 2020-12-15
    rlang::exec(enrichplot::emapplot_cluster, !!!dots) +
      ggtitle(x@ontology) +
      theme(text = element_text(size = 10))
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
  #----------------------------------------------------------------------------#
  #-- Check arguments
  dots <- rlang::list2(...)
  args <- rlang::list2(
    fold_change   = NULL, # foldChange
    layout        = "nicely", #  'star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl' or 'lgl'
    show_category = 12,
    cex_category  = .8,
    cex_line      = .4,
  )
  args <- purrr::list_modify(args, x = x, foldChange = args$fold_change, !!!dots)
  #----------------------------------------------------------------------------#
  #-- update arguments
  #-- enrichResult, gseaResult, compareClusterResult
  if(!is_go_result(x)) {
    warning(glue::glue(
      "x is {class(x)}, expect 'enrichResult'"
    ))
    return(NULL)
  }
  #----------------------------------------------------------------------------#
  p_list <- lapply(seq_len(nrow(x)), function(i){
    args_local <- purrr::list_modify(
      args, geneSetID = i, title = x$Description[i], pvalue_table = TRUE
    )
    # rlang::exec(enrichplot::gseaplot2, !!!args_local)
    rlang::exec(enrichplot::gseaplot, !!!args_local)
  })
  names(p_list) <- seq_len(nrow(x)) # paste0("gsea.", seq_len(nrow(x)))
  p_list
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
  #----------------------------------------------------------------------------#
  #-- Check arguments
  dots <- rlang::list2(...)
  args <- rlang::list2(
    fold_change   = NULL, # foldChange
    layout        = "nicely", #  'star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl' or 'lgl'
    show_category = 12
  )
  dots <- purrr::list_modify(args, !!!dots)
  #-- convert arguments
  dots <- purrr::list_modify(
    dots,
    x            = x,
    foldChange   = dots$fold_change,
    showCategory = dots$show_category,
    font.size    = dots$font_size
  )
  #----------------------------------------------------------------------------#
  #-- update arguments
  #-- enrichResult, gseaResult, compareClusterResult
  if(!is_go_result(x)) {
    warning(glue::glue(
      "x is {class(x)}, expect 'enrichResult'"
    ))
    return(NULL)
  }
  #----------------------------------------------------------------------------#
  rlang::exec(enrichplot::heatplot, !!!dots) +
    ggtitle(x@ontology) +
    theme(text = element_text(size = 10))
}




#' create plots
#' @param x object of enrichGO
#'
#' updated: 2020-12-15, pairwise_termsim(ego)
#'
#' @export
go_treeplot <- function(x, ...) {
  #----------------------------------------------------------------------------#
  #-- Check arguments
  dots <- rlang::list2(...)
  args <- rlang::list2(
    show_category = 12, # showCategory
    fold_change   = NULL, # foldChange
    color_by      = "p.adjust", # pvalue, p.adjust or qvalue, or custome
    n_words       = 4, # nWords
    n_clusters   = 5, # nCluster
    hclust_method = "ward.D", # "ward.D2", "single", "average", "median", "complete"
    layout        = "nicely", #  'star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl' or 'lgl'
    font_size     = 4, # fontsize
    text_width    = 30 # label_format
  )
  dots <- purrr::list_modify(args, !!!dots)
  #-- to env
  for(name in names(dots)) {
    assign(name, dots[[name]])
  }
  #-- update arguments
  dots <- purrr::list_modify(
    dots,
    x            = x,
    showCategory = dots$show_category,
    color        = color_by, #
    nWords       = n_words,
    nCluster     = n_clusters,
    fontsize     = font_size,
    label_format = text_width,
    foldChange   = fold_change,
  )
  #----------------------------------------------------------------------------#
  #-- enrichResult, gseaResult, compareClusterResult
  if(!is_go_result(x)) {
    warning(glue::glue(
      "x is {class(x)}, expect 'enrichResult'"
    ))
    return(NULL)
  }
  #----------------------------------------------------------------------------#
  if(nrow(x) > 1) {
    dots$x <- enrichplot::pairwise_termsim(x) # update x
    rlang::exec(enrichplot::treeplot, !!!dots) +
      ggtitle(x@ontology) +
      theme(text = element_text(size = 10))
  }
}



#' wego plot
#'
#' @param x list of objects of groupGO, enrichGO
#'
#' @return
go_wego_plot <- function(x, ...) {
  #----------------------------------------------------------------------------#
  #-- Check arguments
  dots <- rlang::list2(...)
  args <- rlang::list2(
    show_category = 12, # showCategory
    fold_change   = NULL, # foldChange
    color_by      = "p.adjust" # pvalue, p.adjust or qvalue, or custome
  )
  dots <- purrr::list_modify(args, !!!dots)
  #-- to env
  for(name in names(dots)) {
    assign(name, dots[[name]])
  }
  #----------------------------------------------------------------------------#
  # #-- enrichResult, gseaResult, compareClusterResult
  # if(!is_go_result(x)) {
  #   warning(glue::glue(
  #     "x is {class(x)}, expect 'enrichResult'"
  #   ))
  #   return(NULL)
  # }
  #----------------------------------------------------------------------------#
  if(is(x, "list") & all(purrr::map_lgl(x, is_go_result))) {
    message("Generating wego plot")
    if(is(x, "gseaResult")) {
      message("to-do: wego plot for GSEA result, not available now, skipped")
    }
  } else {
    message(glue::glue("x is {class(x)}, expect list of GO Result"))
    return(NULL)
  }
  # add ont to results
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
