# functions for GO/KEGG plots

#' create plots
#' @param x list of objects of groupGO, enrichGO
#'
#' @export
go_barplot <- function(x, foldChange = NULL, text_width = 40) {
  if(! is_go(x)) {
    warning("groupGOResult, enrichResult expected")
    return(NULL)
  }

  # plot
  barplot(x, dorp = TRUE, showCategory = 12, order = TRUE) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = text_width)) +
    ylab("Number of genes") +
    ggtitle(x@ontology) +
    theme(plot.title = element_text(hjust = 0.5))
}


#' create plots
#' @param x object of enrichGO
#'
#' @export
go_dotplot <- function(x, foldChange = NULL, text_width = 40) {
  if(! is_go(x)) {
    warning("groupGOResult, enrichResult expected")
    return(NULL)
  }

  # plto
  clusterProfiler::dotplot(x, showCategory = 12) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = text_width)) +
    ylab("Gene Ratio") +
    ggtitle(x@ontology) +
    theme(axis.text.y = element_text(size = 10),
          plot.title = element_text(hjust = 0.5))
}


#' create plots
#' @param x object of enrichGO
#'
#' @export
go_cnetplot <- function(x, foldChange = NULL, text_width = 40) {
  if(! is_go(x)) {
    warning("groupGOResult, enrichResult expected")
    return(NULL)
  }

  # plot
  # clusterProfiler::cnetplot(x, foldChange = foldChange, categorySize = "pvalue") +
  clusterProfiler::cnetplot(x, foldChange = foldChange, categorySize = "pvalue") +
    ggtitle(x@ontology) +
    theme(text = element_text(size = 8),
          plot.title = element_text(hjust = .5))
}


#' create plots
#' @param x object of enrichGO
#'
#' @export
go_emapplot <- function(x, foldChange = NULL, text_width = 40) {
  if(! is_go(x)) {
    warning("groupGOResult, enrichResult expected")
    return(NULL)
  }

  # error !!!!
  # Error in graph_to_tree(graph, mode = direction) : Graph must be directed
  if(nrow(x) < 2) {
    return(NULL)
  }

  # plot
  clusterProfiler::emapplot(x, pie_scale = 1.5, layout = "kk") +
    ggtitle(x@ontology) +
    theme(text = element_text(size = 8))
}


#' create plots
#' @param x object of gseaResult
#'
#' @export
go_gsea_plot <- function(x, foldChange = NULL) {
  if(class(x) == "gseaResult") {
    message("Generating GSEA plot")
  } else {
    return(NULL)
  }

  # plot
  p_list <- lapply(seq_len(nrow(x)), function(i){
    enrichplot::gseaplot2(x, geneSetID = i, title = x$Description[i],
                          pvalue_table = TRUE)
  })
  names(p_list) <- seq_len(nrow(x)) # paste0("gsea.", seq_len(nrow(x)))
  p_list
}



#' create plots
#' @param x object of enrichGO
#'
#' @export
go_heatplot <- function(x, foldChange = NULL, text_width = 40) {
  if(! is_go(x)) {
    warning("groupGOResult, enrichResult expected")
    return(NULL)
  }

  # plot
  enrichplot::heatplot(x, foldChange = foldChange) +
    ggtitle(x@ontology) +
    theme(text = element_text(size = 8))
}


#' wego plot
#'
#' @param x list of objects of groupGO, enrichGO
#'
#' @return
go_wego_plot <- function(x) {
  if(! is_go(x, recursive = TRUE)) {
    warning("groupGO, enrichGO expected")
    return(NULL)
  }

  if(is(x[[1]], "gseaResult")) {
    message("Warning, to-do, wego for gseaResult")
    return(NULL)
  }

  # prepare data
  x_list <- lapply(x, function(d){
    d@result %>%
      dplyr::mutate(ontology = d@ontology)
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
