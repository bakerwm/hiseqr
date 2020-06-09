

#' make plots for publication
#'
#' @param x Path to the xls file of DEseq2 output, eg, function
#'   DESeq2_for_featureCounts, result "transcripts_deseq2.csv.fix.xls"
#'
#' @param outdir Directory to save the pdf file, contains
#'   scatter, MA, Volcano plots
#'
#' @param save2df A logical value, whether or not save the plot to PDF file,
#'   default: TRUE
#'
#'
#' @import ggplot2
#' @import ggrepel
#'
#' @export
DESeq2_publish_plot <- function(x, outdir,
                                labels = NULL, # c("nxf2", "piwi", "CG9754"),
                                x_name = NULL,
                                y_name = NULL,
                                fc_cutoff = 2,
                                pval_cutoff = 0.05,
                                save2pdf = TRUE) {
  stopifnot(file.exists(x))
  df <- readr::read_delim(x, "\t", col_types = readr::cols())

  # rename id FBgn0000001_te
  # convert to FB00000001 and te
  # gene/name
  if ("Gene" %in% colnames(df)) {
    df$id <- df$Gene
  } else if ("id" %in% colnames(df)) {
    df$id <- df$id
  } else {
    warning("gene name not fount, Gene|id, expected")
    return(NULL)
  }

  # for TE ids
  # whether FB***._
  # id_tag <- sum(grepl("FBgn.*_", df$id)) / length(df$id)
  #
  # if(ceiling(id_tag) == 1) {
  #   df <- df %>%
  #     tidyr::separate(id, into = c("name", "id"), sep = "_")
  #   df$name <- NULL
  # }

  names(df) <- gsub("[^A-Za-z0-9]", "_", names(df)) # remove unsupport characters
  x_name <- colnames(df)[2]
  y_name <- colnames(df)[3]

  # make plots
  p1 <- DESeq2_de_scatter(df, x_name, y_name,
                          show_sig = TRUE,
                          labels = labels,
                          fc_cutoff = fc_cutoff,
                          pval_cutoff = pval_cutoff,
                          max_labels = 5
  )

  p2 <- DESeq2_de_ma(df,
                     show_sig = TRUE,
                     labels = labels,
                     fc_cutoff = fc_cutoff,
                     pval_cutoff = pval_cutoff,
                     max_labels = 5,
                     show_volcano = FALSE
  )

  p3 <- DESeq2_de_ma(df,
                     show_sig = TRUE,
                     labels = labels,
                     fc_cutoff = fc_cutoff,
                     pval_cutoff = pval_cutoff,
                     max_labels = 5,
                     show_volcano = TRUE
  )

  # add legend
  figure_legend <- "Figure. Differentially expression analysis. (a-b) Comparasion of ene expression is shown as rpm from two conditions.
  Dashed lines indicate two fold change. a. scatter plot, b. MA plot.
  (c) Volcano plot showing enrichment values and corresponding significance levels."

  pg1 <- cowplot::plot_grid(p1, p2, p3,
                            align = "hv",
                            ncol = 2,
                            labels = "auto")
  pg2 <- cowplot::add_sub(pg1, figure_legend, x = 0, hjust = 0, size = 10)
  # cowplot::ggdraw(pg2)

  # save components
  obj_list <- list(
    p1 = p1,
    p2 = p2,
    p3 = p3,
    data = df,
    labels = labels,
    x_name = x_name,
    y_name = y_name,
    fc_cutoff = fc_cutoff,
    pval_cutoff = pval_cutoff
  )

  obj_file <- file.path(outdir, "publish_plot_data.rds")
  saveRDS(obj_list, file = obj_file)

  # save plot
  if (isTRUE(save2pdf)) {
    if (!dir.exists(outdir)) {
      dir.create(outdir, showWarnings = FALSE, recursive = TRUE, mode = "0740")
    }
    rpt_fname <- paste0("DESeq2.", x_name, ".vs.", y_name, ".publish.pdf")
    rpt_file <- file.path(outdir, rpt_fname)
    pdf(rpt_file, width = 8, height = 6, paper = "a4r")
    print(cowplot::ggdraw(pg2))
    dev.off()
  } else {
    return(pg2)
  }
}


#' generate scatter plot for DESeq2 analysis
#'
#' @param data A data frame with gene expression,
#' @param x_name Name in column names of data frame, present on x axis
#' @param y_name Same as x_name, but present on y axis
#' @param show_sig A logical or int value, whether or not to show
#'   significantly changes genes, default: FALSE
#' @param labels A set of gene ids to show in plot, default: NULL
#' @param pval_cutoff A value to select significantly changes genes, default: 0.05
#' @param max_labels A integer value, how many genes to be labeled on plot,
#'   default: 3
#'
#' @import ggplot2
#' @import ggrepeal
#'
#' @export
#'
DESeq2_de_scatter <- function(data, x_name = NULL, y_name = NULL,
                              show_sig = TRUE,
                              labels = NULL,
                              fc_cutoff = 2,
                              pval_cutoff = 0.05,
                              max_labels = 5) {
  stopifnot(is.data.frame(data))

  # required columns
  required_cols <- c("id", "log2FoldChange", "padj")
  stopifnot(all(required_cols %in% colnames(data)))

  # determine labels
  # top DE genes
  sig_labels <- data %>%
    dplyr::filter(abs(log2FoldChange) >= log2(fc_cutoff) &
                    padj <= pval_cutoff) %>%
    dplyr::arrange(padj) %>%
    dplyr::pull(id) %>%
    head(max_labels)

  if (is.null(labels)) {
    labels <- sig_labels
  } else {
    max_labels <- ifelse(length(labels) > max_labels, max_labels, length(labels))
    labels <- labels[1:max_labels]
  }

  p <- plot_scatter(data,
                    x_name = x_name,
                    y_name = y_name,
                    labels = labels,
                    show_sig = show_sig,
                    fc_cutoff = fc_cutoff,
                    pval_cutoff = pval_cutoff
  )

  return(p)
}


#' create MA plot
#'
#' @param data A data frame with gene expression,
#' @param show_sig A logical or int value, whether or not to show
#'   significantly changes genes, default: FALSE
#' @param labels A set of gene ids to show in plot, default: NULL
#' @param pval_cutoff A value to select significantly changes genes, default: 0.05
#' @param max_labels A integer value, how many genes to be labeled on plot,
#'   default: 3
#'
#' @import ggplot2
#' @import ggrepeal
#'
#' @export
#'
DESeq2_de_ma <- function(data, show_sig = FALSE,
                         labels = NULL,
                         fc_cutoff = 2,
                         pval_cutoff = 0.05,
                         max_labels = 3,
                         show_volcano = FALSE) {
  stopifnot(is.data.frame(data))

  # required columns
  required_cols <- c("id", "baseMean", "log2FoldChange", "padj")
  stopifnot(all(required_cols %in% colnames(data)))

  # determine labels
  # top DE genes
  sig_labels <- data %>%
    dplyr::filter(abs(log2FoldChange) >= log2(fc_cutoff) &
                    padj <= pval_cutoff) %>%
    dplyr::arrange(padj) %>%
    dplyr::pull(id) %>%
    head(max_labels)

  if (is.null(labels)) {
    labels <- sig_labels
  } else {
    labels <- labels[1:max_labels]
  }

  # plot
  if (isTRUE(show_volcano)) {
    plot_func <- plot_volcano
  } else {
    plot_func <- plot_ma
  }

  p <- plot_func(
    data, labels, show_sig,
    fc_cutoff, pval_cutoff
  )

  return(p)
}




#' generate scatter plot for DE analysis
#'
#' @param data A data frame for the full version of data, required
#' @param labels string or list, ids to lable on plot
#' @param x_name column name to show on x-axis, default: column-2
#' @param y_name column name to show on y-axis, default: column-3
#' @param xLabel title for x axis, default: x_name
#' @param yLabel title for y axis, default: y_name
#' @param show_sig boolen, whether or not show significant genes in colors
#' @param fc_cutoff numeric, cutoff for foldchange
#' @param pval_cutoff numeric, cutoff for pvalue
#'
#' @import ggplot2
#' @import ggrepel
#'
#' @export
#'
plot_scatter <- function(data, x_name = NULL, y_name = NULL,
                         xLabel = NULL, yLabel = NULL,
                         labels = NULL, show_sig = TRUE,
                         fc_cutoff = 2, pval_cutoff = 0.05,
                         show_volcano = NULL) {
  # data_sig is required
  stopifnot(inherits(data, "data.frame"))

  # convert num
  fc_cutoff <- as.numeric(fc_cutoff)
  pval_cutoff <- as.numeric(pval_cutoff)

  # show sig
  if (isTRUE(show_sig)) {
    stopifnot(all(c("id", "log2FoldChange", "pvalue") %in% names(data)))
    # plot_df <- deseq_sig(data, fc = fc_cutoff, pval = pval_cutoff, type = "all")
    plot_df <- deseq_sig(data, fc = fc_cutoff, pval = pval_cutoff)
  } else {
    plot_df$sig <- "gene"
  }

  # determine x,y
  x_name <- ifelse(is.null(x_name), colnames(data)[2], x_name)
  y_name <- ifelse(is.null(y_name), colnames(data)[3], y_name)
  stopifnot(all(c(x_name, y_name) %in% colnames(data)))

  # labels on axis
  xLabel <- ifelse(is.null(xLabel), x_name, xLabel)
  yLabel <- ifelse(is.null(yLabel), y_name, yLabel)

  # axis title
  xTitle <- bquote(.(xLabel) ~ "[" * log["10"] ~ "rpm]")
  yTitle <- bquote(.(yLabel) ~ "[" * log["10"] ~ "rpm]")

  # pre-process
  if (nrow(plot_df) == 0) {
    return(NULL)
  }

  # function
  to_log10 <- function(data, x, y) {
    stopifnot(is.data.frame(data))
    stopifnot(all(c("id", x, y) %in% names(data)))
    # select columns
    # df <- dplyr::select_(data, .dots = c("id", x, y)) %>%
    #   as.data.frame()
    df <- dplyr::select(data, id, x, y) %>%
      as.data.frame()
    # remove old row names
    rownames(df) <- NULL
    df <- tibble::column_to_rownames(df, "id") %>% log10()
    df <- df[complete.cases(df), ]
    df <- tibble::rownames_to_column(df, "id")
    return(df)
  }

  df_log10 <- to_log10(plot_df, x_name, y_name)

  # combine data
  plot_df2 <- plot_df %>% dplyr::select(-c(x_name, y_name))
  plot_df3 <- merge(df_log10, plot_df2, by = "id")

  if (nrow(plot_df3) == 0) {
    return(NULL)
  }

  # basic plot
  p <- ggplot(plot_df3, aes_string(x_name, y_name, color = "sig")) +
    geom_abline(slope = 1, intercept = 0, color = "grey10") +
    geom_abline(
      slope = 1, intercept = log10(2),
      color = "grey30", linetype = 2
    ) +
    geom_abline(
      slope = 1, intercept = -log10(2),
      color = "grey30", linetype = 2
    ) +
    geom_point(size = .4)

  # axis
  p <- p +
    scale_x_continuous(
      limits = c(0, 6),
      breaks = seq(0, 6, 1),
      labels = seq(0, 6, 1),
      expand = c(0, 0, 0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, 6),
      breaks = seq(0, 6, 1),
      labels = seq(0, 6, 1),
      expand = c(0, 0, 0, 0)
    )

  # add colors
  sigColors <- c("green4", "grey70", "red", "grey30")
  names(sigColors) <- c("down", "not", "up", "gene")
  sigLevels <- sort(unique(plot_df3$sig))
  hitColors <- sigColors[sigLevels]
  p <- p + scale_color_manual(values = hitColors)

  # show labels
  labelDf <- plot_df3 %>% dplyr::filter(id %in% labels)
  if (nrow(labelDf) > 0) {
    p <- p +
      ggrepel::geom_text_repel(
        mapping = aes_string(x_name, y_name, color = "sig"),
        data = labelDf,
        label = labelDf$id,
        nudge_x = 1.2,
        nudge_y = .02,
        point.padding = .4,
        box.padding = .4,
        segment.size = .6,
        segment.color = "grey60",
        direction = "both"
      )
  }

  # highlight points
  ptDf <- plot_df3 %>% dplyr::filter(id %in% show_volcano)
  if (nrow(ptDf) > 0) {
    p <- p +
      geom_point(aes_string(x_name, y_name, color = "sig"),
                 ptDf,
                 shape = 7, size = 2
      )
  }

  # theme
  p <- p +
    xlab(xTitle) +
    ylab(yTitle) +
    labs(color = "Significant") +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = .7),
      plot.title = element_text(color = "black", hjust = .5, size = 14),
      panel.grid = element_blank(),
      axis.line = element_blank(),
      axis.ticks.length = unit(.2, "cm"),
      axis.ticks = element_line(color = "black", size = .5),
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(color = "black", size = 12)
    )

  return(p)
}


#' generate volcano plot for DE analysis
#'
#' @param data A data frame for the full version of data, required
#' @param labels string or list, ids to lable on plot
#' @param show_sig boolen, whether or not show significant genes in colors
#' @param fc_cutoff numeric, cutoff for foldchange
#' @param pval_cutoff numeric, cutoff for pvalue
#'
#' @import ggplot2
#' @import ggrepel
#'
#' @export
#'
plot_ma <- function(data, labels = NULL, show_sig = TRUE,
                    fc_cutoff = 2, pval_cutoff = 0.05) {
  # data_sig is required
  stopifnot(is.data.frame(data))

  stopifnot(all(c("id", "baseMean", "log2FoldChange", "pvalue") %in% colnames(data)))
  if ("padj" %in% colnames(data)) {
    data$pvalue <- data$padj
  }

  # convert num
  fc_cutoff <- as.numeric(fc_cutoff)
  pval_cutoff <- as.numeric(pval_cutoff)

  ## add significant tags
  plot_df <- deseq_sig(data, fc = fc_cutoff, pval = pval_cutoff, type = "all") # sig genes
  plot_df <- dplyr::rename(plot_df, logFC = log2FoldChange, pvalue = pvalue)

  # split data
  if (isFALSE(show_sig) | !"sig" %in% colnames(plot_df)) {
    plot_df$sig <- "gene"
  }

  # log10
  plot_df$baseMean <- log10(plot_df$baseMean + 1)

  # axis title
  xTitle <- expression(paste(log["10"], "(baseMean)"))
  yTitle <- expression(paste(log["2"], "(Fold-change)"))

  # basic plot
  p <- ggplot(plot_df, aes(baseMean, logFC, color = sig)) +
    geom_hline(yintercept = 0, color = "grey10", size = .7) +
    geom_point(size = .4)

  # axis
  p <- p + scale_y_continuous(limits = c(-4, 4))

  # add colors
  sigColors <- c("green4", "grey70", "red", "grey30")
  names(sigColors) <- c("down", "not", "up", "gene")
  sigLevels <- sort(unique(plot_df$sig))
  hitColors <- sigColors[sigLevels]
  p <- p + scale_color_manual(values = hitColors)

  # show labels
  labelDf <- plot_df %>% dplyr::filter(id %in% labels)
  if (nrow(labelDf) > 0) {
    p <- p +
      ggrepel::geom_text_repel(
        mapping = aes(baseMean, logFC, color = sig),
        data = labelDf,
        label = labelDf$id,
        nudge_x = 1.2,
        nudge_y = .02,
        point.padding = .4,
        box.padding = .4,
        segment.size = .6,
        segment.color = "grey60",
        direction = "both"
      )
  }

  # theme
  p <- p +
    xlab(xTitle) +
    ylab(yTitle) +
    labs(color = "Significant") +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = .7),
      plot.title = element_text(color = "black", hjust = .5, size = 14),
      panel.grid = element_blank(),
      axis.line = element_blank(),
      axis.ticks.length = unit(.2, "cm"),
      axis.ticks = element_line(color = "black", size = .5),
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(color = "black", size = 12)
    )

  return(p)
}



#' generate volcano plot for DE analysis
#'
#' @param data A data frame for the full version of data
#' @param labels list, ids to label in plot
#' @param show_sig boolen, whether or not show significant genes in colors
#' @param fc_cutoff numeric, cutoff for foldchange
#' @param pval_cutoff numeric, cutoff for pvalue
#'
#' @import ggplot2
#' @import ggrepel
#'
#' @export
#'
plot_volcano <- function(data, labels = NULL, show_sig = TRUE,
                         fc_cutoff = 2, pval_cutoff = 0.05) {
  # data_sig not required
  stopifnot(is.data.frame(data))
  stopifnot(all(c("id", "log2FoldChange", "pvalue") %in% colnames(data)))
  if ("padj" %in% colnames(data)) {
    data$pvalue <- data$padj
  }

  # convert num
  fc_cutoff <- as.numeric(fc_cutoff)
  pval_cutoff <- as.numeric(pval_cutoff)

  ## add significant tags
  plot_df <- deseq_sig(data, fc = fc_cutoff, pval = pval_cutoff, type = "all") # sig genes
  # plot_df <- dplyr::rename(plot_df, logFC = log2FoldChange, pvalue = pvalue)

  # split data
  if (isFALSE(show_sig) | !"sig" %in% colnames(plot_df)) {
    plot_df$sig <- "gene"
  }

  # axis title
  xtitle <- expression(paste(log["2"], "(mean ratio of Treatment/Control)"))
  ytitle <- expression(paste(-log["10"], "(", italic("P"), " value)"))

  # convert pval to -log10()
  plot_df$logPval <- -log10(plot_df$pvalue)

  # basic plot
  p <- ggplot(plot_df, aes(log2FoldChange, logPval, color = sig)) +
    geom_vline(xintercept = 0, size = .5, color = "black") +
    geom_vline(
      xintercept = c(-log2(fc_cutoff), log2(fc_cutoff)),
      size = .5,
      linetype = 2,
      color = "grey60"
    ) +
    geom_hline(
      yintercept = -log10(pval_cutoff),
      size = .5,
      linetype = 2,
      color = "grey60"
    )

  # add points
  p <- p + geom_point(size = .5)

  # add colors
  sig_colors <- c("green4", "grey70", "red", "grey30")
  names(sig_colors) <- c("down", "not", "up", "gene")
  sig_levels <- sort(unique(plot_df$sig))
  hit_colors <- sig_colors[sig_levels]
  p <- p + scale_color_manual(values = hit_colors)

  # show labels
  data_label <- plot_df %>% dplyr::filter(id %in% labels)
  if (nrow(data_label) > 0) {
    p <- p +
      ggrepel::geom_text_repel(
        mapping = aes(log2FoldChange, logPval, color = sig),
        data = data_label,
        label = data_label$id,
        nudge_x = 1.2,
        nudge_y = .02,
        point.padding = .4,
        box.padding = .4,
        segment.size = .6,
        segment.color = "grey60",
        direction = "both"
      )
  }

  # theme
  p <- p +
    xlab(xtitle) +
    ylab(ytitle) +
    labs(color = "Significant") +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", size = .5),
      axis.ticks = element_line(color = "black", size = .5),
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(color = "black", size = 12),
      axis.ticks.length = unit(.2, "cm")
    )

  return(p)
}
