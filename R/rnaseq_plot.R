

#' make plots for publication
#'
#' @param x Path to the xls file of DEseq2 output, eg, function
#'   DESeq2_for_featureCounts, result "transcripts_deseq2.csv.fix.xls"
#'
#' @param pathPdf Directory to save the pdf file, contains
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
DESeq2_publish_plot <- function(x, pathPdf,
                                labels = NULL, #c("nxf2", "piwi", "CG9754"),
                                xName = NULL,
                                yName = NULL,
                                fcCutoff = 2,
                                pvalCutoff = 0.05,
                                save2pdf = TRUE) {
  stopifnot(file.exists(x))
  df <- deseqCsvReader(x)

  # rename id FBgn0000001_te
  # convert to FB00000001 and te
  id_tag <- sum(grepl("FBgn.*_", df$id)) / length(df$id)
  if(ceiling(id_tag) == 1) {
    df <- df %>%
      tidyr::separate(id, into = c("name", "id"), sep = "_")
    df$name <- NULL
  }

  names(df) <- gsub("[^A-Za-z0-9]", "_", names(df)) # remove unsupport characters
  xName <- colnames(df)[2]
  yName <- colnames(df)[3]

  # make plots
  p1 <- DESeq2_de_scatter(df, xName, yName,
                          showSig = TRUE,
                          labels = labels,
                          fcCutoff = fcCutoff,
                          pvalCutoff = pvalCutoff,
                          maxLabels = 5)

  p2 <- DESeq2_de_ma(df, showSig = TRUE,
                     labels = labels,
                     fcCutoff = fcCutoff,
                     pvalCutoff = pvalCutoff,
                     maxLabels = 5,
                     showVolcano = FALSE)

  p3 <- DESeq2_de_ma(df, showSig = TRUE,
                     labels = labels,
                     fcCutoff = fcCutoff,
                     pvalCutoff = pvalCutoff,
                     maxLabels = 5,
                     showVolcano = TRUE)

  # add legend
  figure_legend <- "Figure. Differentially expression analysis. (a-b) Comparasion of ene expression is shown as rpm from two conditions.
  Dashed lines indicate two fold change. a. scatter plot, b. MA plot.
  (c) Volcano plot showing enrichment values and corresponding significance levels."

  pg1 <- cowplot::plot_grid(p1, p2, p3, align = "hv",
                            ncol = 2, labels = "auto")
  pg2 <- cowplot::add_sub(pg1, figure_legend, x = 0, hjust = 0, size = 10)
  # cowplot::ggdraw(pg2)

  # save plot
  if (isTRUE(save2pdf)) {
    if (! dir.exists(pathPdf)) {
      dir.create(pathPdf, showWarnings = FALSE, recursive = TRUE, mode = "0740")
    }
    rpt_fname <- paste0("DESeq2.", xName, ".vs.", yName, ".publish.pdf")
    rpt_file  <- file.path(pathPdf, rpt_fname)
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
#' @param xName Name in column names of data frame, present on x axis
#' @param yName Same as xName, but present on y axis
#' @param showSig A logical or int value, whether or not to show
#'   significantly changes genes, default: FALSE
#' @param labels A set of gene ids to show in plot, default: NULL
#' @param pvalCutoff A value to select significantly changes genes, default: 0.05
#' @param maxLabels A integer value, how many genes to be labeled on plot,
#'   default: 3
#'
#' @import ggplot2
#' @import ggrepeal
#'
#' @export
#'
DESeq2_de_scatter <- function(data, xName = NULL, yName = NULL,
                              showSig = TRUE,
                              labels = NULL,
                              fcCutoff = 2, pvalCutoff = 0.05,
                              maxLabels = 5) {
  stopifnot(is.data.frame(data))

  # required columns
  reqCols <- c("id", "log2FoldChange", "padj")
  stopifnot(all(reqCols %in% colnames(data)))

  # determine labels
  # top DE genes
  sigLabels <- data %>%
    dplyr::filter(abs(log2FoldChange) >= log2(fcCutoff) &
                    padj <= pvalCutoff) %>%
    dplyr::arrange(padj) %>%
    dplyr::pull(id) %>%
    head(maxLabels)

  if (is.null(labels)) {
    labels <- sigLabels
  } else {
    maxLabels <- ifelse(length(labels) > maxLabels, maxLabels, length(labels))
    labels <- labels[1:maxLabels]
  }

  p <- plot_scatter(data,
                    xName      = xName,
                    yName      = yName,
                    labels     = labels,
                    showSig    = showSig,
                    fcCutoff   = fcCutoff,
                    pvalCutoff = pvalCutoff)

  return(p)
}


#' create MA plot
#'
#' @param data A data frame with gene expression,
#' @param showSig A logical or int value, whether or not to show
#'   significantly changes genes, default: FALSE
#' @param labels A set of gene ids to show in plot, default: NULL
#' @param pvalCutoff A value to select significantly changes genes, default: 0.05
#' @param maxLabels A integer value, how many genes to be labeled on plot,
#'   default: 3
#'
#' @import ggplot2
#' @import ggrepeal
#'
#' @export
#'
DESeq2_de_ma <- function(data, showSig = FALSE, labels = NULL,
                         fcCutoff = 2, pvalCutoff = 0.05, maxLabels = 3,
                         showVolcano = FALSE) {
  stopifnot(is.data.frame(data))

  # required columns
  reqCols <- c("id", "baseMean", "log2FoldChange", "padj")
  stopifnot(all(reqCols %in% colnames(data)))

  # determine labels
  # top DE genes
  sigLabels <- data %>%
    dplyr::filter(abs(log2FoldChange) >= log2(fcCutoff) &
                    padj <= pvalCutoff) %>%
    dplyr::arrange(padj) %>%
    dplyr::pull(id) %>%
    head(maxLabels)

  if (is.null(labels)) {
    labels <- sigLabels
  } else {
    labels <- labels[1:maxLabels]
  }

  # plot
  if (isTRUE(showVolcano)) {
    plot_func <- plot_volcano
  } else {
    plot_func <- plot_ma
  }

  p <- plot_func(data, labels, showSig,
                 fcCutoff, pvalCutoff)

  return(p)
}



#' add sig label to data.frame
#'
#' @param data data.frame
#' @param type all, both, up, down
#' @param fcCutoff numeric, cutoff for foldchange
#' @param pvalCutoff numeric, cutoff for pvalue
#' @param type string, all, both, up, down, not, default: all
#'
#' @import dplyr
#'
#' @export
#'
sig_genes <- function(data, fcCutoff = 2, pvalCutoff = 0.05, type = "all") {
  stopifnot(inherits(data, "data.frame"))

  stopifnot(all(c("id", "log2FoldChange", "pvalue") %in% colnames(data)))
  if ("padj" %in% colnames(data)) {
    data$pvalCheck <- data$padj
  } else {
    data$pvalCheck <- data$pvalue
  }

  type   <- tolower(type) # lower case
  log2FC <- log2(as.numeric(fcCutoff))

  # columns required
  columns1 <- c("id", "log2FoldChange", "pvalue")

  # split group
  df1 <- dplyr::filter(data, log2FoldChange >= log2FC & pvalCheck <= pvalCutoff)
  df2 <- dplyr::filter(data, log2FoldChange <= -log2FC & pvalCheck <= pvalCutoff)
  df3 <- dplyr::filter(data, ! id %in% c(df1$id, df2$id))

  if (nrow(df1) > 0 ) {
    df1$sig <- "up"
  }

  if (nrow(df2) > 0) {
    df2$sig <- "down"
  }

  if (nrow(df3) > 0) {
    df3$sig <- "not"
  }

  # up/down
  if (type == "up") {
    df <- df1
    # df <- df %>% dplyr::filter(sig == "up")
  } else if (type == "down") {
    df <- df2
    # df <- df %>% dplyr::filter(sig == "down")
  } else if (type == "both") {
    df <- dplyr::bind_cols(df1, df2)
    # df <- df %>% dplyr::filter(sig %in% c("up", "down"))
  } else if (type == "not") {
    df <- df3
    # df <- df %>% dplyr::filter(sig %in% c("not"))
  } else {
    df <- dplyr::bind_rows(df1, df2, df3)
  }

  # remove pvalCheck column
  df <- dplyr::select(df, - pvalCheck)

  # return
  return(df)
}



#' generate scatter plot for DE analysis
#'
#' @param data A data frame for the full version of data, required
#' @param labels string or list, ids to lable on plot
#' @param xName column name to show on x-axis, default: column-2
#' @param yName column name to show on y-axis, default: column-3
#' @param xLabel title for x axis, default: xName
#' @param yLabel title for y axis, default: yName
#' @param showSig boolen, whether or not show significant genes in colors
#' @param fcCutoff numeric, cutoff for foldchange
#' @param pvalCutoff numeric, cutoff for pvalue
#'
#' @import ggplot2
#' @import ggrepel
#'
#' @export
#'
plot_scatter <- function(data, xName = NULL, yName = NULL,
                         xLabel = NULL, yLabel = NULL,
                         labels = NULL, showSig = TRUE,
                         fcCutoff = 2, pvalCutoff = 0.05,
                         addPoints = NULL) {
  # data_sig is required
  stopifnot(inherits(data, "data.frame"))

  # convert num
  fcCutoff   <- as.numeric(fcCutoff)
  pvalCutoff <- as.numeric(pvalCutoff)

  # show sig
  if (isTRUE(showSig)) {
    stopifnot(all(c("id", "log2FoldChange", "pvalue") %in% names(data)))
    plotDf <- sig_genes(data, fc = fcCutoff, pval = pvalCutoff, type = "all")
  } else {
    plotDf$sig <- "gene"
  }

  # determine x,y
  xName <- ifelse(is.null(xName), colnames(data)[2], xName)
  yName <- ifelse(is.null(yName), colnames(data)[3], yName)
  stopifnot(all(c(xName, yName) %in% colnames(data)))

  # labels on axis
  xLabel <- ifelse(is.null(xLabel), xName, xLabel)
  yLabel <- ifelse(is.null(yLabel), yName, yLabel)

  # axis title
  xTitle <- bquote(.(xLabel)~"["*log["10"]~"rpm]")
  yTitle <- bquote(.(yLabel)~"["*log["10"]~"rpm]")

  # pre-process
  if (nrow(plotDf) == 0) {
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

  plotDfLog10 <- to_log10(plotDf, xName, yName)

  # combine data
  plotDf2 <- plotDf %>% dplyr::select(-c(xName, yName))
  plotDf3 <- merge(plotDfLog10, plotDf2, by = "id")

  if (nrow(plotDf3) == 0) {
    return(NULL)
  }

  # basic plot
  p <- ggplot(plotDf3, aes_string(xName, yName, color = "sig")) +
    geom_abline(slope = 1, intercept = 0, color = "grey10") +
    geom_abline(slope = 1, intercept = log10(2),
                color = "grey30", linetype = 2) +
    geom_abline(slope = 1, intercept = -log10(2),
                color = "grey30", linetype = 2) +
    geom_point(size = .4)

  # axis
  p <- p +
    scale_x_continuous(limits = c(0, 6),
                       breaks = seq(0, 6, 1),
                       labels = seq(0, 6, 1),
                       expand = c(0, 0, 0, 0)) +
    scale_y_continuous(limits = c(0, 6),
                       breaks = seq(0, 6, 1),
                       labels = seq(0, 6, 1),
                       expand = c(0, 0, 0, 0))

  # add colors
  sigColors <- c("green4", "grey70", "red", "grey30")
  names(sigColors) <- c("down", "not", "up", "gene")
  sigLevels <- sort(unique(plotDf3$sig))
  hitColors <- sigColors[sigLevels]
  p <- p + scale_color_manual(values = hitColors)

  # show labels
  labelDf <- plotDf3 %>% dplyr::filter(id %in% labels)
  if (nrow(labelDf) > 0) {
    p <- p +
      geom_text_repel(
        mapping       = aes_string(xName, yName, color = "sig"),
        data          = labelDf,
        label         = labelDf$id,
        nudge_x       = 1.2,
        nudge_y       = .02,
        point.padding = .4,
        box.padding   = .4,
        segment.size  = .6,
        segment.color = "grey60",
        direction     = "both")
  }

  # highlight points
  ptDf <- plotDf3 %>% dplyr::filter(id %in% addPoints)
  if (nrow(ptDf) > 0) {
    p <- p +
      geom_point(aes_string(xName, yName, color = "sig"),
                 ptDf,
                 shape = 7, size = 2)
  }

  # theme
  p <- p +
    xlab(xTitle) +
    ylab(yTitle) +
    labs(color = "Significant") +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = .7),
          plot.title   = element_text(color = "black", hjust = .5, size = 14),
          panel.grid   = element_blank(),
          axis.line    = element_blank(),
          axis.ticks.length = unit(.2, "cm"),
          axis.ticks   = element_line(color = "black", size = .5),
          axis.text    = element_text(color = "black", size = 10),
          axis.title   = element_text(color = "black", size = 12))

  return(p)
}


#' generate volcano plot for DE analysis
#'
#' @param data A data frame for the full version of data, required
#' @param labels string or list, ids to lable on plot
#' @param showSig boolen, whether or not show significant genes in colors
#' @param fcCutoff numeric, cutoff for foldchange
#' @param pvalCutoff numeric, cutoff for pvalue
#'
#' @import ggplot2
#' @import ggrepel
#'
#' @export
#'
plot_ma <- function(data, labels = NULL, showSig = TRUE,
                    fcCutoff = 2, pvalCutoff = 0.05) {
  # data_sig is required
  stopifnot(is.data.frame(data))

  stopifnot(all(c("id", "baseMean", "log2FoldChange", "pvalue") %in% colnames(data)))
  if ("padj" %in% colnames(data)) {
    data$pvalue <- data$padj
  }

  # convert num
  fcCutoff   <- as.numeric(fcCutoff)
  pvalCutoff <- as.numeric(pvalCutoff)

  ## add significant tags
  plotDf <- sig_genes(data, fc = fcCutoff, pval = pvalCutoff, type = "all") # sig genes
  plotDf <- dplyr::rename(plotDf, logFC = log2FoldChange, pvalue = pvalue)

  # split data
  if(isFALSE(showSig) | ! "sig" %in% colnames(plotDf)) {
    plotDf$sig <- "gene"
  }

  # log10
  plotDf$baseMean <- log10(plotDf$baseMean + 1)

  # axis title
  xTitle <- expression(paste(log["10"], "(baseMean)"))
  yTitle <- expression(paste(log["2"], "(Fold-change)"))

  # basic plot
  p <- ggplot(plotDf, aes(baseMean, logFC, color = sig)) +
    geom_hline(yintercept = 0, color = "grey10", size = .7) +
    geom_point(size = .4)

  # axis
  p <- p + scale_y_continuous(limits = c(-4, 4))

  # add colors
  sigColors <- c("green4", "grey70", "red", "grey30")
  names(sigColors) <- c("down", "not", "up", "gene")
  sigLevels <- sort(unique(plotDf$sig))
  hitColors <- sigColors[sigLevels]
  p <- p + scale_color_manual(values = hitColors)

  # show labels
  labelDf <- plotDf %>% dplyr::filter(id %in% labels)
  if (nrow(labelDf) > 0) {
    p <- p +
      geom_text_repel(
        mapping       = aes(baseMean, logFC, color = sig),
        data          = labelDf,
        label         = labelDf$id,
        nudge_x       = 1.2,
        nudge_y       = .02,
        point.padding = .4,
        box.padding   = .4,
        segment.size  = .6,
        segment.color = "grey60",
        direction     = "both")
  }

  # theme
  p <- p +
    xlab(xTitle) +
    ylab(yTitle) +
    labs(color = "Significant") +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = .7),
          plot.title   = element_text(color = "black", hjust = .5, size = 14),
          panel.grid   = element_blank(),
          axis.line    = element_blank(),
          axis.ticks.length = unit(.2, "cm"),
          axis.ticks   = element_line(color = "black", size = .5),
          axis.text    = element_text(color = "black", size = 10),
          axis.title   = element_text(color = "black", size = 12))

  return(p)
}



#' generate volcano plot for DE analysis
#'
#' @param data A data frame for the full version of data
#' @param labels list, ids to label in plot
#' @param showSig boolen, whether or not show significant genes in colors
#' @param fcCutoff numeric, cutoff for foldchange
#' @param pvalCutoff numeric, cutoff for pvalue
#'
#' @import ggplot2
#' @import ggrepel
#'
#' @export
#'
plot_volcano <- function(data, labels = NULL, showSig = TRUE,
                         fcCutoff = 2, pvalCutoff = 0.05) {
  # data_sig not required
  stopifnot(is.data.frame(data))
  stopifnot(all(c("id", "log2FoldChange", "pvalue") %in% colnames(data)))
  if ("padj" %in% colnames(data)) {
    data$pvalue <- data$padj
  }

  # convert num
  fcCutoff   <- as.numeric(fcCutoff)
  pvalCutoff <- as.numeric(pvalCutoff)

  ## add significant tags
  plotDf <- sig_genes(data, fc = fcCutoff, pval = pvalCutoff, type = "all") # sig genes
  # plotDf <- dplyr::rename(plotDf, logFC = log2FoldChange, pvalue = pvalue)

  # split data
  if(isFALSE(showSig) | ! "sig" %in% colnames(plotDf)) {
    plotDf$sig <- "gene"
  }

  # axis title
  xtitle <- expression(paste(log["2"], "(mean ratio of Treatment/Control)"))
  ytitle <- expression(paste(-log["10"], "(", italic("P"), " value)"))

  # convert pval to -log10()
  plotDf$logPval <- -log10(plotDf$pvalue)

  # basic plot
  p <- ggplot(plotDf, aes(log2FoldChange, logPval, color = sig)) +
    geom_vline(xintercept = 0, size = .5, color = "black") +
    geom_vline(xintercept = c(-log2(fcCutoff), log2(fcCutoff)),
               size = .5,
               linetype = 2,
               color = "grey60") +
    geom_hline(yintercept = -log10(pvalCutoff),
               size = .5,
               linetype = 2,
               color = "grey60")

  # add points
  p <- p + geom_point(size = .5)

  # add colors
  sig_colors <- c("green4", "grey70", "red", "grey30")
  names(sig_colors) <- c("down", "not", "up", "gene")
  sig_levels <- sort(unique(plotDf$sig))
  hit_colors <- sig_colors[sig_levels]
  p <- p + scale_color_manual(values = hit_colors)

  # show labels
  data_label <- plotDf %>% dplyr::filter(id %in% labels)
  if (nrow(data_label) > 0) {
    p <- p + geom_text_repel(
      mapping       = aes(log2FoldChange, logPval, color = sig),
      data          = data_label,
      label         = data_label$id,
      nudge_x       = 1.2,
      nudge_y       = .02,
      point.padding = .4,
      box.padding   = .4,
      segment.size  = .6,
      segment.color = "grey60",
      direction     = "both")
  }

  # theme
  p <- p +
    xlab(xtitle) +
    ylab(ytitle) +
    labs(color = "Significant") +
    theme_bw()   +
    theme(panel.border = element_blank(),
          panel.grid   = element_blank(),
          axis.line    = element_line(color = "black", size = .5),
          axis.ticks   = element_line(color = "black", size = .5),
          axis.text    = element_text(color = "black", size = 10),
          axis.title   = element_text(color = "black", size = 12),
          axis.ticks.length = unit(.2, "cm"))

  return(p)
}

