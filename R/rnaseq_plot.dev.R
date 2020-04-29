#'
#' ##----------------------------------------------------------------------------##
#' ## Extra plots
#' ## publish quality plots
#' ##
#' ## quality control plots...
#' ##
#' ## import: transcripts_deseq2.csv file
#' ##
#'
#'
#' ## 1. scatter plot 1
#'
#' #' scatter plot for DESeq
#' #'
#' #' Description:
#' #' 1. dots, scater plot: x=control, y=treated
#' #' 2. lines, 2-fold change, dashed
#' #' 3. abline,
#' #' 4. log10(score + 0.5) on axis
#' #' 5. colors, up/red, down/green, not/grey
#' #' 6. add numbers n=
#' #' 7. add extra_labels c()
#' #' 8. themes()
#' #'
#' #' make scatter plot
#' #'
#' #' @param csv_file file the deseq.csv file, output of deseqNode()
#' #' @param x_label show labels on x-axis, default: auto
#' #' @param y_label show labels on y-axis, default: auto
#' #' @param sig_in_color bool, whether show significant dots in colors, default: TRUE
#' #' @param sig_in_text vector, list of ids to show, labes in plot, default: NULL
#' #' @param fc_cutoff float, foldchange cutoff, default: 2
#' #' @param pvalue_cutoff float, default: 0.05
#' #' @param save2pdf bool or file,
#' #'
#' #'
#' #' @import ggplot2
#' #' @import dplyr
#' #'
#' #' @export
#' deseqScatterPlot1 <- function(csv_file,
#'                               x_label       = NULL,
#'                               y_label       = NULL,
#'                               fc_cutoff     = 2,
#'                               pvalue_cutoff = 0.05,
#'                               save2pdf      = FALSE){
#'   # convert to data.frame
#'   stopifnot(file.exists(csv_file))
#'   df <- deseqCsvReader(csv_file)
#'
#'   # convert to log10 (log10(v + 0.5))
#'   smp_ctl <- colnames(df)[2]
#'   smp_exp <- colnames(df)[3]
#'
#'   # convert to log10
#'   df2 <- df %>%
#'     dplyr::mutate_(x = smp_ctl,
#'                    y = smp_exp) %>%
#'     dplyr::mutate(x = log10(x + 0.5),
#'                   y = log10(y + 0.5))
#'
#'   # add sig
#'   df2 <- deseqSig(df2)
#'   # # add sig: up, down, not
#'   # df2 <- df2 %>%
#'   #   dplyr::mutate(
#'   #     sig = ifelse(padj <= pvalue_cutoff,
#'   #                  ifelse(log2FoldChange >= log2(fc_cutoff),
#'   #                         "up",
#'   #                         ifelse(log2FoldChange <= -log2(fc_cutoff),
#'   #                                "down",
#'   #                                "not")),
#'   #                  "not")) %>%
#'   #   dplyr::mutate(sig = factor(sig, levels = c("up", "not", "down")))
#'
#'   # prepare axis labels
#'   x_label <- ifelse(is.null(x_label), smp_ctl, x_label)
#'   y_label <- ifelse(is.null(y_label), smp_exp, y_label)
#'
#'   p <- scatterPlot1(df2, x_label, y_label)
#'   return(p)
#' }
#'
#'
#'
#' #' scatter plot for DESeq
#' #'
#' #' Description:
#' #' 1. dots, scater plot: x=control, y=treated
#' #' 2. lines, 2-fold change, dashed
#' #' 3. abline,
#' #' 4. log10(score + 0.5) on axis
#' #' 5. colors, up/red, down/green, not/grey
#' #' 6. add numbers n=
#' #' 7. add extra_labels c()
#' #' 8. themes()
#' #'
#' #' @param data
#' #' @param x_label
#' #' @param y_label
#' #'
#' #' @import ggplot2
#' #' @import scales
#' #' @import dplyr
#' #'
#' #' @export
#' scatterPlot1 <- function(data,
#'                          x_label      = NULL,
#'                          y_label      = NULL,
#'                          axis_max     = 5,
#'                          axis_unit    = 1,
#'                          save2pdf     = TRUE,
#'                          sig_label_show    = 10,
#'                          extra_point_text  = NULL,
#'                          extra_point_size  = 1,
#'                          extra_point_color = "red"){
#'   # check data
#'   stopifnot(is.data.frame(data))
#'   stopifnot(all(c("id", "x", "y", "sig") %in% colnames(data)))
#'   stopifnot(is.character(data$id))
#'   stopifnot(is.numeric(data$x) & is.numeric(data$y))
#'
#'   # axis title
#'   x_title <- bquote(.(x_label)~"["*log["10"]~"FPKM + 0.5]")
#'   y_title <- bquote(.(y_label)~"["*log["10"]~"FPKM + 0.5]")
#'
#'   # sig labels
#'   df_sig <- data %>%
#'     dplyr::filter(sig %in% c("up", "down")) %>%
#'     dplyr::arrange(padj) %>%
#'     head(sig_label_show)
#'
#'   # prepare data
#'   #
#'   # basic part
#'   # make basic plot
#'   p <- ggplot(data, aes(x, y, color = sig)) +
#'     geom_point() +
#'     scale_color_manual(values = c("red2", "grey50", "green2", "grey50")) +
#'     scale_size_manual(values = c(1, 0.5, 1))
#'
#'   # add title
#'   p <- p +
#'     xlab(x_title) + ylab(y_title)
#'
#'   # add sig labels
#'   p <- p +
#'     geom_text_repel(
#'       data          = df_sig,
#'       label         = df_sig$id,
#'       size          = 4,
#'       box.padding   = .5,
#'       segment.size  = 0.5,
#'       segment.color = "grey50",
#'       direction     = "both")
#'
#'   # extra labels
#'   if(! is.null(extra_point_text)) {
#'     df_extra <- dplyr::filter(data, id %in% extra_point_text)
#'     # add text
#'     p <- p +
#'       geom_text_repel(
#'         data          = df_extra,
#'         label         = df_extra$id,
#'         size          = 4,
#'         box.padding   = .5,
#'         segment.size  = 0.5,
#'         segment.color = "grey50",
#'         direction     = "both")
#'     # add point
#'     p <- p +
#'       geom_point(data  = df_extra,
#'                  color = extra_point_color,
#'                  size  = extra_point_size)
#'   }
#'
#'   # axis breaks
#'   p <- p +
#'     scale_x_continuous(limits = c(0, axis_max),
#'                        breaks = seq(0, axis_max, axis_unit),
#'                        labels = seq(0, axis_max, axis_unit),
#'                        expand = c(0, 0)) +
#'     scale_y_continuous(limits = c(0, axis_max),
#'                        breaks = seq(0, axis_max, axis_unit),
#'                        labels = seq(0, axis_max, axis_unit),
#'                        expand = c(0, 0))
#'
#'   # add ablines
#'   p <- p +
#'     geom_abline(slope = 1, intercept = 0, color = "black", size = .5) +
#'     geom_abline(slope = 1, intercept = log10(2), color = "black",
#'                 linetype = 2, size = .5) +
#'     geom_abline(slope = 1, intercept = -log10(2), color = "black",
#'                 linetype = 2, size = .5)
#'
#'   # add theme
#'   p <- p +
#'     theme_linedraw() +
#'     theme(panel.border = element_rect(color = "black", fill = NA, size = .7),
#'           plot.title   = element_text(color = "black", hjust = .5, size = 12),
#'           panel.grid   = element_blank(),
#'           axis.line    = element_blank(),
#'           axis.ticks.length = unit(.15, "cm"),
#'           axis.ticks   = element_line(color = "black", size = .5),
#'           axis.text    = element_text(color = "black", size = 10),
#'           axis.title   = element_text(color = "black", size = 12))
#'
#'   # whether save to pdf
#'   if (is.character(save2pdf)) {
#'     pdf_dir <- dirname(savd2pdf)
#'     if (! dir.exists(pdf_dir)) {
#'       warning("directory of pdf file not exists")
#'     } else {
#'       pdf(rpt_file, width = 4, height = 4, paper = "a4")
#'       print(p)
#'       dev.off()
#'       print(paste("save plot to file:", save2pdf))
#'     }
#'   }
#'
#'   return(p)
#' }
#'
#'
#'
#' #' MA plot
#' #'
#' #' Description:
#' #' 1. dots, scater plot: x=control, y=treated
#' #' 2. lines, 2-fold change, dashed
#' #' 3. abline,
#' #' 4. log10(score + 0.5) on axis
#' #' 5. colors, up/red, down/green, not/grey
#' #' 6. add numbers n=
#' #' 7. add extra_labels c()
#' #' 8. themes()
#' #' make scatter plot
#' #'
#' #' @param csv_file file the deseq.csv file, output of deseqNode()
#' #' @param x_label show labels on x-axis, default: auto
#' #' @param y_label show labels on y-axis, default: auto
#' #' @param sig_in_color bool, whether show significant dots in colors, default: TRUE
#' #' @param sig_in_text vector, list of ids to show, labes in plot, default: NULL
#' #' @param fc_cutoff float, foldchange cutoff, default: 2
#' #' @param pvalue_cutoff float, default: 0.05
#' #' @param save2pdf bool or file,
#' #'
#' #'
#' #' @import ggplot2
#' #' @import dplyr
#' #'
#' #' @export
#' #'
#' deseqScatterPlot1 <- function(csv_file,
#'                               x_label       = NULL,
#'                               y_label       = NULL,
#'                               fc_cutoff     = 2,
#'                               pvalue_cutoff = 0.05,
#'                               save2pdf      = FALSE){
#'   # convert to data.frame
#'   stopifnot(file.exists(csv_file))
#'   df <- deseqCsvReader(csv_file)
#'
#'   # convert to log10 (log10(v + 0.5))
#'   smp_ctl <- colnames(df)[2]
#'   smp_exp <- colnames(df)[3]
#'
#'   # convert to log10
#'   df2 <- df %>%
#'     dplyr::mutate_(x = smp_ctl,
#'                    y = smp_exp) %>%
#'     dplyr::mutate(x = log10(x + 0.5),
#'                   y = log10(y + 0.5))
#'
#'   # add sig
#'   df2 <- deseqSig(df2)
#'
#'   # prepare axis labels
#'   x_label <- ifelse(is.null(x_label), smp_ctl, x_label)
#'   y_label <- ifelse(is.null(y_label), smp_exp, y_label)
#'
#'   p <- scatterPlot1(df2, x_label, y_label)
#'   return(p)
#' }
#'
#'
#'
#' #' generate volcano plot for DE analysis
#' #'
#' #' @param data A data frame for the full version of data, required
#' #' @param labels string or list, ids to lable on plot
#' #' @param showSig boolen, whether or not show significant genes in colors
#' #' @param fcCutoff numeric, cutoff for foldchange
#' #' @param pvalCutoff numeric, cutoff for pvalue
#' #'
#' #' @import ggplot2
#' #' @import ggrepel
#' #'
#' #' @export
#' #'
#' maPlot <- function(data,
#'                    labels       = NULL,
#'                    showSig      = TRUE,
#'                    fcCutoff     = 2,
#'                    pvalueCutoff = 0.05) {
#'   # data_sig is required
#'   stopifnot(is.data.frame(data))
#'
#'   stopifnot(all(c("id", "baseMean", "log2FoldChange", "pvalue") %in% colnames(data)))
#'   if ("padj" %in% colnames(data)) {
#'     data$pvalue <- data$padj
#'   }
#'
#'   # convert num
#'   fcCutoff   <- as.numeric(fcCutoff)
#'   pvalCutoff <- as.numeric(pvalueCutoff)
#'
#'   ## add significant tags
#'   df_sig <- deseqSig(data, fc = fcCutoff, pval = pvalCutoff, type = "all") # sig genes
#'   df_sig <- dplyr::rename(df_sig, logFC = log2FoldChange, pvalue = pvalue)
#'
#'   # split data
#'   if(isFALSE(showSig) | ! "sig" %in% colnames(plotDf)) {
#'     df_sig$sig <- "gene"
#'   }
#'
#'   # log10
#'   df_sig$baseMean <- log10(plotDf$baseMean + 1)
#'
#'   # axis title
#'   xTitle <- expression(paste(log["10"], "(baseMean)"))
#'   yTitle <- expression(paste(log["2"], "(Fold-change)"))
#'
#'   # basic plot
#'   p <- ggplot(df_sig, aes(baseMean, logFC, color = sig)) +
#'     geom_hline(yintercept = 0, color = "grey10", size = .7) +
#'     geom_point(size = .4)
#'
#'   # axis
#'   p <- p + scale_y_continuous(limits = c(-4, 4))
#'
#'   # add colors
#'   sigColors <- c("green4", "grey70", "red", "grey30")
#'   names(sigColors) <- c("down", "not", "up", "gene")
#'   sigLevels <- sort(unique(plotDf$sig))
#'   hitColors <- sigColors[sigLevels]
#'   p <- p + scale_color_manual(values = hitColors)
#'
#'   # show labels
#'   df_label <- df_sig %>% dplyr::filter(id %in% labels)
#'   if (nrow(labelDf) > 0) {
#'     p <- p +
#'       geom_text_repel(
#'         mapping       = aes(baseMean, logFC, color = sig),
#'         data          = labelDf,
#'         label         = labelDf$id,
#'         nudge_x       = 1.2,
#'         nudge_y       = .02,
#'         point.padding = .4,
#'         box.padding   = .4,
#'         segment.size  = .6,
#'         segment.color = "grey60",
#'         direction     = "both")
#'   }
#'
#'   # theme
#'   p <- p +
#'     xlab(xTitle) +
#'     ylab(yTitle) +
#'     labs(color = "Significant") +
#'     theme_classic() +
#'     theme(panel.border = element_rect(color = "black", fill = NA, size = .7),
#'           plot.title   = element_text(color = "black", hjust = .5, size = 14),
#'           panel.grid   = element_blank(),
#'           axis.line    = element_blank(),
#'           axis.ticks.length = unit(.2, "cm"),
#'           axis.ticks   = element_line(color = "black", size = .5),
#'           axis.text    = element_text(color = "black", size = 10),
#'           axis.title   = element_text(color = "black", size = 12))
#'
#'   return(p)
#' }
#'
#'
#'
#'
#' #' generate volcano plot for DE analysis
#' #'
#' #' @param data A data frame for the full version of data
#' #' @param labels list, ids to label in plot
#' #' @param showSig boolen, whether or not show significant genes in colors
#' #' @param fcCutoff numeric, cutoff for foldchange
#' #' @param pvalCutoff numeric, cutoff for pvalue
#' #'
#' #' @import ggplot2
#' #' @import ggrepel
#' #'
#' #' @export
#' #'
#' volcanoPlot <- function(data,
#'                         labels     = NULL,
#'                         showSig    = TRUE,
#'                         fcCutoff   = 2,
#'                         pvalCutoff = 0.05) {
#'   # data_sig not required
#'   stopifnot(is.data.frame(data))
#'   stopifnot(all(c("id", "log2FoldChange", "pvalue") %in% colnames(data)))
#'   if ("padj" %in% colnames(data)) {
#'     data$pvalue <- data$padj
#'   }
#'
#'   # convert num
#'   fcCutoff   <- as.numeric(fcCutoff)
#'   pvalCutoff <- as.numeric(pvalCutoff)
#'
#'   ## add significant tags
#'   plotDf <- sig_genes(data, fc = fcCutoff, pval = pvalCutoff, type = "all") # sig genes
#'   # plotDf <- dplyr::rename(plotDf, logFC = log2FoldChange, pvalue = pvalue)
#'
#'   # split data
#'   if(isFALSE(showSig) | ! "sig" %in% colnames(plotDf)) {
#'     plotDf$sig <- "gene"
#'   }
#'
#'   # axis title
#'   xtitle <- expression(paste(log["2"], "(mean ratio of Treatment/Control)"))
#'   ytitle <- expression(paste(-log["10"], "(", italic("P"), " value)"))
#'
#'   # convert pval to -log10()
#'   plotDf$logPval <- -log10(plotDf$pvalue)
#'
#'   # basic plot
#'   p <- ggplot(plotDf, aes(log2FoldChange, logPval, color = sig)) +
#'     geom_vline(xintercept = 0, size = .5, color = "black") +
#'     geom_vline(xintercept = c(-log2(fcCutoff), log2(fcCutoff)),
#'                size = .5,
#'                linetype = 2,
#'                color = "grey60") +
#'     geom_hline(yintercept = -log10(pvalCutoff),
#'                size = .5,
#'                linetype = 2,
#'                color = "grey60")
#'
#'   # add points
#'   p <- p + geom_point(size = .5)
#'
#'   # add colors
#'   sig_colors <- c("green4", "grey70", "red", "grey30")
#'   names(sig_colors) <- c("down", "not", "up", "gene")
#'   sig_levels <- sort(unique(plotDf$sig))
#'   hit_colors <- sig_colors[sig_levels]
#'   p <- p + scale_color_manual(values = hit_colors)
#'
#'   # show labels
#'   data_label <- plotDf %>% dplyr::filter(id %in% labels)
#'   if (nrow(data_label) > 0) {
#'     p <- p + geom_text_repel(
#'       mapping       = aes(log2FoldChange, logPval, color = sig),
#'       data          = data_label,
#'       label         = data_label$id,
#'       nudge_x       = 1.2,
#'       nudge_y       = .02,
#'       point.padding = .4,
#'       box.padding   = .4,
#'       segment.size  = .6,
#'       segment.color = "grey60",
#'       direction     = "both")
#'   }
#'
#'   # theme
#'   p <- p +
#'     xlab(xtitle) +
#'     ylab(ytitle) +
#'     labs(color = "Significant") +
#'     theme_bw()   +
#'     theme(panel.border = element_blank(),
#'           panel.grid   = element_blank(),
#'           axis.line    = element_line(color = "black", size = .5),
#'           axis.ticks   = element_line(color = "black", size = .5),
#'           axis.text    = element_text(color = "black", size = 10),
#'           axis.title   = element_text(color = "black", size = 12),
#'           axis.ticks.length = unit(.2, "cm"))
#'
#'   return(p)
#' }
#'
#'
#'
#'
#'
