#' Functions for DESeq2, dds quality control
#'
#' Reading files, directories, pipeline, ...
#' Processing data
#' Prepare for plot
#' Plotting
#'
#' to-to: simpify arguments by `...`
#'
#' @name deseq_qc




#------------------------------------------------------------------------------#
# quality check for dds object

#' @describeIn deseq_qc quality control for dds object
#'
#' @param x character file `norm_table.fix.csv`, or `DESeqDataSet`
#' also support the following data types: (be careful)
#' `data.frame`, `.rds`#'
#' @param outdir character saving the results
#' @param log2fc_limits numeric setting the range of log2FoldChange on plots
#' default: NULL, by `scales::breaks_extend()`
#' @param transform bool compute tranformed data, `vst()`, `rlog()`
#' default: FALSE
#' @param n_max integer number of genes to display, defualt: 16
#' @param fc numeric cutoff for foldchange, default: 1, ignore foldchange
#' @param pvalue numeric cutoff for padj, default: 0.1, the main criteria
#' @param p_adjust bool use p-adjust value instead
#' @param overwrite bool overwrite exists file, default: FALSE
#'
#' @import DESeq2
#'
#' @return list of ggplots
#'
#' @export
deseq_qc <- function(x, outdir = NULL, log2fc_limits = c(-4, 4),
                     transform = TRUE, n_max = 20,
                     fc = 1, pvalue = 0.1, p_adjust = TRUE,
                     label_list = NULL, label_max = 8,
                     overwrite = FALSE) {
  #-- qc: dds
  qc_dds <- deseq_qc_dds(x, outdir, transform, return_dds = TRUE)
  qc_res <- deseq_qc_res(x, outdir, log2fc_limits, fc, pvalue, p_adjust,
                         overwrite = overwrite)
  if(inherits(qc_dds, "DESeqDataSet")) {
    qc1 <- list(
      counts   = deseq_qc_counts(qc_dds, outdir, transform, n_max = 12,
                                 overwrite = overwrite),
      mean_sd  = deseq_qc_mean_sd(qc_dds, outdir, transform, overwrite),
      top_gene = deseq_qc_top_gene(qc_dds, outdir, transform, n_max, overwrite),
      dist     = deseq_qc_dist(qc_dds, outdir, transform, overwrite),
      pca      = deseq_qc_pca(qc_dds, outdir, transform, overwrite))
  } else {
    qc1 <- list(counts = NULL, mean_sd = NULL, top_gene = NULL, dist = NULL,
                pca = NULL)
  }
  #-- qc: res
  qc_res <- deseq_qc_res(x, outdir, log2fc_limits, fc, pvalue, p_adjust,
                         overwrite = overwrite)
  if(inherits(qc_res, "data.frame")) {
    qc2 <- list(
      # ma      = deseq_qc_ma(x, outdir, log2fc_limits, fc, pvalue, p_adjust,
      #                       label_list, label_max, overwrite),
      volcano = deseq_qc_volcano(x, outdir, log2fc_limits, fc, pvalue, p_adjust,
                                 label_list, label_max, overwrite),
      # scatter = deseq_qc_scatter(x, outdir, fc, pvalue, p_adjust,
      #                            label_list, label_max, overwrite)
    )
    print(label_list)
  } else {
    qc2 <- list(ma = NULL, volcano = NULL, scatter = NULL)
  }
  #-- return
  c(qc1, qc2)
}



#------------------------------------------------------------------------------#
# loading data for deseq_qc()

#' @describeIn deseq_qc_dds
#' check input data for deseq_qc
#'
#' @param x `DESeqDataSet`
#' @param outdir character saving the results, default: NULL
#' @param transform bool compute tranformed data, `vst()`, `rlog()`
#' default: FALSE
#' @param return_dds bool return `DESeqDataSet` instead
#' @param overwrite bool overwrite exists file, default: FALSE
#'
#' @return list of standard*, vst, rlog objects
#'
#' @export
deseq_qc_dds <- function(x = NULL, outdir = NULL, transform = TRUE,
                         return_dds = FALSE, overwrite = FALSE) {
  #-- Check: x
  dds <- NULL
  if(inherits(x, "DESeqDataSet")) {
    dds <- x
  } else if(inherits(x, "character")) {
    if(file.exists(x) & endsWith(x, ".rds")) {
      dds <- readRDS(x) # update
    }
  }
  if(inherits(outdir, "character")) {
    fs <- file.path(outdir, "deseq_dds.rds")
    if(file.exists(fs)) {
      dds <- readRDS(fs) # update
    }
  }
  #-- Check: dds
  if(!inherits(dds, "DESeqDataSet")) {
    warning(glue::glue(
      "x is {class(dds)}, expect `DESeqDataSet`; ",
      "specify either `x=` or `outdir=`, or both "))
    return(NULL)
  }
  qc_dds_rds <- file.path(outdir, "deseq_qc_dds.rds")
  out <- NULL
  if(isTRUE(transform)) {
    if(inherits(qc_dds_rds, "character")) {
      if(file.exists(qc_dds_rds)) {
        message(glue::glue("loading deseq_qc_dds from file: {qc_dds_rds}"))
        out <- readRDS(qc_dds_rds) # from file
      }
    }
    if(!inherits(out, "list")) {
      # dds <- DESeq2::DESeq(dds) # check
      message("use `outdir=<chr>`, as `vst()` or `rlog()` takes long time")
      if(nrow(dds) > 1000) {
        # 1000 genes? see https://support.bioconductor.org/p/98634/#98637
        out <- list(
          standard = DESeq2::normTransform(dds),
          vst      = DESeq2::vst(dds, blind = FALSE),
          rlog     = DESeq2::rlog(dds, blind = FALSE)
        )
      } else {
        out <- list(standard = DESeq2::normTransform(dds))#
      }
    }
  } else {
    out <- list(standard = DESeq2::normTransform(dds)) # direct output
  }
  #-- save to file
  if(inherits(qc_dds_rds, "character")) {
    if(check_path(outdir)) {
      if(!file.exists(qc_dds_rds) | overwrite) {
        saveRDS(out, qc_dds_rds)
      }
    }
  }
  #-- return
  if(isTRUE(return_dds)) {
    dds
  } else {
    out
  }
}


#' @describeIn deseq_qc_res
#' get data.frame from `norm_table.fix.xls`, or `DESeqDataSet`,
#' for `dds`, see `import_featurecounts()` and `prep_hiseq_deseq()`
#' be careful using the following compatiable data types:
#' `data.frame`, `deseq_dds.rds`,
#'
#' return: data.frame
#' for scatter plot, require mean values of wt/mut
#' x: log10(wt + 1), names(df)[2]
#' y: log10(mut + 1), names(df)[3]
#'
#' assign `sig` column based on log2FoldChange, padj (or pvalue)
#' see: `get_sig_name(..., return_dataframe = TRUE)`
#' sig: for colors, ["up":"red", "not":"grey50", "down":"blue"]
#'   up: padj < pval & log2fc >= log2(fc)
#'   down: padj < pval & log2fc <= -log2(fc)
#'   not: is.na(padj) | padj >= pval
#'
#' ext: for shapes, ["up":2, "dot":20, "down":6]
#'   change dot shapes, based on the range of log2fc_limits
#'   default limits = `scales::breaks_ext(n=5)(log2FoldChange)`
#'   up: log2fc > max(log2fc_limits)
#'   down: log2fc < min(log2fc_limits)
#'   not: log2fc >= min(log2fc_limits) & log2fc <= max(log2fc_limits)
#'
#' @param x character path to the file `norm_table.fix.xls`
#' also support data types: `DESeqDataSet`, ...
#' @param outdir character saving the results
#' @param log2fc_limits numeric setting the range of log2FoldChange on plots
#' default: NULL, by `scales::breaks_extend()`
#' @param fc numeric cutoff for foldchange, default: 1, ignore foldchange
#' @param pvalue numeric cutoff for padj, default: 0.1, the main criteria
#' @param p_adjust bool use p-adjust value instead
#' @param overwrite bool overwrite exists file, default: FALSE
#'
#' @improt ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom patchwork wrap_plots plot_annotation
#'
#' @return ggplot
#'
#' @export
deseq_qc_res <- function(x, outdir = NULL, log2fc_limits = NULL,
                         fc = 1, pvalue = 0.1, p_adjust = TRUE,
                         overwrite = FALSE) {
  #-- Check: loading dds, DESeqDataSet obj
  dds <- NULL # init
  if(inherits(x, "data.frame")) {
    dds <- x #
  } else if(inherits(x, "DESeqDataSet")) {
    dds <- x # overwrite
  } else if(inherits(x, "character")) {
    if(file.exists(x)) {
      if(endsWith(x, ".rds")) {
        dds <- readRDS(x)
      } else if(endsWith(x, ".csv")) {
        dds <- read.csv(x)
      } else {
        message(glue::glue("unknown x: {x}"))
      }
    }
  }
  #-- Check: loading from outdir, priority: outdir > x
  if(inherits(outdir, "character")) {
    fs <- file.path(outdir, "norm_table.fix.xls") # dds
    if(file.exists(fs)) {
      dds <- read.csv(fs) # fixed data
    }
  }
  #-- Check: dds
  if(inherits(dds, "data.frame")) {
    df <- dds
  } else if(inherits(dds, "DESeqDataSet")) {
    df <- deseq_mean(dds, outdir) # see: deseq_utils.R
  } else {
    warning(glue::glue(
      "x is {class(x)}, expect `DESeqDataSet` or `deseq_dds.rds` file, ",
      "specify either `x=` or `outdir=`, or both "))
    return(NULL)
  }
  #-- Check: df
  if(inherits(df, "data.frame")) {
    rc <- c("baseMean", "log2FoldChange", "padj")
    if(!all(rc %in% names(df))) {
      rc_str <- paste(rc, collapse = ", ")
      warning(glue::glue("missing columns: [{rc_str}]"))
      return(NULL)
    }
  } else {
    warning(glue::glue("failed loading data from x={x}, outdir={outdir}"))
    return(NULL)
  }
  #-- Check: wt, mut
  if(inherits(dds, "DESeqDataSet")) {
    coldata <- SummarizedExperiment::colData(dds)
    wt  <- levels(coldata$condition)[1]
    mut <- levels(coldata$condition)[2]
  } else {
    # gene_id, wt, mut, ...
    wt  <- names(df)[2]
    mut <- names(df)[3]
  }
  if(!all(c(wt, mut) %in% names(df))) {
    warning(glue::glue("missing reuired columns: {wt}, {mut}"))
    return(NULL)
  }
  #-- Check: limits for x, y
  if(inherits(log2fc_limits, "numeric")) {
    breaks <- log2fc_limits
  } else {
    breaks <- scales::breaks_extended(n = 5)(df$log2FoldChange)
  }
  #-- run: sig for colors, up, not, down
  # in case, pvalue exists in data.frame
  if(!rlang::has_name(df, "sig")) {
    df <- get_sig_name(df, fc, pvalue, p_adjust, return_dataframe = TRUE)
  }
  #-- run: ext for shapes, up, dot, down
  # fix lgo2fc, by shapes
  out <- df %>%
    dplyr::mutate(
      log10pval = -log10(padj),
      log10basemean = log10(baseMean + 1),
      ext = ifelse(
        is.na(log2FoldChange), "dot", ifelse(
          log2FoldChange > max(breaks), "up", ifelse(
            log2FoldChange < min(breaks), "down", "dot")))) %>%
    dplyr::mutate(log2fc = ifelse(
      ext == "up", max(breaks), ifelse(
        ext == "down", min(breaks), log2FoldChange)))
  #-- save to file
  qc_res_rds <- file.path(outdir, "deseq_qc_res.rds")
  if(inherits(qc_res_rds, "character")) {
    if(check_path(outdir)) {
      if(!file.exists(qc_res_rds) | overwrite) {
        saveRDS(out, qc_res_rds)
      }
    }
  }
  #-- return
  out
}


#------------------------------------------------------------------------------#
# main: qc plots

#' @describeIn deseq_qc_counts
#' check counts for each genes
#'
#' @param x DESeqDataSet
#' @param outdir character saving the results
#' @param transform bool compute tranformed data, `vst()`, `rlog()`
#' default: FALSE
#' @param n_max integer number of genes to display, defualt: 16
#' @param res `DESeqResults` object, choose top genes by `padj`, default: NULL
#' @param overwrite bool overwrite exists file, default: FALSE
#'
#' @return ggplot
#'
#' @export
deseq_qc_counts <- function(x, outdir = NULL, transform = TRUE, n_max = 16,
                            res = NULL, overwrite = FALSE) {
  #-- Check: load pre-data
  qc_data <- deseq_qc_dds(x, outdir, transform)
  dds <- x
  if(inherits(qc_data, "NULL")) {
    return(NULL)
  }
  #-- Check: coldata
  coldata <- colData(dds) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("group")
  vars <- setNames(coldata$condition, nm = coldata$group)
  #-- Check: res, for top sig genes
  if(inherits(res, "DESeqResults")) {
    top_gene <- rownames(res)[order(res$padj)[1:n_max]]
  } else {
    top_gene <- order(rowMeans(DESeq2::counts(dds, normalized = TRUE)),
                      decreasing = TRUE)[1:n_max]
  }
  #-- run: plot
  df <- DESeq2::counts(dds, normalized = TRUE)[top_gene, ] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(-gene, names_to = "group", values_to = "count") %>%
    dplyr::mutate(condition = dplyr::recode(group, !!!vars))
  # plot
  p <- ggplot(df, aes(condition, count, color = condition)) +
    geom_point(position = position_jitter(width = 0.1, height = 0), size = 1) +
    facet_wrap(~gene, ncol = 4) +
    theme_bw() +
    theme(legend.position = "none")
  # save plot
  if(inherits(outdir, "character")) {
    p_png <- file.path(outdir, "deseq_qc_counts.png")
    if(check_path(outdir)) {
      if(!file.exists(p_png) | overwrite) {
        ggsave(p_png, p, width = 8, height = 8, dpi = 200, units = "in")
      }
    }
  }
  p # return data
}


#' @describeIn deseq_qc_mean_sd
#' check mean standard diviation (SD), after transformation by vst, vlog
#'
#' @param x DESeqDataSet
#' @param outdir character saving the results
#' @param transform bool compute tranformed data, `vst()`, `rlog()`
#' default: FALSE
#' @param overwrite bool overwrite exists file, default: FALSE
#'
#' @importFrom vsn meanSdPlot
#' @importFrom patchwork wrap_plots plot_annotation
#'
#' @return ggplot
#'
#' @export
deseq_qc_mean_sd <- function(x, outdir = NULL, transform = TRUE,
                             overwrite = FALSE) {
  qc_data <- deseq_qc_dds(x, outdir, transform)
  if(is.null(qc_data)) {
    return(NULL)
  }
  p_list <- lapply(names(qc_data), function(i) {
    tf   <- vsn::meanSdPlot(assay(qc_data[[i]]), plot = FALSE)
    if(rlang::has_name(tf, "gg")) {
      tf$gg +
        ggtitle(i) +
        theme_bw()
    }
  })
  # combine plots
  if(length(p_list) > 1) {
    p <- patchwork::wrap_plots(p_list, ncol = 2, guides = "collect") +
      patchwork::plot_annotation(title = "Effects of transformations")
    w <- 8 # fig
    h <- 8 # fig
  } else {
    p <- p_list[[1]] + ggtitle("Effects of transformations")
    w <- 6
    h <- 6
  }
  # save plot
  if(inherits(outdir, "character")) {
    p_png <- file.path(outdir, "deseq_qc_counts.png")
    if(check_path(outdir)) {
      if(!file.exists(p_png) | overwrite) {
        ggsave(p_png, p, width = w, height = w, dpi = 200, units = "in")
      }
    }
  }
  p # return data
}


#' @describeIn deseq_qc_top_gene
#' Check values for top genes
#'
#' @param x `DEseqDataSet`, see `deseq_qc_dds(x=)`
#' @param outdir character saving the results
#' @param transform bool compute tranformed data, `vst()`, `rlog()`
#' default: FALSE
#' @param n_max integer number of genes to display, defualt: 30
#' @param overwrite bool overwrite exists file, default: FALSE
#'
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplotify as.ggplot
#' @importFrom patchwork wrap_plots plot_annotation
#'
#' @return ggplot
#'
#' @export
deseq_qc_top_gene <- function(x, outdir = NULL, transform = TRUE,
                              n_max = 30, overwrite = FALSE) {
  qc_data <- deseq_qc_dds(x, outdir, transform)
  dds <- x
  # dds     <- deseq_qc_dds(x, outdir, transform, return_dds = TRUE)
  if(is.null(qc_data)) {
    return(NULL)
  }
  # top genes
  top_gene <- order(rowMeans(counts(dds, normalized = TRUE)),
                    decreasing = TRUE)[1:n_max]
  dsn <- as.data.frame(colData(dds)[,c("condition")])
  # colors and breaks
  val <- unlist(lapply(qc_data, function(i) {
    assay(i)[top_gene, ]
  }))
  breaks <- seq(min(val), max(val), length.out = 255)
  colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdYlBu")))(255)
  # plot
  p_list <- lapply(names(qc_data), function(i) {
    assay(qc_data[[i]])[top_gene, ] %>%
      pheatmap::pheatmap(breaks = breaks, color = colors,
                         cluster_rows = FALSE, show_rownames = FALSE,
                         cluster_cols = FALSE, border_color = NA,
                         silent = TRUE) %>%
      ggplotify::as.ggplot() +
      ggtitle(i)
  })
  # combine plots
  if(length(p_list) > 1) {
    p <- patchwork::wrap_plots(p_list, nrow = 1, guides = "collect") +
      patchwork::plot_annotation(title = "Top ranked genes")
    w <- 6 # fig
    h <- 5 # fig
  } else {
    p <- p_list[[1]] + ggtitle("Top ranked genes")
    w <- 3
    h <- 5
  }
  # save plot
  if(inherits(outdir, "character")) {
    p_png <- file.path(outdir, "deseq_qc_top_gene.png")
    if(check_path(outdir)) {
      if(!file.exists(p_png) | overwrite) {
        ggsave(p_png, p, width = w, height = w, dpi = 200, units = "in")
      }
    }
  }
  p # return data
}


#' @describeIn deseq_qc_dist
#' check distance between samples
#'
#' @param x `DEseqDataSet`, see `deseq_qc_dds(x=)`
#' @param outdir character saving the results
#' @param transform bool compute tranformed data, `vst()`, `rlog()`
#' default: FALSE
#' @param overwrite bool overwrite exists file, default: FALSE
#'
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplotify as.ggplot
#' @importFrom patchwork wrap_plots plot_annotation
#'
#' @return ggplot
#'
#' @export
deseq_qc_dist <- function(x, outdir = NULL, transform = TRUE,
                          overwrite = FALSE) {
  qc_data <- deseq_qc_dds(x, outdir, transform)
  if(is.null(qc_data)) {
    return(NULL)
  }
  colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255)
  p_list <- lapply(names(qc_data), function(i) {
    sample_dist <- dist(t(assay(qc_data[[i]])))
    as.matrix(sample_dist) %>%
      pheatmap::pheatmap(
        clustering_distance_cols = sample_dist,
        clustering_distance_rows = sample_dist,
        color = colors,
        silent = TRUE) %>%
      ggplotify::as.ggplot() +
      ggtitle(i)
  })
  # combine plots
  if(length(p_list) > 1) {
    p <- patchwork::wrap_plots(p_list, nrow = 2, guides = "collect") +
      patchwork::plot_annotation(title = "Distance between samples")
    w <- 8 # fig
    h <- 8 # fig
  } else {
    p <- p_list[[1]] + ggtitle("Distance between samples")
    w <- 5
    h <- 5
  }
  # save plot
  if(inherits(outdir, "character")) {
    p_png <- file.path(outdir, "deseq_qc_dist.png")
    if(check_path(outdir)) {
      if(!file.exists(p_png) | overwrite) {
        ggsave(p_png, p, width = w, height = w, dpi = 200, units = "in")
      }
    }
  }
  p # return data
}


#' @describeIn deseq_qc_pca
#' check PCA plot for dds
#'
#' @param x `DEseqDataSet`, see `deseq_qc_dds(x=)`
#' @param outdir character saving the results
#' @param transform bool compute tranformed data, `vst()`, `rlog()`
#' default: FALSE
#' @param overwrite bool overwrite exists file, default: FALSE
#'
#' @improt ggplot2
#' @importFrom BiocGenerics plotPCA
#' @importFrom ggrepel geom_text_repel
#' @importFrom patchwork wrap_plots plot_annotation
#'
#' @return ggplot
#'
#' @export
deseq_qc_pca <- function(x, outdir = NULL, transform = FALSE,
                         overwrite = FALSE) {
  qc_data <- deseq_qc_dds(x, outdir, transform)
  if(is.null(qc_data)) {
    return(NULL)
  }
  p_list <- lapply(names(qc_data), function(i) {
    pca_data <- plotPCA(qc_data[[i]], returnData = TRUE)
    var <- round(attr(pca_data, "percentVar")*100, 1)
    var_text <- paste0(c("PC1", "PC2"), ": ", var, "% variance")
    # for batch
    if(! rlang::has_name(pca_data, "batch")) {
      pca_data$batch <- pca_data$condition
    }
    ggplot(pca_data, aes(PC1, PC2, color = condition, shape = batch)) +
      geom_point() +
      ggrepel::geom_text_repel(aes(label = name)) +
      xlab(var_text[1]) +
      ylab(var_text[2]) +
      ggtitle(i) +
      theme_bw()
  })
  # combine plots
  if(length(p_list) > 1) {
    p <- patchwork::wrap_plots(p_list, nrow = 2, guides = "collect") +
      patchwork::plot_annotation(title = "PCA plot")
    w <- 7 # fig
    h <- 6 # fig
  } else {
    p <- p_list[[1]] + ggtitle("PCA plot")
    w <- 4
    h <- 4
  }
  # save plot
  if(inherits(outdir, "character")) {
    p_png <- file.path(outdir, "deseq_qc_pca.png")
    if(check_path(outdir)) {
      if(!file.exists(p_png) | overwrite) {
        ggsave(p_png, p, width = w, height = w, dpi = 200, units = "in")
      }
    }
  }
  p # return data
}


#' @describeIn deseq_qc_ma
#' check MA plot for res
#'
#' @param x see `deseq_qc_res(x=)`
#' @param outdir character saving the results
#' @param ylim numeric, limits for log2FoldChange, default: NULL,
#' calculate by `scales::breaks_ext(n=5)(df$log2FoldChange)`
#' @param fc numeric cutoff for foldchange, default: 1, ignore foldchange
#' @param pvalue numeric cutoff for padj, default: 0.1, the main criteria
#' @param p_adjust bool use p-adjust value instead
#' @param overwrite bool overwrite exists file, default: FALSE
#'
#' @improt ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom patchwork wrap_plots plot_annotation
#'
#' @return ggplot
#'
#' @export
deseq_qc_ma <- function(x, outdir = NULL, ylim = NULL,
                        fc = 1, pvalue = 0.1, p_adjust = TRUE,
                        label_list = NULL, label_max = 8, overwrite = FALSE) {
  #-- Check: arguments
  df <- deseq_qc_res(x, outdir, ylim, fc, pvalue, p_adjust)
  if(!inherits(df, "data.frame")) {
    warning("failed reading data from `res`")
    return(NULL)
  }
  #-- Check: limits for y
  fclim <- df$log2FoldChange
  fclim <- fclim[!is.na(fclim)]
  if(inherits(ylim, "numeric")) {
    ylim <- c(max(min(ylim), min(fclim)),
              min(max(ylim), max(fclim)))
  } else {
    ylim <- c(min(fclim), max(fclim))
  }
  # foldchange, at least [-2, 2]
  if(min(ylim) > -2) ylim[1] = -2
  if(max(ylim) < 2)  ylim[2] = 2
  breaks <- scales::breaks_extended(n = 5)(ylim)
  #-- run: title
  title <- glue::glue("criteria: foldChange >= {fc}, pvalue < {pvalue}")
  #-- plot
  p <- df %>%
    ggplot(aes(log10basemean, log2FoldChange, color = sig, shape = ext)) +
    geom_hline(yintercept = c(-1, 1), size = .5, color = "grey30", linetype = 2) +
    geom_hline(yintercept = 0, size = .5, color = "grey30") +
    geom_point(size = .7) +
    scale_y_continuous(name   = "log2 fold change",
                       breaks = breaks,
                       limits = c(min(breaks), max(breaks))) +
    scale_color_manual(values = c("up"   = "red",
                                  "not"  = "grey60",
                                  "down" = "blue")) +
    scale_shape_manual(values = c("up"   = 2,
                                  "dot"  = 20,
                                  "down" = 6)) +
    xlab("log10 mean of normalized counts") +
    ggtitle(title) +
    theme_bw() +
    theme(legend.position = "none")
  p <- deseq_qc_add_sig_label(p, lable_list, label_max)
  # save plot
  if(inherits(outdir, "character")) {
    p_png <- file.path(outdir, "deseq_qc_ma.png")
    if(check_path(outdir)) {
      if(!file.exists(p_png) | overwrite) {
        ggsave(p_png, p, width = 4, height = 3, dpi = 200, units = "in")
      }
    }
  }
  p # return data
}



#' @describeIn deseq_qc_volcano
#' check volcano plot for res
#'
#' @param x see `deseq_qc_res(x=)`
#' @param outdir character saving the results
#' @param xlim numeric, limits for log2FoldChange, default: NULL,
#' calculate by `scales::breaks_ext(n=5)(df$log2FoldChange)`
#' @param fc numeric cutoff for foldchange, default: 1, ignore foldchange
#' @param pvalue numeric cutoff for padj, default: 0.1, the main criteria
#' @param p_adjust bool use p-adjust value instead
#' @param overwrite bool overwrite exists file, default: FALSE
#'
#' @improt ggplot2
#' @importFrom BiocGenerics plotPCA
#' @importFrom ggrepel geom_text_repel
#' @importFrom patchwork wrap_plots plot_annotation
#'
#' @return ggplot
#'
#' @export
deseq_qc_volcano <- function(x, outdir = NULL, xlim = NULL,
                             fc = 1, pvalue = 0.1, p_adjust = TRUE,
                             label_list = NULL, label_max = 8,
                             overwrite = FALSE) {
  #-- Check: input
  df <- deseq_qc_res(x, outdir, xlim, fc, pvalue, p_adjust)
  if(!inherits(df, "data.frame")) {
    warning("failed reading data from `res`")
    return(NULL)
  }
  #-- Check: limits for x
  fclim <- df$log2FoldChange
  fclim <- fclim[!is.na(fclim)]
  if(inherits(xlim, "numeric")) {
    xlim <- c(max(min(xlim), min(fclim)),
              min(max(xlim), max(fclim)))
  } else {
    xlim <- c(min(fclim), max(fclim))
  }
  # foldchange, at least [-2, 2]
  if(min(xlim) > -2) xlim[1] = -2
  if(max(xlim) < 2)  xlim[2] = 2
  breaks <- scales::breaks_extended(n = 5)(xlim)
  #-- run: title
  title <- glue::glue("criteria: foldChange >= {fc}, pvalue < {pvalue}")
  #-- plot:
  p <- df %>%
    ggplot(aes(log2fc, log10pval, color = sig, shape = ext)) +
    geom_vline(xintercept = c(-1, 1), size = .5, color = "grey30", linetype = 2) +
    geom_vline(xintercept = 0, size = .5, color = "grey30") +
    geom_point(size = .7) +
    scale_x_continuous(name   = "log2 fold change",
                       breaks = breaks,
                       limits = c(min(breaks), max(breaks))) +
    scale_color_manual(values = c("up"   = "red",
                                  "not"  = "grey60",
                                  "down" = "blue")) +
    scale_shape_manual(values = c("up"   = 2,
                                  "dot"  = 20,
                                  "down" = 6)) +
    ylab("-log10 pvalue adjust") +
    ggtitle(title) +
    theme_bw() +
    theme(legend.position = "none")
  p <- deseq_qc_add_sig_label(p, lable_list, label_max)
  # save plot
  if(inherits(outdir, "character")) {
    p_png <- file.path(outdir, "deseq_qc_volcano.png")
    if(check_path(outdir)) {
      if(!file.exists(p_png) | overwrite) {
        ggsave(p_png, p, width = 4, height = 4, dpi = 200, units = "in")
      }
    }
  }
  p # return data
}


#' @describeIn deseq_qc_scatter
#' for scatter plot
#' x: log10(wt + 1)
#' y: log10(mut + 1)
#'
#' @param x see `deseq_qc_res(x=)`
#' @param outdir character saving the results
#' @param fc numeric cutoff for foldchange, default: 1, ignore foldchange
#' @param pvalue numeric cutoff for padj, default: 0.1, the main criteria
#' @param p_adjust bool use p-adjust value instead
#' @param density_points bool use `stat_density_2d()`, for
#'  large number dots
#' @param overwrite bool overwrite exists file, default: FALSE
#'
#' @improt ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom patchwork wrap_plots plot_annotation
#'
#' @return ggplot
#'
#' @export
deseq_qc_scatter <- function(x, outdir = NULL,
                             fc = 1, pvalue = 0.1, p_adjust = TRUE,
                             density_points = FALSE,
                             label_list = NULL, label_max = 8,
                             overwrite = FALSE) {
  #-- Check: input
  df <- deseq_qc_res(x, outdir, NULL, fc, pvalue, p_adjust)
  if(!inherits(df, "data.frame")) {
    warning("failed reading data from `res`")
    return(NULL)
  }
  #-- run: extract wt, mut
  wt  <- names(df)[2]
  mut <- names(df)[3]
  df <- df %>%
    dplyr::mutate(wt  = log10(!!as.name(wt) + 1),
                  mut = log10(!!as.name(mut) + 1))
  #-- run: plot
  breaks <- scales::breaks_extended(n = 5)(c(df$wt, df$mut))
  if(isTRUE(density_points)) {
    p1 <- df %>%
      ggplot(aes(wt, mut, color = sig)) +
      stat_density_2d(
        aes(fill = ..density..),
        data = dplyr::filter(df, sig == "not"),
        geom = "raster", contour = FALSE)
  } else {
    p1 <- df %>%
      ggplot(aes(wt, mut, color = sig)) +
      geom_point(size = .5)
  }
  #-- run: title
  title <- glue::glue("criteria: foldChange >= {fc}, pvalue < {pvalue}")
  p <- p1 +
    scale_color_manual(values = c("up"   = "red",
                                  "not"  = "grey60",
                                  "down" = "blue")) +
    scale_fill_gradient(low = "white", high = "black") +
    geom_abline(intercept = 0, slope = 1, linetype = 1, color = "grey30") +
    geom_abline(intercept = c(log10(2), -log10(2)), slope = 1, linetype = 2,
                color = "grey50") +
    geom_point(data = dplyr::filter(df, sig %in% c("up", "down"))) +
    scale_x_continuous(breaks = breaks, limits = c(min(breaks), max(breaks)),
                       name = glue::glue("log10 count of {wt}")) +
    scale_y_continuous(breaks = breaks, limits = c(min(breaks), max(breaks)),
                       name = glue::glue("log10 count of {mut}")) +
    ggtitle(title) +
    theme_bw() +
    theme(panel.grid = element_blank())
  # add sig labels
  p <- deseq_qc_add_sig_label(p, lable_list, label_max)
  # save plot
  if(inherits(outdir, "character")) {
    p_png <- file.path(outdir, "deseq_qc_scatter.png")
    if(check_path(outdir)) {
      if(!file.exists(p_png) | overwrite) {
        ggsave(p_png, p, width = 5, height = 4, dpi = 200, units = "in")
      }
    }
  }
  p # return data
}



#' @describeIn deseq_res_for_plot2
#' prepare res for plot
#' sig: for colors, ["up", "not", "down"]
#'   up: padj < pval & log2fc >= log2(fc)
#'   down: padj < pval & log2fc <= -log2(fc)
#'   not: is.na(padj) | padj >= pval
#'
#' ext: for shapes, ["up":2, "dot":20, "down":6]
#'   change dot shapes, based on the range of log2fc_limits
#'   default limits = `scales::breaks_ext(n=5)(log2FoldChange)`
#'   up: log2fc > max(log2fc_limits)
#'   down: log2fc < min(log2fc_limits)
#'   not: log2fc >= min(log2fc_limits) & log2fc <= max(log2fc_limits)
#'
#' @param res DESeqResults, data.frame or csv file
#' @param log2fc_limits numeric, limits for log2FoldChange, default: NULL,
#' calculate by `scales::breaks_ext(n=5)(df$log2FoldChange)`
#' @param fc numeric cutoff for foldchange, default: 1, ignore foldchange
#' @param pvalue numeric cutoff for padj, default: 0.1, the main criteria
#' @param p_adjust bool use p-adjust value instead
#' @importFrom scales breaks_ext
#'
#' @return data.frame
#'
#' @export
deseq_res_for_plot2 <- function(res, log2fc_limits = NULL,
                                fc = 1, pvalue = 0.1, p_adjust = p_adjust) {
  #-- Check: arguments
  if(inherits(res, "DESeqResults")) {
    df <- as.data.frame(res)
  } else if(inherits(res, "data.frame")) {
    df <- res
  } else if(inherits(res, "character")) {
    if(file.exists(res) & endsWith(res, ".csv")) {
      df <- read.csv(res)
    } else {
      df <- data.frame()
    }
  } else {
    warning(glue::glue("unknown res, expect `DESeqResults` or `data.frame`, ",
                       "get: {class(res)}"))
    return(NULL)
  }
  #-- Check: required columns
  rc <- c("baseMean", "log2FoldChange", "padj")
  if(!all(rc %in% names(df))) {
    rc_str <- paste(rc, collapse = ", ")
    warning(glue::glue("missing columns: [{rc_str}]"))
    return(NULL)
  }
  #-- Check: limits for x, y
  if(inherits(log2fc_limits, "numeric")) {
    breaks <- log2fc_limits
  } else {
    breaks <- scales::breaks_extended(n = 5)(df$log2FoldChange)
  }
  #-- run: sig for colors, up, not, down
  # in case, pvalue exists in data.frame
  pval <- pvalue
  df <- df %>%
    dplyr::mutate(
      log10pval = -log10(padj),
      log10basemean = log10(baseMean + 1),
      sig = ifelse(
        is.na(padj), "not", ifelse(
          padj < pval, ifelse(
            log2FoldChange >= log2(fc), "up", ifelse(
              log2FoldChange <= -log2(fc), "down", "not")), "not")))
  #-- run: ext for shapes, up, dot, down
  # fix lgo2fc, by shapes
  df %>%
    dplyr::mutate(
      ext = ifelse(
        is.na(log2FoldChange), "dot", ifelse(
          log2FoldChange > max(breaks), "up", ifelse(
            log2FoldChange < min(breaks), "down", "dot")))) %>%
    dplyr::mutate(log2fc = ifelse(
      ext == "up", max(breaks), ifelse(
        ext == "down", min(breaks), log2FoldChange)))
}




#' @describeIn deseq_qc_add_sig_label
#' add sig labels to plot, based on `sig` column and
#' is designed for:
#' `deseq_qc_ma()`, `deseq_qc_volcano()`, `deseq_qc_scatter()`
#'
#' @param x ggplot
#' @param label_list character gene names
#' priority: symbol > gene_id
#'
#' @return ggplot
#'
#' @export
deseq_qc_add_sig_label <- function(x, lable_list = NULL, label_max = 8) {
  if(!inherits(x, "ggplot")) {
    message(glue::glue("x is {class(x)}, expect `ggplot`"))
    return(x)
  }
  # retrieve the `data` from `ggplot`
  if("data" %in% names(x)) {
    df <- x$data
    rc <- c("gene_id", "sig", "padj")
    if(all(rc %in% names(df))) {
      #-- run: choose label_list genes
      t1 <- df$gene_id %in% label_list
      if("SYMBOL" %in% names(df)) {
        t2 <- df$SYMBOL %in% label_list
        t1 <- t1 | t2
      }
      df1 <- df[t1, ]
      #-- run: show top sig genes
      if(nrow(df1) < label_max) {
        df2 <- df %>%
          dplyr::arrange(padj) %>%
          head(label_max - nrow(df1))
        df1 <- rbind(df1, df2)
      }
      #-- run: show label
      if(nrow(df1) > 0) {
        label_col <- ifelse("SYMBOL" %in% names(df), "SYMBOL", "gene_id")
        x +
          ggrepel::geom_text_repel(
            mapping = aes(label = !!as.name(label_col)),
            data    = df1,
            # color              = "grey10",
            size               = 4,
            force              = .2,
            direction          = "both",
            point.padding      = .2,
            max.overlaps = Inf,
            box.padding  = 0.5,
            # max.overlaps       = 30,
            # box.padding        = .1,
            segment.color      = "grey20",
            # segment.size       = .4,
            # min.segment.length = 0
          ) +
          geom_point(
            data = df1, size = 2, shape = 20
          )
      }
    } else {
      rc_str <- paste(rc, collapse = ", ")
      message(glue::glue("missing columns: {rc_str}"))
    }
  }
}
























