#' Functions for deseq quality control, utils
#'
#' deseq_qc_dds(), parsing dds data from deseq dir
#' deseq_qc_res(), parsing res data from deseq dir
#' deseq_qc_counts()
#' deseq_qc_mean_sd()
#' deseq_qc_top_gene()
#' deseq_qc_dist()
#' deseq_qc_pca()
#' deseq_qc_ma()
#' deseq_qc_volcano()
#' deseq_qc_scatter()
#' deseq_qc_add_sig_label() add labels for ma, volcano, scatter
#'
#' to-to: simpify arguments by `...`
#'
#' @name deseq_qc_plot

#------------------------------------------------------------------------------#
# main: qc plots

#' @describeIn deseq_qc_counts
#'
#' x could be "dds", deseq_deseq2, `norm_table.fix.csv`
#' output: data.frame (merge)
#'
#' @param x DESeqDataSet
#' @param outdir character saving the results
#' @param transform_method bool compute tranformed data, `vst()`, `rlog()`
#' default: standard
#' @param n_max integer number of genes to display, defualt: 16
#' @param overwrite bool overwrite exists file, default: FALSE
#'
#' @return ggplot
#'
#' @export
deseq_qc_counts <- function(x, ...) {
  message("run 'deseq_qc_counts()' ...")
  #----------------------------------------------------------------------------#
  #-- Check: default values
  dots <- rlang::list2(...)
  args <- rlang::list2(
    outdir = NULL,
    n_max = 16,
    transform_method = "standard", # vst, rlog
    overwrite  = FALSE,
    readable   = TRUE
  )
  #-- update dots, for child functions
  args <- purrr::list_modify(args, !!!dots)
  for(name in names(args)) {
    if(rlang::is_empty(name)) next
    assign(name, args[[name]])
  }
  #----------------------------------------------------------------------------#
  #-- n_max
  if(inherits(n_max, "numeric")) {
    n_max <- round(n_max)
    if(n_max < 1 | n_max > 30) {
      n_max <- 16
    }
  } else {
    n_max <- 16
  }
  #-- transform
  if(inherits(transform_method, "character")) {
    if(!transform_method %in% c("standard", "vst", "rlog")) {
      warning(glue::glue(
        "transform_method is {transform_method}, expect ",
        "'standard', 'vst', 'rlog'; set to 'standard' "
      ))
      transform_method <- "standard"
    }
  } else {
    warning(glue::glue(
      "transform_method is {class(transform_method)}, expect 'character'",
      "'standard', 'vst', 'rlog'; set to 'standard' "
    ))
    transform_method <- "standard"
  }
  #-- x
  df <- deseq_mean(x, outdir = outdir, ...)
  if(inherits(df, "data.frame")) {
    message("x is 'data.frame', deseq_qc_counts()")
  } else if(is_hiseq_dir(x, "deseq_deseq2")) {
    df <- deseq_qc_res(x, ...) #
  } else {
    warning(glue::glue(
      "'x' is {class(x)}, expect 'DESeqDataSet' or 'deseq_deseq2'"
    ))
    return(NULL)
  }
  #-- check required columns: gene_id, baseMean
  df2 <- NULL
  if(inherits(df, "data.frame")) {
    wt  <- names(df)[2]
    mut <- names(df)[3]
    wt_names  <- names(df)[startsWith(names(df), wt)][-1] # remove wt
    mut_names <- names(df)[startsWith(names(df), mut)][-1] # remove mut
    df_names  <- data.frame(
      sample    = c(wt_names, mut_names),
      condition = c(rep(wt, length(wt_names)), rep(mut, length(mut_names)))
    )
    if(any(grepl("baseMean", names(df)))) {
      bidx <- grep("baseMean", names(df))[1]
      df1  <- df[, c(1, 4:(bidx-1))] %>%
        tidyr::pivot_longer(-"gene_id", names_to = "group", values_to = "count")
      df2 <- dplyr::left_join(df1, df_names, by = c("group" = "sample"))
    }
  }
  #-- plotting
  if(inherits(df2, "data.frame")) {
    # top n counts: row mean / pvalue
    .col_pvalue <- names(df)[names(df) %in% c("pvalue", "padj")]
    if(length(.col_pvalue) > 0) {
      top_gene <- dplyr::arrange(df, !!.col_pvalue[1]) %>%
        dplyr::pull(gene_id) %>%
        head(n_max)
    } else {
      top_gene <- dplyr::group_by(df2, gene_id) %>%
        dplyr::summarise(mean_count = mean(count)) %>%
        dplyr::arrange(desc(mean_count)) %>%
        dplyr::pull(gene_id) %>%
        head(n_max)
    }
    # make plot
    title <- glue::glue("counts, transform: {transform_method}")
    dplyr::filter(df2, gene_id %in% top_gene) %>%
      ggplot(aes(condition, count, color = condition)) +
      geom_point(position = position_jitter(width = 0.1, height = 0), size = 1) +
      facet_wrap(~gene_id, ncol = 4) +
      ggtitle(title) +
      theme_bw() +
      theme(legend.position = "none")
  }
}


#' @describeIn deseq_qc_mean_sd
#' check mean standard diviation (SD), after transformation by vst, vlog
#'
#' @param x DESeqDataSet
#' @param outdir character saving the results
#' @param transform_method bool compute tranformed data, `vst()`, `rlog()`
#' default: standard
#' @param overwrite bool overwrite exists file, default: FALSE
#'
#' @importFrom vsn meanSdPlot
#' @importFrom patchwork wrap_plots plot_annotation
#'
#' @return ggplot
#'
#' @export
deseq_qc_mean_sd <- function(x, ...) {
  message("run 'deseq_qc_mean_sd()' ...")
  #----------------------------------------------------------------------------#
  #-- Check: default values
  dots <- rlang::list2(...)
  args <- rlang::list2(
    outdir = NULL,
    transform_method = "standard", # vst, rlog
    overwrite  = FALSE,
    readable   = TRUE
  )
  #-- update dots, for child functions
  args <- purrr::list_modify(args, !!!dots)
  for(name in names(args)) {
    if(rlang::is_empty(name)) next
    assign(name, args[[name]])
  }
  #----------------------------------------------------------------------------#
  #-- transform
  if(inherits(transform_method, "character")) {
    if(!transform_method %in% c("standard", "vst", "rlog")) {
      warning(glue::glue(
        "transform_method is {transform_method}, expect ",
        "'standard', 'vst', 'rlog'; set to 'standard' "
      ))
      transform_method <- "standard"
    }
  } else {
    warning(glue::glue(
      "transform_method is {class(transform_method)}, expect 'character'",
      "'standard', 'vst', 'rlog'; set to 'standard' "
    ))
    transform_method <- "standard"
  }
  #-- data
  dds_trans <- deseq_qc_dds(x, outdir = outdir, return_data = "dds_trans")
  #-- plot
  title <- glue::glue("tranform: {transform_method}")
  if(inherits(dds_trans, "list")) {
    trans <- dds_trans[[transform_method]]
    if(inherits(trans, "DESeqTransform")) {
      ms <- vsn::meanSdPlot(assay(trans), plot = FALSE) # rank, sd, px, py, gg
      if("gg" %in% names(ms)) {
        ms[["gg"]] +
          ggtitle(title) +
          theme_bw()
      }
    }
  }
}


#' @describeIn deseq_qc_top_gene
#' Check values for top genes
#'
#' @param x `DESeqDataSet`, see `deseq_qc_dds(x=)`
#' @param outdir character saving the results
#' @param transform_method bool compute tranformed data, `vst()`, `rlog()`
#' default: standard
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
deseq_qc_top_gene <- function(x, ...) {
  message("run 'deseq_qc_top_gene()' ...")
  #----------------------------------------------------------------------------#
  #-- Check: default values
  dots <- rlang::list2(...)
  args <- rlang::list2(
    outdir = NULL,
    n_max = 20,
    transform_method = "standard", # vst, rlog
    overwrite  = FALSE,
    readable   = TRUE
  )
  #-- update dots, for child functions
  args <- purrr::list_modify(args, !!!dots)
  for(name in names(args)) {
    if(rlang::is_empty(name)) next
    assign(name, args[[name]])
  }
  #----------------------------------------------------------------------------#
  #-- n_max
  if(inherits(n_max, "numeric")) {
    n_max <- round(n_max)
    if(n_max < 1 | n_max > 1000) {
      n_max <- 20
    }
  } else {
    n_max <- 20
  }
  #-- transform
  if(inherits(transform_method, "character")) {
    if(!transform_method %in% c("standard", "vst", "rlog")) {
      warning(glue::glue(
        "transform_method is {transform_method}, expect ",
        "'standard', 'vst', 'rlog'; set to 'standard' "
      ))
      transform_method <- "standard"
    }
  } else {
    warning(glue::glue(
      "transform_method is {class(transform_method)}, expect 'character'",
      "'standard', 'vst', 'rlog'; set to 'standard' "
    ))
    transform_method <- "standard"
  }
  #-- data
  dds_trans <- deseq_qc_dds(x, outdir = outdir, return_data = "dds_trans")
  colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdYlBu")))(255)
  #-- plot
  title <- glue::glue("tranform: {transform_method}")
  if(inherits(dds_trans, "list")) {
    trans <- dds_trans[[transform_method]]
    if(inherits(trans, "DESeqTransform")) {
      top_gene <- order(rowMeans(assay(trans)), decreasing = TRUE)[1:n_max]
      assay(trans)[top_gene, ] %>%
        pheatmap::pheatmap(
          color         = colors,
          cluster_rows  = FALSE,
          show_rownames = FALSE,
          cluster_cols  = FALSE,
          border_color  = NA,
          silent        = TRUE
        ) %>%
        ggplotify::as.ggplot() +
        ggtitle(title)
    }
  }
}


#' @describeIn deseq_qc_dist
#' check distance between samples
#'
#' @param x DESeqDataSet
#' @param outdir character saving the results
#' @param transform_method bool compute tranformed data, `vst()`, `rlog()`
#' default: standard
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
deseq_qc_dist <- function(x, ...) {
  message("run 'deseq_qc_dist()' ...")
  #----------------------------------------------------------------------------#
  #-- Check: default values
  dots <- rlang::list2(...)
  args <- rlang::list2(
    outdir = NULL,
    transform_method = "standard", # vst, rlog
    overwrite  = FALSE,
    readable   = TRUE
  )
  #-- update dots, for child functions
  args <- purrr::list_modify(args, !!!dots)
  for(name in names(args)) {
    if(rlang::is_empty(name)) next
    assign(name, args[[name]])
  }
  #----------------------------------------------------------------------------#
  #-- transform
  if(inherits(transform_method, "character")) {
    if(!transform_method %in% c("standard", "vst", "rlog")) {
      warning(glue::glue(
        "transform_method is {transform_method}, expect ",
        "'standard', 'vst', 'rlog'; set to 'standard' "
      ))
      transform_method <- "standard"
    }
  } else {
    warning(glue::glue(
      "transform_method is {class(transform_method)}, expect 'character'",
      "'standard', 'vst', 'rlog'; set to 'standard' "
    ))
    transform_method <- "standard"
  }
  #-- data
  dds_trans <- deseq_qc_dds(x, outdir = outdir, return_data = "dds_trans")
  colors    <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(255)
  title     <- glue::glue("tranform: {transform_method}")
  #-- plot
  if(inherits(dds_trans, "list")) {
    trans <- dds_trans[[transform_method]]
    if(inherits(trans, "DESeqTransform")) {
      sample_dist <- dist(t(assay(trans)))
      as.matrix(sample_dist) %>%
        pheatmap::pheatmap(
          color  = colors,
          silent = TRUE,
          clustering_distance_cols = sample_dist,
          clustering_distance_rows = sample_dist) %>%
        ggplotify::as.ggplot() +
        ggtitle(title)
    }
  }
}


#' @describeIn deseq_qc_pca
#' check PCA plot for dds
#'
#' @param x `DESeqDataSet`, see `deseq_qc_dds(x=)`
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
deseq_qc_pca <- function(x, ...) {
  message("run 'deseq_qc_pca()' ...")
  #----------------------------------------------------------------------------#
  #-- Check: default values
  dots <- rlang::list2(...)
  args <- rlang::list2(
    outdir = NULL,
    transform_method = "standard", # vst, rlog
    overwrite  = FALSE,
    readable   = TRUE
  )
  #-- update dots, for child functions
  args <- purrr::list_modify(args, !!!dots)
  for(name in names(args)) {
    if(rlang::is_empty(name)) next
    assign(name, args[[name]])
  }
  #----------------------------------------------------------------------------#
  #-- transform
  if(inherits(transform_method, "character")) {
    if(!transform_method %in% c("standard", "vst", "rlog")) {
      warning(glue::glue(
        "transform_method is {transform_method}, expect ",
        "'standard', 'vst', 'rlog'; set to 'standard' "
      ))
      transform_method <- "standard"
    }
  } else {
    warning(glue::glue(
      "transform_method is {class(transform_method)}, expect 'character'",
      "'standard', 'vst', 'rlog'; set to 'standard' "
    ))
    transform_method <- "standard"
  }
  #-- data
  dds_trans <- deseq_qc_dds(x, outdir = outdir, return_data = "dds_trans")
  colors    <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255)
  title     <- glue::glue("tranform: {transform_method}")
  #-- plot
  if(inherits(dds_trans, "list")) {
    trans <- dds_trans[[transform_method]]
    if(inherits(trans, "DESeqTransform")) {
      if("batch" %in% names(colData(trans))) {
        intgroup = c("condition", "batch")
      } else {
        intgroup = c("condition")
      }
      pca_data <- plotPCA(trans, intgroup = intgroup, returnData = TRUE)
      var <- round(attr(pca_data, "percentVar")*100, 1)
      var_text <- paste0(c("PC1", "PC2"), ": ", var, "% variance")
      # for batch
      if(!"batch" %in% names(pca_data)) {
        pca_data$batch <- pca_data$condition
      }
      ggplot(pca_data, aes(PC1, PC2, color = condition, shape = batch)) +
        geom_point(size = 2) +
        ggrepel::geom_text_repel(aes(label = name)) +
        xlab(var_text[1]) +
        ylab(var_text[2]) +
        ggtitle(title) +
        theme_bw()
    }
  }
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
deseq_qc_ma <- function(x, ...) {
  message("run 'deseq_qc_ma()' ...")
  #----------------------------------------------------------------------------#
  #-- Check: default values
  dots <- rlang::list2(...)
  args <- rlang::list2(
    outdir     = NULL,
    fc         = 2,
    pvalue     = 0.05,
    p_adjust   = TRUE,
    ylim       = c(-2, 2),
    label_list = NULL,
    label_max  = 8,
    overwrite  = FALSE,
    readable   = TRUE,
    # .col_label = "gene_id",
    add_sig    = FALSE,
    shrink_method = "standard" # apeglm, ashr, normal
  )
  #-- update dots, for child functions
  args <- purrr::list_modify(args, !!!dots)
  args$log2fc_limits <- args$ylim
  #-- update global
  for(name in names(args)) {
    if(rlang::is_empty(name)) next
    assign(name, args[[name]])
  }
  #----------------------------------------------------------------------------#
  #-- run: load data.frame
  df <- deseq_qc_res(x, !!!args)
  #-- Check columns
  rc <- c("log10basemean", "log2fc", "ext", "gene_id", "log2FoldChange")
  if(inherits(df, "data.frame")) {
    if(!all(rc %in% names(df))) {
      rc_str <- paste(rc, collapse = ", ")
      warning(glue::glue(
        "missing columns: [{rc_str}], check argument; 'x'={x}"
      ))
      return(NULL)
    }
  } else {
    warning(glue::glue(
      "could not find `norm_table.fix.csv` file, check 'x'"
    ))
    return(NULL)
  }
  #----------------------------------------------------------------------------#
  #-- check: log2 limits
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
  ylim <- ylim * 1.1 # extend ylim by 10%
  breaks <- scales::breaks_extended(n = 5)(ylim)
  #----------------------------------------------------------------------------#
  #-- add 'sig', force, re-run
  df <- get_sig_name(df, return_dataframe = TRUE, force = TRUE, !!!args)
  title <- glue::glue("criteria: foldChange >= {fc}, pvalue < {pvalue}")
  #----------------------------------------------------------------------------#
  #-- plot
  p <- df %>%
    ggplot(aes(log10basemean, log2fc, color = sig, shape = ext)) +
    geom_hline(yintercept = c(-1, 1), size = .5, color = "grey30", linetype = 2) +
    geom_hline(yintercept = 0, size = .5, color = "grey30") +
    geom_point(size = .4) +
    scale_y_continuous(name   = "log2 fold change",
                       breaks = breaks,
                       limits = ylim) +
    # limits = c(min(breaks), max(breaks))) +
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
  #-- add sig labels
  if(add_sig) {
    p <- deseq_qc_add_sig_label(p)
  }
  #-- return
  p
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
deseq_qc_volcano <- function(x, ...) {
  message("run 'deseq_qc_volcano()' ...")
  #----------------------------------------------------------------------------#
  #-- Check: default values
  dots <- rlang::list2(...)
  args <- rlang::list2(
    outdir = NULL,
    fc = 2,
    pvalue     = 0.05,
    p_adjust   = TRUE,
    xlim       = c(-2, 2),
    label_list = NULL,
    label_max  = 8,
    overwrite  = FALSE,
    readable   = TRUE,
    # .col_label = "gene_id",
    add_sig    = FALSE,
    shrink_method = "standard"
  )
  #-- update dots, for child functions
  args <- purrr::list_modify(args, !!!dots) # update log2fc_limit
  args$log2fc_limits <- args$xlim
  #-- update global
  for(name in names(args)) {
    if(rlang::is_empty(name)) next
    assign(name, args[[name]])
  }
  #----------------------------------------------------------------------------#
  #-- Check: arguments: update log2fc
  df <- deseq_qc_res(x, !!!args)
  #-- Check columns
  rc <- c("log10pval", "log2fc", "ext", "gene_id", "log2FoldChange")
  if(inherits(df, "data.frame")) {
    if(!all(rc %in% names(df))) {
      rc_str <- paste(rc, collapse = ", ")
      warning(glue::glue(
        "missing columns: [{rc_str}], check argument; 'x'={x}"
      ))
      return(NULL)
    }
  } else {
    warning(glue::glue(
      "could not find `norm_table.fix.csv` file, check 'x'"
    ))
    return(NULL)
  }
  #-- y-limits;
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
  xlim <- xlim * 1.1 # extend ylim by 10%
  breaks <- scales::breaks_extended(n = 5)(xlim)
  #----------------------------------------------------------------------------#
  #-- update 'sig', force, re-run
  df <- get_sig_name(df, return_dataframe = TRUE, force = TRUE, !!!args)
  title <- glue::glue("criteria: foldChange >= {fc}, pvalue < {pvalue}")
  #----------------------------------------------------------------------------#
  #-- plot:
  p <- df %>%
    ggplot(aes(log2fc, log10pval, color = sig, shape = ext)) +
    geom_vline(xintercept = c(-1, 1), size = .5, color = "grey30", linetype = 2) +
    geom_vline(xintercept = 0, size = .5, color = "grey30") +
    geom_point(size = .5, alpha = 0.5) +
    scale_x_continuous(name   = "log2 fold change",
                       breaks = breaks,
                       limits = xlim) +
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
  #-- add sig labels
  if(add_sig) {
    p <- deseq_qc_add_sig_label(p)
  }
  #-- return
  p
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
deseq_qc_scatter <- function(x, ...) {
  message("run 'deseq_qc_scatter()' ...")
  #----------------------------------------------------------------------------#
  #-- Check: default values
  dots <- rlang::list2(...)
  args <- rlang::list2(
    outdir = NULL,
    fc = 2,
    pvalue     = 0.05,
    p_adjust   = TRUE,
    label_list = NULL,
    label_max  = 8,
    overwrite  = FALSE,
    readable   = TRUE,
    density_point = FALSE,
    # .col_label = "gene_id",
    color_by   = "sig", # sig, tissue
    add_sig    = FALSE,
    shrink_method = "standard"
  )
  #-- update dots, for child functions
  args <- purrr::list_modify(args, !!!dots)
  #-- update global
  for(name in names(args)) {
    if(rlang::is_empty(name)) next
    assign(name, args[[name]])
  }
  #----------------------------------------------------------------------------#
  #-- Check: arguments: update log2fc
  df <- deseq_qc_res(x, !!!args)
  #-- Check columns
  rc <- c("log2FoldChange", "log2fc", "ext", "gene_id", "log2FoldChange")
  if(inherits(df, "data.frame")) {
    if(!all(rc %in% names(df))) {
      rc_str <- paste(rc, collapse = ", ")
      warning(glue::glue(
        "missing columns: [{rc_str}], check argument; 'x'={x}"
      ))
      return(NULL)
    }
  } else {
    warning(glue::glue(
      "could not find `norm_table.fix.csv` file, check 'x'"
    ))
    return(NULL)
  }
  #----------------------------------------------------------------------------#
  #-- breaks
  wt  <- names(df)[2]
  mut <- names(df)[3]
  df <- df %>%
    dplyr::mutate(wt  = log10(!!as.name(wt) + 1),
                  mut = log10(!!as.name(mut) + 1))
  #-- run: determine limits, breaks
  breaks  <- scales::breaks_extended(n = 5)(c(df$wt, df$mut))
  xlimits <- c(min(c(df$wt, df$mut)), max(c(df$wt, df$mut)))
  ylimits <- xlimits
  #-- y-limits;
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
  ylim <- ylim * 1.1 # extend ylim by 10%
  breaks <- scales::breaks_extended(n = 5)(ylim)
  #----------------------------------------------------------------------------#
  #-- update 'sig', force, re-run
  df <- get_sig_name(df, return_dataframe = TRUE, force = TRUE, !!!args)
  title <- glue::glue("criteria: foldChange >= {fc}, pvalue < {pvalue}")
  #-- run: plot
  if(isTRUE(density_point)) {
    p1 <- df %>%
      # # ggplot(aes(wt, mut, color = sig)) +
      # ggplot(aes_string("wt", "mut", color = color_by)) +
      ggplot(aes(wt, mut, color = !!color_by)) +
      stat_density_2d(
        aes(fill = ..density..),
        data = dplyr::filter(df, sig == "not"),
        geom = "raster", contour = FALSE)
  } else {
    p1 <- df %>%
      # ggplot(aes(wt, mut, color = sig)) +
      # ggplot(aes_string("wt", "mut", color = color_by)) +
      ggplot(aes(wt, mut, color = !!color_by)) +
      geom_point(size = .4, alpha = .5)
  }
  #----------------------------------------------------------------------------#
  #-- plot
  # add colors
  cc <- list(
    sig = c("up"   = "red", "not"  = "grey60", "down" = "blue"),
    tissue = c("germline" = "#ff2e16", "intermediate" = "#ffeb3d",
               "other" = "#85878b", "soma"  = "#2ca748")
  )
  p1 <- p1 +
    scale_color_manual(values = c("up"   = "red",
                                  "not"  = "grey60",
                                  "down" = "blue")) +
    # scale_color_manual(values = cc[[color_by]]) +
    scale_fill_gradient(low = "white", high = "black") +
    geom_abline(intercept = 0, slope = 1, linetype = 1, color = "grey30") +
    geom_abline(intercept = c(log10(2), -log10(2)), slope = 1, linetype = 2,
                color = "grey50") +
    geom_point(data = dplyr::filter(df, sig %in% c("up", "down")),
               size = .6) +
    scale_x_continuous(breaks = breaks, limits = xlimits,
                       name = glue::glue("log10 count of {wt}")) +
    scale_y_continuous(breaks = breaks, limits = ylimits,
                       name = glue::glue("log10 count of {mut}")) +
    ggtitle(title) +
    theme_bw() +
    theme(panel.grid = element_blank())
  #-- add sig labels
  if(add_sig) {
    p1 <- deseq_qc_add_sig_label(p1)
  }
  #-- return
  p1
}


#' @describeIn deseq_qc_scatter2
#' for scatter plot, colored by specific column
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
deseq_qc_scatter2 <- function(x, ...) {
  message("run 'deseq_qc_scatter2()' ...")
  #----------------------------------------------------------------------------#
  #-- Check: default values
  args <- rlang::list2(
    outdir = NULL,
    fc = 2,
    pvalue     = 0.05,
    p_adjust   = TRUE,
    label_list = NULL,
    label_max  = 8,
    overwrite  = FALSE,
    readable   = TRUE,
    density_point = FALSE,
    # .col_label = "gene_id",
    color_by   = "sig", # sig, tissue
    add_sig    = FALSE,
    shrink_method = "standard"
  )
  args <- purrr::list_modify(args, ...)
  #-- update global
  for(name in names(args)) {
    if(rlang::is_empty(name)) next
    assign(name, args[[name]])
  }
  #----------------------------------------------------------------------------#
  #-- Check: arguments: update log2fc
  df <- deseq_qc_res(x, !!!args) # loading data
  #-- Check columns
  rc <- c("gene_id", "log2FoldChange", "log2fc", "ext")
  if(inherits(df, "data.frame")) {
    if(!all(rc %in% names(df))) {
      rc_str <- paste(rc, collapse = ", ")
      warning(glue::glue(
        "missing columns: [{rc_str}], check argument; 'x'={x}"
      ))
      return(NULL)
    }
  } else {
    warning(glue::glue(
      "could not find `norm_table.fix.csv` file, check 'x'"
    ))
    return(NULL)
  }
  if(! color_by %in% names(df)) {
    color_by <- "sig"
    message(glue::glue(
      "column [{color_by}] not found, use [sig] instead"
    ))
  }
  #----------------------------------------------------------------------------#
  #-- convert log10 scale
  wt  <- names(df)[2]
  mut <- names(df)[3]
  df <- df %>%
    dplyr::mutate(wt  = log10(!!as.name(wt) + 1),
                  mut = log10(!!as.name(mut) + 1))
  #-- run: determine limits, breaks
  breaks  <- scales::breaks_extended(n = 5)(c(df$wt, df$mut))
  xlimits <- c(min(c(df$wt, df$mut)), max(c(df$wt, df$mut)))
  ylimits <- xlimits
  #-- y-limits;
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
  ylim <- ylim * 1.1 # extend ylim by 10%
  breaks <- scales::breaks_extended(n = 5)(ylim)
  #----------------------------------------------------------------------------#
  #-- update 'sig', force, re-run
  df <- get_sig_name(df, return_dataframe = TRUE, force = TRUE, !!!args)
  title <- glue::glue("criteria: foldChange >= {fc}, pvalue < {pvalue}")
  #-- run: plot
  p1 <- df %>%
    ggplot(aes(wt, mut, color = !!as.name(color_by)))
  if(isTRUE(density_point)) {
    p2 <- p1 +
      stat_density_2d(
        aes(fill = ..density..),
        data = dplyr::filter(df, sig == "not"),
        geom = "raster", contour = FALSE)
  } else {
    p2 <- p1 +
      geom_point(size = 1, alpha = .8)
  }
  #----------------------------------------------------------------------------#
  #-- plot
  # add colors
  cc <- list(
    sig = c("up" = "red", "not" = "grey60", "down" = "blue"),
    tissue = c("germline" = "#ff2e16", "intermediate" = "#ffeb3d",
               "other" = "#85878b", "soma"  = "#2ca748")
  )
  p2 <- p2 +
    scale_color_manual(values = cc[[color_by]]) +
    scale_fill_gradient(low = "white", high = "black") +
    geom_abline(intercept = 0, slope = 1, linetype = 1, color = "grey30") +
    geom_abline(intercept = c(log10(2), -log10(2)), slope = 1, linetype = 2,
                color = "grey50") +
    geom_point(data = dplyr::filter(df, sig %in% c("up", "down")),
               size = .6) +
    scale_x_continuous(breaks = breaks, limits = xlimits,
                       name = glue::glue("log10 count of {wt}")) +
    scale_y_continuous(breaks = breaks, limits = ylimits,
                       name = glue::glue("log10 count of {mut}")) +
    ggtitle(title) +
    theme_bw() +
    theme(panel.grid = element_blank())
  #-- add sig labels
  if(add_sig) {
    p2 <- deseq_qc_add_sig_label(p2)
  }
  #-- return
  p2
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
deseq_qc_add_sig_label <- function(x, ...) {
  #----------------------------------------------------------------------------#
  #-- Check: arguments
  dots <- rlang::list2(...)
  args <- rlang::list2(
    label_list = NULL,
    label_max  = 8,
    .col_sig   = "sig",
    .col_label = "auto" # gene_id, SYMBOL, Gene, ...
  )
  #-- update dots, for child functions
  args <- purrr::list_modify(args, !!!dots)
  #-- update args,
  for(name in names(args)) {
    if(rlang::is_empty(name)) next
    assign(name, args[[name]])
  }
  #----------------------------------------------------------------------------#
  #-- Check: ggplot input
  if(!inherits(x, "ggplot")) {
    message(glue::glue("x is {class(x)}, expect `ggplot`"))
    return(x)
  } else {
    # retrieve the `data` from `ggplot`
    df <- x[["data"]]
    #-- for sig column
    if(!all(.col_sig %in% names(df))) {
      message(glue::glue("missing 'sig' columns: {.col_sig}"))
      return(x)
    }
    #--------------------------------------------------------------------------#
    #-- coordinates for sig text
    #-- sig text
    sig_table <- table(df[[.col_sig]])
    sig_text  <- paste(
      paste0(names(sig_table), ": ", sig_table), collapse = "\n"
    )
    #-- sig coordinates c(0.1, 0.9)
    x_col <- rlang::as_name(x[["mapping"]][["x"]])
    y_col <- rlang::as_name(x[["mapping"]][["y"]])
    x_breaks <- scales::breaks_extended(n = 5)(df[[x_col]])
    y_breaks <- scales::breaks_extended(n = 5)(df[[y_col]])
    df_sig_text <- setNames(data.frame(matrix(ncol = ncol(df), nrow = 1)),
                                     nm = names(df)) %>%
      dplyr::mutate(
        !!x_col := min(x_breaks),
        !!y_col := max(y_breaks),
        label   = sig_text
      )
    # df_sig_text <- setNames(data.frame(matrix(ncol = ncol(df), nrow = 1)),
    #         nm = names(df)) %>%
    #   dplyr::mutate(
    #     !!x_col := max(df[[x_col]], na.rm = TRUE) * 0.15,
    #     !!y_col := max(df[[y_col]], na.rm = TRUE) * 0.95,
    #     label    = sig_text
    #   )
    #--------------------------------------------------------------------------#
    #-- for label
    rc2 <- c("SYMBOL", "gene_id", "gene", "Gene", "id", .col_label)
    rc2 <- rc2[rc2 %in% names(df)]
    if(length(rc2) > 0) {
      .col_label <- rc2[1]
    } else {
      message(glue::glue("missing label column, eg: 'SYMBOL', 'gene_id'"))
      return(x)
    }
    #-- subset data.frame; label_list
    r_idx <- df %>%
      dplyr::select(all_of(.col_label)) %>%
      dplyr::mutate(across(everything(), function(i) i %in% label_list)) %>%
      rowMeans()
    df_label <- df[r_idx > 0, ]
    #-- add more
    n_left <- label_max - nrow(df_label)
    if(n_left > 0) {
      df_ex <- filt_sig_gene(df, type = "sig") %>% head(n_left)
      df_label <- rbind(df_label, df_ex)
    }
    # fix .col_label for NA
    df_label <- df_label %>%
      dplyr::mutate(!!.col_label := ifelse(
        is.na(!!sym(.col_label)) | !!sym(.col_label) == "NA",
        gene_id, !!sym(.col_label)))
    #--------------------------------------------------------------------------#
    #-- add labels
    if(nrow(df_label) > 0) {
      x +
        ggrepel::geom_text_repel(
          mapping = aes(label = .data[[.col_label]]),
          data    = df_label,
          color          = "grey10",
          size           = 4,
          force          = .2,
          direction      = "both",
          point.padding  = .2,
          max.overlaps   = Inf,
          # max.overlaps   = 30,
          box.padding    = 0.5, # additional padding around each text label
          segment.color  = "grey20",
          # segment.size   = .4,
          min.segment.length = 0, # draw all line segments
          max.time = 1, max.iter = 1e5 # stop after 1 second, or after 100,000 iterations
        ) +
        geom_point(
          data = df_label, size = .6, shape = 20
        ) +
        ggrepel::geom_text_repel(
          mapping = aes(label = label),
          data    = df_sig_text,
          hjust   = 0.5,
          color   = "grey50",
          size    = 2.5,
          max.overlaps = Inf
        )
    } else {
      return(x)
    }
  }
}

