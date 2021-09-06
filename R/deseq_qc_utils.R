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
#' @name deseq_qc_utils

#------------------------------------------------------------------------------#
# loading data for deseq_qc()

#' @describeIn deseq_qc_dds
#'
#' parse dds_trans data, from file: `deseq_res.rds`
#' see `run_deseq(dds, outdir, ...)`
#'
#' @param x DESeqDataSet object or outdir (deseq_deseq2)
#'
#' @return list of standard*, vst, rlog objects
#'
#' @export
deseq_qc_dds <- function(x = NULL, ...) {
  #----------------------------------------------------------------------------#
  #-- Check: default values
  dots <- rlang::list2(...)
  args <- rlang::list2(
    return_data = "dds", # dds, dds_trans, res, res_lfc
    outdir      = NULL,
    overwrite   = FALSE
  )
  #-- update dots, for child functions
  dots_args <- lapply(names(args), function(i) {
    if(!i %in% names(dots)) {
      args[i]
    }
  })
  dots <- c(dots, unlist(dots_args, recursive = FALSE, use.names = TRUE))
  #-- update global
  for(name in names(dots)) {
    assign(name, dots[[name]])
  }
  #----------------------------------------------------------------------------#
  #-- return_data
  if(inherits(return_data, "logical")) {
    if(return_data) return_data <- "dds" # default
  } else if(inherits(return_data, "character")) {
    rc <- c("dds", "dds_trans", "res", "res_lfc")
    if(!return_data[1] %in% rc) {
      rc_str <- paste(rc, collapse = ", ")
      warning(glue::glue(
        "'return_data' is {return_data}, expect ['TRUE', {rc_str}], \n",
        "use 'dds' instead"
      ))
      return_data <- "dds"
    }
  }
  #-- Check: x, DESeqDataSeq
  res <- NULL
  if(is_hiseq_dir(x, "deseq_deseq2")) {
    deseq_res_rds <- list_hiseq_file(x, "deseq_res_rds", "deseq_deseq2")
    if(inherits(deseq_res_rds, "character")) {
      if(file.exists(deseq_res_rds)) {
        res <- readRDS(deseq_res_rds)
      }
    }
  } else if(inherits(x, "DESeqDataSet")) {
    # message("run `vst()` and `rlog()` might takes too long, \n",
    #         "set `outdir` to deseq_deseq2 directory, will speed up")
    res <- run_deseq_res(x, outdir = outdir, shrink = FALSE) #
  } else {
    warning(glue::glue(
      "x is '{class(x)}, expect 'DESeqDataSet', deseq_deseq2 dir"
    ))
  }
  #-- output
  if(inherits(res, "list")) {
    res[[return_data]]
  }
}


#' @describeIn deseq_qc_res
#'
#' @example
#' ma, volcano: log2FolcChange, pvalue
#' scatter: wt vs mut
#'
#'
#' get data.frame from `norm_table.fix.xls`,
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
deseq_qc_res <- function(x, ...) {
  #----------------------------------------------------------------------------#
  dots <- rlang::list2(...)
  args <- rlang::list2(
    shrink_method = "standard", # standard, normal, apeglm, ashr
    transform_method = "standard", # standard, vst, rlog
    log2fc_limits    = c(-2, 2),
    fc = 2,
    pvalue = 0.05,
    p_adjust = TRUE,
    overwrite = FALSE,
    .col_sig  = "sig",
    readable  = TRUE,
    genome    = NULL
  )
  #-- update dots, for child functions
  dots_args <- lapply(names(args), function(i) {
    if(!i %in% names(dots)) {
      args[i]
    }
  })
  dots <- c(dots, unlist(dots_args, recursive = FALSE, use.names = TRUE))
  #-- update global
  for(name in names(dots)) {
    assign(name, dots[[name]])
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
  #-- shrink
  if(inherits(shrink_method, "character")) {
    if(!shrink_method %in% c("standard", "normal", "apeglm", "ashr")) {
      warning(glue::glue(
        "shrink_method is {shrink_method}, expect ",
        "'standard', 'normal', 'apeglm', 'ashr'; set to 'standard' "
      ))
      shrink_method <- "standard"
    }
  } else {
    warning(glue::glue(
      "shrink_method is {class(shrink_method)}, expect ",
      "'standard', 'normal', 'apeglm', 'ashr'; set to 'standard' "
    ))
    shrink_method <- "standard"
  }
  #-- log2fc_limits, c(-2, 2)
  if(inherits(log2fc_limits, "numeric")) {
    if(length(log2fc_limits) == 2
       & all(log2fc_limits >= -10)
       & all(log2fc_limits <= 10)
    ) {
      class(log2fc_limits)
    } else {
      warning(glue::glue(
        "log2fc_limits is '{log2fc_limits}', not valid; set to 'c(-2, 2)'"
      ))
      log2fc_limits <- c(-2, 2)
    }
  } else {
    warning(glue::glue(
      "'log2fc_limits' is {class(log2fc_limits)}, expect 'numeric', ",
      "set to 'c(-2, 2)'"
    ))
    log2fc_limits <- c(-2, 2)
  }
  #-- fc, pvalue
  if(inherits(fc, "numeric") & inherits(pvalue, "numeric")) {
    if(fc <= 0 | pvalue <= 0 | pvalue > 1) {
      warning(glue::glue(
        "'fc' is {fc}, expect (0-Inf), greater than 0; ",
        "'pvalue' is {pvalue}`, expect (0-1)")
      )
      fc <- 2
      pvalue <- 0.05
    }
  } else {
    warning(glue::glue(
      "`fc` is {class(fc)}, expect `numeric`, ",
      "`pvalue` is {class(pvalue)}, expect `numeric`"))
    fc <- 2
    pvalue <- 0.05
  }
  if(!inherits(p_adjust, "logical")) {
    warning(glue::glue(
      "`p_value` is {class(pvalue)}, expect logical, ",
      "set to 'TRUE'"
    ))
    p_adjust <- TRUE
  }
  #-- overwrite
  if(!inherits(overwrite, "logical")) {
    message(glue::glue(
      "overwrite is {class(overwrite)}, expect logical; ",
      "set to 'FALSE'"
    ))
    overwrite <- FALSE
  }
  #----------------------------------------------------------------------------#
  #-- run: loading data from x
  df <- NULL
  if(inherits(x, "data.frame")) {
    message(glue::glue(
      "ignore 'shrink_method={shrink_method}'"
    ))
    df <- x
  } else if(is_hiseq_dir(x, hiseq_type = "deseq_deseq2")) {
    # re-create table; see `run_deseq_res()`
    # 1. dds trans; log2, normalized
    dt_list <- deseq_qc_dds(x, return_data = "dds_trans") #log2
    dt <- dt_list[[transform_method]] # dds_trans
    if(is.null(dt)) {
      warning("unknown x, `deseq_qc_res()` failed; `dds_trans` not found")
      return(NULL)
    }
    # 2. results table, fc, pvalue
    res <- deseq_qc_dds(x, return_data = "res_lfc") # see shrink_method
    dr <- res[[shrink_method]] # dds_res
    if(is.null(dr)) {
      warning("unknown x, `deseq_qc_res()` failed; `res_lfc` not found")
      return(NULL)
    }
    # 3. combine tables
    dt_count <- 2^assay(dt)
    dr <- as.data.frame(dr)
    df <- merge(dt_count, dr, by = "row.names" )
    colnames(df)[1] <- "gene_id"
    # 4. mean
    df <- deseq_mean(df) # add mean
  } else if(inherits(x, "character")) {
    if(file.exists(x) & endsWith(x, ".fix.csv")) {
      message(glue::glue(
        "ignore 'shrink_method={shrink_method}'"
      ))
      df <- read.csv(x) # for norm_table.fix.csv
    }
  } else {
    df <- NULL
  }
  #----------------------------------------------------------------------------#
  #-- Check: required columns
  if(inherits(df, "data.frame")) {
    rc <- c("baseMean", "log2FoldChange", "padj")
    if(!all(rc %in% names(df))) {
      rc_str <- paste(rc, collapse = ", ")
      warning(glue::glue("missing columns: [{rc_str}]"))
      return(NULL)
    }
  } else {
    warning(glue::glue("failed loading data from x={x}"))
    return(NULL)
  }
  #-- for scatter:
  # gene_id, wt, mut, for scatter
  wt  <- names(df)[2]
  mut <- names(df)[3]
  #-- breaks for x, y axis
  breaks <- scales::breaks_extended(n = 5)(log2fc_limits)
  #-- run: add sig labels
  if(.col_sig %in% names(df) & !overwrite) {
    message(glue::glue(".col_sig = '{.col_sig}' exists"))
  } else {
    df2 <- get_sig_name(df, return_dataframe = TRUE, !!!dots)
  }
  #----------------------------------------------------------------------------#
  #-- run: readable
  if(readable) {
    if(is_hiseq_dir(x, "deseq_deseq2")) {
      gene_table <- list_hiseq_file(x, "gene_readable_csv", "deseq_deseq2")
      Lgenome <- list_hiseq_file(x, "genome", "deseq_deseq2")
    } else {
      gene_table <- NULL
    }
    df <- set_readable(df, genome = genome, gene_table = gene_table)
  }
  #-- run: update log2fc, ext; for point-shapes
  # for shapes, up, dot, down
  # update log2fc, check outlier by log2fc_limits
  df %>%
    dplyr::mutate(
      log10pval     = -log10(padj),
      log10basemean = log10(baseMean + 1),
      ext = ifelse(
        is.na(log2FoldChange), "dot", ifelse(
          log2FoldChange > max(breaks), "up", ifelse(
            log2FoldChange < min(breaks), "down", "dot")))) %>%
    dplyr::mutate(log2fc = ifelse(
      ext == "up", max(breaks), ifelse(
        ext == "down", min(breaks), log2FoldChange)))
}



#------------------------------------------------------------------------------#
# deprecated functions

#' #' @describeIn deseq_qc_top_gene
#' #' Check values for top genes
#' #'
#' #' @param x `DEseqDataSet`, see `deseq_qc_dds(x=)`
#' #' @param outdir character saving the results
#' #' @param transform bool compute tranformed data, `vst()`, `rlog()`
#' #' default: FALSE
#' #' @param n_max integer number of genes to display, defualt: 30
#' #' @param overwrite bool overwrite exists file, default: FALSE
#' #'
#' #' @importFrom pheatmap pheatmap
#' #' @importFrom RColorBrewer brewer.pal
#' #' @importFrom ggplotify as.ggplot
#' #' @importFrom patchwork wrap_plots plot_annotation
#' #'
#' #' @return ggplot
#' #'
#' #' @export
#' deseq_qc_top_gene <- function(x, outdir = NULL, ...) {
#'   #-- Check: default values
#'   n_max <- 30
#'   transform <- TRUE
#'   overwrite <- FALSE
#'   #-- Check: input
#'   dots <- rlang::list2(...)
#'   for(name in names(dots)) {
#'     assign(name, dots[[name]])
#'   }
#'   #-- arguments
#'   qc_data <- deseq_qc_dds(x, outdir, transform)
#'   if(!inherits(x, "DESeqDataSet")) {
#'     dds <- deseq_qc_dds(x, outdir, transform, return_dds = TRUE)
#'   } else {
#'     dds <- x
#'   }
#'   if(inherits(qc_data, "NULL")) {
#'     return(NULL)
#'   }
#'   if(!inherits(dds, "DESeqDataSet")) {
#'     return(NULL)
#'   }
#'   # top genes
#'   top_gene <- order(rowMeans(counts(dds, normalized = TRUE)),
#'                     decreasing = TRUE)[1:n_max]
#'   dsn <- as.data.frame(colData(dds)[,c("condition")])
#'   # colors and breaks
#'   val <- unlist(lapply(qc_data, function(i) {
#'     assay(i)[top_gene, ]
#'   }))
#'   breaks <- seq(min(val), max(val), length.out = 255)
#'   colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdYlBu")))(255)
#'   # plot
#'   p_list <- lapply(names(qc_data), function(i) {
#'     assay(qc_data[[i]])[top_gene, ] %>%
#'       pheatmap::pheatmap(breaks = breaks, color = colors,
#'                          cluster_rows = FALSE, show_rownames = FALSE,
#'                          cluster_cols = FALSE, border_color = NA,
#'                          silent = TRUE) %>%
#'       ggplotify::as.ggplot() +
#'       ggtitle(i)
#'   })
#'   # combine plots
#'   if(length(p_list) > 1) {
#'     p <- patchwork::wrap_plots(p_list, nrow = 1, guides = "collect") +
#'       patchwork::plot_annotation(title = "Top ranked genes")
#'     w <- 6 # fig
#'     h <- 5 # fig
#'   } else {
#'     p <- p_list[[1]] + ggtitle("Top ranked genes")
#'     w <- 3
#'     h <- 5
#'   }
#'   # save plot
#'   if(inherits(outdir, "character")) {
#'     p_png <- file.path(outdir, "deseq_qc_top_gene.png")
#'     if(check_path(outdir)) {
#'       if(!file.exists(p_png) | overwrite) {
#'         ggsave(p_png, p, width = w, height = w, dpi = 200, units = "in")
#'       }
#'     }
#'   }
#'   p # return data
#' }
#'
#'
#' #' @describeIn deseq_qc_dist
#' #' check distance between samples
#' #'
#' #' @param x `DEseqDataSet`, see `deseq_qc_dds(x=)`
#' #' @param outdir character saving the results
#' #' @param transform bool compute tranformed data, `vst()`, `rlog()`
#' #' default: FALSE
#' #' @param overwrite bool overwrite exists file, default: FALSE
#' #'
#' #' @importFrom pheatmap pheatmap
#' #' @importFrom RColorBrewer brewer.pal
#' #' @importFrom ggplotify as.ggplot
#' #' @importFrom patchwork wrap_plots plot_annotation
#' #'
#' #' @return ggplot
#' #'
#' #' @export
#' deseq_qc_dist <- function(x, outdir = NULL, ...) {
#'   #-- Check: default values
#'   transform <- TRUE
#'   overwrite <- FALSE
#'   #-- Check: input
#'   dots <- rlang::list2(...)
#'   for(name in names(dots)) {
#'     assign(name, dots[[name]])
#'   }
#'   #-- Check: arguments
#'   qc_data <- deseq_qc_dds(x, outdir, transform)
#'   # if(!inherits(x, "DESeqDataSet")) {
#'   #   dds <- deseq_qc_dds(x, outdir, transform, return_dds = TRUE)
#'   # } else {
#'   #   dds <- x
#'   # }
#'   if(inherits(qc_data, "NULL")) {
#'     return(NULL)
#'   }
#'   # if(!inherits(dds, "DESeqDataSet")) {
#'   #   return(NULL)
#'   # }
#'   colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255)
#'   p_list <- lapply(names(qc_data), function(i) {
#'     sample_dist <- dist(t(assay(qc_data[[i]])))
#'     as.matrix(sample_dist) %>%
#'       pheatmap::pheatmap(
#'         clustering_distance_cols = sample_dist,
#'         clustering_distance_rows = sample_dist,
#'         color = colors,
#'         silent = TRUE) %>%
#'       ggplotify::as.ggplot() +
#'       ggtitle(i)
#'   })
#'   # combine plots
#'   if(length(p_list) > 1) {
#'     p <- patchwork::wrap_plots(p_list, nrow = 2, guides = "collect") +
#'       patchwork::plot_annotation(title = "Distance between samples")
#'     w <- 8 # fig
#'     h <- 8 # fig
#'   } else {
#'     p <- p_list[[1]] + ggtitle("Distance between samples")
#'     w <- 5
#'     h <- 5
#'   }
#'   # save plot
#'   if(inherits(outdir, "character")) {
#'     p_png <- file.path(outdir, "deseq_qc_dist.png")
#'     if(check_path(outdir)) {
#'       if(!file.exists(p_png) | overwrite) {
#'         ggsave(p_png, p, width = w, height = w, dpi = 200, units = "in")
#'       }
#'     }
#'   }
#'   p # return data
#' }
#'
#'
#' #' @describeIn deseq_qc_pca
#' #' check PCA plot for dds
#' #'
#' #' @param x `DEseqDataSet`, see `deseq_qc_dds(x=)`
#' #' @param outdir character saving the results
#' #' @param transform bool compute tranformed data, `vst()`, `rlog()`
#' #' default: FALSE
#' #' @param overwrite bool overwrite exists file, default: FALSE
#' #'
#' #' @improt ggplot2
#' #' @importFrom BiocGenerics plotPCA
#' #' @importFrom ggrepel geom_text_repel
#' #' @importFrom patchwork wrap_plots plot_annotation
#' #'
#' #' @return ggplot
#' #'
#' #' @export
#' deseq_qc_pca <- function(x, outdir = NULL, ...) {
#'   #-- Check: default values
#'   transform <- TRUE
#'   overwrite <- FALSE
#'   #-- Check: input
#'   dots <- rlang::list2(...)
#'   for(name in names(dots)) {
#'     assign(name, dots[[name]])
#'   }
#'   #-- arguments
#'   qc_data <- deseq_qc_dds(x, outdir, transform)
#'   # if(!inherits(x, "DESeqDataSet")) {
#'   #   dds <- deseq_qc_dds(x, outdir, transform, return_dds = TRUE)
#'   # } else {
#'   #   dds <- x
#'   # }
#'   if(inherits(qc_data, "NULL")) {
#'     return(NULL)
#'   }
#'   # if(!inherits(dds, "DESeqDataSet")) {
#'   #   return(NULL)
#'   # }
#'   p_list <- lapply(names(qc_data), function(i) {
#'     pca_data <- plotPCA(qc_data[[i]], returnData = TRUE)
#'     var <- round(attr(pca_data, "percentVar")*100, 1)
#'     var_text <- paste0(c("PC1", "PC2"), ": ", var, "% variance")
#'     # for batch
#'     if(! rlang::has_name(pca_data, "batch")) {
#'       pca_data$batch <- pca_data$condition
#'     }
#'     ggplot(pca_data, aes(PC1, PC2, color = condition, shape = batch)) +
#'       geom_point() +
#'       ggrepel::geom_text_repel(aes(label = name)) +
#'       xlab(var_text[1]) +
#'       ylab(var_text[2]) +
#'       ggtitle(i) +
#'       theme_bw()
#'   })
#'   # combine plots
#'   if(length(p_list) > 1) {
#'     p <- patchwork::wrap_plots(p_list, nrow = 2, guides = "collect") +
#'       patchwork::plot_annotation(title = "PCA plot")
#'     w <- 7 # fig
#'     h <- 6 # fig
#'   } else {
#'     p <- p_list[[1]] + ggtitle("PCA plot")
#'     w <- 4
#'     h <- 4
#'   }
#'   # save plot
#'   if(inherits(outdir, "character")) {
#'     p_png <- file.path(outdir, "deseq_qc_pca.png")
#'     if(check_path(outdir)) {
#'       if(!file.exists(p_png) | overwrite) {
#'         ggsave(p_png, p, width = w, height = w, dpi = 200, units = "in")
#'       }
#'     }
#'   }
#'   p # return data
#' }
#'
#'
#' #' @describeIn deseq_qc_ma
#' #' check MA plot for res
#' #'
#' #' @param x see `deseq_qc_res(x=)`
#' #' @param outdir character saving the results
#' #' @param ylim numeric, limits for log2FoldChange, default: NULL,
#' #' calculate by `scales::breaks_ext(n=5)(df$log2FoldChange)`
#' #' @param fc numeric cutoff for foldchange, default: 1, ignore foldchange
#' #' @param pvalue numeric cutoff for padj, default: 0.1, the main criteria
#' #' @param p_adjust bool use p-adjust value instead
#' #' @param overwrite bool overwrite exists file, default: FALSE
#' #'
#' #' @improt ggplot2
#' #' @importFrom ggrepel geom_text_repel
#' #' @importFrom patchwork wrap_plots plot_annotation
#' #'
#' #' @return ggplot
#' #'
#' #' @export
#' deseq_qc_ma <- function(x, outdir = "NULL", ...) {
#'   #-- Check: default values
#'   fc <- 2
#'   pvalue     <- 0.05
#'   p_adjust   <- TRUE
#'   ylim       <- c(-2, 2)
#'   label_list <- NULL
#'   label_max  <- 8
#'   overwrite  <- FALSE
#'   #-- Check: input
#'   dots <- rlang::list2(...)
#'   for(name in names(dots)) {
#'     assign(name, dots[[name]])
#'   }
#'   #-- Check: arguments: update log2fc
#'   df <- deseq_qc_res(x, outdir, log2fc_limits = ylim, ...)
#'   if(!inherits(df, "data.frame")) {
#'     warning("failed reading data from `res`")
#'     return(NULL)
#'   }
#'   #-- Check: limits for y, at least c(-2, 2)
#'   fclim <- df$log2FoldChange
#'   fclim <- fclim[!is.na(fclim)]
#'   if(inherits(ylim, "numeric")) {
#'     ylim <- c(max(min(ylim), min(fclim)),
#'               min(max(ylim), max(fclim)))
#'   } else {
#'     ylim <- c(min(fclim), max(fclim))
#'   }
#'   # foldchange, at least [-2, 2]
#'   if(min(ylim) > -2) ylim[1] = -2
#'   if(max(ylim) < 2)  ylim[2] = 2
#'   breaks <- scales::breaks_extended(n = 5)(ylim)
#'   # extend ylim by 10%
#'   ylim[1] <- ylim[1] * 1.1
#'   ylim[2] <- ylim[2] * 1.1
#'   #-- update data.frame by `ylim` ?
#'   df <- deseq_qc_res(x, outdir, log2fc_limits = ylim, ...)
#'   #-- run: title
#'   title <- glue::glue("criteria: foldChange >= {fc}, pvalue < {pvalue}")
#'   #-- plot
#'   p <- df %>%
#'     ggplot(aes(log10basemean, log2fc, color = sig, shape = ext)) +
#'     geom_hline(yintercept = c(-1, 1), size = .5, color = "grey30", linetype = 2) +
#'     geom_hline(yintercept = 0, size = .5, color = "grey30") +
#'     geom_point(size = .4) +
#'     scale_y_continuous(name   = "log2 fold change",
#'                        breaks = breaks,
#'                        limits = ylim) +
#'     # limits = c(min(breaks), max(breaks))) +
#'     scale_color_manual(values = c("up"   = "red",
#'                                   "not"  = "grey60",
#'                                   "down" = "blue")) +
#'     scale_shape_manual(values = c("up"   = 2,
#'                                   "dot"  = 20,
#'                                   "down" = 6)) +
#'     xlab("log10 mean of normalized counts") +
#'     ggtitle(title) +
#'     theme_bw() +
#'     theme(legend.position = "none")
#'   p <- deseq_qc_add_sig_label(p, ...)
#'   # save plot
#'   if(inherits(outdir, "character")) {
#'     p_png <- file.path(outdir, "deseq_qc_ma.png")
#'     if(check_path(outdir)) {
#'       if(!file.exists(p_png) | overwrite) {
#'         ggsave(p_png, p, width = 5, height = 4, dpi = 200, units = "in")
#'       }
#'     }
#'   }
#'   p # return data
#' }
#'
#'
#'
#' #' @describeIn deseq_qc_volcano
#' #' check volcano plot for res
#' #'
#' #' @param x see `deseq_qc_res(x=)`
#' #' @param outdir character saving the results
#' #' @param xlim numeric, limits for log2FoldChange, default: NULL,
#' #' calculate by `scales::breaks_ext(n=5)(df$log2FoldChange)`
#' #' @param fc numeric cutoff for foldchange, default: 1, ignore foldchange
#' #' @param pvalue numeric cutoff for padj, default: 0.1, the main criteria
#' #' @param p_adjust bool use p-adjust value instead
#' #' @param overwrite bool overwrite exists file, default: FALSE
#' #'
#' #' @improt ggplot2
#' #' @importFrom BiocGenerics plotPCA
#' #' @importFrom ggrepel geom_text_repel
#' #' @importFrom patchwork wrap_plots plot_annotation
#' #'
#' #' @return ggplot
#' #'
#' #' @export
#' deseq_qc_volcano <- function(x, outdir = NULL, ...) {
#'   #-- Check: default values
#'   fc <- 2
#'   pvalue     <- 0.05
#'   p_adjust   <- TRUE
#'   xlim       <- c(-2, 2)
#'   label_list <- NULL
#'   label_max  <- 8
#'   overwrite  <- FALSE
#'   #-- Check: input
#'   dots <- rlang::list2(...)
#'   for(name in names(dots)) {
#'     assign(name, dots[[name]])
#'   }
#'   #-- Check: input: update log2fc
#'   df <- deseq_qc_res(x, outdir, log2fc_limits = xlim, ...)
#'   if(!inherits(df, "data.frame")) {
#'     warning("failed reading data from `res`")
#'     return(NULL)
#'   }
#'   #-- Check: limits for y, at least c(-2, 2)
#'   fclim <- df$log2FoldChange
#'   fclim <- fclim[!is.na(fclim)]
#'   if(inherits(xlim, "numeric")) {
#'     xlim <- c(max(min(xlim), min(fclim)),
#'               min(max(xlim), max(fclim)))
#'   } else {
#'     xlim <- c(min(fclim), max(fclim))
#'   }
#'   # foldchange, at least [-2, 2]
#'   if(min(xlim) > -2) xlim[1] = -2
#'   if(max(xlim) < 2)  xlim[2] = 2
#'   breaks <- scales::breaks_extended(n = 5)(xlim)
#'   # extend ylim by 10%
#'   xlim[1] <- xlim[1] * 1.1
#'   xlim[2] <- xlim[2] * 1.1
#'   #-- update data.frame by `xlim` ?
#'   df <- deseq_qc_res(x, outdir, log2fc_limits = xlim, ...)
#'   #-- run: title
#'   title <- glue::glue("criteria: foldChange >= {fc}, pvalue < {pvalue}")
#'   #-- plot:
#'   p <- df %>%
#'     ggplot(aes(log2fc, log10pval, color = sig, shape = ext)) +
#'     geom_vline(xintercept = c(-1, 1), size = .5, color = "grey30", linetype = 2) +
#'     geom_vline(xintercept = 0, size = .5, color = "grey30") +
#'     geom_point(size = .5, alpha = 0.5) +
#'     scale_x_continuous(name   = "log2 fold change",
#'                        breaks = breaks,
#'                        limits = xlim) +
#'     scale_color_manual(values = c("up"   = "red",
#'                                   "not"  = "grey60",
#'                                   "down" = "blue")) +
#'     scale_shape_manual(values = c("up"   = 2,
#'                                   "dot"  = 20,
#'                                   "down" = 6)) +
#'     ylab("-log10 pvalue adjust") +
#'     ggtitle(title) +
#'     theme_bw() +
#'     theme(legend.position = "none")
#'   p <- deseq_qc_add_sig_label(p, ...)
#'   # save plot
#'   if(inherits(outdir, "character")) {
#'     p_png <- file.path(outdir, "deseq_qc_volcano.png")
#'     if(check_path(outdir)) {
#'       if(!file.exists(p_png) | overwrite) {
#'         ggsave(p_png, p, width = 4, height = 4, dpi = 200, units = "in")
#'       }
#'     }
#'   }
#'   p # return data
#' }
#'
#'
#' #' @describeIn deseq_qc_scatter
#' #' for scatter plot
#' #' x: log10(wt + 1)
#' #' y: log10(mut + 1)
#' #'
#' #' @param x see `deseq_qc_res(x=)`
#' #' @param outdir character saving the results
#' #' @param fc numeric cutoff for foldchange, default: 1, ignore foldchange
#' #' @param pvalue numeric cutoff for padj, default: 0.1, the main criteria
#' #' @param p_adjust bool use p-adjust value instead
#' #' @param density_points bool use `stat_density_2d()`, for
#' #'  large number dots
#' #' @param overwrite bool overwrite exists file, default: FALSE
#' #'
#' #' @improt ggplot2
#' #' @importFrom ggrepel geom_text_repel
#' #' @importFrom patchwork wrap_plots plot_annotation
#' #'
#' #' @return ggplot
#' #'
#' #' @export
#' deseq_qc_scatter <- function(x, outdir = NULL, ...) {
#'   #-- Check: default values
#'   fc <- 1
#'   pvalue     <- 0.1
#'   p_adjust   <- TRUE
#'   label_list <- NULL
#'   label_max  <- 8
#'   overwrite  <- FALSE
#'   density_points <- FALSE
#'   #-- Check: input
#'   dots <- rlang::list2(...)
#'   for(name in names(dots)) {
#'     assign(name, dots[[name]])
#'   }
#'   #-- Check: input
#'   df <- deseq_qc_res(x, outdir, ...)
#'   if(!inherits(df, "data.frame")) {
#'     warning("failed reading data from `res`")
#'     return(NULL)
#'   }
#'   #-- run: extract wt, mut
#'   wt  <- names(df)[2]
#'   mut <- names(df)[3]
#'   df <- df %>%
#'     dplyr::mutate(wt  = log10(!!as.name(wt) + 1),
#'                   mut = log10(!!as.name(mut) + 1))
#'   #-- run: determine limits, breaks
#'   breaks  <- scales::breaks_extended(n = 5)(c(df$wt, df$mut))
#'   xlimits <- c(min(c(df$wt, df$mut)), max(c(df$wt, df$mut)))
#'   ylimits <- xlimits
#'   #-- run: plot
#'   if(isTRUE(density_points)) {
#'     p1 <- df %>%
#'       ggplot(aes(wt, mut, color = sig)) +
#'       stat_density_2d(
#'         aes(fill = ..density..),
#'         data = dplyr::filter(df, sig == "not"),
#'         geom = "raster", contour = FALSE)
#'   } else {
#'     p1 <- df %>%
#'       ggplot(aes(wt, mut, color = sig)) +
#'       geom_point(size = .4, alpha = .5) #+
#'     # scale_size_manual(values = c("up" = .6, "not" = .3, "down" = .6))
#'   }
#'   #-- run: title
#'   title <- glue::glue("criteria: foldChange >= {fc}, pvalue < {pvalue}")
#'   p <- p1 +
#'     scale_color_manual(values = c("up"   = "red",
#'                                   "not"  = "grey60",
#'                                   "down" = "blue")) +
#'     scale_fill_gradient(low = "white", high = "black") +
#'     geom_abline(intercept = 0, slope = 1, linetype = 1, color = "grey30") +
#'     geom_abline(intercept = c(log10(2), -log10(2)), slope = 1, linetype = 2,
#'                 color = "grey50") +
#'     geom_point(data = dplyr::filter(df, sig %in% c("up", "down")),
#'                size = .6) +
#'     scale_x_continuous(breaks = breaks, limits = xlimits,
#'                        name = glue::glue("log10 count of {wt}")) +
#'     scale_y_continuous(breaks = breaks, limits = ylimits,
#'                        name = glue::glue("log10 count of {mut}")) +
#'     ggtitle(title) +
#'     theme_bw() +
#'     theme(panel.grid = element_blank())
#'   # add sig labels
#'   p <- deseq_qc_add_sig_label(p, ...)
#'   # save plot
#'   if(inherits(outdir, "character")) {
#'     p_png <- file.path(outdir, "deseq_qc_scatter.png")
#'     if(check_path(outdir)) {
#'       if(!file.exists(p_png) | overwrite) {
#'         ggsave(p_png, p, width = 5, height = 4, dpi = 200, units = "in")
#'       }
#'     }
#'   }
#'   p # return data
#' }
#'
#'
#'
#' #' @describeIn deseq_res_for_plot2
#' #' prepare res for plot
#' #' sig: for colors, ["up", "not", "down"]
#' #'   up: padj < pval & log2fc >= log2(fc)
#' #'   down: padj < pval & log2fc <= -log2(fc)
#' #'   not: is.na(padj) | padj >= pval
#' #'
#' #' ext: for shapes, ["up":2, "dot":20, "down":6]
#' #'   change dot shapes, based on the range of log2fc_limits
#' #'   default limits = `scales::breaks_ext(n=5)(log2FoldChange)`
#' #'   up: log2fc > max(log2fc_limits)
#' #'   down: log2fc < min(log2fc_limits)
#' #'   not: log2fc >= min(log2fc_limits) & log2fc <= max(log2fc_limits)
#' #'
#' #' @param res DESeqResults, data.frame or csv file
#' #' @param log2fc_limits numeric, limits for log2FoldChange, default: NULL,
#' #' calculate by `scales::breaks_ext(n=5)(df$log2FoldChange)`
#' #' @param fc numeric cutoff for foldchange, default: 1, ignore foldchange
#' #' @param pvalue numeric cutoff for padj, default: 0.1, the main criteria
#' #' @param p_adjust bool use p-adjust value instead
#' #' @importFrom scales breaks_ext
#' #'
#' #' @return data.frame
#' #'
#' #' @export
#' deseq_res_for_plot2 <- function(res, log2fc_limits = NULL,
#'                                 fc = 1, pvalue = 0.1, p_adjust = p_adjust) {
#'   #-- Check: arguments
#'   if(inherits(res, "DESeqResults")) {
#'     df <- as.data.frame(res)
#'   } else if(inherits(res, "data.frame")) {
#'     df <- res
#'   } else if(inherits(res, "character")) {
#'     if(file.exists(res) & endsWith(res, ".csv")) {
#'       df <- read.csv(res)
#'     } else {
#'       df <- data.frame()
#'     }
#'   } else {
#'     warning(glue::glue("unknown res, expect `DESeqResults` or `data.frame`, ",
#'                        "get: {class(res)}"))
#'     return(NULL)
#'   }
#'   #-- Check: required columns
#'   rc <- c("baseMean", "log2FoldChange", "padj")
#'   if(!all(rc %in% names(df))) {
#'     rc_str <- paste(rc, collapse = ", ")
#'     warning(glue::glue("missing columns: [{rc_str}]"))
#'     return(NULL)
#'   }
#'   #-- Check: limits for x, y
#'   if(inherits(log2fc_limits, "numeric")) {
#'     breaks <- log2fc_limits
#'   } else {
#'     breaks <- scales::breaks_extended(n = 5)(df$log2FoldChange)
#'   }
#'   #-- run: sig for colors, up, not, down
#'   # in case, pvalue exists in data.frame
#'   pval <- pvalue
#'   df <- df %>%
#'     dplyr::mutate(
#'       log10pval = -log10(padj),
#'       log10basemean = log10(baseMean + 1),
#'       sig = ifelse(
#'         is.na(padj), "not", ifelse(
#'           padj < pval, ifelse(
#'             log2FoldChange >= log2(fc), "up", ifelse(
#'               log2FoldChange <= -log2(fc), "down", "not")), "not")))
#'   #-- run: ext for shapes, up, dot, down
#'   # fix lgo2fc, by shapes
#'   df %>%
#'     dplyr::mutate(
#'       ext = ifelse(
#'         is.na(log2FoldChange), "dot", ifelse(
#'           log2FoldChange > max(breaks), "up", ifelse(
#'             log2FoldChange < min(breaks), "down", "dot")))) %>%
#'     dplyr::mutate(log2fc = ifelse(
#'       ext == "up", max(breaks), ifelse(
#'         ext == "down", min(breaks), log2FoldChange)))
#' }
#'
#'
#'
#'
#' #' @describeIn deseq_qc_add_sig_label
#' #' add sig labels to plot, based on `sig` column and
#' #' is designed for:
#' #' `deseq_qc_ma()`, `deseq_qc_volcano()`, `deseq_qc_scatter()`
#' #'
#' #' @param x ggplot
#' #' @param label_list character gene names
#' #' priority: symbol > gene_id
#' #'
#' #' @return ggplot
#' #'
#' #' @export
#' # deseq_qc_add_sig_label <- function(x, label_list = NULL, label_max = 8) {
#' deseq_qc_add_sig_label <- function(x, ...) {
#'   label_list <- NULL
#'   label_max  <- 8
#'   dots <- rlang::list2(...)
#'   for(name in names(dots)) {
#'     assign(name, dots[[name]])
#'   }
#'   if(!inherits(x, "ggplot")) {
#'     message(glue::glue("x is {class(x)}, expect `ggplot`"))
#'     return(x)
#'   }
#'   # retrieve the `data` from `ggplot`
#'   if("data" %in% names(x)) {
#'     df <- x$data
#'     rc <- c("gene_id", "sig", "padj")
#'     if(all(rc %in% names(df))) {
#'       #-- run: choose label_list genes
#'       t1 <- df$gene_id %in% label_list
#'       if("SYMBOL" %in% names(df)) {
#'         t2 <- df$SYMBOL %in% label_list
#'         t1 <- t1 | t2
#'       }
#'       df1 <- df[t1, ]
#'       #-- run: show top sig genes
#'       if(nrow(df1) < label_max) {
#'         df2 <- df %>%
#'           dplyr::arrange(padj) %>%
#'           head(label_max - nrow(df1))
#'         df1 <- rbind(df1, df2)
#'       }
#'       #-- run: show label
#'       if(nrow(df1) > 0) {
#'         label_col <- ifelse("SYMBOL" %in% names(df), "SYMBOL", "gene_id")
#'         x +
#'           ggrepel::geom_text_repel(
#'             mapping = aes(label = !!as.name(label_col)),
#'             data    = df1,
#'             color          = "grey10",
#'             size           = 4,
#'             force          = .2,
#'             direction      = "both",
#'             point.padding  = .2,
#'             max.overlaps   = Inf,
#'             # max.overlaps   = 30,
#'             box.padding    = 0.5, # additional padding around each text label
#'             segment.color  = "grey20",
#'             # segment.size   = .4,
#'             min.segment.length = 0, # draw all line segments
#'             max.time = 1, max.iter = 1e5 # stop after 1 second, or after 100,000 iterations
#'           ) +
#'           geom_point(
#'             data = df1, size = .6, shape = 20
#'           )
#'       }
#'     } else {
#'       rc_str <- paste(rc, collapse = ", ")
#'       message(glue::glue("missing columns: {rc_str}"))
#'     }
#'   }
#' }






# # to-do: deprecated
# deseq_qc_res <- function(x, ...) {
#   #-- default values
#   shrink_method <- NULL # null, normal, apeglm, ashr
#   log2fc_limits <- c(-2, 2)
#   fc        <- 2
#   pvalue    <- 0.05
#   p_adjust  <- TRUE
#   overwrite <- FALSE
#   dots <- rlang::list2(...)
#   for(name in names(dots)) {
#     assign(name, dots[[name]])
#   }
#   #-- Check: exists file
#   qc_res_rds <- file.path(outdir, "deseq_qc_res.rds")
#   if(inherits(qc_res_rds, "character")) {
#     if(file.exists(qc_res_rds) & !overwrite) {
#       message(glue::glue("loading deseq_qc_res from file: {qc_res_rds}"))
#       return(readRDS(qc_res_rds))
#     }
#   }
#   #-- Check: x, DESeqDataSet obj
#   dds <- NULL # init
#   if(inherits(x, "DESeqDataSet")) { # dds
#     dds <- x #
#   } else if(inherits(x, "character")) {
#     if(file.exists(x)) {
#       if(endsWith(x, ".csv")) { # norm_table.fix.csv
#         dds <- read.csv(x)
#       } else if(endsWith(x, ".rds")) { # deseq_dds.rds
#         dds <- readRDS(x)
#       } else {
#         message(glue::glue("unknown x: {x}"))
#       }
#     }
#   } else if(inherits(x, "data.frame")) {
#     dds <- x # ?! watching
#   }
#   #-- Check: outdir, looking for `norm_table.fix.csv`; priority first
#   if(inherits(outdir, "character")) {
#     fs <- file.path(outdir, "norm_table.fix.xls") # dds
#     message(glue::glue("will load `table` from file \n",
#                        "you can change `outdir` to `NULL` or a new directory; ",
#                        "and re-run this function again; \n",
#                        "file: {fs}"))
#     if(file.exists(fs)) {
#       dds <- read.csv(fs) # overwrite
#     }
#   }
#   #-- Check: update table
#   if(inherits(dds, "data.frame")) {
#     df <- dds
#   } else if(inherits(dds, "DESeqDataSet")) {
#     df <- deseq_mean(dds, outdir) # see: deseq_utils.R
#   } else {
#     warning(glue::glue(
#       "x is {class(x)}, expect `DESeqDataSet` or `deseq_dds.rds` file, ",
#       "specify either `x=` or `outdir=`, or both "))
#     return(NULL)
#   }
#   #-- Check: table, required columns
#   if(inherits(df, "data.frame")) {
#     rc <- c("baseMean", "log2FoldChange", "padj")
#     if(!all(rc %in% names(df))) {
#       rc_str <- paste(rc, collapse = ", ")
#       warning(glue::glue("missing columns: [{rc_str}]"))
#       return(NULL)
#     }
#   } else {
#     warning(glue::glue("failed loading data from x={x}, outdir={outdir}"))
#     return(NULL)
#   }
#   #-- Check: for scatter plot: wt, mut
#   if(inherits(dds, "DESeqDataSet")) {
#     coldata <- SummarizedExperiment::colData(dds)
#     wt  <- levels(coldata$condition)[1]
#     mut <- levels(coldata$condition)[2]
#   } else {
#     # gene_id, wt, mut, ...
#     wt  <- names(df)[2]
#     mut <- names(df)[3]
#   }
#   if(!all(c(wt, mut) %in% names(df))) {
#     warning(glue::glue("missing reuired columns: {wt}, {mut}"))
#     return(NULL)
#   }
#   #-- Check: update log2fc, for ma, volcano
#   if(inherits(log2fc_limits, "numeric")) {
#     breaks <- scales::breaks_extended(n = 5)(log2fc_limits)
#   } else {
#     breaks <- scales::breaks_extended(n = 5)(df$log2FoldChange)
#   }
#   #-- run: sig for colors, up, not, down
#   # in case, pvalue exists in data.frame
#   if(!rlang::has_name(df, "sig")) {
#     # df <- get_sig_name(df, fc, pvalue, p_adjust, return_dataframe = TRUE)
#     df <- get_sig_name(df, return_dataframe = TRUE, ...)
#   }
#   #-- run: ext
#   # for shapes, up, dot, down
#   # update log2fc, check outlier by log2fc_limits
#   out <- df %>%
#     dplyr::mutate(
#       log10pval     = -log10(padj),
#       log10basemean = log10(baseMean + 1),
#       ext = ifelse(
#         is.na(log2FoldChange), "dot", ifelse(
#           log2FoldChange > max(breaks), "up", ifelse(
#             log2FoldChange < min(breaks), "down", "dot")))) %>%
#     dplyr::mutate(log2fc = ifelse(
#       ext == "up", max(breaks), ifelse(
#         ext == "down", min(breaks), log2FoldChange)))
#   #-- save to file
#   qc_res_rds <- file.path(outdir, "deseq_qc_res.rds")
#   if(inherits(qc_res_rds, "character")) {
#     if(check_path(outdir)) {
#       if(!file.exists(qc_res_rds) | overwrite) {
#         saveRDS(out, qc_res_rds)
#       }
#     }
#   }
#   #-- return
#   out
# }



#' #' @describeIn deseq_res_for_plot2
#' #' prepare res for plot
#' #' sig: for colors, ["up", "not", "down"]
#' #'   up: padj < pval & log2fc >= log2(fc)
#' #'   down: padj < pval & log2fc <= -log2(fc)
#' #'   not: is.na(padj) | padj >= pval
#' #'
#' #' ext: for shapes, ["up":2, "dot":20, "down":6]
#' #'   change dot shapes, based on the range of log2fc_limits
#' #'   default limits = `scales::breaks_ext(n=5)(log2FoldChange)`
#' #'   up: log2fc > max(log2fc_limits)
#' #'   down: log2fc < min(log2fc_limits)
#' #'   not: log2fc >= min(log2fc_limits) & log2fc <= max(log2fc_limits)
#' #'
#' #' @param res DESeqResults, data.frame or csv file
#' #' @param log2fc_limits numeric, limits for log2FoldChange, default: NULL,
#' #' calculate by `scales::breaks_ext(n=5)(df$log2FoldChange)`
#' #' @param fc numeric cutoff for foldchange, default: 1, ignore foldchange
#' #' @param pvalue numeric cutoff for padj, default: 0.1, the main criteria
#' #' @param p_adjust bool use p-adjust value instead
#' #' @importFrom scales breaks_ext
#' #'
#' #' @return data.frame
#' #'
#' #' @export
#' deseq_res_for_plot2 <- function(res, log2fc_limits = NULL,
#'                                 fc = 1, pvalue = 0.1, p_adjust = p_adjust) {
#'   #-- Check: arguments
#'   if(inherits(res, "DESeqResults")) {
#'     df <- as.data.frame(res)
#'   } else if(inherits(res, "data.frame")) {
#'     df <- res
#'   } else if(inherits(res, "character")) {
#'     if(file.exists(res) & endsWith(res, ".csv")) {
#'       df <- read.csv(res)
#'     } else {
#'       df <- data.frame()
#'     }
#'   } else {
#'     warning(glue::glue("unknown res, expect `DESeqResults` or `data.frame`, ",
#'                        "get: {class(res)}"))
#'     return(NULL)
#'   }
#'   #-- Check: required columns
#'   rc <- c("baseMean", "log2FoldChange", "padj")
#'   if(!all(rc %in% names(df))) {
#'     rc_str <- paste(rc, collapse = ", ")
#'     warning(glue::glue("missing columns: [{rc_str}]"))
#'     return(NULL)
#'   }
#'   #-- Check: limits for x, y
#'   if(inherits(log2fc_limits, "numeric")) {
#'     breaks <- log2fc_limits
#'   } else {
#'     breaks <- scales::breaks_extended(n = 5)(df$log2FoldChange)
#'   }
#'   #-- run: sig for colors, up, not, down
#'   # in case, pvalue exists in data.frame
#'   pval <- pvalue
#'   df <- df %>%
#'     dplyr::mutate(
#'       log10pval = -log10(padj),
#'       log10basemean = log10(baseMean + 1),
#'       sig = ifelse(
#'         is.na(padj), "not", ifelse(
#'           padj < pval, ifelse(
#'             log2FoldChange >= log2(fc), "up", ifelse(
#'               log2FoldChange <= -log2(fc), "down", "not")), "not")))
#'   #-- run: ext for shapes, up, dot, down
#'   # fix lgo2fc, by shapes
#'   df %>%
#'     dplyr::mutate(
#'       ext = ifelse(
#'         is.na(log2FoldChange), "dot", ifelse(
#'           log2FoldChange > max(breaks), "up", ifelse(
#'             log2FoldChange < min(breaks), "down", "dot")))) %>%
#'     dplyr::mutate(log2fc = ifelse(
#'       ext == "up", max(breaks), ifelse(
#'         ext == "down", min(breaks), log2FoldChange)))
#' }
#'
#'
#'

