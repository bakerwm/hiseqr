#' Functions for DESeq2, dds quality control
#'
#' deseq_qc_counts()
#' deseq_qc_mean_sd()
#' deseq_qc_top_gene()
#' deseq_qc_dist()
#' deseq_qc_pca()
#' deseq_qc_ma()
#' deseq_qc_volcano()
#' deseq_qc_scatter()
#'
#' to-to:
#' - simpify arguments by `...`
#' - distance between samples (bam cor)
#'
#' @name deseq_qc


#' @describeIn deseq_qc quality control for dds object
#'
#' @param x character file `norm_table.fix.csv`, or `DESeqDataSet`
#' also support the following data types: (be careful)
#' `data.frame`, `.rds`#'
#'
#' generating the following plots:
#'
#' deseq_qc_counts()
#' deseq_qc_mean_sd() # standard deviation
#' deseq_qc_top_gene()
#' deseq_qc_dist() # correlation
#' deseq_qc_pca()  # correlation
#' deseq_qc_ma()
#' deseq_qc_volcano()
#' deseq_qc_scatter()
#'
#' @param outdir character saving the results
#' @param log2fc_limits numeric setting the range of log2FoldChange on plots
#' default: NULL, by `scales::breaks_extend()`
#' @param transform logical compute tranformed data, `vst()`, `rlog()`
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
deseq_qc <- function(x, ...) {
  #----------------------------------------------------------------------------#
  #-- Check: default values
  dots <- rlang::list2(...)
  args <- rlang::list2(
    n_max = 20,
    fc = 2,
    pvalue     = 0.05,
    p_adjust   = TRUE,
    ylim       = c(-2, 2),
    label_list = NULL,
    label_max  = 8,
    overwrite  = FALSE,
    readable   = TRUE,
    transform  = TRUE,
    density_points = FALSE,
    # .col_label = "gene_id",
    shrink_method = "standard" # apeglm, ashr, normal
  )
  #-- update, for child functions
  args <- purrr::list_modify(args, !!!dots)
  for(name in names(args)) {
    if(rlang::is_empty(name)) next
    assign(name, args[[name]])
  }
  #----------------------------------------------------------------------------#
  if(!is_hiseq_dir(x, "deseq_deseq2")) {
    message(glue::glue(
      "'x' is {class(x)}, expect 'deseq_deseq2' directory"
    ))
    return(NULL)
  }
  dds_rds <- list_hiseq_file(x, "deseq_dds_rds", hiseq_type = "deseq_deseq2")
  if(file.exists(dds_rds)) {
    dds       <- readRDS(dds_rds)
    dds_trans <- deseq_qc_dds(dds, outdir = x, return_data = "dds_trans")
    res_lfc   <- deseq_qc_dds(dds, outdir = x, return_data = "res_lfc")
    # 1. counts
    png1 <- file.path(x, "deseq_qc_counts.png")
    p1 <- deseq_qc_counts(dds, !!!args)
    ggsave(png1, p1, width = 6, height = 6)
    # 2. mean sd; transformd
    tmp <- lapply(names(dds_trans), function(s) {
      png2 <- file.path(x, paste0("deseq_qc_mean_sd.", s, ".png"))
      p2 <- deseq_qc_mean_sd(dds, transform_method = s, outdir = x, !!!args)
      ggsave(png2, p2, width = 4, height = 3)
    })
    # 3. top genes
    tmp <- lapply(names(dds_trans), function(s) {
      png3 <- file.path(x, paste0("deseq_qc_top_gene.", s, ".png"))
      p3 <- deseq_qc_top_gene(dds, transform_method = s, outdir = x, !!!args)
      ggsave(png3, p3, width = 3, height = 5)
    })
    # 4. dist
    tmp <- lapply(names(dds_trans), function(s) {
      png4 <- file.path(x, paste0("deseq_qc_dist.", s, ".png"))
      p4 <- deseq_qc_dist(dds, transform_method = s, outdir = x, !!!args)
      ggsave(png4, p4, width = 5, height = 4.5)
    })
    # 5. pca
    tmp <- lapply(names(dds_trans), function(s) {
      png5 <- file.path(x, paste0("deseq_qc_pca.", s, ".png"))
      p5 <- deseq_qc_pca(dds, transform_method = s, outdir = x, !!!args)
      ggsave(png5, p5, width = 6, height = 4)
    })
    # for sig plots
    # 6. ma
    tmp <- lapply(names(res_lfc), function(s) {
      png6 <- file.path(x, paste0("deseq_qc_ma.", s, ".png"))
      p6 <- deseq_qc_ma(x, !!!args, shrink_method = s)
      p6 <- deseq_qc_add_sig_label(p6)
      ggsave(png6, p6, width = 5, height = 4)
    })
    # 7. volcano
    tmp <- lapply(names(res_lfc), function(s) {
      png7 <- file.path(x, paste0("deseq_qc_volcano.", s, ".png"))
      p7 <- deseq_qc_volcano(x, !!!args, shrink_method = s)
      p7 <- deseq_qc_add_sig_label(p7)
      ggsave(png7, p7, width = 5, height = 5)
    })
    # 8. ma
    tmp <- lapply(names(res_lfc), function(s) {
      png8 <- file.path(x, paste0("deseq_qc_scatter.", s, ".png"))
      p8 <- deseq_qc_scatter(x, !!!args, shrink_method = s)
      p8 <- deseq_qc_add_sig_label(p8)
      ggsave(png8, p8, width = 5, height = 4)
    })
  }
}



