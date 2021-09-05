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
deseq_qc <- function(x, outdir = NULL, ...) {
  #-- default values
  log2fc_limits = c(-2, 2)
  transform <- TRUE
  n_max <- 20
  fc <- 2
  pvalue <- 0.05
  p_adjust   <- TRUE
  label_list <- NULL
  label_max  <- 8
  overwrite  <- FALSE
  density_points <- FALSE
  #-- arguments
  dots <- rlang::list2(...)
  for(name in names(dots)) {
    assign(name, dots[[name]])
  }
  #-- qc: dds
  # qc_dds <- deseq_qc_dds(x, outdir = outdir, return_dds = TRUE)
  # qc_res <- deseq_qc_res(x, outdir = outdir, ...)
  if(inherits(qc_dds, "DESeqDataSet")) {
    qc1 <- list(
      counts   = deseq_qc_counts(x, ...),
      mean_sd  = deseq_qc_mean_sd(x, ...),
      top_gene = deseq_qc_top_gene(x, ...),
      dist     = deseq_qc_dist(x, ...),
      pca      = deseq_qc_pca(x, ...))
  } else {
    qc1 <- list(counts = NULL, mean_sd = NULL, top_gene = NULL, dist = NULL,
                pca = NULL)
  }
  #-- qc: res
  qc_res <- deseq_qc_res(x, ...)
  if(inherits(qc_res, "data.frame")) {
    qc2 <- list(
      ma      = deseq_qc_ma(x, ...),
      volcano = deseq_qc_volcano(x, ...),
      scatter = deseq_qc_scatter(x, ...)
    )
    print(label_list)
  } else {
    qc2 <- list(ma = NULL, volcano = NULL, scatter = NULL)
  }
  #-- return
  c(qc1, qc2)
}



