#' Functions for DESeq2 down-stream analysis
#'
#' Reading files, directories, pipeline, ...
#' Processing data
#' Prepare for plot
#' Plotting
#'
#' @name deseq


#' @describeIn deseq
#'
#' @param dds DESeqDataSet
#' @param outdir character saving the results
#' @param strandness character could be "sens", "anti", default "sens"
#' @param shrink character use `lfcShrink` function to calculate shrunken LFC
#'        could be ["apeglm", "ashr", "normal"], default: "apeglm"
#' @param overwrite bool overwrite exists file, default: FALSE
#'
#' @import DESeq2
#' @import apeglm
#' @import ggplot2
#'
#' @return DESeqResults, res (shrinked)
#'
#' @export
deseq <- function(dds, outdir = NULL, shrink = "apeglm", cpu = 4,
                  fc = 1, pvalue = 0.1, p_adjust = TRUE,
                  readable = TRUE, genome = NULL, overwrite = FALSE) {
  if(is.null(outdir)) {
    outdir <- tempdir()
  }
  dd <- run_deseq(dds, outdir, shrink, cpu, overwrite)
  if(is.null(dd)) {
    warning("`deseq()` failed, see above messages")
    return(NULL)
  }
  if(shrink %in% c("apeglm", "ashr", "normal")) {
    res <- dd$reslfc
  } else {
    res <- dd$res
  }
  if(check_path(outdir)) {
    #-- Check: run set_readable()?
    .to_readable <- function(x) {
      b1 <- isTRUE(readable) & inherits(genome, "character")
      b2 <- isTRUE(overwrite)
      b3 <- !file.exists(x)
      b1 & (b2 | b3)
    }
    #-- run: saving norm table
    norm_table <- file.path(outdir, "norm_table.csv")
    df1a <- DESeq2::counts(dd$dds, normalized = TRUE) # normalized counts
    df1  <- merge(as.data.frame(df1a), as.data.frame(res), by = "row.names")
    colnames(df1)[1] <- "gene_id"
    # if(.to_readable(norm_table)) {
    #   df1 <- set_readable(df1, genome)
    # }
    write.csv(df1, norm_table, quote = TRUE, row.names = FALSE)
    #-- run: saving norm table, fix
    norm_fix_table <- file.path(outdir, "norm_table.fix.csv")
    df2 <- deseq_mean(dds, outdir)
    if(.to_readable(norm_fix_table)) {
      df2 <- set_readable(df2, genome)
    }
    write.csv(df2, norm_fix_table, quote = TRUE, row.names = FALSE)
    #-- run: fpkm
    fpkm_table <- file.path(outdir, "fpkm_table.csv")
    if("basepairs" %in% names(mcols(dd$dds))) {
      df3a <- DESeq2::fpkm(dd$dds)
      df3  <- merge(as.data.frame(df3a), as.data.frame(res), by = "row.names")
      if(.to_readable(fpkm_table)) {
        df3 <- set_readable(df3, genome)
      }
      write.csv(df3, fpkm_table, quote = TRUE, row.names = FALSE)
    }
    #-- run: quality-control, require outdir
    tmp <- deseq_qc(norm_fix_table, outdir, n_max = 20,
                    fc = fc, pvalue = pvalue, p_adjust = p_adjust,
                    overwrite = overwrite)
  }
  res
}


#' @describeIn run_deseq
#' run DESeq2::results() for dds
#'
#' @param dds DESeqDataSet
#' @param outdir character saving the results
#' @param shrink character use `lfcShrink` function to calculate shrunken LFC
#'        could be ["apeglm", "ashr", "normal"], default: "apeglm"
#' @param cpu integer, number of CPU to run in parallel, default: 4
#' @param overwrite bool overwrite exists file, default: FALSE
#'
#' @export
run_deseq <- function(dds, outdir = NULL, shrink = "apeglm", cpu = 4,
                      overwrite = FALSE) {
  #-- Check: pre-data
  if(inherits(outdir, "character")) {
    deseq_data_rds <- file.path(outdir, "deseq_data.rds")
    if(file.exists(deseq_data_rds) & !isTRUE(overwrite)) {
      message(glue::glue("loading deseq data from: {deseq_data_rds}"))
      return(readRDS(deseq_data_rds))
    }
  }
  #-- Check: arguments
  if(!inherits(dds, "DESeqDataSet")) {
    warning(glue::glue("expect `DESeqDataSet`, got {class(dds)}"))
    return(NULL)
  }
  coldata <- colData(dds)
  if(! rlang::has_name(coldata, "condition")) {
    warning("`condition` not found, check `colData(dds)`")
    return(NULL)
  }
  #-- run: DESeq analysis
  # levels wt mut #
  # dds <- DESeq2::DESeq(dds)
  wt  <- levels(coldata$condition)[1] #
  mut <- levels(coldata$condition)[2] #
  coef <- deseq_sanitize_str(paste0("condition_", mut, "_vs_", wt))
  BiocParallel::register(BiocParallel::MulticoreParam(cpu))
  res <- DESeq2::results(dds, contrast = c("condition", mut, wt),
                         parallel = TRUE)
  res <- res[order(res$padj), ] # Order by adjusted p-value
  # updated values
  if(shrink %in% c("normal", "apeglm", "ashr")) {
    message(glue::glue("shrink log2 fold changes by: {shrink}"))
    reslfc <- DESeq2::lfcShrink(dds, coef = coef, type = shrink,
                                parallel = TRUE)
  } else {
    reslfc <- NULL # skipped
  }
  out <- list(
    dds    = dds,
    res    = res,
    reslfc = reslfc,
    wt     = wt,
    mut    = mut
  )
  if(inherits(outdir, "character")) {
    deseq_data_rds <- file.path(outdir, "deseq_data.rds")
    if(!file.exists(deseq_data_rds) | isTRUE(overwrite)) {
      check_path(outdir)
      message(glue::glue("setting tempdir outdir = {outdir}"))
      saveRDS(out, deseq_data_rds)
    }
  }
  out
}























