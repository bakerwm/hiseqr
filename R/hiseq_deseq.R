#' Functions for DESeq2 down-stream analysis
#'
#' Reading files, directories, pipeline, ...
#' Processing data
#' Prepare for plot
#' Plotting
#'
#' @import prep_hiseq_desqe deseq
#'
#'
#' @name hiseq_deseq




#' @describeIn hiseq_deseq
#'
#' @param dds DESeqDataSet
#' @param outdir character saving the results
#' @param strandness character could be "sens", "anti", default "sens"
#' @param fix_batch bool fix batch effect, default: TRUE
#' @param shrink character use `lfcShrink` function to calculate shrunken LFC
#'        could be ["apeglm", "ashr", "normal"], default: "apeglm"
#' @param cpu integer, number of CPU to run in parallel, default: 4
#' @param overwrite bool overwrite exists file, default: FALSE
#'
#' @import DESeq2
#' @import apeglm
#' @import ggplot2
#'
#' @return DESeqResults, res (shrinked)
#'
#' @export
hiseq_deseq <- function(x, outdir = NULL, shrink = "apeglm", fix_batch = TRUE,
                        fc = 1, pvalue = 0.1, p_adjust = TRUE,
                        strandness = "sens", cpu = 4, readable = TRUE,
                        genome = NULL, overwrite = FALSE) {
  #-- Check: args
  if(!is_hiseq_dir(x, "_rx")) {
    warning(glue::glue("not a `rnaseq_rx` dir: {x}"))
    return(NULL)
  }
  if(!inherits(outdir, "character")) {
    outdir <- list_hiseq_file(x, "deseq_dir", "_rx")
  }
  #-- Check: pre-compute data
  dds_rds <- file.path(outdir, "deseq_dds.rds")
  if(file.exists(dds_rds)) {
      dds <- readRDS(dds_rds)
  } else {
    dds <- hiseq_prep_deseq(x, strandness, fix_batch)
    if((check_path(outdir) & !file.exists(dds_rds)) | overwrite) {
      saveRDS(dds, dds_rds)
    }
  }
  #-- run: deseq
  genome <- list_hiseq_file(x, "genome", "rx")
  if(inherits(dds, "DESeqDataSet")) {
    res <- deseq(dds, outdir, shrink, cpu, fc = fc, pvalue = pvalue,
                 p_adjust = p_adjust, genome = genome, overwrite = overwrite)
  } else {
    warning("`deseq()` failed")
    return(NULL)
  }
  #-- run: transcripts_deseq2.csv
  norm_table <- file.path(outdir, "norm_table.csv")
  res_df <- read.csv(norm_table)
  res_csv <- file.path(outdir, "transcripts_deseq2.csv")
  write.csv(res_df, res_csv, quote = TRUE, row.names = FALSE)
  #-- run: fix mean
  norm_fix_table <- file.path(outdir, "norm_table.fix.csv")
  res_fix_df <- read.csv(norm_fix_table)
  fix_csv <- file.path(outdir, "transcripts_deseq2.fix.csv")
  write.csv(res_fix_df, fix_csv, quote = TRUE, row.names = FALSE)
  # return
  res
}





















