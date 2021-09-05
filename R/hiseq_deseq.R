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
hiseq_deseq <- function(x, outdir = NULL, ...) {
  #-- Check: default values
  fix_batch  <- TRUE   # for DESeq(), design: ~ condition + batch
  shrink     <- "apeglm"   # shrink log2fc, "apeglm", "ashr", "normal"
  n_max      <- 20
  cpu        <- 4
  fc         <- 2
  pvalue     <- 0.1
  p_adjust   <- TRUE   # for DESeq2 `padj`
  genome     <- NULL   # character
  strandness <- "sens" # "sens", "anti"
  readable   <- TRUE   # add SYMBOL, ENTREZID, by gene_id
  overwrite  <- FALSE
  transform  <- TRUE   # dds transformation, vst(), vlog()
  label_list <- NULL   # for ma,volcano,scatter
  label_max  <- 8      # for ma,volcano,scatter
  density_points <- FALSE  # for scatter plot
  log2fc_limits  <- c(-2, 2)   # for ma,volcano,scatter
  #-- Check: arguments
  dots <- rlang::list2(...)
  for(name in names(dots)) {
    assign(name, dots[[name]])
  }
  #-- Check: args
  if(!is_hiseq_dir(x, "_rx")) {
    warning(glue::glue("not a `rnaseq_rx` dir: {x}"))
    return(NULL)
  }
  if(!inherits(outdir, "character")) {
    outdir <- list_hiseq_file(x, "deseq_dir", "_rx")
  }
  outdir <- normalizePath(outdir) # absolute path
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
  if(inherits(dds, "DESeqDataSet")) {
    genome <- list_hiseq_file(x, "genome", "rx") # add genome
    dots$genome <- genome
    res <- deseq(dds, outdir, !!!dots)
  } else {
    warning("`deseq()` failed")
    return(NULL)
  }
  #-- run: save config
  name_csv <- file.path(outdir, "smp_name.csv")
  name_rds <- file.path(outdir, "smp_name.rds")
  coldata  <- colData(dds)
  name_df  <- data.frame(
    smp_name  = coldata$smp_name,
    label     = rownames(coldata),
    condition = coldata$condition
  )
  saveRDS(name_df, name_rds)
  write.csv(name_df, name_csv, row.names = FALSE)
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



