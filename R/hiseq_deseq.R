#' Functions for DESeq2 down-stream analysis
#'
#' Reading files, directories, pipeline, ...
#' Processing data
#' Prepare for plot
#' Plotting
#'
#' @import prep_hiseq_desqe deseq
#'
#' @name hiseq_deseq


#' @describeIn hiseq_deseq
#'
#' @param x character path to the 'rnaseq_rx' directory
#' @param outdir character saving the results
#' @param strandness character could be "sens", "anti", default "sens"
#' @param fix_batch bool fix batch effect, default: TRUE
#' @param shrink logical `shrink` LFC by ["apeglm", "ashr", "normal"],
#'   default: TRUE
#' @param transform logical transform dds by `vst()`, `rlog()`, default: TRUE
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
hiseq_deseq <- function(x, ...) {
  #----------------------------------------------------------------------------#
  #-- Check: args
  if(!is_hiseq_dir(x, "_rx")) {
    warning(glue::glue("not a `rnaseq_rx` dir: {x}"))
    return(NULL)
  }
  #-- Check: default values
  dots <- rlang::list2(...)
  args <- rlang::list2(
    outdir     = NULL, # default, tempdir()
    strandness = "sens", # "sens", "anti"
    fix_batch  = TRUE, # for DESeq(), design: ~ condition + batch
    shrink     = TRUE, # lfcShrink(), normal, apeglm, ashr
    transform  = TRUE, # dds transformation, vst(), vlog()
    n_max      = 20,   # for top_gene, counts,
    cpu        = 2,    # for DESeq2::DESeq(), in parallel
    fc         = 2,    # for sig genes; fc + pvalue
    pvalue     = 0.05, # for sig genes; fc + pvalue
    p_adjust   = TRUE, # for sig genes; fc + pvalue
    genome     = NULL, # character
    readable   = TRUE, # add SYMBOL, ENTREZID, by gene_id
    label_list = NULL,  # for ma,volcano,scatter
    label_max  = 8,     # for ma,volcano,scatter
    density_points = FALSE,   # for scatter plot
    log2fc_limits  = c(-2, 2), # for ma,volcano,scatter
    overwrite  = FALSE
  )
  args <- purrr::list_modify(args, !!!dots)
  #-- update: outdir
  if(!inherits(args$outdir, "character")) {
    args$outdir <- list_hiseq_file(x, "deseq_dir", "_rx")
  }
  args$outdir <- normalizePath(args$outdir) # absolute path
  # #-- update args, for child functions
  # dots_args <- lapply(names(args), function(i) {
  #   if(!i %in% names(dots)) {
  #     args[i]
  #   }
  # })
  # dots <- c(dots, unlist(dots_args, recursive = FALSE, use.names = TRUE))
  #-- update global
  for(name in names(args)) {
    assign(name, args[[name]])
  }
  #----------------------------------------------------------------------------#
  #-- Check: pre-compute data
  dds_rds <- file.path(outdir, "deseq_dds.rds")
  if(file.exists(dds_rds)) {
      dds <- readRDS(dds_rds)
  } else {
    #--------------------------------------------------------------------------#
    #-- check: aligner: salmon, STAR
    aligner <- list_hiseq_file(x, "aligner", "_rx")
    if(aligner == "salmon") {
      dds <- hiseq_prep_deseq4salmon(x, fix_batch)
    } else {
      dds <- hiseq_prep_deseq4fc(x, strandness, fix_batch)
    }
    if((check_path(outdir) & !file.exists(dds_rds)) | overwrite) {
      saveRDS(dds, dds_rds)
    }
  }
  #----------------------------------------------------------------------------#
  #-- run: deseq
  if(inherits(dds, "DESeqDataSet")) {
    genome <- list_hiseq_file(x, "genome", "rx") # add genome
    args$genome <- genome
    res <- deseq(dds, !!!args)
  } else {
    warning("`deseq()` failed")
    return(NULL)
  }
  #----------------------------------------------------------------------------#
  #-- run: save config
  # sample names, sanitized_str
  name_csv <- file.path(outdir, "smp_name.csv")
  name_rds <- file.path(outdir, "smp_name.rds")
  coldata  <- colData(dds)
  name_df  <- data.frame(
    smp_name  = coldata$smp_name,
    label     = rownames(coldata),
    condition = coldata$condition
  )
  # saveRDS(name_df, name_rds)
  write.csv(name_df, name_csv, row.names = FALSE)
  #-- run: transcripts_deseq2.csv
  norm_table <- file.path(outdir, "norm_table.csv")
  res_df  <- read.csv(norm_table)
  res_csv <- file.path(outdir, "transcripts_deseq2.csv")
  write.csv(res_df, res_csv, quote = TRUE, row.names = FALSE)
  #-- run: fix mean
  norm_fix_table <- file.path(outdir, "norm_table.fix.csv")
  res_fix_df <- read.csv(norm_fix_table)
  fix_csv    <- file.path(outdir, "transcripts_deseq2.fix.csv")
  write.csv(res_fix_df, fix_csv, quote = TRUE, row.names = FALSE)
  #----------------------------------------------------------------------------#
  # return
  res
}



