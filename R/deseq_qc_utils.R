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

#------------------------------------------------------------------------------#
# loading data for deseq_qc()

#' deseq_qc_dds
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
  #-- update, for child functions
  args <- purrr::list_modify(args, !!!dots)
  for(name in names(args)) {
    if(rlang::is_empty(name)) next
    assign(name, args[[name]])
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
  #----------------------------------------------------------------------------#
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
    res <- run_deseq_res(x, outdir = outdir, shrink = TRUE) #
  } else if(dir.exists(x)) {
    dds_rds <- file.path(x, "DESeq2_dds.rds")
    if(file.exists(dds_rds)) {
      dds <- readRDS(dds_rds)
      res <- run_deseq_res(dds, outdir = outdir, shrink = TRUE)
    }
  } else {
    warning(glue::glue(
      "x is '{class(x)}, expect 'DESeqDataSet', deseq_deseq2 dir"
    ))
  }
  #----------------------------------------------------------------------------#
  #-- output
  if(inherits(res, "list")) {
    res[[return_data]]
  }
}


#' deseq_qc_res
#'
#' @example
#'   ma, volcano: log2FolcChange, pvalue
#'   scatter: wt vs mut
#' @param x character path to the file `norm_table.fix.xls`
#'   also support data types: `DESeqDataSet`, ...
#' @param log2fc_limits numeric setting the range of log2FoldChange on plots
#'   default: NULL, by `scales::breaks_extend()`
#' @param fc numeric cutoff for foldchange, default: 1, ignore foldchange
#' @param pvalue numeric cutoff for padj, default: 0.1, the main criteria
#' @param p_adjust bool use p-adjust value instead
#' @param overwrite bool overwrite exists file, default: FALSE
#'
#'
#' @return ggplot
#'
#' @export
deseq_qc_res <- function(x, ...) {
  #----------------------------------------------------------------------------#
  dots <- rlang::list2(...)
  args <- rlang::list2(
    shrink_method    = "standard", # standard, normal, apeglm, ashr
    transform_method = "standard", # standard, vst, rlog
    log2fc_limits    = c(-2, 2),
    fc        = 2,
    pvalue    = 0.05,
    p_adjust  = TRUE,
    overwrite = FALSE,
    .col_sig  = "sig",
    readable  = TRUE,
    genome    = NULL
  )
  #-- update, for child functions
  args <- purrr::list_modify(args, !!!dots)
  #-- update global
  for(name in names(args)) {
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
      # pass
      log2fc_limits <- as.numeric(log2fc_limits)
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
  .col_pvalue <- ifelse(p_adjust, "padj", "pvalue") #!!!! pvalue
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
  # } else if(dir.exists(x)) {
  #   # 1. dds trans; log2, normalized
  #   dt_list <- deseq_qc_dds(x, return_data = "dds_trans") #log2
  #   dt <- dt_list[[transform_method]] # dds_trans
  #   if(is.null(dt)) {
  #     warning("unknown x, `deseq_qc_res()` failed; `dds_trans` not found")
  #     return(NULL)
  #   }
  #   # 2. results table, fc, pvalue
  #   res <- deseq_qc_dds(x, return_data = "res_lfc")
  #   dr <- res[[shrink_method]] # dds_res
  #   if(is.null(dr)) {
  #     warning("unknown x, `deseq_qc_res()` failed; `res_lfc` not found")
  #     return(NULL)
  #   }
  #   # 3. combine tables
  #   dt_count <- 2^assay(dt)
  #   dr <- as.data.frame(dr)
  #   df <- merge(dt_count, dr, by = "row.names" )
  #   colnames(df)[1] <- "gene_id"
  #   # 4. mean
  #   df <- deseq_mean(df) # add mean
  } else if(dir.exists(x) || is_hiseq_dir(x, hiseq_type = "deseq_deseq2")) {
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
    message(glue::glue(
      "ignore 'shrink_method={shrink_method}'"
    ))
    if(file.exists(file.path(x, "transcripts_deseq2.fix.csv"))) {
      df <- read.csv(file.path(x, "transcripts_deseq2.fix.csv"))
    } else if(file.exists(x) & endsWith(x, ".fix.csv")) {
      df <- read.csv(x) # for norm_table.fix.csv
    }
  } else {
    df <- NULL
  }
  #----------------------------------------------------------------------------#
  #-- Check: required columns
  if(inherits(df, "data.frame")) {
    rc <- c("baseMean", "log2FoldChange", .col_pvalue)
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
    df2 <- get_sig_name(df, return_dataframe = TRUE, !!!args)
  }
  #----------------------------------------------------------------------------#
  #-- run: readable
  if(readable) {
    if(is_hiseq_dir(x, "deseq_deseq2")) {
      gene_table <- list_hiseq_file(x, "gene_readable_csv", "deseq_deseq2")
      genome <- list_hiseq_file(x, "genome", "deseq_deseq2")
    } else {
      gene_table <- NULL
    }
    df <- set_readable(df, genome = genome, gene_table = gene_table)
  }
  #----------------------------------------------------------------------------#
  #-- add 'sig', force, re-run
  df <- get_sig_name(df, return_dataframe = TRUE, force = TRUE, !!!args)
  #----------------------------------------------------------------------------#
  #-- run: update log2fc, ext; for point-shapes
  # for shapes, up, dot, down
  # update log2fc, check outlier by log2fc_limits
  df %>%
    dplyr::mutate(
      log10pval     = -log10(!!sym(.col_pvalue)),
      log10basemean = log10(baseMean + 1),
      ext = ifelse(
        is.na(log2FoldChange), "dot", ifelse(
          log2FoldChange > max(breaks), "up", ifelse(
            log2FoldChange < min(breaks), "down", "dot")))) %>%
    dplyr::mutate(log2fc = ifelse(
      ext == "up", max(breaks), ifelse(
        ext == "down", min(breaks), log2FoldChange)))
}


