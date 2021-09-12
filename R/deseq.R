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
#' @param fix_batch bool fix batch effect, from replicates, suffix; default: TRUE
#' @param shrink logical `shrink` LFC by ["apeglm", "ashr", "normal"],
#'   default: TRUE
#' @param transform logical transform dds by `vst()`, `rlog()`, default: TRUE
#' @param overwrite bool overwrite exists file, default: FALSE
#' @param readable add `ENTREZID`, `SYMBOL`, based on gene_id, require: `genome`
#' @param genome character name of the organism, could be dm6, fruitfly
#' @param fc numeric cutoff for foldchange, default: 2, ignore foldchange
#' @param pvalue numeric cutoff for padj, default: 0.05, the main criteria
#' @param p_adjust bool use p-adjust value instead
#' @param cpu integer, number of CPU to run in parallel, default: 2
#'
#' @import DESeq2
#' @import apeglm
#' @import ggplot2
#'
#' @return DESeqResults, res (shrinked)
#'
#' @export
deseq <- function(dds, ...) {
  message(">>> run 'deseq()' ...")
  #----------------------------------------------------------------------------#
  #-- Check: default values
  dots <- rlang::list2(...)
  args <- rlang::list2(
    outdir     = NULL,
    cpu        = 2,
    n_max      = 20,
    fc         = 2,
    pvalue     = 0.05,
    p_adjust   = TRUE,
    label_list = NULL,
    label_max  = 8,
    overwrite  = FALSE,
    readable   = TRUE,
    shrink     = TRUE,
    transform  = TRUE,
    fix_batch  = TRUE,
    density_points = FALSE
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
  #-- Check: arguments
  if(!inherits(outdir, "character")) {
    outdir <- tempdir()
  }
  check_path(outdir) # create
  outdir <- normalizePath(outdir)
  #-- config: for hiseqr package
  config_list <- rlang::list2(
    hiseq_type        = 'deseq_deseq2',
    outdir            = outdir,
    fix_batch         = fix_batch,
    shrink            = shrink,
    genome            = genome,
    cpu               = cpu,
    fc                = fc,
    pvalue            = pvalue,
    p_adjust          = p_adjust,
    config_yaml       = file.path(outdir, "config.yaml"),
    deseq_dds_rds     = file.path(outdir, "deseq_dds.rds"),
    deseq_res_rds     = file.path(outdir, "deseq_res.rds"),
    deseq_qc_dds_rds  = file.path(outdir, "deseq_qc_dds.rds"), # dds?
    deseq_qc_res_rds  = file.path(outdir, "deseq_qc_res.rds"), # res
    fpkm_csv          = file.path(outdir, "fpkm_table.csv"),
    norm_csv          = file.path(outdir, "norm_table.csv"),
    norm_fix_csv      = file.path(outdir, "norm_table.fix.csv"),
    smp_name_csv      = file.path(outdir, "smp_name.csv"),
    gene_readable_csv = file.path(outdir, "gene_readable.csv")
  )
  #----------------------------------------------------------------------------#
  #-- run: main
  saveRDS(dds, config_list$deseq_dds_rds) # deseq_dds.rds
  dd <- run_deseq_res(dds, !!!dots) # deseq_res.rds, dds, res, res_lfc
  if(!inherits(dd, "list")) {
    warning("`run_deseq_des()` failed, see above messages")
    return(NULL)
  }
  #----------------------------------------------------------------------------#
  #-- run: save config, and qc
  if(check_path(outdir)) {
    #--------------------------------------------------------------------------#
    #-- run: save config
    yaml::write_yaml(config_list, config_list$config_yaml)
    #--------------------------------------------------------------------------#
    #-- run: gene_readable , SYMBOL, ENTREZ
    if(inherits(genome, "character")
       & isTRUE(readable)
       & !file.exists(config_list$norm_fix_csv)
    ) {
      df4 <- data.frame(
        gene_id = rownames(dd$res),
        stringsAsFactors = FALSE
      )
      df4 <- set_readable(df4,
                          genome = genome,
                          gene_table = config_list$gene_readable_csv)
      #------------------------------------------------------------------------#
      # update: TE, piRC list; unique;
      # most freq prefix for gene_id
      prefix_table <- table(substr(df4[["gene_id"]], 1, 4))
      if(max(prefix_table) / sum(prefix_table) > 0.9) {
        prefix_top <- names(prefix_table[prefix_table == max(prefix_table)])
        # update SYMBOL for TE,piRC
        df4 <- df4 %>%
          dplyr::mutate(SYMBOL = ifelse(startsWith(gene_id, prefix_top), SYMBOL, gene_id))
      }
      write.csv(df4, config_list$gene_readable_csv, row.names = FALSE,
                quote = TRUE)
    }
    #--------------------------------------------------------------------------#
    #-- run: saving norm table
    df1a <- DESeq2::counts(dds, normalized = TRUE) # normalized count
    df1  <- merge(as.data.frame(df1a), as.data.frame(dd$res), by = "row.names")
    colnames(df1)[1] <- "gene_id"
    write.csv(df1, config_list$norm_csv, quote = TRUE, row.names = FALSE)
    #--------------------------------------------------------------------------#
    #-- run: saving norm table, fix
    df2 <- deseq_mean(dds, outdir = outdir)
    if(inherits(genome, "character")
       & isTRUE(readable)
    ) {
      message("aaa")
      df2 <- set_readable(df2,
                          genome = genome,
                          gene_table = config_list$gene_readable_csv)
    }
    write.csv(df2, config_list$norm_fix_csv, quote = TRUE, row.names = FALSE)
    #--------------------------------------------------------------------------#
    #-- run: fpkm
    if("basepairs" %in% names(mcols(dds))) {
      df3a <- DESeq2::fpkm(dds)
      df3 <- merge(as.data.frame(df3a), as.data.frame(dd$res), by = "row.names")
      if(inherits(genome, "character")
         & is_valid_organism(genome)
         & isTRUE(readable)
         & !file.exists(config_list$fpkm_csv)
      ) {
        df3 <- set_readable(df3,
                            genome = genome,
                            gene_table = config_list$gene_readable_csv)
      }
      write.csv(df3, config_list$fpkm_csv, quote = TRUE, row.names = FALSE)
    }
    #--------------------------------------------------------------------------#
    #-- run: quality-control, require outdir
    tmp <- deseq_qc(outdir, !!!dots)
  }
  #----------------------------------------------------------------------------#
  dd$res # original res
}


#' @describeIn run_deseq_res
#' run DESeq2::results() for dds
#' vst(), rlog() for dds
#' results() for res
#' lfcShrink() for res
#'
#' @param dds DESeqDataSet
#' @param outdir character saving deseq_res to file: deseq_res.rds
#' @param shrink logical `shrink` LFC by ["apeglm", "ashr", "normal"],
#'   default: TRUE
#' @param transform logical transform dds by `vst()`, `rlog()`, default: TRUE
#' @param cpu integer, number of CPU to run in parallel, default: 2
#' @param overwrite bool overwrite exists file, default: FALSE
#'
#' @return list(dds=, dds_trans = list(), res=, res_lfc=list())
#'
#' @export
run_deseq_res <- function(dds, ...) {
  #----------------------------------------------------------------------------#
  #-- Check: default values
  dots <- rlang::list2(...)
  args <- rlang::list2(
    outdir     = NULL,
    cpu        = 2,
    overwrite  = FALSE,
    readable   = TRUE,
    shrink     = TRUE,
    transform  = TRUE,
    fix_batch  = TRUE
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
  #-- dds
  if(!inherits(dds, "DESeqDataSet")) {
    warning(glue::glue("`dds` is {class(dds)}, expect `DESeqDataSet`"))
    return(NULL)
  }
  #-- cpu
  if(inherits(cpu, "numeric")) {
    if(cpu < 0 | cpu > 8) {
      message(glue::glue("illegal 'cpu' = {cpu}, use 2 instead"))
      cpu <- 2
    }
  } else {
    message(glue::glue("illegal 'cpu' = {cpu}, use 2 instead"))
    cpu <- 2
  }
  cpu <- round(cpu)
  #-- shrink
  if(!inherits(shrink, "logical")) {
    shrink <- TRUE
  }
  #-- transform
  if(!inherits(transform, "logical")) {
    transform <- TRUE
  }
  #-- overwrite
  if(!inherits(overwrite, "logical")) {
    overwrite <- FALSE
  }
  #-- Check: config
  deseq_res_rds <- list_hiseq_file(outdir, "deseq_res_rds", "deseq_deseq2")
  if(inherits(deseq_res_rds, "character")) {
    if(file.exists(deseq_res_rds) & !overwrite) {
      message(glue::glue("loading deseq_res data from: {deseq_res_rds}"))
      return(readRDS(deseq_res_rds))
    }
  }
  #----------------------------------------------------------------------------#
  #-- run: main
  coldata <- colData(dds)
  if(!rlang::has_name(coldata, "condition")) {
    warning("required column `condition` not found, check `colData(dds)`")
    return(NULL)
  }
  #-- run: DESeq analysis
  BiocParallel::register(BiocParallel::MulticoreParam(cpu))
  #----------------------------------------------------------------------------#
  #-- run: `vst()`, `vlog()`, norm, `DESeqTransform`
  if(isTRUE(transform)) {
    if(nrow(dds) > 1000) {
      # why > 1000 genes? see https://support.bioconductor.org/p/98634/#98637
      dds_trans <- list(
        standard = DESeq2::normTransform(dds),
        vst      = DESeq2::vst(dds, blind = FALSE),
        rlog     = DESeq2::rlog(dds, blind = FALSE)
      )
    } else {
      dds_trans <- list(standard = DESeq2::normTransform(dds))#
    }
  } else {
    dds_trans <- list(standard = DESeq2::normTransform(dds))#
  }
  # dds <- DESeq2::DESeq(dds)
  wt   <- levels(coldata$condition)[1] #
  mut  <- levels(coldata$condition)[2] #
  coef <- deseq_sanitize_str(paste0("condition_", mut, "_vs_", wt))
  res  <- DESeq2::results(
    dds,
    contrast = c("condition", mut, wt), # b vs a
    parallel = TRUE
  )
  res <- res[order(res$padj), ] # Order by adjusted p-value
  #----------------------------------------------------------------------------#
  #-- run: shrink, option
  if(shrink) {
    res_lfc <- sapply(c("normal", "apeglm", "ashr"), function(s) {
      message(glue::glue("using '{s}` for LFC shrinkage"))
      DESeq2::lfcShrink(
        dds,
        coef     = coef,
        type     = s,
        parallel = TRUE
      )
    })
  } else {
    res_lfc <- list()
  }
  # update res_lfc
  res_lfc[["standard"]] <- res
  out <- list(
    dds       = dds,
    dds_trans = dds_trans,
    res       = res,
    res_lfc   = res_lfc
  )
  #-- output file
  if(!inherits(deseq_res_rds, "character")) {
    deseq_res_rds <- file.path(outdir, "deseq_res.rds")
  }
  #----------------------------------------------------------------------------#
  #-- save
  if(inherits(outdir, "character")) {
    if(check_path(outdir)) {
      if(!file.exists(deseq_res_rds) | overwrite) {
        message(glue::glue("saving `deseq_res` to file: {deseq_res_rds}"))
        saveRDS(out, deseq_res_rds)
      }
    }
  }
  #----------------------------------------------------------------------------#
  out
}

