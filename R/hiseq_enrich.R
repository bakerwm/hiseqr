#' Functions for enrich analysis
#'
#' Path/DB: GO, KEGG (DO, MSigDb, Reactome, MeSH)
#' - enrich: over-representation analysis
#' - GSEA: Gene Set Enrichment analysis
#'
#' including:
#' - enrich_go
#' - enrich_kegg
#' - gsea_go
#' - gsea_kegg
#' - group_go
#' ...
#'
#' Adapt the functions from clusterProfiler, enrichplot, ...
#'
#' `enrichGO()`, `enrichKEGG()`, `gseGO`, `gseKEGG`, `groupGO`, ...
#'
#' @name hiseq_enrich


#' @describeIn hiseq_enrich
#' enrich analysis for rnaseq_rx directory
#'
#' @param x character path to the rnaseq_rx directory
#'
#' @example
#' hiseq_enrich(x)
#' hiseq_enrich(x, sig_type = "up", outdir = "up")
#'
#' @return results
#'
#' @export
hiseq_enrich<- function(x, ...) {
  message(">>> run hiseq_enrich()")
  #----------------------------------------------------------------------------#
  #-- Check: arguments
  #-- scan shrink method, sig type (sig, up, down)
  # for(shrink in c("ashr", "standard", "apeglm", "normal")) {
  for(shrink in c("apeglm", "standard")) {
    for(sig in c("sig", "up", "down")) {
      message(glue::glue(
        ">>> run hiseq_enrich() for '{shrink}'-'{sig}'"
      ))
      dots <- rlang::list2(...) # init
      args <- hiseq_prep_enrich(
        x,
        !!!dots,
        shrink_method = shrink,
        sig_type = sig
      ) # valid args
      if(is.null(args)) {
        message("hiseq_enrich() skipped for '{shrink}'-'{sig}'")
        next
      }
      #-- to env
      for(name in names(args)) {
        assign(name, args[[name]])
      }
      #------------------------------------------------------------------------#
      #-- update subdir
      args$outdir <- file.path(args$outdir, shrink, sig) # update directory
      check_path(args$outdir)
      #-- save arguments to file
      args_rds <- file.path(args$outdir, "args.rds")
      saveRDS(dots, args_rds)
      #-- run
      if(inherits(gene_list, "character")) {
        go(gene_list, organism, !!!args)
        kegg(gene_list, organism, !!!args)
      }
    }
  }
}


#' @describeIn prep_enrich Prepare data for Enrich analysis
#'
#' @description
#'
#' @param x path to the directory of rnaseq_rx
#'
#' @import readr
#' @import configr
#' @import dplyr
#'
#' @export
hiseq_prep_enrich <- function(x, ...) {
  #----------------------------------------------------------------------------#
  #-- Check: arguments
  dots <- rlang::list2(...)
  args <- rlang::list2(
    outdir        = NULL,
    sig_type      = "sig", # up, down, not, sig(up+down), all(up+down+not)
    shrink_method = "standard", # standard, apeglm, ashr, normal
    fc            = 2,
    pvalue        = 0.05,
    p_adjust      = TRUE,
    genome        = NULL,
    organism      = NULL,
    fold_change   = NULL,
    keytype       = NULL,
    orgdb         = NULL,
    show_category = 12,
    level         = 2,   # for GO group, level
    text_width    = 40,
    pval_cutoff   = 0.05, # try to return enrich results for all
    qval_cutoff   = 0.05, # see pval_cutoff
    readable      = TRUE # for enrich, readable
  )
  dots <- purrr::list_modify(args, !!!dots) # from arguments
  for(name in names(dots)) {
    assign(name, dots[[name]])
  }
  #----------------------------------------------------------------------------#
  #-- Check: args
  if(!is_hiseq_dir(x, "rnaseq_rx")) {
    warning(glue::glue("x is {x}, expect 'rnaseq_rx'"))
    return(NULL)
  }
  #-- Check: deseq data
  deseq_dir <- list_hiseq_file(x, "deseq_dir", "rx")
  if(!is_hiseq_dir(deseq_dir, "deseq_deseq2")) {
    warning(glue::glue(
      "'deseq' data not found, please run 'hiseq_deseq(x, ...)' again"
    ))
    return(NULL)
  }
  #-- organism
  genome   <- list_hiseq_file(x, "genome", "rx")
  organism <- get_organism_name(genome)
  #-- outdir
  if(!inherits(outdir, "character")) {
    enrich_dir <- list_hiseq_file(x, "enrich_dir", "rx")
    outdir <- ifelse(inherits(enrich_dir, "character"), enrich_dir, getwd())
  }
  check_path(outdir)
  #-- orgdb
  orgdb <- get_orgdb(genome)
  #----------------------------------------------------------------------------#
  #-- for: gene_list
  #-- pvalue, padjust
  if(inherits(fc, "numeric") & inherits(pvalue, "numeric")) {
    if(fc <= 0 | pvalue <= 0 | pvalue > 1) {
      warning(glue::glue(
        "'fc' is {fc}, expect (0-Inf), greater than 0; ",
        "'pvalue' is {pvalue}`, expect (0-1)"
      ))
      fc <- 2
      pvalue <- 0.05
    }
  } else {
    warning(glue::glue(
      "`fc` is {class(fc)}, expect `numeric`, ",
      "`pvalue` is {class(pvalue)}, expect `numeric`"
    ))
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
  #----------------------------------------------------------------------------#
  #-- for: gene_list
  #-- update dots, for deseq.res table
  dots <- purrr::list_modify(
    dots,
    outdir   = outdir,
    genome   = genome,
    organism = organism,
    orgdb    = orgdb,
    fc       = fc,
    pvalue   = pvalue,
    p_adjust = p_adjust
  )
  df <- deseq_qc_res(deseq_dir, !!!dots)
  if(inherits(df, "data.frame")) {
    #-- save data.frame to file
    gene_table <- file.path(outdir, "gene_table.csv")
    write.csv(df, gene_table, row.names = FALSE, quote = TRUE)
    #--------------------------------------------------------------------------#
    #-- filter by criterias
    df2 <- filt_sig_gene(df, type = dots$sig_type)
    if(nrow(df2) == 0) {
      message(glue::glue(
        "no '{sig_type}' for shrink:'{shrink_method}'"
      ))
      return(NULL)
    }
    gene_list <- as.character(df2[["gene_id"]])
    fold_change <- setNames(
      df2[["log2FoldChange"]], nm = df2[["gene_id"]]
    )
    fold_change <- sort(fold_change, decreasing = TRUE)
  }
  #-- guess: keytype
  if(!inherits(keytype, "character")) {
    keytype <- guess_keytype(gene_list, genome)
  }
  #----------------------------------------------------------------------------#
  #-- update dots, again
  purrr::list_modify(
    dots,
    gene_list   = gene_list,
    fold_change = fold_change,
    keytype     = keytype
  )
}


