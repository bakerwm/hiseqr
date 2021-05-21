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
#' @name enrich





#' @describeIn enrich_hub A port for enrich analysis
#'
#'
#' @param x path to the RNAseqRx directory: a.vs.b
#' @param ... extra arguments for demseq2_main function
#'
#' extra arguments:
#' - overwrite    : FALSE
#' - text_width   : 40
#' - pval_cutoff  : 0.9, conflict with deseq2_main()
#' - qval_cutoff  : 0.9
#' - level        : 2
#' - readable     : bool
#' - layout       : nicely, kk
#' - show_category: 12
#'
#'
#' @import dplyr
#' @import DESeq2
#'
#' @example
#' library(clusterProfiler)
#' library(DOSE)
#' data(geneList)
#' dots <- list(
#'     gene_list   = geneList[geneList > 1] %>% names,
#'     organism    = "hg38",
#'     outdir      = "demo",
#'     keytype     = "ENTREZID",
#'     fold_change = geneList
#'   )
#' do.call(enrich_hub, dots)
#'
#'
#'
#' @export
enrich_hub <- function(gene_list, organism, outdir = NULL, ...) {
  outdir   <- ifelse(is(outdir, "character"), outdir, getwd()) # working dir
  dots     <- rlang::list2(...)
  arg_vars <- list(
    orgdb         = NULL,
    keytype       = NULL,
    fold_change   = NULL,
    overwrite     = FALSE,
    level         = 2,
    readable      = TRUE,
    show_category = 12,
    text_width    = 40,
    pval_cutoff   = 0.9,
    qval_cutoff   = 0.9
  )
  arg_vars   <- purrr::list_modify(arg_vars, !!!dots)
  arg_vars   <- purrr::list_modify(arg_vars,
                                   gene_list = gene_list,
                                   organism  = organism,
                                   outdir    = outdir)
  #--Update args
  edata      <- do.call(prep_enrich, arg_vars)
  arg_vars   <- purrr::list_modify(arg_vars, !!!edata)
  if(is(edata, "list")) {
    #--GO: group, enrich, gsea analysis ...
    go_res   <- do.call(run_go, arg_vars)
    #--KEGG: enrich, gsea ...
    kegg_res <- do.call(run_kegg, arg_vars)
  }
}






#' @describeIn prep_enrich Prepare data for Enrich analysis
#'
#' @description
#'
#' @param x path to the directory of RNAseqRx
#'
#' @import readr
#' @import configr
#' @import dplyr
#'
#'
#'
#'
#'
#' @export
prep_enrich <- function(gene_list, organism, outdir = NULL, ...) {
  #--Check: args
  if(! is(gene_list, "character")) {
    warning(paste0("`gene_list` require characters, failed, [",
                   class(gene_list), "]"))
    return(NULL)
  }
  organism2 <- get_organism_name(organism)
  if(! is(organism2, "character")) {
    warning(paste0("`organism` invalid, [", organism, "]"))
    return(NULL)
  }
  organism <- organism2 # update name
  outdir <- ifelse(is(outdir, "character"), outdir, getwd()) # working dir
  #--Check: args, dots
  dots <- rlang::list2(...)
  arg_vars <- rlang::list2(
    fold_change = NULL,
    keytype     = NULL
  )
  arg_vars <- purrr::list_modify(arg_vars, !!!dots)
  arg_vars <- purrr::list_modify(arg_vars,
                                 gene_list = gene_list,
                                 organism  = organism,
                                 outdir    = outdir
  )
  #--Check: args, keytype
  if(! is(arg_vars$keytype, "character")) {
    arg_vars$keytype <- guess_keytype(gene_list, organism)
  }
  arg_vars$orgdb <- get_orgdb(arg_vars$organism)
  #--Output:
  # gsea_gene <- prep_gsea_input(gene_list, arg_vars$fold_change, organism)
  # kegg_code <- get_kegg_code(organism)
  # kegg_gene <- to_kegg_gene_id(gene_list, organism, arg_vars$keytype, simplify = TRUE)
  # # kegg_keytype <- ifelse(organism == "Drosophila melanogaster",
  # #                        "ncbi-geneid", "kegg")
  # kegg_keytype <- get_kegg_keytype(organism) # for OrgDb,
  # kegg_keyType <- get_kegg_keyType(organism) # for clusterProfiler
  list(
    gene_list    = gene_list,
    organism     = organism,
    outdir       = outdir,
    orgdb        = arg_vars$orgdb,
    keytype      = arg_vars$keytype
  )
}




