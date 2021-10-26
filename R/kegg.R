#' Functions for KEGG analysis
#'
#' overrepresent analysis
#' GSEA analysis
#'
#'
#' enrich kegg
#' 1. csv files, results
#' 2. barplot
#' 3. wego plot,
#' 4. dotplot
#' 5. netplot
#' 6. ...
#'
#' gsea kegg
#' 1. csv files, results
#' 2. barplot
#' 3. dotplot
#' 4. netplot
#' 5. ...
#'
#' @name kegg



#' @describeIn run_kegg Main port for KEGG analysis
#'
#' enrich: over-representation analysis
#' gsea: gsea
#'
#'
#' organize KEGG output
#' generate tables, plots
#' @param gene_list vector of gene list
#' @param organism name of genome, eg: Homo sapiens
#' @param outdir string
#' @param ... extra argument
#'
#' fold_change, pval_cutoff, qval_cutoff
#'
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
#' do.call(run_kegg, dots)
#'
#'
#' @export
kegg <- function(gene_list, organism, ...) {
  #----------------------------------------------------------------------------#
  dots <- rlang::list2(...)
  args <- prep_kegg(gene_list, organism, !!!dots) # update arguments
  #-- to env
  for(name in names(args)) {
    assign(name, args[[name]])
  }
  #----------------------------------------------------------------------------#
  if(!is_valid_kegg_input(!!!args)) {
    message("kegg() skipped, invalid arguments, check above message")
    return(NULL)
  }
  #-- outdir, updat
  outdir <- file.path(outdir, "kegg_enrich") # update: outdir
  check_path(outdir)
  #----------------------------------------------------------------------------#
  if(is_valid_kegg_input(!!!args)) {
    tryCatch(
      {
        kegg_enrich(gene_list, organism, !!!args)
        kegg_gsea(gene_list, organism, !!!args)
      },
      error = function(cond) {
        warning(">>> kegg() failed")
      }
    )
  } else {
    warning(glue::glue(
      ">>> kegg() failed, args not valid"
    ))
  }
}


#' @describeIn prep_kegg
#' Prepare data for kegg analysis
#'
#' 1. Guess the keytype of the gene
#' 2. Convert gene to entrezid
#' 3. add fold_change to gene
#' 4. extract orgdb
#'
#' !important: skip orgdb in arg_vars
#' use saveDb() and loadDb() to save/read OrgDb from file
#' Because the object is a reference to a sqlite data base.
#'
#' genes, orgsnism, orgdb, keytype, gsea_gene
#'
#' genes with fold_change, or other ranking
#' @param gene_list gene names
#' @param organism string, eg: Homo sapiens
#' @param ... extra argument
#'
#' @param fold_change numeric, fold change
#' @param keytype character
#'
#' @export
prep_kegg <- function(gene_list, organism, ...) {
  #----------------------------------------------------------------------------#
  #-- Check: args
  dots <- rlang::list2(...)
  args <- rlang::list2(
    outdir       = NULL,
    orgdb        = NULL,
    keytype      = NULL,
    fold_change  = NULL,
    overwrite    = FALSE,
    level        = 2,
    readable     = TRUE,
    # gsea_gene    = NULL,
    kegg_code    = NULL,
    kegg_gene    = NULL,
    kegg_gsea_gene = NULL,
    kegg_keytype = NULL,
    kegg_keyType = NULL,
    pval_cutoff  = 0.05,
    qval_cutoff  = 0.05
  )
  dots <- purrr::list_modify(args, !!!dots) # args, overwrite by dots
  #-- force, update positional args
  dots <- purrr::list_modify(
    dots,
    gene_list = gene_list,
    organism  = organism
  )
  #-- to env
  for(name in names(dots)) {
    assign(name, dots[[name]])
  }
  #----------------------------------------------------------------------------#
  #-- outdir
  if(!inherits(outdir, "character")) {
    outdir <- getwd() # current directory
  }
  #-- organism, "dm6" -> "Drosophila melanogaster"
  organism <- get_organism_name(organism)
  if(!inherits(orgdb, "OrgDb")) {
    orgdb <- get_orgdb(organism)
  }
  if(!is_valid_keytype(keytype, organism = organism)) {
    keytype <- guess_keytype(gene_list, organism)
  }
  kegg_code <- get_kegg_code(organism)
  kegg_gene <- to_kegg_gene_id(
    gene_list = gene_list,
    organism  = organism,
    keytype   = keytype,
    simplified  = TRUE
  )
  #-- keytypes for kegg
  kegg_keytype <- ifelse(
    organism == "Drosophila melanogaster",
    "ncbi-geneid", "kegg"
  )
  kegg_keytype <- get_kegg_keytype(organism) # for OrgDb
  kegg_keyType <- get_kegg_keyType(organism) # for clusterProfiler
  #----------------------------------------------------------------------------#
  #-- update: dots, return
  purrr::list_modify(
    dots,
    # gsea_gene    = gsea_gene,
    outdir       = outdir,
    organism     = organism,
    orgdb        = orgdb,
    keytype      = keytype,
    kegg_code    = kegg_code,
    kegg_gene    = kegg_gene,
    kegg_gsea_gene = kegg_gsea_gene,
    kegg_keytype = kegg_keytype,
    kegg_keyType = kegg_keyType
  )
}


#' @describeIn check_go_input Check the input for GO analysis
#'
#' Required:
#' - gene_list, chr, characters (group, enrich)
#' - outdir, chr, character
#' - organism, chr, character
#' - orgdb, OrgDb, AnnotationDbi, one of organism, orgdb required
#'
#' optional:
#' - for_gsea, logical, gene_list, num, sorted, named (GSEA)
#'
#' @export
is_valid_kegg_input <- function(...) {
  #----------------------------------------------------------------------------#
  #-- Check: args
  dots <- rlang::list2(...)
  args <- rlang::list2(
    gene_list    = NULL,
    organism     = NULL,
    outdir       = NULL,
    for_gsea     = FALSE,
    orgdb        = NULL,
    keytype      = NULL,
    fold_change  = NULL,
    level        = 2,
    gsea_gene    = NULL,
    kegg_code    = NULL,
    kegg_gene    = NULL,
    kegg_gsea_gene = NULL,
    kegg_keytype = NULL,
    overwrite    = FALSE,
    readable     = TRUE
  )
  dots <- purrr::list_modify(args, !!!dots) # args, overwrite by dots
  #-- to env
  for(name in names(dots)) {
    assign(name, dots[[name]])
  }
  #----------------------------------------------------------------------------#
  #  d_str <- paste(as.character(names(dots)), collapse = ",")
  # message(glue::glue("<<< dots for kegg_gsea: kegg_gene: {dots$kegg_gene};  {d_str}"))
  #-- Check: types
  f_gene_list <- inherits(gene_list, "character") & length(gene_list) > 0
  f_outdir    <- inherits(outdir, "character") & check_path(outdir)
  f_organism  <- is_valid_organism(organism)
  f_for_gsea  <- inherits(for_gsea, "logical")
  f_orgdb     <- inherits(orgdb, "OrgDb")
  f_keytype   <- is_valid_keytype(keytype, organism = organism)
  f_level     <- inherits(level, "numeric")
  f_fold_change <- is.null(fold_change) | (inherits(fold_change, "numeric") & inherits(names(fold_change), "character"))
  # f_gsea_gene <- is.null(gsea_gene) | (inherits(gsea_gene, "numeric") & inherits(names(gsea_gene), "character"))
  f_kegg_code <- is.null(kegg_code) | inherits(kegg_code, "character")
  f_kegg_gene <- is.null(kegg_gene) | inherits(kegg_gene, "character")
  f_kegg_gsea_gene <- is.null(kegg_gsea_gene) | (inherits(kegg_gsea_gene, "numeric") & inherits(names(kegg_gsea_gene), "character"))
  f_kegg_keytype <- is.null(kegg_keytype) | inherits(kegg_keytype, "character")
  f_kegg_keyType <- is.null(kegg_keyType) | inherits(kegg_keyType, "character")
  f_overwrite <- inherits(overwrite, "logical")
  f_readable  <- inherits(overwrite, "logical")
  # "gsea_gene"
  args_list <- c(
    "gene_list", "outdir", "organism", "for_gsea", "orgdb", "keytype",
    "fold_change", "level",  "kegg_code", "kegg_gene",
    "kegg_gsea_gene", "kegg_keytype", "kegg_keyType", "overwrite", "readable"
  )
  f_class <- unlist(lapply(args_list, function(i) class(get(i))))
  f_res   <- unlist(lapply(args_list, function(i) get(paste0("f_", i))))
  msg <- paste(
    paste0(f_res, ": '", args_list, "' is ", f_class),
    collapse = "\n"
  )
  if(!all(f_res)) {
    warning(glue::glue(msg))
  }
  all(f_res)
}

