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
  #-- Check: args
  dots <- rlang::list2(...)
  dots <- prep_kegg(gene_list, organism, !!!dots) # update arguments
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
  #-- Check: valid args
  if(!is_valid_kegg_input(!!!dots)) {
    message("kegg() skipped, invalid arguments, check above message")
    return(NULL)
  }
  #-- outdir, updat
  outdir <- file.path(outdir, "kegg_enrich") # update: outdir
  check_path(outdir)
  #----------------------------------------------------------------------------#
  if(is_valid_kegg_input(!!!dots)) {
    kegg_enrich(gene_list, organism, !!!dots)
    kegg_gsea(gene_list, organism, !!!dots)
  } else {
    warning(glue::glue(
      "kegg() failed, check above message"
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
    gsea_gene    = NULL,
    kegg_code    = NULL,
    kegg_gene    = NULL,
    kegg_gsea_gene = NULL,
    kegg_keytype = NULL,
    pval_cutoff  = 0.9,
    qval_cutoff  = 0.9
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
    simplify = TRUE
  )
  # gsea_gene <- prep_gsea_input(
  #   gene_list   = gene_list,
  #   fold_change = fold_change,
  #   organism = organism
  # )
  kegg_keytype <- ifelse(
    organism == "Drosophila melanogaster",
    "ncbi-geneid", "kegg"
  )
  # kegg_gsea_gene <- prep_kegg_gsea_input(
  #   gene_list,
  #   fold_change,
  #   organism,
  #   keytype = keytype)
  #-- keytype
  kegg_keytype <- get_kegg_keytype(organism) # for OrgDb
  kegg_keyType <- get_kegg_keyType(organism) # for clusterProfiler
  #----------------------------------------------------------------------------#
  #-- update: dots, return
  purrr::list_modify(
    dots,
    outdir       = outdir,
    organism     = organism,
    orgdb        = orgdb,
    keytype      = keytype,
    kegg_code    = kegg_code,
    # kegg_gene    = kegg_gene,
    gsea_gene    = gsea_gene,
    kegg_keytype = kegg_keytype
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
  #-- Check: types
  f_gene_list <- inherits(gene_list, "character") & length(gene_list) > 0
  f_outdir    <- inherits(outdir, "character") & check_path(outdir)
  f_organism  <- is_valid_organism(organism)
  f_for_gsea  <- inherits(for_gsea, "logical")
  f_orgdb     <- inherits(orgdb, "OrgDb")
  f_keytype   <- is_valid_keytype(keytype, organism = organism)
  f_level     <- inherits(level, "numeric")
  f_fold_change <- is.null(fold_change) | (inherits(fold_change, "numeric") & inherits(names(fold_change), "character"))
  f_gsea_gene <- is.null(gsea_gene) | (inherits(gsea_gene, "numeric") & inherits(names(gsea_gene), "character"))
  f_kegg_code <- is.null(kegg_code) | inherits(kegg_code, "character")
  f_kegg_gene <- is.null(kegg_gene) | inherits(kegg_gene, "character")
  f_kegg_gsea_gene <- is.null(kegg_gsea_gene) | (inherits(kegg_gsea_gene, "numeric") & inherits(names(kegg_gsea_gene), "character"))
  f_kegg_keytype <- is.null(kegg_keytype) | inherits(kegg_keytype, "character")
  f_kegg_keyType <- is.null(kegg_keyType) | inherits(kegg_keyType, "character")
  f_overwrite <- inherits(overwrite, "logical")
  f_readable  <- inherits(overwrite, "logical")
  args_list <- c(
    "gene_list", "outdir", "organism", "for_gsea", "orgdb", "keytype",
    "fold_change", "level", "gsea_gene", "kegg_code", "kegg_gene",
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



#' #' create plots
#' #' @param gene_list object of enrichGO
#' #' @param ... passing extra arguments
#' #' parent function: get_go_plots(),
#' #'
#' #' fold_change
#' #' text_width
#' #'
#' #'
#' #' @export
#' enrich_kegg_plots <- function(x, ...) {
#'   #--Default values: END
#'   if(class(x) == "enrichResult") {
#'     if(nrow(x)) {
#'       list(barplot  = go_barplot(x, ...),
#'            dotplot  = go_dotplot(x, ...),
#'            cnetplot = go_cnetplot(x, ...),
#'            emapplot = go_emapplot(x, ...),
#'            # emapplot_cluster = go_emapplot_cluster(x, ...),
#'            heatplot = go_heatplot(x, ...))
#'     } else {
#'       warning("`x` is enrichGOResult, but contains 0 rows data")
#'       NULL
#'     }
#'   } else {
#'     warning("`x` not enrichGOResult")
#'     NULL
#'   }
#' }
#'
#'
#'
#'
#' #' create plots
#' #' @param x object of gseaGO
#' #'
#' #' @export
#' gsea_kegg_plots <- function(x, fold_change = NULL, ...) {
#'   if(class(x) == "gseaResult") {
#'     if(nrow(x)) {
#'       list(gseaplot = go_gsea_plot(x)) # go_gsea_plot
#'     } else {
#'       warning("`x` is gseaResult, but contains 0 rows data")
#'       NULL
#'     }
#'   } else {
#'     warning("`x` not gseaResult")
#'     NULL
#'   }
#' }








#--Deprecated functions --------------------------------------------------------

#' #' run enrichKEGG()
#' #' @param gene_list gene names
#' #' @param organism kegg_code, eg: dme
#' #'
#' #' @export
#' kegg_enrich <- function(gene_list, organism,
#'                         pval_cutoff = 0.05,
#'                         qval_cutoff = 0.05) {
#'   # choose keytype
#'   keytype <- ifelse(organism %in% c("dme"), "ncbi-geneid", "kegg")
#'
#'   # enrich
#'   print("Running enrichKEGG")
#'   kk <- clusterProfiler::enrichKEGG(gene          = gene_list,
#'                                     organism      = organism,
#'                                     keyType       = keytype,
#'                                     pval_cutoff  = pval_cutoff,
#'                                     pAdjustMethod = "BH",
#'                                     qval_cutoff  = qval_cutoff)
#'
#'   # readable
#'   # convert kegg_code to scientific_name
#'   sci_name <- clusterProfiler::search_kegg_organism(organism, by = "kegg_code")$scientific_name
#'   orgdb <- get_orgdb(sci_name)
#'
#'   if(is(kk, "enrichResult")) {
#'     kk <- clusterProfiler::setReadable(kk, orgdb, "ENTREZID")
#'   }
#'   kk
#' }
#'







#' #' save barplot of go
#' #'
#' #' @param x list of objects of groupGO, enrichGO, KEGG, ...
#' #' @param outdir path, to save the results
#' #'
#' #' @export
#' get_kegg_plots <- function(x, outdir, fold_change = NULL) {
#'   # x should be go obj
#'   if(is_go(x, recursive = TRUE)) {
#'     message("Generating plots for KEGG analysis")
#'   } else {
#'     warning("input failed, enrichResult, gseaResult expected")
#'     return(NULL)
#'   }
#'
#'   # prepare plot_list
#'   if(is_go(x)) {
#'     # single object
#'     plot_func <- switch(class(x),
#'                          "enrichResult"  = go_enrich_plots, # yes, it is correct
#'                          "gseaResult"    = go_gsea_plots) # yes, it is correct
#'     plot_func(x, fold_change)
#'   } else if(is.list(x)) {
#'     # list of obj
#'     p_list <- lapply(x, function(i){
#'       get_kegg_plots(i, outdir, fold_change)
#'     })
#'     names(p_list) <- names(x)
#'     p_list # return
#'   }
#' }




# run_kegg <- function(gene_list, organism, outdir, ...) {
#   # force input, ENTREZID #
#   input_gene  <- go_input(gene_list, organism, fold_change)
#   orgdb       <- input_gene$orgdb # global
#   gene_list   <- input_gene$gene
#   keytype     <- input_gene$keytype
#   fold_change <- input_gene$fold_change
#   gsea_gene   <- input_gene$gsea_gene
#   kegg_code   <- input_gene$kegg_code
#
#   ## run KEGG analysis
#   kegg_rds <- file.path(outdir, "kegg_data.rds")
#   if(file.exists(kegg_rds)) {
#     kegg_obj <- readRDS(kegg_rds)
#   } else {
#     kegg_obj <- list(enrich = list(
#       kegg = kegg_enrich(
#         gene_list, kegg_code,
#         pval_cutoff, qval_cutoff)),
#       gsea   = list(
#         kegg = kegg_gsea(
#           gene_list, kegg_code, fold_change)))
#     # save to obj
#     saveRDS(kegg_obj, file = kegg_rds)
#   }
#
#   ## Generate KEGG plots
#   ## 1. enrich
#   enrich_dir <- file.path(outdir, "kegg_enrich")
#   table1     <- save_table(kegg_obj$enrich, enrich_dir)
#   plot1      <- get_go_plots(kegg_obj$enrich, enrich_dir, fold_change)
#   save_plot(plot1, enrich_dir)
#
#   ## 2. GSEA
#   gsea_dir <- file.path(outdir, "kegg_gsea")
#   table2   <- save_table(kegg_obj$gsea, gsea_dir)
#   plot2    <- get_go_plots(kegg_obj$gsea, gsea_dir, fold_change)
#   save_plot(plot2, gsea_dir)
#
#   # save obj
#   plot_obj <- list(enrich = plot1,
#                    gsea   = plot2)
#   plot_rds <- file.path(outdir, "kegg_plots.rds")
#   saveRDS(plot_obj, file = plot_rds)
# }
#
#

#' #' run gseKEGG()
#' #' @param gene_list decreasing sorted numeric vector
#' #' @param organism kegg_code, eg: dme
#' #'
#' #' @export
#' gsea_kegg <- function(gene_list, organism, fold_change,
#'                       pval_cutoff = 0.05) {
#'   # input
#'   gene_list_fc <- gsea_input(gene_list, fold_change)
#'   if(is.null(gene_list_fc)) {
#'     return(NULL)
#'   }
#'
#'   # choose keytype
#'   keytype <- ifelse(organism %in% c("dme"), "ncbi-geneid", "kegg")
#'
#'   # GSEA
#'   print("Running gseKEGG") #!!!!
#'   gsea <- clusterProfiler::gseKEGG(geneList      = gene_list_fc,
#'                                    organism      = organism,
#'                                    keyType       = keytype,
#'                                    minGSSize     = 120,
#'                                    pval_cutoff  = pval_cutoff,
#'                                    pAdjustMethod = "BH",
#'                                    verbose       = FALSE)
#'
#'   # readable
#'   # convert kegg_code to scientific_name
#'   sci_name <- clusterProfiler::search_kegg_organism(organism, by = "kegg_code")$scientific_name
#'   orgdb <- get_orgdb(sci_name)
#'
#'   if(is(gsea, "gseaResult")) {
#'     gsea <- clusterProfiler::setReadable(gsea,
#'                                          OrgDb   = orgdb,
#'                                          keyType = "ENTREZID")
#'   }
#'   gsea
#' }

