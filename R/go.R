#' Functions for GO analysis
#'
#' organize GO output
#' generate tables, plots
#'
#' go_group
#' 1. csv files, results
#' 2. barplot,
#' 3. wego plot,
#'
#' go_enrich
#' 1. csv files, results
#' 2. barplot
#' 3. wego plot,
#' 4. dotplot
#' 5. netplot
#' 6. ...
#'
#' @name go



#' @describeIn go
#' Main GO analysis function for gene list
#'
#' @param gene_list character Genes for GO analysis
#' @param organism character Name of the organism, eg: dm3, fruitfly
#' @param outdir string Path to the directory, saving go_group results
#' ..., pass the following arguments from parent function
#'
#' or else, set default values:
#' orgdb    = NULL
#' keytype  = NULL
#' fold_change = NULL
#' level    = 2
#' readable = TRUE
#' overwrite = FALSE
#'
#' @param orgdb OrgDb from AnnotationDbi, overwrite `organism`
#' @param keytype character Name of the keytype for input, default: NULL
#' @param level int GO levels, default: 2
#' @param readable bool Args for groupGO function, convert to gene symbol
#'
#' ... pass
#' text_width
#' organism, and so on, to children functions
#'
#' go_group(), keytype, level, readable
#' go_enrich(), keytype, readable, pval_cutoff, qval_cutoff
#' go_gsea(), keytype, pval_cutoff
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
#' do.call(run_go, dots)
#'
#' @export
go <- function(gene_list, organism, ...) {
  message("run go() ...")
  #----------------------------------------------------------------------------#
  #-- Check: args
  dots <- rlang::list2(...)
  dots <- prep_go(gene_list, organism, !!!dots) # update dots
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
  #-- init args
  if(is_valid_go_input(!!!dots)) {
    go_group(gene_list, organism, !!!dots)
    go_enrich(gene_list, organism, !!!dots)
    go_gsea(go_gsea, !!!dots)
  }
}


#' @describeIn prep_go
#' Prepare data for GO analysis
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
prep_go <- function(gene_list, organism, ...) {
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
    kegg_keytype = NULL
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
  #-- fixed name
  kegg_keytype <- ifelse(
    organism == "Drosophila melanogaster",
    "ncbi-geneid", "kegg"
  )
  #----------------------------------------------------------------------------#
  #-- update: dots, return
  purrr::list_modify(
    dots,
    outdir       = outdir,
    organism     = organism,
    orgdb        = orgdb,
    keytype      = keytype,
    kegg_code    = kegg_code,
    kegg_gene    = kegg_gene,
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
is_valid_go_input <- function(...) {
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
  f_kegg_keytype <- is.null(kegg_keytype) | inherits(kegg_keytype, "character")
  f_overwrite <- inherits(overwrite, "logical")
  f_readable  <- inherits(overwrite, "logical")
  args_list <- c(
    "gene_list", "outdir", "organism", "for_gsea", "orgdb", "keytype",
    "fold_change", "level", "gsea_gene", "kegg_code", "kegg_gene",
    "kegg_keytype", "overwrite", "readable")
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


#' save obj to outdir
#'
#' @param x list of go results
#' @param outdir string, path to save the plots
#' @param name prefix for the plots
#'
#' @return
save_go_table <- function(x, outdir, name = NULL) {
  if(! dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE, mode = "0755")
  }
  if(is_go_result(x)) {
    if(is(name, "character")) {
      fname <- paste0(name[1], ".csv")
    } else {
      fname <- basename(tempfile("go_data.", fileext = ".csv"))
    }
    table_csv <- file.path(outdir, fname)
    message(paste0("Saving GO results to csv: ", fname))
    if(file.exists(table_csv)) {
      message(paste0("file exists, skipped: ", table_csv))
    } else {
      write.csv(x@result, table_csv, quote = TRUE, row.names = FALSE)
    }
    table_csv
  } else if(is(x, "list") & all(purrr::map_lgl(x, is_go_result))) {
    lapply(names(x), function(i) {
      x_i <- x[[i]]
      prefix <- paste(c(name, i), collapse = ".")
      save_go_table(x_i, outdir, prefix)
    })
  } else {
    warning("`x` expect GO Results, or list of GO Results, failed")
    NULL
  }
}


#' save plots to outdir
#'
#' @param x list of go plots
#' @param outdir string, path to save the plots
#' @param name prefix for the plots
#'
#' @import ggplot2
#'
#' @return
save_go_plot <- function(x, outdir, name = NULL) {
  if(! dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE, mode = "0755")
  }
  if(is(x, "gg")) {
    if(is(name, "character")) {
      fname <- paste0(name[1], ".png")
    } else {
      fname <- basename(tempfile("go_plot.", fileext = ".png"))
    }
    plot_file <- file.path(outdir, fname)
    message(glue::glue("Saving GO plot to file: {fname}"))
    if(file.exists(plot_file)) {
      message(paste0("file exists, skipped: ", basename(plot_file)))
    } else {
      tryCatch(
        {
          ggplot2::ggsave(plot_file, plot = x, width = 8, height = 8, dpi = 200)
        },
        error=function(cond) {
          warning(glue::glue("ggsave() failed, {class(x)}"))
          return(NULL)
        }
      )
    }
    plot_file
  } else if(is(x, "list")) {
    lapply(names(x), function(i) {
      x_i <- x[[i]]
      prefix <- paste(c(name, i), collapse = ".")
      save_go_plot(x_i, outdir, prefix)
    })
  } else {
    warning("`x` expect ggplot2, or list of ggplot2, failed")
    NULL
  }
}



