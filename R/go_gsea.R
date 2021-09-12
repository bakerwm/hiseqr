#' Functions for GO analysis
#'
#' organize GO output
#' generate tables, plots
#'
#' group_go
#' 1. csv files, results
#' 2. barplot,
#' 3. wego plot,
#'
#' enrich_go
#' 1. csv files, results
#' 2. barplot
#' 3. wego plot,
#' 4. dotplot
#' 5. netplot
#' 6. ...
#'
#' @name go_gsea



#' @describeIn go_gsea
#' gseGO analysis from clusterProfiler
#'
#' for clusterProfiler::groupGO, and enrichGO, gseGO
#'
#' @param gene_list numeric
#' @param organism character Name of the organism, eg: dm3, fruitfly
#' @param fold_change numeric, with names
#' @param ... pass arguments: , orgdb, keytype, pval_cutoff, qval_cutoff,
#' readable, with default values:
#' orgdb = NULL
#' keytype = NULL
#' pval_cutoff = 0.9
#'
#' @param orgdb OrgDb from AnnotationDbi, overwrite `organism`
#' @param keytype character Name of the keytype for input, default: NULL
#' @param pval_cutoff float Cutoff for p-value, default: 0.9
#'
#'
#' @import clusterProfiler
#' @example gseGO(gene =, OrgDb = , ont = )
#'
#' @export
go_gsea <- function(gene_list, organism, ...) {
  message(">>> run go_gsea()")
  #----------------------------------------------------------------------------#
  #-- Check: args
  dots <- rlang::list2(...)
  # dots <- purrr::list_modify(
  #   dots,
  #   gene_list = gene_list,
  #   organism  = organism
  # )
  dots <- prep_go_gsea(gene_list, organism, !!!dots) # update arguments
  #-- to env
  for(name in names(dots)) {
    assign(name, dots[[name]])
  }
  #----------------------------------------------------------------------------#
  #-- init args
  #-- Check: valid args
  if(is.null(dots) || !is_valid_go_input(!!!dots)) {
    message("go_gsea() skipped, invalid arguments, check above message")
    return(NULL)
  }
  if(!inherits(gsea_gene, "numeric")
     | !inherits(names(gsea_gene), "character")) {
    message("gesa_go() skipped, invalied 'gsea_gene'")
    return(NULL)
  }
  #-- outdir, updat
  outdir <- file.path(outdir, "go_gsea") # update: outdir
  check_path(outdir)
  #----------------------------------------------------------------------------#
  #-- Ontology:
  onts <- c("BP", "CC", "MF")
  #-- run
  go_data <- sapply(onts, function(ont) {
    message(glue::glue(">>> run gseaGO() for {ont}"))
    go_ont_data_rds <- file.path(
      outdir,
      glue::glue("go_gsea.data.{ont}.rds")
    )
    go_ont_plot_rds <- file.path(
      outdir,
      glue::glue("go_gsea.plot.{ont}.rds")
    )
    if(file.exists(go_ont_data_rds) & !overwrite) {
      go_ont_data <- readRDS(go_ont_data_rds)
    } else {
      go_ont_data <- tryCatch(
        {
          go_ont_data <- clusterProfiler::gseGO(
            gene         = gsea_gene,
            keyType      = keytype,
            OrgDb        = orgdb,
            ont          = ont,
            # nPerm        = 1000,
            minGSSize    = 120,
            # maxGSSize    = 500,
            pvalueCutoff = pval_cutoff,
            verbose      = FALSE)
          #--Save to rds
          saveRDS(go_ont_data, file = go_ont_data_rds)
          return(go_ont_data)
        },
        error=function(cond) {
          warning("gseaGO() failed")
          return(NULL)
        }
      )
    }
    #--------------------------------------------------------------------------#
    #-- run: go plots
    if(file.exists(go_ont_plot_rds) & !overwrite) {
      go_ont_plot <- readRDS(go_ont_plot_rds)
    } else if(inherits(go_ont_data, "enrichResult")) {
      go_ont_plot <- go_gsea_plots(go_ont_data, !!!dots)# !!!! group_go_plot
      saveRDS(go_ont_plot, file = go_ont_plot_rds)
    } else {
      warning("go_gsea() failed")
      go_ont_plot <- NULL
    }
    #--------------------------------------------------------------------------#
    #-- run: save to png files
    prefix <- glue::glue("go_gsea.plot.{ont}")
    save_go_plot(go_ont_plot, outdir, prefix)
    #-- run: save to table
    prefix <- glue::glue("go_gsea.data.{ont}")
    save_go_table(go_ont_data, outdir, prefix)
    #--Return:
    go_ont_data
  }, USE.NAMES = TRUE)
  #--GO plotting: wego plot: to-do
  # wego <- go_wego_plot(go_data)
  # save_go_plot(wego, outdir, "group_go.plot.wego")
}


#' @describeIn gsea_input Prepare data for GSEA analysis
#' require sorted values (fold_change, ...), with names (gene)
#'
#'
#' @param gene_list character, gene name
#' @param organism character, name of the genome
#' @param fold_change numeric, named
#'
#' @export
prep_go_gsea <- function(gene_list, organism, ...) {
  # (gene_list, fold_change, organism = NULL, orgdb = NULL) {
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
    kegg_keytype = NULL,
    .cutoff      = 0.8  # for guessing, te,piRC names
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
  if(!is(gene_list, "character") | length(gene_list) == 0) {
    warning("`gene_list` expect characters, failed")
    return(NULL)
  }
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
  #-- fold_change: numeric, named
  if(inherits(fold_change, "numeric")
     & inherits(names(fold_change), "character")) {
    message("fold_change is numeric, with names")
  } else {
    warning(glue::glue(
      "fold_change is {class(fold_change)}, expect named-numeric"
    ))
    return(NULL)
  }
  # gene_list_fc <- fold_change[gene_list]
  # gene_list_fc <- purrr::discard(gene_list_fc, is.na)
  #----------------------------------------------------------------------------#
  #-- prepare gsea_gene
  #-- round-1
  #-- assign 'fold_change' value to 'gene_list'
  fc1  <- fold_change[gene_list] # found
  fc0  <- length(gene_list) - length(fc1) # not found
  pct1 <- round(length(fc1) / length(gene_list) * 100, 2) # divided by zero ?!
  message(glue::glue(
    "{length(fc1)} of {length(gene_list)} ",
    "({pct1}%) genes found with fold_change, ",
    "{length(fc0)} genes not."
  ))
  if(pct1 > 0) {
    gsea_gene <- fc1
  } else {
    #----------------------------------------------------------------------------#
    #-- round-2
    #-- in case, gene_list and names(fold_change) not match
    #-- so guess both ketypes of gene_list and fold_change
    fc_keytype <- guess_keytype(names(fold_change), organism)
    g_str  <- paste(gene_list[1:3], collapse = ", ")
    fc_str <- paste(names(fold_change)[1:3], collapse = ", ")
    message(glue::glue(
      "gene_list keytype is '{keytype}': {g_str} ..., \n",
      "names(fold_change) keytype is '{fc_keytype}': {fc_str} ..., \n",
      "the keytypes not match, try converting the names ..."
    ))
    #-- check keytypes of fold_change
    if(is_valid_keytype(keytype, organism = organism)
       & is_valid_keytype(fc_keytype, organism = organism)) {
      trans_table <- convert_id(
        gene_list,
        from_keytype = keytype,
        to_keytype   = fc_keytype,
        organism     = organism,
        rm_na        = FALSE
      )
      #-- check, genes exists or not
      tt2  <- setNames(trans_table[[1]], nm = trans_table[[2]]) # convert
      fc2  <- fold_change[trans_table[[2]]] # valid fold_change
      fc0  <- length(gene_list) - length(fc2)
      pct2 <- round(length(fc2) / length(gene_list) * 100, 2)
      message(glue::glue(
        "{length(fc2)} of {length(gene_list)} ",
        "({pct2}%) genes found with fold_change, ",
        "{length(fc1)} genes not."
      ))
      if(pct2 > 0) {
        gsea_gene <- setNames(fc2, nm = tt2[names(fc2)]) # convert names
      }
    } else {
      message(glue::glue(
        "invalid ketypes: gene_list is {keytype}, fold_change is {fc_keytype}"
      ))
    }
  }
  gsea_gene <- sort(gsea_gene, decreasing = TRUE)
  #----------------------------------------------------------------------------#
  #-- return
  purrr::list_modify(
    dots,
    outdir    = outdir,
    organism  = organism,
    orgdb     = orgdb,
    keytype   = keytype,
    gsea_gene = gsea_gene,
    kegg_gene = kegg_gene
  )
}


#' create plots
#' @param x object of gseaGO
#'
#' @export
go_gsea_plots <- function(x, fold_change = NULL, ...) {
  if(class(x) == "gseaResult") {
    if(nrow(x)) {
      list(gseaplot = go_gsea_plot(x))
    } else {
      warning("`x` is gseaResult, but contains 0 rows data")
      NULL
    }
  } else {
    warning("`x` not gseaResult")
    NULL
  }
}

