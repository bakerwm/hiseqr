#' Functions for KEGG analysis
#'
#' overrepresent analysis
#' GSEA analysis
#'
#' gsea kegg
#' 1. csv files, results
#' 2. barplot
#' 3. dotplot
#' 4. netplot
#' 5. ...
#'
#' @name kegg


#' @describeIn kegg_gsea
#' run gseKEGG()
#'
#' @param gene_list decreasing sorted numeric vector
#' @param organism kegg_code, eg: dme
#'
#' @export
kegg_gsea <- function(gene_list, organism, ...) {
  message(">>> run kegg_gsea()")
  #----------------------------------------------------------------------------#
  dots <- rlang::list2(...)
  dots <- prep_kegg_gsea(gene_list, organism, !!!dots) # update arguments
  #-- to env
  for(name in names(dots)) {
    assign(name, dots[[name]])
  }
  #----------------------------------------------------------------------------#
  d_str <- paste(as.character(names(dots)), collapse = ",")
  message(glue::glue("dots for kegg_gsea: kegg_gene: {dots$kegg_gene}; {d_str}"))
  #----------------------------------------------------------------------------#
  # must have kegg_gene (NULL ?)
  # if(inherits(kegg_gene, "character")
  #    && length(kegg_gene) > 1) {
  #   message("kegg_gene is valid")
  # } else {
  #   message("kegg_enrich() skipped, invalid kegg_gene")
  #   return(NULL)
  # }
  if(inherits(kegg_gene, "character")
     && length(kegg_gene) > 1
     && inherits(kegg_gsea_gene, "numeric")
     && inherits(names(kegg_gsea_gene), "character")) {
    message("kegg_gene and kegg_gsea_gene are valid ")
  } else {
    message(glue::glue(
      "kegg_gsea() failed, invalied 'kegg_gene', 'kegg_gsea_gene', \n",
      "kegg_gene is {class(kegg_gene)}, expect character, \n",
      "kegg_gsea_gene is {class(kegg_gsea_gene)}, expect numeric, named"
    ))
    return(NULL)
  }
  #----------------------------------------------------------------------------#
  # other args must be qualified
  if(!is_valid_kegg_input(!!!dots)) {
    message("kegg_gsea() skipped, invalid arguments, check above message")
    return(NULL)
  }
  #-- outdir, updat
  outdir <- file.path(outdir, "kegg_gsea") # update: outdir
  check_path(outdir)
  #----------------------------------------------------------------------------#
  #-- run kegg
  kegg_data_rds <- file.path(outdir, glue::glue("kegg_gsea.data.rds"))
  kegg_plot_rds <- file.path(outdir, glue::glue("kegg_gsea.plot.rds"))
  if(file.exists(kegg_data_rds) & !overwrite) {
    kegg_data <- readRDS(kegg_data_rds)
  } else {
    kegg_data <- tryCatch(
      {
        kegg_data <- clusterProfiler::gseKEGG(
          geneList      = kegg_gsea_gene,
          organism      = kegg_code,
          keyType       = kegg_keyType,
          minGSSize     = 120,
          pvalueCutoff  = pval_cutoff,
          pAdjustMethod = "BH",
          verbose       = FALSE
        )
        # readable
        if(inherits(kegg_data, "gseaResult") & inherits(orgdb, "OrgDb")) {
          kegg_data <- clusterProfiler::setReadable(
            x       = kegg_data,
            OrgDb   = orgdb,
            keyType = keytype
          )
        }
        saveRDS(kegg_data, file = kegg_data_rds)
        return(kegg_data)
      },
      error=function(cond) {
        warning("enrich_kegg() failed")
        return(NULL)
      }
    )
  }
  #-- run, plot
  if(file.exists(kegg_plot_rds) & !overwrite) {
    kegg_plot <- readRDS(kegg_plot_rds)
  } else {
    if(inherits(kegg_data, "gseaResult")) {
      kegg_plot <- go_gsea_plot(kegg_data, !!!dots) # see go...plot
      saveRDS(kegg_plot, file = kegg_plot_rds)
    } else {
      kegg_plot <- NULL
    }
  }
  #-- run: save to png files
  prefix <- gsub(".rds$", "", basename(kegg_plot_rds))
  save_go_plot(kegg_plot, outdir, prefix)
  #-- run: save to table
  prefix <- gsub(".rds$", "", basename(kegg_data_rds))
  save_go_table(kegg_data, outdir, prefix)
  #-- Return: data
  kegg_data
}


#' @describeIn prep_kegg_gsea
#' Prepare data for KEGG GSEA analysis
#'
#' require sorted values (fold_change, ...), with names (gene)
#' require: ENTREZID (FLYBASECG for fruitfly)
#'
#' @param gene_list gene name
#' @param fold_change numeric
#'
#' @export
prep_kegg_gsea <- function(gene_list, organism, ...) {
  #----------------------------------------------------------------------------#
  #-- Check: args
  dots <- rlang::list2(...)
  dots <- prep_kegg(gene_list, organism, !!!dots)#
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
    .cutoff      = 0.8  # for guessing, te,piRC names
  )
  dots <- purrr::list_modify(args, !!!dots) # args, overwrite by dots
  #-- force, update positional args
  # dots <- purrr::list_modify(
  #   dots,
  #   gene_list = gene_list,
  #   organism  = organism
  # )
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
  #----------------------------------------------------------------------------#
  #-- prepare kegg_gsea_gene
  #-- assign 'fold_change' to 'gene_list'
  #-- round-1, keytype == fc_keytype
  fc1  <- fold_change[gene_list] # found
  fc0  <- length(gene_list) - length(fc1) # not found
  pct1 <- round(length(fc1) / length(gene_list) * 100, 2) # divided by zero ?!
  if(pct1 > 0) {
    message(glue::glue(
      "{length(fc1)} of {length(gene_list)} ",
      "({pct1}%) genes found with fold_change, ",
      "{length(fc0)} genes not."
    ))
    kegg_gsea_gene <- fc1
    # kegg_gsea_gene <- sort(fc1, decreasing = TRUE) #
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
       & is_valid_keytype(fc_keytype, organism = organism)
      ) {
      trans_table <- convert_id(
        gene_list,
        from_keytype = keytype,
        to_keytype   = fc_keytype,
        organism     = organism,
        rm_na        = TRUE
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
        kegg_gsea_gene <- setNames(fc2, nm = tt2[names(fc2)]) # convert names
        # kegg_gsea_gene <- sort(kegg_gsea_gene, decreasing = TRUE) #
      }
    } else {
      message(glue::glue(
        "invalid ketypes: ",
        "gene_list is '{keytype}', ",
        "fold_change is '{fc_keytype}'"
      ))
    }
    #--------------------------------------------------------------------------#
    # check kegg_gsea_gene, keytype -> kegg_keytype
    kegg_keytype <- get_kegg_keytype(organism)
    if(inherits(kegg_gsea_gene, "numeric")
       & inherits(names(kegg_gsea_gene), "character")
       & !keytype == kegg_keytype
       & inherits(kegg_gsea_gene, "numeric")
       & inherits(names(kegg_gsea_gene), "character")
      ) {
      # convert table
      trans_table <- to_kegg_gene_id(
        gene_list    = names(kegg_gsea_gene),
        organism     = organism,
        keytype      = keytype,
        return_table = TRUE,
        rm_na        = TRUE
      )
      tt3  <- setNames(trans_table[[2]], nm = trans_table[[1]]) #
      fc3  <- setNames(
        kegg_gse_gene,
        nm = tt3[names(kegg_gse_gene)]
      )
      fc0  <- length(kegg_gsea_gene) - length(fc3)
      pct3 <- round(length(fc3) / length(kegg_gsea_gene) * 100, 2)
      message(glue::glue(
        "{length(fc3)} of {length(kegg_gse_gene)} ",
        "({pct3}%) genes found with fold_change, ",
        "{length(fc0)} genes not."
      ))
      if(pct3 > 0) {
        kegg_gsea_gene <- fc3[!is.na(names(fc3))] # remove na
        # names(kegg_gsea_gene) <- tt[names(kegg_gsea_gene)]
        # kegg_gsea_gene <- kegg_gsea_gene[!is.na(names(kegg_gsea_gene))]
        # kegg_gsea_gene <- sort(kegg_gsea_gene, decreasing = TRUE)
      }
    }
  }
  #-- sort
  kegg_gsea_gene <- sort(kegg_gsea_gene, decreasing = TRUE)
  #----------------------------------------------------------------------------#
  #-- return
  purrr::list_modify(
    dots,
    outdir    = outdir,
    organism  = organism,
    orgdb     = orgdb,
    keytype   = keytype,
    kegg_gsea_gene = kegg_gsea_gene
  )
}


#' #' @export
#' kegg_gsea_plot <- function(x, ...) {
#'   go_gsea_plot(x, ...)
#' }


