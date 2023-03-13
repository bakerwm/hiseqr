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


#' kegg_enrich
#'
#' @param gene_list character Genes for KEGG analysis
#' @param organism character Name of the organism, eg: dm3, fruitfly
#' @param ... pass arguments: , orgdb, keytype, pval_cutoff, qval_cutoff,
#' readable, with default values:
#' keytype = NULL
#' pval_cutoff = 0.9
#' qval_cutoff = 0.9
#' readable = TRUE
#'
#' @param keytype character Name of the keytype for input, default: NULL
#' @param pval_cutoff float Cutoff for p-value, default: 0.9
#' @param qval_cutoff float Cutoff for q-value, default: 0.9
#' @param readable bool Args for enrichKEGG function, convert to gene symbol
#'
#' pval_cutoff, qval_cutoff, set to 0.9, in order to return
#' enrich results anyway, filter in downstream analysis
#'
#'
#' @example enrich_kegg(gene = , organism = , keytype = ,
#' readable = TRUE, pval_cutoff = 0.05, qval_cutoff = 0.05, ...)
#'
#' @export
kegg_enrich <- function(gene_list, organism, ...) {
  message(glue::glue(">>> run kegg_enrich()"))
  #----------------------------------------------------------------------------#
  dots <- rlang::list2(...)
  dots <- prep_kegg(gene_list, organism, !!!dots) # update arguments
  #-- to env
  for(name in names(dots)) {
    assign(name, dots[[name]])
  }
  #----------------------------------------------------------------------------#
  #-- Check: valid args
  if(!exists("kegg_gene", inherits = FALSE)) {
    warning("kegg_gene not exists")
    return(NULL)
  }
  # warning(glue::glue("<<< 2. kegg_gene: {dots$kegg_gene}"))
  if(is_valid_kegg_input(!!!dots) & inherits(kegg_gene, "character") & length(kegg_gene) > 1) {
    message("kegg_gene is valid")
  } else {
    message("kegg_enrich() skipped, invalid arguments, check above message")
    return(NULL)
  }
  #-- outdir, updat
  outdir <- file.path(outdir, "kegg_enrich") # update: outdir
  check_path(outdir)
  #----------------------------------------------------------------------------#
  #-- run kegg
  kegg_data_rds <- file.path(outdir, glue::glue("kegg_enrich.data.rds"))
  kegg_plot_rds <- file.path(outdir,glue::glue("kegg_enrich.plot.rds"))
  if(!file.exists(kegg_data_rds) | !overwrite) {
      tmp <- tryCatch(
        {
          kk <- clusterProfiler::enrichKEGG(
            gene          = kegg_gene,
            organism      = kegg_code,
            keyType       = kegg_keyType,
            pAdjustMethod = "BH",
            pvalueCutoff  = pval_cutoff,
            qvalueCutoff  = qval_cutoff
          )
          # readable
          if(inherits(kk, "enrichResult") & inherits(orgdb, "OrgDb")) {
            kk <- clusterProfiler::setReadable(
              x       = kk,
              OrgDb   = orgdb,
              keyType = kegg_keytype
            )
          }
          saveRDS(kk, file = kegg_data_rds)
        },
        error=function(cond) {
          warning("enrich_kegg() failed")
          return(NULL)
        }
      )
  }
  if(file.exists(kegg_data_rds)) {
    kegg_data <- readRDS(kegg_data_rds)
  } else {
    kegg_data <- NULL
  }
  if(!file.exists(kegg_plot_rds) | overwrite) {
    kegg_plot <- go_enrich_plot(kegg_data, !!!dots) # go plots
    saveRDS(kegg_plot, file = kegg_plot_rds)
  }
  if(file.exists(kegg_plot_rds)) {
    kegg_plot <- readRDS(kegg_plot_rds)
  } else {
    kegg_plot <- NULL
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

