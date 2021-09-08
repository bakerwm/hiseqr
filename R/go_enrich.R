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
#' @name go_gsea



#' @describeIn go_enrich enrichGO analysis from clusterprofiler
#'
#' for clusterProfiler::enrichGO
#'
#' @param gene_list character Genes for GO analysis
#' @param organism character Name of the organism, eg: dm3, fruitfly
#' @param ... pass arguments: , orgdb, keytype, pval_cutoff, qval_cutoff,
#' readable, with default values:
#' orgdb = NULL
#' keytype = NULL
#' pval_cutoff = 0.9
#' qval_cutoff = 0.9
#' readable = TRUE
#'
#' @param orgdb OrgDb from AnnotationDbi, overwrite `organism`
#' @param keytype character Name of the keytype for input, default: NULL
#' @param pval_cutoff float Cutoff for p-value, default: 0.9
#' @param qval_cutoff float Cutoff for q-value, default: 0.9
#' @param readable bool Args for groupGO function, convert to gene symbol
#'
#' @description pval_cutoff, qval_cutoff, set to 0.9, in order to return
#' enrich results anyway, filter in downstream analysis
#'
#' @import clusterProfiler
#' @import clusterProfiler.dplyr
#' @import stringr
#' @import cowplot
#'
#' @example go_enrich(gene = , organism = , keytype = ,
#' readable = TRUE, pval_cutoff = 0.05, qval_cutoff = 0.05, ...)
#'
#' @export
go_enrich <- function(gene_list, organism, ...) {
  #----------------------------------------------------------------------------#
  #-- Check: args
  dots <- rlang::list2(...)
  dots <- purrr::list_modify(
    dots,
    gene_list = gene_list,
    organism  = organism
  )
  dots <- prep_go(gene_list, organism, !!!dots) # update arguments
  #-- to env
  for(name in names(dots)) {
    assign(name, dots[[name]])
  }
  #----------------------------------------------------------------------------#
  #-- Check: valid args
  if(!is_valid_go_input(!!!dots)) {
    message("go_group() skipped, invalid arguments, check above message")
    return(NULL)
  }
  #-- outdir, updat
  outdir <- file.path(outdir, "go_enrich") # update: outdir
  check_path(outdir)
  #----------------------------------------------------------------------------#
  onts <- c("BP", "CC", "MF")
  #--GO analysis
  go_data <- sapply(onts, function(ont) {
    message(glue::glue("run enrichGO() for {ont}"))
    go_ont_data_rds <- file.path(
      outdir,
      glue::glue("go_enrich.data.{ont}.rds")
    )
    go_ont_plot_rds <- file.path(
      outdir,
      glue::glue("go_enrich.plot.{ont}.rds")
    )
    if(file.exists(go_ont_data_rds) & !overwrite) {
      go_ont_data <- readRDS(go_ont_data_rds)
    } else {
      go_ont_data <- tryCatch(
        {
          go_ont_data <- clusterProfiler::enrichGO(
            gene          = gene_list,
            OrgDb         = orgdb,
            ont           = ont,
            keyType       = keytype,
            pvalueCutoff  = pval_cutoff,
            qvalueCutoff  = qval_cutoff,
            pAdjustMethod = "BH",
            readable      = readable)
          #-------------------------------------------#
          # simplify results, redundant
          if(class(go_ont_data) == "enrichResult") {
            go_ont_data <- clusterProfiler::simplify(
              x = go_ont_data,
              cutoff = 0.7,
              by = "p.adjust",
              select_fun = min
            ) # redundant
            #--Save to rds
            saveRDS(go_ont_data, file = go_ont_data_rds)
          }
          return(go_ont_data)
        },
        error=function(cond) {
          warning("enrichGO() failed")
          return(NULL)
        }
      )
    }
    #--------------------------------------------------------------------------#
    #-- run: go plots
    if(file.exists(go_ont_plot_rds) & !overwrite) {
      go_ont_plot <- readRDS(go_ont_plot_rds)
    } else if(inherits(go_ont_data, "enrichResult")) {
      go_ont_plot <- go_enrich_polt(go_ont_data)# !!!! go_group_plot
      saveRDS(go_ont_plot, file = go_ont_plot_rds)
    } else {
      warning("go_group() failed")
      go_ont_plot <- NULL
    }
    #--------------------------------------------------------------------------#
    #-- run: save to png files
    prefix <- glue::glue("go_enrich.plot.{ont}")
    save_go_plot(go_ont_plot, outdir, prefix)
    #-- run: save to table
    prefix <- glue::glue("go_enrich.data.{ont}")
    save_go_table(go_ont_data, outdir, prefix)
    #--Return:
    go_ont_data
  }, USE.NAMES = TRUE)
  #-- GO plotting: wego plot
  wego <- go_wego_plot(go_data)
  save_go_plot(wego, outdir, "go_enrich.plot.wego")
}



#' create plots
#' @param gene_list object of enrichGO
#' @param ... passing extra arguments
#' parent function: get_go_plots(),
#'
#' fold_change
#' text_width
#'
#'
#' @export
go_enrich_plot <- function(x, ...) {
  #--Default values: BEGIN
  dots <- rlang::list2(...)
  args <- list(
    fold_change = NULL,
    text_width  = 40
  )
  dots <- purrr::list_modify(args, !!!dots)
  #--Default values: END
  if(class(x) == "enrichResult") {
    if(nrow(x)) {
      list(barplot  = go_barplot(x, !!!dots),
           dotplot  = go_dotplot(x, !!!dots),
           cnetplot = go_cnetplot(x, !!!dots),
           emapplot = go_emapplot(x, !!!dots),
           # emapplot_cluster = go_emapplot_cluster(x, ...),
           heatplot = go_heatplot(x, !!!dots))
    } else {
      warning("`x` is enrichGOResult, but contains 0 rows data")
      NULL
    }
  } else {
    warning("`x` not enrichGOResult")
    NULL
  }
}
