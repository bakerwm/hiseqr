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
#' @name go_group



#' @describeIn go_group groupGO analysis using clusterProfiler
#'
#' mission:
#' - perform groupGO analysis
#' - save groupGOResult to local file: go_group_data_BP.rds # BP, CC, MF
#' - make plots (barplot, wego)
#' - save plots (ggplot) to local file: go_group_plot_BP.rds # BP, CC, MF
#'
#'
#' @param gene_list character Genes for GO analysis
#' @param organism character Name of the organism, eg: dm3, fruitfly
#' @param outdir string Path to the directory, saving go_group results
#' ..., pass the following arguments from parent function
#'
#' or else, set default values:
#' orgdb    = NULL
#' keytype  = NULL
#' level    = 2
#' readable = TRUE
#'
#' @param orgdb OrgDb from AnnotationDbi, overwrite `organism`
#' @param keytype character Name of the keytype for input, default: NULL
#' @param level int GO levels, default: 2
#' @param readable bool Args for groupGO function, convert to gene symbol
#'
#' @import clusterProfiler
#' @import clusterProfiler.dplyr
#' @import stringr
#' @import cowplot
#'
#' @example go_group(gene =, OrgDb = , keyType = , level = , ont = , readable = )
#'
#' @export
go_group <- function(gene_list, organism, ...) {
  #----------------------------------------------------------------------------#
  #-- Check: args
  dots <- rlang::list2(...)
  dots <- prep_go(gene_list, organism, !!!dots) # update arguments
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
  if(!is_valid_go_input(!!!dots)) {
    message("go_group() skipped, invalid arguments, check above message")
    return(NULL)
  }
  #-- outdir, update
  outdir <- file.path(outdir, "go_group") # update: outdir
  check_path(outdir)
  #----------------------------------------------------------------------------#
  #-- run
  onts <- c("BP", "CC", "MF")
  go_data <- sapply(onts, function(ont) {
    message(glue::glue("run groupGO() for {ont}"))
    go_ont_data_rds <- file.path(
      outdir,
      glue::glue("go_group.data.{ont}.rds")
    )
    go_ont_plot_rds <- file.path(
      outdir,
      glue::glue("go_group.plot.{ont}.rds")
    )
    #--------------------------------------------------------------------------#
    #-- run: go_group
    if(file.exists(go_ont_data_rds) & !overwrite) {
      go_ont_data <- readRDS(go_ont_data_rds)
    } else {
      go_ont_data <- tryCatch(
        {
          gd <- clusterProfiler::groupGO(
            gene     = gene_list,
            OrgDb    = orgdb,
            keyType  = keytype,
            ont      = ont,
            level    = level,
            readable = readable)
          saveRDS(gd, file = go_ont_data_rds)
          return(gd)
        },
        error=function(cond) {
          warning("groupGO() failed")
          return(NULL)
        }
      )
    }
    #--------------------------------------------------------------------------#
    #-- run: go plots
    if(file.exists(go_ont_plot_rds) & !overwrite) {
      go_ont_plot <- readRDS(go_ont_plot_rds)
    } else if(inherits(go_ont_data, "groupGOResult")) {
      go_ont_plot <- go_group_plot(go_ont_data, !!!dots)# !!!! go_group_plot
      saveRDS(go_ont_plot, file = go_ont_plot_rds)
    } else {
      warning("go_group() failed")
      go_ont_plot <- NULL
    }
    #--------------------------------------------------------------------------#
    #-- run: save to png files
    prefix <- glue::glue("go_group.plot.{ont}")
    save_go_plot(go_ont_plot, outdir, prefix)
    #-- run: save to table
    prefix <- glue::glue("go_group.data.{ont}")
    save_go_table(go_ont_data, outdir, prefix)
    #--Return:
    go_ont_data
  }, USE.NAMES = TRUE)
  #--GO plotting: wego plot
  wego <- go_wego_plot(go_data)
  save_go_plot(wego, outdir = outdir, name = "go_group.plot.wego")
}


#' create plots
#' @param x object of enrichGO
#' @param ... passing extra arguments
#' parent function: get_go_plots(),
#'
#' text_width
#'
#' @export
go_group_plot <- function(x, ...) {
  #--Default values: BEGIN
  dots <- rlang::list2(...)
  args <- list(
    fold_change = NULL,
    text_width  = 40
  )
  dots <- purrr::list_modify(args, !!!dots)
  if(inherits(x, "groupGOResult")) {
    if(nrow(x)) {
      list(barplot = go_barplot(x, !!!dots))
    } else {
      warning("`x` is groupGOResult data, but contains 0 rows data")
      NULL
    }
  } else {
    warning("`x` not groupGOResult")
    NULL
  }
}


