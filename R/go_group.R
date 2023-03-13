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



#' go_group
#' 
#' @param gene_list character Genes for GO analysis
#' @param organism character Name of the organism, eg: dm3, fruitfly
#' @param outdir string Path to the directory, saving go_group results
#' @param ... pass the following arguments from parent function
#'
#' @param orgdb OrgDb from AnnotationDbi, overwrite `organism`
#' @param keytype character Name of the keytype for input, default: NULL
#' @param level int GO levels, default: 2
#' @param readable bool Args for groupGO function, convert to gene symbol
#'
#'
#' @example
#'    go_group(gene =, OrgDb = , keyType = , level = , ont = , readable = )
#'
#' @export
go_group <- function(gene_list, organism, ...) {
  message(">>> run go_group()")
  #----------------------------------------------------------------------------#
  #-- Check: args
  dots <- rlang::list2(...)
  args <- prep_go(gene_list, organism, !!!dots) # update arguments
  # dots <- purrr::list_modify(
  #   dots,
  #   gene_list = gene_list,
  #   organism  = organism
  # )
  #-- to env
  for(name in names(args)) {
    assign(name, args[[name]])
  }
  #----------------------------------------------------------------------------#
  #-- Check: valid args
  if(is.null(args) || !is_valid_go_input(!!!args)) {
    message(">>> go_group() skipped, invalid arguments, check above message")
    return(NULL)
  }
  #-- outdir, update
  outdir <- file.path(outdir, "go_group") # update: outdir
  check_path(outdir)
  #----------------------------------------------------------------------------#
  #-- run
  onts <- c("BP", "CC", "MF")
  go_data <- sapply(onts, function(ont) {
    message(glue::glue(">>> run groupGO() for {ont}"))
    go_ont_data_rds <- file.path(outdir, glue::glue("go_group.data.{ont}.rds"))
    go_ont_plot_rds <- file.path(outdir, glue::glue("go_group.plot.{ont}.rds"))
    #--------------------------------------------------------------------------#
    #-- run: go_group
    if(!file.exists(go_ont_data_rds) | overwrite) {
      tmp <- tryCatch(
        {
          ggo <- clusterProfiler::groupGO(
            gene     = gene_list,
            OrgDb    = orgdb,
            keyType  = keytype,
            ont      = ont,
            level    = level,
            readable = readable
          )
          saveRDS(ggo, file = go_ont_data_rds)
        },
        error=function(cond) {
          warning("groupGO() failed")
          return(NULL)
        }
      )
    }
    if(file.exists(go_ont_data_rds)) {
      go_ont_data <- readRDS(go_ont_data_rds)
    } else {
      go_ont_data <- NULL
    }
    #--------------------------------------------------------------------------#
    #-- run: go plots
    if(file.exists(go_ont_plot_rds) & ! overwrite) {
      go_ont_plot <-  readRDS(go_ont_plot_rds)
    } else {
      go_ont_plot <- go_group_plot(go_ont_data, !!!dots) # !!!! go_group_plot
      saveRDS(go_ont_plot, file = go_ont_plot_rds)
    }
    #--------------------------------------------------------------------------#
    #-- run: save to png files
    prefix <- glue::glue("go_group.plot.{ont}")
    message("saving to png files")
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


#' go_group_plot
#' 
#' @param x object of enrichGO
#' @param ... passing extra arguments
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
  if(is_go_result(x, "groupGOResult")) {
    if(nrow(x) > 0) {
      list(barplot = go_barplot(x, !!!dots))
    } else {
      warning("`x` is groupGOResult data, but contains 0 rows data")
    }
  } else {
    warning("`x` not groupGOResult")
  }
}


