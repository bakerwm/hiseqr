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
#' @name go



#' @describeIn run_go Main GO analysis function for gene list
#'
#' @param gene_list character Genes for GO analysis
#' @param organism character Name of the organism, eg: dm3, fruitfly
#' @param outdir string Path to the directory, saving group_go results
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
#' group_go(), keytype, level, readable
#' enrich_go(), keytype, readable, pval_cutoff, qval_cutoff
#' gsea_go(), keytype, pval_cutoff
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
#'
#'
#' @export
run_go <- function(gene_list, organism, outdir = NULL, ...) {
  outdir <- ifelse(is(outdir, "character"), outdir, getwd()) # working dir
  dots   <- list(...)
  arg_vars <- list(
    orgdb       = NULL,
    keytype     = NULL,
    fold_change = NULL,
    overwrite   = FALSE,
    level       = 2,
    readable    = TRUE
  )
  arg_vars <- purrr::list_modify(arg_vars, !!!dots)
  arg_vars <- purrr::list_modify(arg_vars,
                                 gene_list = gene_list,
                                 organism  = organism,
                                 outdir    = outdir) #@@@ arg_vars updated
  arg_valid <- do.call(is_valid_go_input, arg_vars)
  if(! arg_valid) {
    stop("Check above errors, prepare correct arguments for GO analysis")
  }
  # return: organism, orgdb, keytype, gsea_gene, kegg_code
  go_input <- do.call(prep_go_input, arg_vars)
  if(is(go_input, "list")) {
    arg_vars <- purrr::list_modify(arg_vars, !!!go_input) #@@@ arg_vars, updated
    do.call(group_go, arg_vars)
    do.call(enrich_go, arg_vars)
    do.call(gsea_go, arg_vars)
  }
}



#' @describeIn group_go groupGO analysis using clusterProfiler
#'
#' mission:
#' - perform groupGO analysis
#' - save groupGOResult to local file: group_go_data_BP.rds # BP, CC, MF
#' - make plots (barplot, wego)
#' - save plots (ggplot) to local file: group_go_plot_BP.rds # BP, CC, MF
#'
#'
#' @param gene_list character Genes for GO analysis
#' @param organism character Name of the organism, eg: dm3, fruitfly
#' @param outdir string Path to the directory, saving group_go results
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
#' @example group_go(gene =, OrgDb = , keyType = , level = , ont = , readable = )
#'
#' @export
group_go <- function(gene_list, organism, outdir = NULL, ...) {
  outdir <- ifelse(is(outdir, "character"), outdir, getwd()) # working dir
  dots   <- list(...)
  arg_vars <- list(
    orgdb     = NULL,
    keytype   = NULL,
    overwrite = FALSE,
    level     = 2,
    readable  = TRUE
  )
  arg_vars <- purrr::list_modify(arg_vars, !!!dots)
  arg_vars <- purrr::list_modify(arg_vars,
                                 gene_list = gene_list,
                                 organism  = organism,
                                 outdir    = outdir)
  arg_valid <- do.call(is_valid_go_input, arg_vars)
  if(! arg_valid) {
    stop("Check above errors, prepare correct arguments for GO analysis")
  }
  #--Ontology:
  arg_vars$outdir <- file.path(arg_vars$outdir, "group_go") # update: outdir
  if(! dir.exists(arg_vars$outdir)) {
    dir.create(arg_vars$outdir, recursive = TRUE, mode = "0755")
  }
  onts <- c("BP", "CC", "MF")
  #--GO analysis:
  go_data <- sapply(onts, function(ont) {
    message(glue::glue("Running groupGO() for {ont}"))
    go_ont_data_rds <- file.path(arg_vars$outdir, paste("group_go", "data", ont, "rds", sep = "."))
    go_ont_plot_rds <- file.path(arg_vars$outdir, paste("group_go", "plot", ont, "rds", sep = "."))
    if(file.exists(go_ont_data_rds) & ! arg_vars$overwrite) {
      go_ont_data <- readRDS(go_ont_data_rds)
    } else {
      go_ont_data <- tryCatch(
        {
          go_ont_data <- clusterProfiler::groupGO(
            gene     = arg_vars$gene_list,
            OrgDb    = arg_vars$orgdb,
            keyType  = arg_vars$keytype,
            ont      = ont,
            level    = arg_vars$level,
            readable = arg_vars$readable)
          saveRDS(go_ont_data, file = go_ont_data_rds)
          return(go_ont_data)
        },
        error=function(cond) {
          warning("groupGO() failed")
          return(NULL)
        }
      )
    }
    #--GO plotting: plot
    if(is.null((go_ont_data))) {
      warning("groupGO() failed, not results")
    } else {
      if(file.exists(go_ont_plot_rds) & ! arg_vars$overwrite) {
        go_ont_plot <- readRDS(go_ont_plot_rds)
      } else {
        go_ont_plot <- group_go_plots(go_ont_data)# !!!! group_go_plot
        saveRDS(go_ont_plot, file = go_ont_plot_rds)
      }
      #--GO plotting: save to png
      prefix <- gsub(".rds$", "", basename(go_ont_plot_rds))
      save_go_plot(go_ont_plot, arg_vars$outdir, prefix)
      #--GO table: save result to csv
      prefix <- gsub(".rds$", "", basename(go_ont_data_rds))
      save_go_table(go_ont_data, arg_vars$outdir, prefix)
    }
    #--Return:
    go_ont_data
  }, USE.NAMES = TRUE)
  #--GO plotting: wego plot
  wego <- go_wego_plot(go_data)
  save_go_plot(wego, outdir = arg_vars$outdir, name = "group_go.plot.wego")
}




#' @describeIn enrich_go enrichGO analysis from clusterprofiler
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
#' @example enrich_go(gene = , organism = , keytype = ,
#' readable = TRUE, pval_cutoff = 0.05, qval_cutoff = 0.05, ...)
#'
#' @export
enrich_go <- function(gene_list, organism, outdir = NULL, ...) {
  outdir <- ifelse(is(outdir, "character"), outdir, getwd()) # working dir
  dots   <- list(...)
  arg_vars <- list(
    orgdb     = NULL,
    keytype   = NULL,
    overwrite = FALSE,
    readable  = TRUE,
    pval_cutoff = 0.9,
    qval_cutoff = 9.9
  )
  arg_vars <- purrr::list_modify(arg_vars, !!!dots)
  arg_vars <- purrr::list_modify(arg_vars,
                                 gene_list = gene_list,
                                 organism  = organism,
                                 outdir    = outdir)
  arg_valid <- do.call(is_valid_go_input, arg_vars)
  if(! arg_valid) {
    stop("Check above errors, prepare correct arguments for GO analysis")
  }
  #--Ontology:
  arg_vars$outdir <- file.path(arg_vars$outdir, "enrich_go") # update: outdir
  if(! dir.exists(arg_vars$outdir)) {
    dir.create(arg_vars$outdir, recursive = TRUE, mode = "0755")
  }
  onts <- c("BP", "CC", "MF")
  #--GO analysis
  go_data <- sapply(onts, function(ont) {
    message(glue::glue("Running enrichGO() for {ont}"))
    go_ont_data_rds <- file.path(arg_vars$outdir, paste("enrich_go", "data", ont, "rds", sep = "."))
    go_ont_plot_rds <- file.path(arg_vars$outdir, paste("enrich_go", "plot", ont, "rds", sep = "."))
    if(file.exists(go_ont_data_rds) & ! arg_vars$overwrite) {
      go_ont_data <- readRDS(go_ont_data_rds)
    } else {
      go_ont_data <- tryCatch(
        {
          ego <- clusterProfiler::enrichGO(
            gene          = arg_vars$gene_list,
            OrgDb         = arg_vars$orgdb,
            ont           = ont,
            keyType       = arg_vars$keytype,
            pvalueCutoff  = arg_vars$pval_cutoff,
            qvalueCutoff  = arg_vars$qval_cutoff,
            pAdjustMethod = "BH",
            readable      = arg_vars$readable)
          if(class(ego) == "enrichResult") {
            go_ont_data <- clusterProfiler::simplify(
              ego, cutoff = 0.7, by = "p.adjust",
              select_fun = min) # redundant
            #--Save to rds
            saveRDS(go_ont_data, file = go_ont_data_rds)
          } else {
            go_ont_data <- ego # skipped
          }
          return(go_ont_data)
        },
        error=function(cond) {
          warning("enrichGO() failed")
          return(NULL)
        }
      )
    }
    #--GO plotting: plot
    if(is.null(go_ont_data)) {
      message("enrichGO() failed, no results")
    } else {
      if(file.exists(go_ont_plot_rds) & ! arg_vars$overwrite) {
        go_ont_plot <- readRDS(go_ont_plot_rds)
      } else {
        go_ont_plot <- enrich_go_plots(go_ont_data) # !!!! enrich_go_plot
        saveRDS(go_ont_plot, file = go_ont_plot_rds)
      }
      #--GO plotting: save to png
      prefix <- gsub(".rds$", "", basename(go_ont_plot_rds))
      save_go_plot(go_ont_plot, arg_vars$outdir, prefix)
      #--GO table: save result to csv
      prefix <- gsub(".rds$", "", basename(go_ont_data_rds))
      save_go_table(go_ont_data, arg_vars$outdir, prefix)
    }
    #--Return: data
    go_ont_data
  }, USE.NAMES = TRUE)
  #--GO plotting: wego plot
  wego <- go_wego_plot(go_data)
  save_go_plot(wego, arg_vars$outdir, "group_go.plot.wego")
}



#' @describeIn gsea_go gseGO analysis from clusterProfiler
#'
#' for clusterProfiler::groupGO, and enrichGO, gseGO
#'
#' @param gsea_gene numeric
#' @param organism character Name of the organism, eg: dm3, fruitfly
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
gsea_go <- function(gsea_gene, organism, outdir = NULL, ...) {
  if(missing(gsea_gene)) {
    warning("`gsea_gene` missing, `gsea_go()` skipped")
    return(NULL)
  }
  outdir <- ifelse(is(outdir, "character"), outdir, getwd()) # working dir
  dots   <- list(...)
  arg_vars <- list(
    orgdb     = NULL,
    keytype   = NULL,
    overwrite = FALSE,
    readable  = TRUE,
    pval_cutoff = 0.9,
    qval_cutoff = 9.9
  )
  arg_vars <- purrr::list_modify(arg_vars, !!!dots)
  arg_vars <- purrr::list_modify(arg_vars,
                                 gsea_gene = gsea_gene,
                                 organism  = organism,
                                 outdir    = outdir)
  arg_valid <- do.call(is_valid_go_input, arg_vars)
  if(! arg_valid) {
    stop("Check above errors, prepare correct arguments for GO analysis")
  }
  #--Ontology:
  arg_vars$outdir <- file.path(arg_vars$outdir, "gsea_go") # update: outdir
  if(! dir.exists(arg_vars$outdir)) {
    dir.create(arg_vars$outdir, recursive = TRUE, mode = "0755")
  }
  onts <- c("BP", "CC", "MF")
  #--GO analysis
  go_data <- sapply(onts, function(ont) {
    message(glue::glue("Running gseaGO() for {ont}"))
    go_ont_data_rds <- file.path(arg_vars$outdir, paste("gsea_go", "data", ont, "rds", sep = "."))
    go_ont_plot_rds <- file.path(arg_vars$outdir, paste("gsea_go", "plot", ont, "rds", sep = "."))
    if(file.exists(go_ont_data_rds) & ! arg_vars$overwrite) {
      go_ont_data <- readRDS(go_ont_data_rds)
    } else {
      go_ont_data <- tryCatch(
        {
          go_ont_data <- clusterProfiler::gseGO(
            gene         = arg_vars$gsea_gene,
            keyType      = arg_vars$keytype,
            OrgDb        = arg_vars$orgdb,
            ont          = ont,
            # nPerm        = 1000,
            minGSSize    = 120,
            # maxGSSize    = 500,
            pvalueCutoff = arg_vars$pval_cutoff,
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
    #--GO plotting: plot
    if(is.null(go_ont_data)) {
      message("gseaGO() failed, no results")
    } else {
      if(file.exists(go_ont_plot_rds) & ! arg_vars$overwrite) {
        go_ont_plot <- readRDS(go_ont_plot_rds)
      } else {
        go_ont_plot <- gsea_go_plots(go_ont_data) # !!!! enrich_go_plot
        saveRDS(go_ont_plot, file = go_ont_plot_rds)
      }
      #--GO plotting: save to png
      prefix <- gsub(".rds$", "", basename(go_ont_plot_rds))
      save_go_plot(go_ont_plot, arg_vars$outdir, prefix)
      #--GO table: save result to csv
      prefix <- gsub(".rds$", "", basename(go_ont_data_rds))
      save_go_table(go_ont_data, arg_vars$outdir, prefix)
    }
    #--Return: data
    go_ont_data
  }, USE.NAMES = TRUE)
  #--GO plotting: wego plot: to-do
  # wego <- go_wego_plot(go_data)
  # save_go_plot(wego, arg_vars$outdir, "group_go.plot.wego")
}





#' create plots
#' @param x object of enrichGO
#' @param ... passing extra arguments
#' parent function: get_go_plots(),
#'
#' text_width
#'
#' @export
group_go_plots <- function(x, ...) {
  if(class(x) == "groupGOResult") {
    if(nrow(x)) {
      list(barplot = go_barplot(x, ...))
    } else {
      warning("`x` is groupGOResult data, but contains 0 rows data")
      NULL
    }
  } else {
    warning("`x` not groupGOResult")
    NULL
  }
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
enrich_go_plots <- function(x, ...) {
  #--Default values: BEGIN
  dots <- list(
    fold_change = NULL,
    text_width  = 40
  )
  dots <- purrr::list_modify(dots, !!!list(...))
  #--Default values: END
  if(class(x) == "enrichResult") {
    if(nrow(x)) {
      list(barplot  = go_barplot(x, ...),
           dotplot  = go_dotplot(x, ...),
           cnetplot = go_cnetplot(x, ...),
           emapplot = go_emapplot(x, ...),
           # emapplot_cluster = go_emapplot_cluster(x, ...),
           heatplot = go_heatplot(x, ...))
    } else {
      warning("`x` is enrichGOResult, but contains 0 rows data")
      NULL
    }
  } else {
    warning("`x` not enrichGOResult")
    NULL
  }
}



#' create plots
#' @param x object of gseaGO
#'
#' @export
gsea_go_plots <- function(x, fold_change = NULL, ...) {
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




#--Deprecated, rewrite, re-construct



# run_go <- function(gene_list, organism = NULL, outdir = NULL, ...) {
#   #--Default values: BEGIN
#   dots <- list(
#     keytype     = NULL,
#     fold_change = NULL,
#     level       = 2,
#     overwrite   = FALSE,
#     readable    = TRUE
#   )
#   dots <- purrr::list_modify(dots, !!!list(...))
#   #--Default values: END
#   if(! is_valid_go_input(gene_list, organism, outdir = outdir, ...)) {
#     stop("Check above errors, prepare correct arguments for GO analysis")
#   }
#   #--Init: dots
#   # orgdb, keytype, gsea_gene, kegg_code
#   go_input <- prep_go_input(gene_list, organism, ...) # fold_change, keytype
#   dots     <- purrr::list_modify(dots, gene_list = gene_list, organism = organism)
#
#   #--Run: go analysis, plotting ...
#   # if(is(go_input, "list")) {
#   #   dots <- purrr::list_modify(dots, !!!go_input)
#   #   ## deprecated: re-analysis for each ont, takes too-much long time.
#   #   # go_data_rds <- file.path(outdir, "go_data.rds")
#   #   # go_plot_rds <- file.path(outdir, "go_plot.rds")
#   #   #--Run GO: groupGO, enrichGO, gseaGO, ...
#   #   if(file.exists(go_data_rds) && ! isTRUE(dots$overwrite)) {
#   #     message(glue::glue("Loading pre-computed GO results: {go_data_rds}"))
#   #     go_data <- readRDS(go_data_rds)
#   #   } else {
#   #     go_data <- list(
#   #       group  = do.call(group_go, dots),
#   #       enrich = do.call(enrich_go, dots),
#   #       gsea   = do.call(gsea_go, dots),
#   #       # group  = group_go(gene_list, organism, ...),
#   #       # enrich = enrich_go(gene_list, organism, ...),
#   #       # gsea   = gsea_go(go_input$gsea_gene, organism, go_input$fold_change, ...),
#   #       organism    = dots$organism,
#   #       gene_list   = dots$gene_list,
#   #       fold_change = dots$fold_change,
#   #       keytype     = dots$keytype
#   #     )
#   #     saveRDS(go_data, file = go_data_rds)
#   #   }
#   #   # #--Run GO: generating plots
#   #   # if(file.exists(go_plot_rds) && ! isTRUE(dots$overwrite)) {
#   #   #   message(glue::glue("Loading pre-computed GO plots: {go_plot_rds}"))
#   #   #   go_plot <- readRDS(go_plot_rds)
#   #   # } else {
#   #   #   #--Plot GO: groupGO ...
#   #   #   group_go_res <- go_data$group
#   #   #   group_go_dir <- file.path(dots$outdir, "group_go")
#   #   #   table1       <- save_go_table(group_go_res, group_go_dir)
#   #   #   plot1        <- get_go_plots(group_go_res, ...) # fold_change,
#   #   #   p1 <- save_go_plot(plot1, group_go_dir, "group_go")
#   #   #
#   #   #   #--Plot: enrichGO ...
#   #   #   enrich_go_res <- go_data$enrich
#   #   #   enrich_go_dir <- file.path(dots$outdir, "enrich_go")
#   #   #   table2        <- save_go_table(enrich_go_res, enrich_go_dir)
#   #   #   plot2         <- get_go_plots(enrich_go_res, ...)
#   #   #   p2 <- save_go_plot(plot2, enrich_go_dir, "enrich_go")
#   #   #
#   #   #   #--Plot: gseaGO ...
#   #   #   gsea_go_res <- go_data$gsea
#   #   #   gsea_go_dir <- file.path(dots$outdir, "gsea_go")
#   #   #   table3      <- save_go_table(gsea_go_res, gsea_go_dir)
#   #   #   plot3       <- get_go_plots(gsea_go_res, ...)
#   #   #   p3 <- save_go_plot(plot3, gsea_go_dir, "gsea_go")
#   #   #
#   #   #   #--Save: save ggplots to rds file
#   #   #   go_plot <- list(group  = plot1,
#   #   #                   enrich = plot2,
#   #   #                   gsea   = plot3)
#   #   #   saveRDS(go_plot, file = go_plot_rds)
#   #   # }
#   #   list(go_data = go_data,
#   #        go_plot = go_plot)
#   # }
# }




# group_go <- function(gene_list, organism, ...) {
#   #--Default values: BEGIN
#   dots <- list(
#     orgdb    = NULL,
#     keytype  = NULL,
#     level    = 2,
#     readable = TRUE
#   )
#   dots <- purrr::list_modify(dots, !!!list(...))
#   #--Default values: END
#   if(is(gene_list, "character")) {
#     # check orgdb
#     if(is(dots$orgdb, "OrgDb")) {
#       dots$orgdb <- dots$orgdb
#     } else if(is(organism, "character")) {
#       dots$orgdb <- get_orgdb(organism)
#     }
#     if(! is(dots$orgdb, "OrgDb")) {
#       stop("`orgdb` and `organism`, not valid")
#     }
#     # check keytype
#     if(is.null(dots$keytype)) {
#       organism <- get_orgdb_metadata(dots$orgdb, "ORGANISM")
#       dots$keytype  <- guess_keytype(gene_list, organism)
#     }
#     if(! is_valid_keytype(dots$keytype, dots$orgdb)) {
#       kt_list <- paste(AnnotationDbi::keytypes(dots$orgdb), collapse = ", ")
#       msg     <- paste("`keytype` not valid: [", dots$keytype, "]",
#                        "choose from: ", kt_list, sep = " ")
#       stop(msg)
#     }
#     # run clusterProfiler::groupGO
#     onts <- c("BP", "CC", "MF")
#     sapply(onts, function(ont) {
#       message(glue::glue("Running groupGO() for {ont}"))
#       clusterProfiler::groupGO(gene     = gene_list,
#                                OrgDb    = dots$orgdb,
#                                keyType  = dots$keytype,
#                                ont      = ont,
#                                level    = dots$level,
#                                readable = dots$readable)
#     }, USE.NAMES = TRUE)
#   } else {
#     stop("`x` not characters")
#   }
# }




# enrich_go <- function(gene_list, organism, ...) {
#   #--Default values: BEGIN
#   dots <- list(
#     orgdb       = NULL,
#     keytype     = NULL,
#     pval_cutoff = 0.9,
#     qval_cutoff = 0.9,
#     readable    = TRUE
#   )
#   dots <- purrr::list_modify(dots, !!!list(...))
#   #--Default values: END
#   if(is(gene_list, "character")) {
#     # check orgdb
#     if(is(dots$orgdb, "OrgDb")) {
#       dots$orgdb <- dots$orgdb
#     } else if(is(organism, "character")) {
#       dots$orgdb <- get_orgdb(organism)
#     }
#     if(! is(dots$orgdb, "OrgDb")) {
#       stop("`orgdb` and `organism`, not valid OrgDb")
#     }
#     # check keytype
#     if(is.null(dots$keytype)) {
#       organism <- get_orgdb_metadata(dots$orgdb, "ORGANISM")
#       dots$keytype  <- guess_keytype(gene_list, organism)
#     }
#     if(! is_valid_keytype(dots$keytype, dots$orgdb)) {
#       kt_list <- paste(AnnotationDbi::keytypes(dots$orgdb), collapse = ", ")
#       msg     <- glue::glue(
#         "`keytype` not valid: [{dots$keytype}], ",
#         "choose from: {kt_list}")
#       stop(msg)
#     }
#     # run clusterProfiler::enrichGO
#     onts <- c("BP", "CC", "MF")
#     sapply(onts, function(ont) {
#       message(glue::glue("Running enrichGO() for {ont}"))
#       ego <- clusterProfiler::enrichGO(gene          = gene_list,
#                                        OrgDb         = dots$orgdb,
#                                        ont           = ont,
#                                        keyType       = dots$keytype,
#                                        pvalueCutoff  = dots$pval_cutoff,
#                                        qvalueCutoff  = dots$qval_cutoff,
#                                        pAdjustMethod = "BH",
#                                        readable      = dots$readable)
#       clusterProfiler::simplify(ego, cutoff = 0.7, by = "p.adjust",
#                                 select_fun = min) # redundant
#     }, USE.NAMES = TRUE)
#   } else {
#     stop("`x` not characters")
#   }
# }


#' get_go_plots <- function(x, ...) {
#'   if(is_go_result(x)) {
#'     plot_func <- switch (class(x),
#'                          "groupGOResult" = group_go_plots,
#'                          "enrichResult"  = enrich_go_plots,
#'                          "gseaResult"    = gsea_go_plots)
#'     plot_func(x, ...)
#'   } else if(is(x, "list") && all(purrr::map(x, is_go_result) %>% unlist)) {
#'     # create wego plot
#'     p_wego <- go_wego_plot(x)
#'     # plots for each ont
#'     p_ont  <- lapply(x, function(i) {
#'       if(is_go_result(i)) {
#'         get_go_plots(i, ...)
#'       }
#'     })
#'     names(p_ont) <- names(x)
#'     list(
#'       wego = p_wego,
#'       ont  = p_ont
#'     )
#'   } else {
#'     warning("`x` not GO results, or list of GO results")
#'     NULL
#'   }
#' }
#'
#'
#'


# gsea_go <- function(gene_list, organism,  ...) {
#   #--Default values: BEGIN
#   dots <- list(
#     orgdb       = NULL,
#     keytype     = NULL,
#     pval_cutoff = 0.9
#   )
#   dots <- purrr::list_modify(dots, !!!list(...))
#   #--Default values: END
#   if(is(gene_list, "numeric") && is(names(gene_list), "character")) {
#     # check orgdb
#     if(is(dots$orgdb, "OrgDb")) {
#       dots$orgdb <- dots$orgdb
#     } else if(is(organism, "character")) {
#       dots$orgdb <- get_orgdb(organism)
#     }
#     if(! is(dots$orgdb, "OrgDb")) {
#       stop("`orgdb` and `organism`, not valid OrgDb")
#     }
#     # check keytype
#     if(is.null(dots$keytype)) {
#       organism <- get_orgdb_metadata(dots$orgdb, "ORGANISM")
#       dots$keytype  <- guess_keytype(names(gene_list), organism)
#     }
#     if(! is_valid_keytype(dots$keytype, dots$orgdb)) {
#       kt_list <- paste(AnnotationDbi::keytypes(dots$orgdb), collapse = ", ")
#       msg     <- glue::glue(
#         "`keytype` not valid: [{dots$keytype}], ",
#         "choose from: {kt_list}")
#       stop(msg)
#     }
#     # check input: sorted numeric, named
#     # run clusterProfiler::gseaGO
#     onts <- c("BP", "CC", "MF")
#     sapply(onts, function(ont) {
#       message(glue::glue("Running gseGO() for {ont}"))
#       clusterProfiler::gseGO(gene         = gene_list,
#                              keyType      = dots$keytype,
#                              OrgDb        = dots$orgdb,
#                              ont          = ont,
#                              # nPerm        = 1000,
#                              minGSSize    = 120,
#                              # maxGSSize    = 500,
#                              pvalueCutoff = dots$pval_cutoff,
#                              verbose      = FALSE)
#     }, USE.NAMES = TRUE)
#
#   } else {
#     stop("`x` not named, numeric values")
#   }
# }
#




