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



# to-do
is_valid_kegg_input <- function(...) { x = 1 }


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
run_kegg <- function(gene_list, organism, outdir = NULL, ...) {
  outdir <- ifelse(is(outdir, "character"), outdir, getwd()) # working dir
  dots   <- list(...)
  arg_vars <- list(
    keytype     = NULL,
    fold_change = NULL,
    overwrite   = FALSE,
    is_kegg_id  = FALSE,
    readable    = TRUE
  )
  arg_vars <- purrr::list_modify(arg_vars, !!!dots)
  arg_vars <- purrr::list_modify(arg_vars,
                                 gene_list = gene_list,
                                 organism  = organism,
                                 outdir    = outdir) #@@@ arg_vars updated
  # return: organism, orgdb, keytype, gsea_gene, kegg_code
  kegg_input <- do.call(prep_kegg_input, arg_vars)
  if(is(kegg_input, "list")) {
    arg_vars <- purrr::list_modify(arg_vars, !!!kegg_input) #@@@ arg_vars, updated
    # arg_vars$gene_list <- arg_vars$kegg_gene
    if(length(arg_vars$kegg_gene) == 0) {
      stop("`gene_list`, `kegg_gene` is empty, not valid kegg_id")
    }
    do.call(enrich_kegg, arg_vars)
    do.call(gsea_kegg, arg_vars)
  }
}




#' @describeIn enrich_kegg enrichKEGG analysis from clusterprofiler
#'
#' for clusterProfiler::enrichKEGG
#'
#' gene_list, the could be: kegg, ncbi-geneid, ncbi-proteinid or uniprot
#'
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
#' @description pval_cutoff, qval_cutoff, set to 0.9, in order to return
#' enrich results anyway, filter in downstream analysis
#'
#' @import clusterProfiler
#' @import clusterProfiler.dplyr
#' @import stringr
#' @import cowplot
#'
#' @example enrich_kegg(gene = , organism = , keytype = ,
#' readable = TRUE, pval_cutoff = 0.05, qval_cutoff = 0.05, ...)
#'
#' @export
enrich_kegg <- function(gene_list, organism, outdir = NULL, ...) {
  outdir <- ifelse(is(outdir, "character"), outdir, getwd()) # working dir
  dots   <- rlang::list2(...)
  arg_vars <- list(
    keytype   = NULL,
    overwrite = FALSE,
    readable  = TRUE,
    pval_cutoff = 0.9,
    qval_cutoff = 0.9
  )
  arg_vars <- purrr::list_modify(arg_vars, !!!dots)
  arg_vars <- purrr::list_modify(arg_vars,
                                 gene_list = gene_list,
                                 organism  = organism,
                                 outdir    = outdir)
  # outdir, enrich
  arg_vars$outdir <- file.path(arg_vars$outdir, "enrich_kegg") # update: outdir
  if(! dir.exists(arg_vars$outdir)) {
    dir.create(arg_vars$outdir, recursive = TRUE, mode = "0755")
  }
  #--KEGG analysis
  message(glue::glue("Running enrich_kegg()..."))
  kegg_data_rds <- file.path(arg_vars$outdir, paste("enrich_kegg", "data", "rds", sep = "."))
  kegg_plot_rds <- file.path(arg_vars$outdir, paste("enrich_kegg", "plot", "rds", sep = "."))
  if(file.exists(kegg_data_rds) & ! arg_vars$overwrite) {
    kegg_data <- readRDS(kegg_data_rds)
  } else {
    kegg_data <- tryCatch(
      {
        kk <- clusterProfiler::enrichKEGG(
          gene          = arg_vars$kegg_gene,
          organism      = arg_vars$kegg_code,
          keyType       = arg_vars$kegg_keyType,
          pAdjustMethod = "BH",
          pvalueCutoff  = arg_vars$pval_cutoff,
          qvalueCutoff  = arg_vars$qval_cutoff)
        # readable
        arg_vars$orgdb <- get_orgdb(arg_vars$organism)
        if(is_go_result(kk)) {
          kegg_data <- clusterProfiler::setReadable(kk, arg_vars$orgdb,
                                                    arg_vars$keytype)
        } else {
          kegg_data <- kk
        }
        #--Save to rds
        saveRDS(kegg_data, file = kegg_data_rds)
        return(kegg_data)
      },
      error=function(cond) {
        warning("enrichKEGG() failed")
        return(NULL)
      }
    )
  }
  #--KEGG plotting: plot
  if(file.exists(kegg_plot_rds) & ! arg_vars$overwrite) {
    kegg_plot <- readRDS(kegg_plot_rds)
  } else {
    if(is_go_result(kegg_data)) {
      kegg_plot <- enrich_kegg_plots(kegg_data) # !!!! enrich_kegg_plot
      saveRDS(kegg_plot, file = kegg_plot_rds)
    } else {
      kegg_plot <- NULL
    }
  }
  #--KEGG plotting: save to png
  prefix <- gsub(".rds$", "", basename(kegg_plot_rds))
  save_go_plot(kegg_plot, arg_vars$outdir, prefix)
  #--KEGG table: save result to csv
  prefix <- gsub(".rds$", "", basename(kegg_data_rds))
  save_go_table(kegg_data, arg_vars$outdir, prefix)
  #--Return: data
  kegg_data
}








#' run gseKEGG()
#' @param gene_list decreasing sorted numeric vector
#' @param organism kegg_code, eg: dme
#'
#' @export
gsea_kegg <- function(gene_list, organism, outdir = NULL, ...) {
  if(missing(gene_list)) {
    warning("`gene_list` missing, `gsea_kegg()` skipped")
    return(NULL)
  }
  outdir <- ifelse(is(outdir, "character"), outdir, getwd()) # working dir
  dots   <- rlang::list2(...)
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
  arg_vars$outdir <- file.path(arg_vars$outdir, "gsea_kegg") # update: outdir
  if(! dir.exists(arg_vars$outdir)) {
    dir.create(arg_vars$outdir, recursive = TRUE, mode = "0755")
  }
  #--Run
  message(glue::glue("Running gsea_kegg()..."))
  kegg_data_rds <- file.path(arg_vars$outdir, paste("gsea_kegg", "data", "rds", sep = "."))
  kegg_plot_rds <- file.path(arg_vars$outdir, paste("gsea_kegg", "plot", "rds", sep = "."))
  if(file.exists(kegg_data_rds) & ! arg_vars$overwrite) {
    kegg_data <- readRDS(kegg_data_rds)
  } else {
    kegg_data <- tryCatch(
      {
        kk <- clusterProfiler::gseKEGG(
          geneList      = arg_vars$kegg_gsea_gene,
          organism      = arg_vars$kegg_code,
          keyType       = arg_vars$kegg_keyType,
          minGSSize     = 120,
          pvalueCutoff   = arg_vars$pval_cutoff,
          pAdjustMethod = "BH",
          verbose       = FALSE)
        # readable
        arg_vars$orgdb <- get_orgdb(arg_vars$organism)
        kegg_data <- clusterProfiler::setReadable(kk, arg_vars$orgdb,
                                                  arg_vars$keytype)
        #--Save to rds
        saveRDS(kegg_data, file = kegg_data_rds)
        return(kegg_data)
      },
      error=function(cond) {
        warning("enrichKEGG() failed")
        return(NULL)
      }
    )
  }
  #--KEGG plotting: plot
  if(file.exists(kegg_plot_rds) & ! arg_vars$overwrite) {
    kegg_plot <- readRDS(kegg_plot_rds)
  } else {
    if(is_go_result(kegg_data)) {
      kegg_plot <- gsea_kegg_plots(kegg_data) # !!!! enrich_kegg_plot
      saveRDS(kegg_plot, file = kegg_plot_rds)
    } else {
      kegg_plot <- NULL
    }
  }
  #--KEGG plotting: save to png
  prefix <- gsub(".rds$", "", basename(kegg_plot_rds))
  save_go_plot(kegg_plot, arg_vars$outdir, prefix)
  #--KEGG table: save result to csv
  prefix <- gsub(".rds$", "", basename(kegg_data_rds))
  save_go_table(kegg_data, arg_vars$outdir, prefix)
  #--Return: data
  kegg_data
}




#' @export
enrich_kegg_plots <- function(...) {
  enrich_go_plots(...)
}


#' @export
gsea_kegg_plots <- function(...) {
  gsea_go_plots(...)
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

