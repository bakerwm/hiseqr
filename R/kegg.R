
#' organize KEGG output
#' generate tables, plots
#' @param gene_list vector of gene list
#' @param organism name of genome, eg: Homo sapiens
#' @param outdir string
#' @param foldChange numeric, with names
#'
#' @export
run_kegg <- function(gene_list, organism, outdir,
                     foldChange   = NULL,
                     pvalueCutoff = 0.5,
                     qvalueCutoff = 0.5) {
  # outdir
  if(! dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  ##------------------------------##
  ## parse args
  # force input, ENTREZID #
  input_gene <- go_input(gene_list, organism,
                         foldChange,  # for enrich only
                         to_entrezid = TRUE)
  orgdb      <- input_gene$orgdb # global
  gene_list  <- input_gene$gene
  keytype    <- input_gene$keytype
  foldChange <- input_gene$foldChange
  gsea_gene  <- input_gene$gsea_gene
  kegg_code  <- input_gene$kegg_code

  ##--------------------------------------------------------------------------##
  ## run KEGG analysis
  kegg_rds <- file.path(outdir, "kegg_data.rds")
  if(file.exists(kegg_rds)) {
    kegg_obj <- readRDS(kegg_rds)
  } else {
    # run
    kegg_obj <- list(enrich = list(kegg = kegg_enrich(gene_list, kegg_code,
                                                      pvalueCutoff, qvalueCutoff)),
                     gsea   = list(kegg = kegg_gsea(gene_list, kegg_code, foldChange)))
    # save to obj
    saveRDS(kegg_obj, file = kegg_rds)
  }

  ##--------------------------------------------------------------------------##
  ## Generate KEGG plots
  ## 1. enrich
  enrich_dir <- file.path(outdir, "kegg_enrich")
  table1     <- save_table(kegg_obj$enrich, enrich_dir)
  plot1      <- get_go_plots(kegg_obj$enrich, enrich_dir, foldChange)
  save_plot(plot1, enrich_dir)

  ## 2. GSEA
  gsea_dir <- file.path(outdir, "kegg_gsea")
  table2   <- save_table(kegg_obj$gsea, gsea_dir)
  plot2    <- get_go_plots(kegg_obj$gsea, gsea_dir, foldChange)
  save_plot(plot2, gsea_dir)

  # save obj
  plot_obj <- list(enrich = plot1,
                   gsea   = plot2)
  plot_rds <- file.path(outdir, "kegg_plots.rds")
  saveRDS(plot_obj, file = plot_rds)
}


#' run enrichKEGG()
#' @param gene_list gene names
#' @param organism kegg_code, eg: dme
#'
#' @export
kegg_enrich <- function(gene_list, organism,
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05) {
  # choose keytype
  keytype <- ifelse(organism %in% c("dme"), "ncbi-geneid", "kegg")

  # enrich
  print("Running enrichKEGG")
  kk <- clusterProfiler::enrichKEGG(gene          = gene_list,
                                    organism      = organism,
                                    keyType       = keytype,
                                    pvalueCutoff  = pvalueCutoff,
                                    pAdjustMethod = "BH",
                                    qvalueCutoff  = qvalueCutoff)

  # readable
  # convert kegg_code to scientific_name
  sci_name <- clusterProfiler::search_kegg_organism(organism, by = "kegg_code")$scientific_name
  orgdb <- get_orgdb(sci_name)

  if(is(kk, "enrichResult")) {
    kk <- clusterProfiler::setReadable(kk, orgdb, "ENTREZID")
  }
  kk
}


#' run gseKEGG()
#' @param gene_list decreasing sorted numeric vector
#' @param organism kegg_code, eg: dme
#'
#' @export
kegg_gsea <- function(gene_list, organism, foldChange,
                      pvalueCutoff = 0.05) {
  # input
  gene_list_fc <- gsea_input(gene_list, foldChange)
  if(is.null(gene_list_fc)) {
    return(NULL)
  }

  # choose keytype
  keytype <- ifelse(organism %in% c("dme"), "ncbi-geneid", "kegg")

  # GSEA
  print("Running gseKEGG") #!!!!
  gsea <- clusterProfiler::gseKEGG(geneList      = gene_list_fc,
                                   organism      = organism,
                                   keyType       = keytype,
                                   minGSSize     = 120,
                                   pvalueCutoff  = pvalueCutoff,
                                   pAdjustMethod = "BH",
                                   verbose       = FALSE)

  # readable
  # convert kegg_code to scientific_name
  sci_name <- clusterProfiler::search_kegg_organism(organism, by = "kegg_code")$scientific_name
  orgdb <- get_orgdb(sci_name)

  if(is(gsea, "gseaResult")) {
    gsea <- clusterProfiler::setReadable(gsea,
                                         OrgDb   = orgdb,
                                         keyType = "ENTREZID")
  }
  gsea
}

#'
#' #' save barplot of go
#' #'
#' #' @param x list of objects of groupGO, enrichGO, KEGG, ...
#' #' @param outdir path, to save the results
#' #'
#' #' @export
#' get_kegg_plots <- function(x, outdir, foldChange = NULL) {
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
#'     plot_func(x, foldChange)
#'   } else if(is.list(x)) {
#'     # list of obj
#'     p_list <- lapply(x, function(i){
#'       get_kegg_plots(i, outdir, foldChange)
#'     })
#'     names(p_list) <- names(x)
#'     p_list # return
#'   }
#' }




