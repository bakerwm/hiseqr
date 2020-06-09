
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
#' @export
run_go <- function(gene_list, organism, outdir,
                   foldChange = NULL,
                   level = 2,
                   readable = TRUE,
                   pvalueCutoff = 0.5,
                   qvalueCutoff = 0.5,
                   to_entrezid = TRUE) {
  # outdir
  if(! dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  ##------------------------------##
  ## parse args
  # force input, ENTREZID #
  input_gene <- go_input(gene_list, organism, foldChange, to_entrezid)
  orgdb      <- input_gene$orgdb
  gene_list  <- input_gene$gene
  keytype    <- input_gene$keytype
  foldChange <- input_gene$foldChange
  gsea_gene  <- input_gene$gsea_gene

  ##--------------------------------------------------------------------------##
  ## run GO analysis
  go_rds <- file.path(outdir, "go_data.rds")
  if(file.exists(go_rds)) {
    message("Loading previous compute GO analysis results")
    go_obj <- readRDS(go_rds)
  } else {
    # run
    go_obj <- list(group = go_group(gene_list, orgdb, "ENTREZID",
                                    level = level, readable = readable),
                   enrich = go_enrich(gene_list, orgdb, "ENTREZID",
                                      pvalueCutoff  = pvalueCutoff,
                                      qvalueCutoff  = qvalueCutoff,
                                      readable      = readable),
                   gsea = go_gsea(gene_list, orgdb, foldChange,
                                  pvalueCutoff = pvalueCutoff),
                   gene_list  = gene_list,
                   foldChange = foldChange,
                   keytype    = keytype,
                   organism   = organism)
    # save obj
    saveRDS(go_obj, file = go_rds)
  }

  ##--------------------------------------------------------------------------##
  ## Generate plots
  # 1 group
  group_dir <- file.path(outdir, "go_group")
  table1    <- save_table(go_obj$group, group_dir)
  plot1     <- get_go_plots(go_obj$group, group_dir, foldChange)
  save_plot(plot1, group_dir)

  # 2 enrich
  enrich_dir <- file.path(outdir, "go_enrich")
  table2     <- save_table(go_obj$enrich, enrich_dir)
  plot2      <- get_go_plots(go_obj$enrich, enrich_dir, foldChange)
  save_plot(plot2, enrich_dir)

  # 3. GSEA
  gsea_dir <- file.path(outdir, "go_gsea")
  table3   <- save_table(go_obj$gsea, gsea_dir)
  plot3    <- get_go_plots(go_obj$gsea, gsea_dir, foldChange)
  save_plot(plot3, gsea_dir)

  # save obj
  plot_obj <- list(group  = plot1,
                   enrich = plot2,
                   gsea   = plot3)
  plot_rds <- file.path(outdir, "go_plots.rds")
  saveRDS(plot_obj, file = plot_rds)
}


#' for clusterProfiler::groupGO, and enrichGO
#'
#' @import clusterProfiler
#' @import clusterProfiler.dplyr
#' @import stringr
#' @import cowplot
#'
#' @example go_group(gene =, OrgDb = , keyType = , level = , ont = , readable = )
#'
#'
#' @export
go_group <- function(gene_list, OrgDb,
                     keytype = "ENTREZID",
                     level = 2,
                     readable = FALSE) {
  # ont
  onts <- setNames(c("BP", "CC", "MF"), nm = c("BP", "CC", "MF"))
  # output
  lapply(onts, function(ont) {
    print(glue::glue("Running groupGO() for {ont}"))
    clusterProfiler::groupGO(gene     = gene_list,
                             OrgDb    = OrgDb,
                             keyType  = keytype,
                             ont      = ont,
                             level    = level,
                             readable = readable)
  })
}


#' go analysis
#'
#' for clusterProfiler::enrichGO
#'
#' @import clusterProfiler
#' @import clusterProfiler.dplyr
#' @import stringr
#' @import cowplot
#'
#' @example go_enrich(gene = , OrgDb = , keyType = , ont = ,
#' pvalueCutoff = , pAdjustMethod = , universe = ,
#' qvalueCutoff = , minGSSize = , maxGSSize = ,
#' readable = , pool = )
#'
#' @export
go_enrich <- function(gene_list, OrgDb,
                      keytype = "ENTREZID",
                      readable = FALSE,
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05) {
  onts <- setNames(c("BP", "CC", "MF"), nm = c("BP", "CC", "MF"))
  # output
  lapply(onts, function(ont) {
    print(glue::glue("Running enrichGO() for {ont}"))
    ego <- clusterProfiler::enrichGO(gene          = gene_list,
                                     OrgDb         = OrgDb,
                                     ont           = ont,
                                     keyType       = keytype,
                                     pvalueCutoff  = pvalueCutoff,
                                     qvalueCutoff  = qvalueCutoff,
                                     pAdjustMethod = "BH",
                                     readable      = readable)
    clusterProfiler::simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min) # redundant
  })
}


#' for clusterProfiler::groupGO, and enrichGO, gseGO
#'
#' @param gene_list numeric, decreasing ordered,
#' @param OrgDb
#' @pvalueCutoff float default 0.05
#'
#' @import clusterProfiler
#' @example gseGO(gene =, OrgDb = , ont = )
#'
#' @export
go_gsea <- function(gene_list, OrgDb,
                    foldChange,
                    pvalueCutoff = 0.05) {
  # input
  gene_list_fc <- gsea_input(gene_list, foldChange)
  if(is.null(gene_list_fc)) {
    return(NULL)
  }
  # ont
  onts <- setNames(c("BP", "CC", "MF"), nm = c("BP", "CC", "MF"))
  # output
  lapply(onts, function(ont) {
    print(glue::glue("Running gseGO() for {ont}"))
    clusterProfiler::gseGO(gene         = gene_list_fc,
                           OrgDb        = OrgDb,
                           ont          = ont,
                           minGSSize    = 120,
                           pvalueCutoff = pvalueCutoff)
  })
}


#' save barplot of go
#'
#' @param x list of objects of groupGO, enrichGO, KEGG, ...
#' @param outdir path, to save the results
#'
#' @export
get_go_plots <- function(x, outdir, foldChange = NULL) {
  # x should be go obj
  if(is_go(x, recursive = TRUE)) {
    message("Generating plots for GO analysis")
  } else {
    warning("input failed, groupGO, enrichGO expected")
    return(NULL)
  }

  # prepare plot_list
  if(is_go(x)) {
    # single object
    plot_func <- switch (class(x),
                         "groupGOResult" = go_group_plots,
                         "enrichResult"  = go_enrich_plots,
                         "gseaResult"    = go_gsea_plots)
    plot_func(x, foldChange)
  } else if(is.list(x)) {
    # list of obj
    # for wego plot
    p_list1 <- list(p_wego = go_wego_plot(x))
    # plots for each ont/
    p_list2 <- lapply(x, function(i){
      get_go_plots(i, outdir, foldChange)
    })
    names(p_list2) <- names(x)
    c(p_list1, p_list2) # return
  } else {
    message(paste0("GO output expected, ", class(x), " got"))
  }
}


#' create plots
#' @param x object of enrichGO
#'
#' @export
go_group_plots <- function(x, foldChange = NULL, text_width = 40) {
  if(class(x) == "groupGOResult") {
    message("Generating plots for GO group")
  } else {
    return(list())
  }

  # empty data
  if(nrow(x) == 0) {
    message("Empty groupGORsult result detected")
    return(list())
  }

  # output
  list(barplot = go_barplot(x))
}


#' create plots
#' @param x object of enrichGO
#'
#' @export
go_enrich_plots <- function(x, foldChange = NULL, text_width = 40) {
  if(class(x) == "enrichResult") {
    message("Generating plots for GO enrichment")
  } else {
    return(list())
  }

  # empty data
  if(nrow(x) == 0) {
    message("Empty enrichResult result detected")
    return(list())
  }

  # output
  list(barplot  = go_barplot(x),
       dotplot  = go_dotplot(x),
       cnetplot = go_cnetplot(x, foldChange),
       emapplot = go_emapplot(x),
       heatplot = go_heatplot(x, foldChange))
}


#' create plots
#' @param x object of gseaGO
#'
#' @export
go_gsea_plots <- function(x, foldChange) {
  if(class(x) == "gseaResult") {
    message("Generating plots for GO GSEA")
  } else {
    return(list())
  }

  # empty data
  if(nrow(x) == 0) {
    message("Empty gseaResult result detected")
    return(list())
  }

  list(gseaplot = go_gsea_plot(x))
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
save_plot <- function(x, outdir, name = NULL) {
  # output
  if(! dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE, mode = "0755")
  }

  # for plot
  if(is(x, "gg")) {
    # save to file
    if(is.character(name)) {
      f <- file.path(outdir, paste0(name[1], ".png"))
    } else {
      f <- tempfile("plot.", outdir, ".png")
    }
    # log
    message(paste0("Saving plot to file: ", basename(f)))
    ggplot2::ggsave(f, plot = x, width = 6, height = 6, dpi = 150)
  } else if(is.list(x)) {
    # level-1:
    tmp <- lapply(names(x), function(i){
      # update name
      prefix <- ifelse(is.null(name), i, paste(name, i, sep = ".")) #
      g <- x[[i]] # object
      save_plot(x[[i]], outdir, prefix)
    })
  } else {
    warning(paste0("ggplot2 or list of ggplot2 obj expected, ",
                   class(x), " got."))
  }
}


#' save obj to outdir
#'
#' @param x list of go plots
#' @param outdir string, path to save the plots
#' @param name prefix for the plots
#'
#' @return
save_table <- function(x, outdir, name = NULL) {
  # output
  if(! dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE, mode = "0755")
  }

  # save table
  if(is_go(x)) {
    # filename
    if(is.character(name)) {
      f <- file.path(outdir, paste0(name[1], ".csv"))
    } else {
      f <- tempfile("table.", outdir, ".csv")
    }
    # log
    message(paste0("Saving table to file: ", basename(f)))
    write.csv(x@result, f, quote = TRUE, row.names = FALSE)
  } else if(is.list(x)) {
    # level-1:
    tmp <- lapply(names(x), function(i){
      # update name
      prefix = ifelse(is.null(name), i, paste(name, i, sep = ".")) #
      g <- x[[i]] # object
      save_table(x[[i]], outdir, prefix)
    })
  } else {
    warning(paste0("GO results obj expected, ",
                   class(x), " got."))
  }
}



