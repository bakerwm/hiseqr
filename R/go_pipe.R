# name: go_pipe
#
# run go analysis, for specific input: deseq_dir |  gene_list
#' run go analysis
#'
#' @param x string, deseq_dir or gene names
#' @param organism string scientific name, eg: Drosophila melanogaster
#' @param outdir string, where to save the results
#' @param foldChange numeric, with gene_names,
#' @param feature string, required for deseq_dir
#'
#' @export
go_pipe <- function(x, organism, outdir,
                    foldChange = NULL,
                    feature = "gene") {
  x_type <- is_hiseq_dir(x[1])

  if(is.null(x_type)) {
    message("## Run GO analysis for gene list...")

    # prepare genes
    if(! is_valid_organism(organism)) {
      stop(paste0("unknown organism: ", organism))
    }

    # run
    hiseqr::run_go(x, organism, outdir, foldChange)
    hiseqr::run_kegg(x, organism, outdir, foldChange)

    # report
    go_report(outdir)

  } else if(x_type == "deseq_single") {
    message("# Run GO analysis for DESeq single dir...")
    x <- normalizePath(x)

    ## config
    args <- deseq_single_dir(x, feature)

    ## number of sig genes
    sig_gene <- get_sig_gene(x, feature)

    ## up/down/not/up+down
    for(s in c("up", "down", "up_and_down")){
      message(paste0("# GO and KEGG analysis for genes: ", s))
      outdir2 <- file.path(args$enrichdir, s)

      # sig genes
      if(s == "up_and_down") {
        sig_gene_df <- rbind(sig_gene[["up"]], sig_gene[["down"]])
      } else {
        sig_gene_df <- sig_gene[[s]]
      }

      if(is.null(sig_gene_df)) {
        warning("significant genes not found, GO skipped ...")
        return(NULL)
      }

      # prepare args
      gene_list    <- sig_gene_df$Gene
      foldChange   <- setNames(sig_gene_df$log2FoldChange, gene_list)
      foldChange   <- sort(foldChange, decreasing = TRUE)
      organism     <- search_species(args$genome) # scientific_name
      pvalueCutoff <- 0.1
      qvalueCutoff <- 0.1

      # GO
      tmp1 <- run_go(gene_list, organism, outdir2,
                     foldChange   = foldChange,
                     level        = 2,
                     readable     = TRUE,
                     pvalueCutoff = pvalueCutoff,
                     qvalueCutoff = qvalueCutoff,
                     to_entrezid  = TRUE)

      # KEGG
      tmp2 <- run_kegg(gene_list, organism, outdir2,
                       foldChange   = foldChange,
                       pvalueCutoff = pvalueCutoff,
                       qvalueCutoff = qvalueCutoff)

      ##---------------------##
      # report
      # GO report for gene list
      go_report(outdir2, feature)
    }

    ##---------------------##
    # report
    go_report(x, feature)
  }
}



#' generate report for directory
#'
#' @param x path to GO output or deseq_dir
#' @param feature string support for deseq_dir, default: gene
#'
#' @export
go_report <- function(x, feature = "gene") {
  x_type <- is_hiseq_dir(x[1])
  outdir <- normalizePath(x)

  if(is_enrich_dir(x)) {
    message("## Run GO report for single dir...")

    ##---------------------##
    # report
    template <- system.file("rnaseq", "go_enrich_report_single.Rmd",
                            package = "hiseqr")
    template_to <- file.path(outdir, basename(template))
    ## copy Rmd
    file.copy(template, template_to, overwrite = TRUE)
    out_html <- file.path(outdir, "go_enrich_report_single.html")

    ## run Rmd
    rmarkdown::render(input       = template,
                      output_file = out_html,
                      params      = list(input_dir  = outdir))

  } else if(is.null(x_type)) {
    warning(paste0("unknown input x: ", x))
  } else if(x_type == "deseq_single") {
    message("# Run GO analysis for DESeq dir...")

    ## args
    args <- deseq_single_dir(x, feature)

    ##---------------------##
    # report
    template <- system.file("rnaseq", "go_enrich_report.Rmd",
                            package = "hiseqr")
    template_to <- file.path(args$report_dir, basename(template))
    ## copy Rmd
    file.copy(template, template_to)
    out_html <- file.path(args$report_dir,
                          paste("go_enrich", "report.html", sep = "."))

    ## run Rmd
    rmarkdown::render(input       = template,
                      output_file = out_html,
                      params      = list(deseq_dir  = x,
                                         feature    = feature))
  } else {
    warning(paste0("unknown input: ", x[1]))
  }
}




