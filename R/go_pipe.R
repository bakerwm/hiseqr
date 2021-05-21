# name: go_pipe
#
# run go analysis, for specific input: deseq_dir |  gene_list
#' run go analysis
#'
#' @param x string, deseq_dir or gene names
#' @param organism string scientific name, eg: Drosophila melanogaster
#' @param outdir string, where to save the results
#' @param fold_change numeric, with gene_names,
#' @param feature string, required for deseq_dir
#'
#' @export
go_pipe <- function(x, organism, outdir,
                    fold_change = NULL,
                    feature = "gene") {
  if(is.character(x)) {
    x <- x[1]
  } else {
    stop("'x' require character, path to RNAseqRx dir")
  }
  if(! is_valid_organism(organism)) {
    stop(paste0("unknown organism: ", organism))
  }
  x_type <- get_hiseq_type(x)
  if(x_type == "hiseq_rx") {
    message("Run GO analysis ...")
    hiseqr::run_go(x, organism, outdir, fold_change)
    hiseqr::run_kegg(x, organism, outdir, fold_change)
    go_report(outdir)
  } else if(x_type == "hiseq_rx") {
    message("# Run GO analysis for DESeq dir...")
    x <- normalizePath(x)
    args <- read_hiseq(x)$args
    # args <- deseq_single_dir(x, feature)
    sig_gene <- get_sig_gene(x, feature)
    for(sig in c("up", "down", "up_and_down")){
      message(paste0("# GO and KEGG analysis for genes: ", sig))
      sig_outdir <- file.path(args$enrichdir, sig)
      if(sig == "up_and_down") {
        sig_gene_df <- rbind(sig_gene[["up"]], sig_gene[["down"]])
      } else {
        sig_gene_df <- sig_gene[[sig]]
      }

      if(is.null(sig_gene_df)) {
        warning("significant genes not found, GO skipped ...")
        return(NULL)
      }

      # prepare args
      gene_list   <- sig_gene_df$Gene
      fold_change <- setNames(sig_gene_df$log2fold_change, gene_list)
      fold_change <- sort(fold_change, decreasing = TRUE)
      # organism    <- search_species(args$genome) # scientific_name
      pval_cutoff <- 0.1
      qval_cutoff <- 0.1

      # GO
      tmp1 <- run_go(gene_list, organism, sig_outdir,
                     fold_change = fold_change,
                     level       = 2,
                     readable    = TRUE,
                     pval_cutoff = pval_cutoff,
                     qval_cutoff = qval_cutoff)

      # KEGG
      tmp2 <- run_kegg(gene_list, organism, sig_outdir,
                       fold_change = fold_change,
                       pval_cutoff = pval_cutoff,
                       qval_cutoff = qval_cutoff)

      ##---------------------##
      # GO report for gene list
      go_report(sig_outdir, feature)
    }

    ##---------------------##
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
  x_type <- read_hiseq(x[1])
  outdir <- normalizePath(x[1])

  if(is_enrich_dir(x)) {
    message("## Run GO report for single dir...")

    ##---------------------##
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




