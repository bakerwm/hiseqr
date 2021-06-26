#' Functions for salmon output
#'
#' Apply DESeq2 for RNAseq DE analysis
#'
#' @import fs
#' @import ggplot2
#' @import ggthemes
#' @import dplyr
#' @import tidyr
#' @import readr
#' @import tibble
#' @import patchwork
#' @import DEseq2
#'
#'
#' @name salmon



#' @describeIn rnaseq_salmon_hub A port for RNAseq analysis, DE analysis
#'
#' @param x path to the RNAseqRx directory, a.vs.b
#' @param ... extra arguments for demseq2_main function#'
#'
#'
#' deseq2_main() extra arguments:
#' - fc_cutoff    : 2
#' - pval_cutoff  : 0.05
#' - readable     : TRUE
#'
#'
#' @param organism character, the name of the organism, eg: dm6, hg38, human,
#' @param readable bool, Add symbol and entrezid to table
#' organism, fc_cutoff, pval_cutoff, readable,
#'
#' @import dplyr
#' @import DESeq2
#'
#' @export
rnaseq_salmon_hub <- function(x, ...) {
  # genome    <- list_hiseq_file(project_dir, "genome", "rx")
  # make sure: aligner: salmon
  if(is_hiseq_dir(x, 'rnaseq_rx')) {
    aligner <- list_hiseq_file(x, "aligner", "rx")
    if(aligner == "salmon") {
      message(glue::glue("RNAseq for salmon output: {x}"))
    } else {
      on.exit(glue::glue("not a salmon output, got {aligner}"))
    }
  } else {
    on.exit(glue::glue("not a rnaseq_rx directory, {x}"))
  }
  # arguments
  wt_quant  <- list_hiseq_file(x, "wt_quant", "rx") # px$wt_quant
  wt_name   <- list_hiseq_file(x, "wt_name", "rx") # pd$args$wt_name
  mut_quant <- list_hiseq_file(x, "mut_quant", "rx") # px$mut_quant
  mut_name  <- list_hiseq_file(x, "mut_name", "rx") # pd$args$mut_name
  deseq_dir <- list_hiseq_file(x, "deseq_dir", "rx")
  genome    <- list_hiseq_file(x, "genome", "rx")
  fix_xls   <- list_hiseq_file(x, "deseq_fix_xls", "rx")
  # check tx2gene
  salmon_index <- list_hiseq_file(x, "salmon_index", "rx")
  if(is(salmon_index, "character")) {
    tx2gene_csv <- file.path(salmon_index, "tx2gene.csv")
    if(! file.exists(tx2gene_csv)) {
      on.exit(glue::glue("tx2gene.csv not found: {tx2gene_csv}"))
    }
  } else {
    on.exit(glue::glue("salmon_index not found, {x}"))
  }
  # prepare data
  sf_list <- c(wt_quant, mut_quant)
  names(sf_list) <- basename(dirname(dirname(sf_list)))
  # design
  samples <- data.frame(
    name      = names(sf_list),
    condition = c(rep(wt_name,  length = length(wt_quant)),
                  rep(mut_name, length = length(mut_quant)))
  ) %>%
    dplyr::mutate(condition = factor(
      condition, levels = c(wt_name, mut_name)))
  # experiment design
  tx2gene <- read.csv(tx2gene_csv)
  # load data
  txi <- tximport::tximport(sf_list, type = "salmon", tx2gene = tx2gene)
  # DEseq2 analysis
  dds_txi <- DESeq2::DESeqDataSetFromTximport(
    txi, colData = samples, design = ~condition
  )
  # run
  hiseqr::deseq2_main(dds_txi, deseq_dir, organism = genome)
  hiseqr::make_publish_plots(x)
  if(run_go == 1) {
    hiseqr::rnaseq_enrich_hub(x, organism = genome)
  }
  rnaseq_report(x, file.path(x, "report"))
}



