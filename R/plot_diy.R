#' plot_diy
#'
#'


#' load te ids from dm6
#'
#' @export
load_te_id <- function(genome = "dm6") {
  f1 <- system.file("te_dm6.rds", package = "hiseqr")
  if(file.exists(f1)) {
    readRDS(f1)
  }
}

#' load te tissue specific annotation
#'
#' Data from this publication:
#' 1.Malone, C. D. et al. Specialized piRNA pathways act in germline and somatic tissues of the Drosophila ovary. Cell vol. 137 522–535 (2009). DOI: [10.1016/j.cell.2009.03.040](https://doi.org/10.1016/j.cell.2009.03.040)
#'
#' @export
load_te_tissue <- function(...) {
  f1 <- system.file("te_dm6_tissue.csv", package = "hiseqr")
  if(file.exists(f1)) {
    read.csv(f1)
  }
}




#' check input is valid rnaseq_rx dir
#'
#' support:
#' _rx
#' _rx/deseq
#'
.is_valid_rnaseq_rx <- function(x) {
  if(inherits(x, "character")) {
    sapply(x, is_hiseq_dir, hiseq_type = "rnaseq_rx")
  } else {
    FALSE
  }
}


#' plot_scatter_te
#'
#' plot RNAseq for TE
#'
#' @export
plot_scatter_te <- function(x, ...) {
  # te, piRC, gene
  args <- list(gene_group = "te", color_by = "sig") # default
  args <- purrr::list_modify(args, ...) # updated
  gene_group <- args$gene_group # global
  color_by   <- args$color_by
  # checka rags
  if(! .is_valid_rnaseq_rx(x)) {
    message(glue::glue("Not a valid rnaseq_rx dir: {x}"))
    return(NULL)
  }
  #-- 1. only support one directory
  if(length(x) > 1) {
    x <- x[1]
    message(glue::glue("More than 1 dir found, choose the first one"))
  }
  #-- loading te ids
  genome <- list_hiseq_file(x, "genome")
  te_id <- load_te_id(genome)
  te_id <- te_id[[gene_group]] # check-point
  #-- loading deseq_data
  deseq_dir <- list_hiseq_file(x, "deseq_dir")
  #-- load deseq2 table
  df1 <- deseq_qc_res(deseq_dir)
  #-- 2. check if te exists, subset
  df2 <- df1[df1$gene_id %in% te_id, ]
  if(nrow(df2) > 0) {
    message(glue::glue("Found {nrow(df2)} of {length(te_id)} TE/genes"))
  } else {
    warning("No TE found")
    return(NULL)
  }
  #-- 3. plotting
  g_siga <- dplyr::filter(df2, ! sig == "not") %>% dplyr::pull(SYMBOL)
  gs     <- table(df2$sig)
  te_tissue <- load_te_tissue() # dm6
  names(te_tissue) <- c("gene_id", "tissue")
  df2 <- dplyr::left_join(df2, te_tissue, by = "gene_id")
  #-- basic
  p1 <- deseq_qc_scatter2(df2, color_by = color_by, add_sig = TRUE) +
    ggtitle(glue::glue(
      "TE (up={gs['up']}, down={gs['down']}, not={gs['not']})"
    )) +
    theme(plot.title = element_text(size = 10))
  #-- return
  p1
}
















