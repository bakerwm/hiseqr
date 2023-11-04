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
#' 1.Malone, C. D. et al. Specialized piRNA pathways act in germline and
#' somatic tissues of the Drosophila ovary. Cell vol. 137 522–535 (2009).
#' DOI: [10.1016/j.cell.2009.03.040](https://doi.org/10.1016/j.cell.2009.03.040)
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
#' @export
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
  # p1 <- deseq_qc_scatter2(df2, color_by = color_by, add_sig = TRUE, ...) +
  args$add_sig <- NULL # remove add_sig from args
  p1 <- deseq_qc_scatter2(df2, add_sig = TRUE, !!!args) +
    ggtitle(glue::glue(
      "TE (up={gs['up']}, down={gs['down']}, not={gs['not']})"
    )) +
    theme(plot.title = element_text(size = 10))
  #-- return
  p1
}





# # heatmap
# library(hiseqr)
# library(dplyr)
# library(ggplot2)
# library(SummarizedExperiment)
#
#
# x <- "~/work/yu_2023/projects/20230309_dlj_ChrRNA_yy249/results/RNAseq_salmon/ChrRNA_ovary_3127_DMSO_24hr.vs.ChrRNA_ovary_3127_5ad_IAA_24hr/"
# genome <- hiseqr::list_hiseq_file(x, "genome")
# gene_group <- "te"
# te_id <- load_te_id(genome)
# te_id <- te_id[[gene_group]] # check-point
# # deseq_dir <- list_hiseq_file(x, "deseq_dir")
# # df1 <- deseq_qc_res(deseq_dir)
# # df2 <- df1[df1$gene_id %in% te_id, ]
# #
# # df3 <- dplyr::select(df2, 1:3) %>%
# #   tibble::column_to_rownames("gene_id") %>%
# #   as.matrix()
# idx <- list_hiseq_file(x, "salmon_index")
# t2g_csv <- file.path(idx, "tx2gene.csv")
# t2g <- readr::read_csv(t2g_csv, show_col_types = FALSE)
#
# f <- list_hiseq_file(x, "wt_dirs")
# sf <- list.files(f, "quant.sf", recursive = TRUE, full.names = TRUE)
# txi <- tximport::tximport(sf[1], type = "salmon", tx2gene = t2g, abundanceCol = "TPM")
# df1 <- cbind(as.data.frame(txi$length),
#              as.data.frame(txi$abundance),
#              as.data.frame(txi$counts))
# df1 <- round(df1, 4)
# quant_dir <- list_hiseq_file(x, "quant_dir")


## WT dirs
read_salmon_rn <- function(x) {
  # locate the info.json
  idx     <- list_hiseq_file(x, "index_list")
  t2g_csv <- file.path(idx[1], "tx2gene.csv")
  t2g     <- readr::read_csv(t2g_csv, show_col_types = FALSE)
  # suppressPackageStartupMessages(library(tximport))
  sf  <- list.files(x, "quant.sf", full.names = TRUE, recursive = TRUE)
  txi <- tximport::tximport(sf, type = "salmon", tx2gene = t2g,
                            abundanceCol = "TPM")
  # gene, length, count, tpm
  df1 <- cbind(as.data.frame(txi$length),
               as.data.frame(txi$abundance),
               as.data.frame(txi$counts))
  df1 <- round(df1, 4)
  colnames(df1) <- c("length", "TPM", "count")
  df1$smp_name <- list_hiseq_file(x, "smp_name")
  tibble::rownames_to_column(df1, "gene_id")
  # df1
}


read_salmon_rx <- function(x) {
  # dirs
  wt_name  <- list_hiseq_file(x, "wt_name")
  wt_dirs  <- list_hiseq_file(x, "wt_dirs")
  mut_name <- list_hiseq_file(x, "mut_name")
  mut_dirs <- list_hiseq_file(x, "mut_dirs")
  # read quant.sf
  df1 <- lapply(c(wt_dirs, mut_dirs), read_salmon_rn) %>%
    dplyr::bind_rows() %>%
    dplyr::select(gene_id, smp_name, TPM) %>%
    tidyr::pivot_wider(names_from = "smp_name", values_from = "TPM")
  # mean
  data.frame(
    gene_id  = df1$gene_id,
    wt  = dplyr::select(df1, starts_with(wt_name)) %>% rowMeans(na.rm = TRUE),
    mut = dplyr::select(df1, starts_with(mut_name)) %>% rowMeans(na.rm = TRUE)
  ) %>%
    dplyr::rename(
      !!wt_name  := wt,
      !!mut_name := mut
    )
}




#' plot heatmap for DESeq TE
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_heatmap_te <- function(x, ...) {
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
  te_id  <- load_te_id(genome)
  te_id  <- te_id[[gene_group]] # check-point
  #-- TPM
  df1 <- read_salmon_rx(x)
  #-- subset
  df2 <- df1[df1$gene_id %in% te_id, ]
  if(nrow(df2) > 0) {
    message(glue::glue("Found {nrow(df2)} of {length(te_id)} TE/genes"))
  } else {
    warning("No TE found")
    return(NULL)
  }
  #-- tissue specific
  # g_siga <- dplyr::filter(df2, ! sig == "not") %>% dplyr::pull(SYMBOL)
  # gs     <- table(df2$sig)
  te_tissue <- load_te_tissue() # dm6
  te_germ   <- setNames(te_tissue$tissue, nm = te_tissue$id)
  df2 <- dplyr::mutate(df2, tissue = te_germ[gene_id])
  #-- simplify names
  rank_smp <- c(list_hiseq_file(x, "wt_name"), list_hiseq_file(x, "mut_name"))
  rank_smp2 <- deseq_sanitize_str(rank_smp, n = 10)
  smp <- setNames(rank_smp2, nm = rank_smp)
  colnames(df2) <- c("gene_id", rank_smp2, "tissue")
  # #-- data
  # ma  <- df2 %>%
  #   dplyr::select(-tissue) %>%
  #   tibble::remove_rownames() %>%
  #   tibble::column_to_rownames("gene_id") %>%
  #   as.matrix()
  # ma <- log10(ma + .5)
  #-- heatmap
  #-- sample names
  # rank data.frame
  rank_col <- colnames(df2)[2]
  rank_id <- dplyr::arrange(df2, !!as.name(rank_col)) %>%
    dplyr::pull(gene_id)
  df3 <- df2 %>%
    tidyr::pivot_longer(
      cols = 2:3,
      names_to = "sample", values_to = "TPM"
    ) %>%
    dplyr::mutate(
      TPM = log10(TPM + .5),
      gene_id = factor(gene_id, levels = rank_id),
      sample  = factor(sample, levels = rank_smp2)
    )
  #-- heatmap, ggplot2::geom_tile()
  p1 <- ggplot(df3, aes(sample, gene_id)) +
    geom_tile(aes(fill = TPM), color = "white") +
    scale_fill_gradient(low = "white", high = "#3556ad") +
    scale_x_discrete(position = "top", expand = expansion(mult = c(0, .5))) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.text.y = element_blank()
    )
  cc <- c("germline" = "#ff2e16", "intermediate" = "#ffeb3d",
          "other" = "#85878b", "soma"  = "#2ca748")
  p_anno <- ggplot(df3, aes(x = 1, y = gene_id)) +
    geom_point(aes(color = tissue), size = 3) +
    scale_x_continuous(
      limits = c(0.9, 1.1),
      position = "top", expand = expansion(mult = c(0, 0))
    ) +
    scale_color_manual(values = cc) +
    theme(
      # legend.position = "None",
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      # axis.line.x = element_line(color = "blue"),
      axis.text = element_text(color = "black")
    )
  p <- patchwork::wrap_plots(p_anno, p1, nrow = 1, guides = "collect") +
    patchwork::plot_layout(widths = c(1, 6))
  #-- save to file
  deseq_dir <- list_hiseq_file(x, "deseq_dir")
  file_ht   <- file.path(deseq_dir, "fig3.publish.heatmap.te.png")
  ggsave(file_ht, p, width = 4, height = 15, units = "in")
}


