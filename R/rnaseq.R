## functions for RNAseq

## required pkgs
## apeglm
## MASS

## ----------------------------------------------------------------------------##
## Run DESeq analysis using DESeq2
##

#' deseq_hub
#' design for hiseq package, directory structure
#'
#' @param x path to the directory of project directory: a.vs.b
#' @param feature character, gene|te|..., default, gene
#'
#' deseq_dir: prj_dir/feature/deseq
#' count_dir: prj_dir/feature/count
#'
deseq_hub <- function(x, feature = "gene", outdir = NULL) {
  hiseq_type <- is_hiseq_dir(x)
  if (!hiseq_type == "deseq_single") {
    stop(glue::glue("DESeq dir expected, {hseq_type} found, {x}"))
  }

  # call for count files
  pk_files <- list_arguments_file(x)
  args <- load_pickle(pk_files[[feature]])

  # reading data
  pdata <- count_to_matrix(
    count_ctl = args$count_ctl,
    count_exp = args$count_exp
  )

  # checkout outdir
  if(is.null(outdir)) {
    outdir <- x # same as :deseq
  }

  deseq_hub2(args$count_ctl, args$count_exp, outdir, args$genome)

  # # for DESeqDataSet
  # dds <- deseq_import_from_matrix(
  #   ma = pdata$ma,
  #   smp_control = pdata$smp_ctl,
  #   smp_treatment = pdata$smp_exp
  # )
  #
  # # add meta data to dds $basepairs
  # df2 <- read_fc(args$count_ctl[1], gene_length = TRUE)
  # fData <- data.frame(basepairs = df2$Length)
  # S4Vectors::mcols(dds) <- S4Vectors::DataFrame(S4Vectors::mcols(dds), fData)
  #
  # # for DEseq2
  # deseq_node(dds,
  #            outdir = args$deseqdir,
  #            pval_cutoff = 0.1,
  #            organism = args$genome,
  #            readable = TRUE
  # )
}


#' deseq_hub2
#' design for count.txt files
#'
#' @param count_ctl path to the count.txt for control
#' @param count_exp path to the count.txt for treatment
#' @param outdir path to the output directory of deseq
#' @param feature character, gene|te|..., default, gene
#'
#' deseq_dir: prj_dir/feature/deseq
#' count_dir: prj_dir/feature/count
#'
deseq_hub2 <- function(count_ctl, count_exp, outdir, genome) {
  # check
  stopifnot(all(file.exists(count_ctl)))
  stopifnot(all(file.exists(count_exp)))

  # import data
  pdata <- count_to_matrix(count_ctl, count_exp)

  # for DESeqDataSet
  dds <- deseq_import_from_matrix(
    ma = pdata$ma,
    smp_control = pdata$smp_ctl,
    smp_treatment = pdata$smp_exp
  )

  # add meta data to dds $basepairs
  df2 <- read_fc(count_ctl[1], gene_length = TRUE)
  fData <- data.frame(basepairs = df2$Length)
  S4Vectors::mcols(dds) <- S4Vectors::DataFrame(S4Vectors::mcols(dds), fData)

  # for DEseq2
  deseq_node(dds,
             outdir = outdir,
             pval_cutoff = 0.1,
             organism = genome,
             readable = TRUE)

  ## Publish quality figures
  cnt_fix1 <- file.path(outdir, "transcripts_deseq2.fix.xls")
  stopifnot(file.exists(cnt_fix1))
  message("Generating publishable plots")
  message(paste0("found DESeq2 ouptut: ", cnt_fix1))
  tmp1 <- DESeq2_publish_plot(cnt_fix1, outdir, save2pdf = TRUE)
}


#' run DESeq2 analysis
#' import: dds
#' output: basic analysis
#' gene-level
#'
#' @param dds count data in matrix
#' @param outdir Directory to save the results
#' @param pval_cutoff Cutoff to filt the records, padj for DESeq2 output,
#'   default: 0.1
#'
#' @import DESeq2
#' @import apeglm
#' @import ggplot2
#' @import pheatmap
#' @import RColorBrewer
#' @import SummarizedExperiment
#'
#' @export
deseq_node <- function(dds,
                       outdir,
                       fc_cutoff = 2,
                       pval_cutoff = 0.05,
                       organism = NULL,
                       readable = TRUE) {
  # check input param
  stopifnot(inherits(dds, "DESeqDataSet"))
  # dir
  if (!dir.exists(outdir)) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE, mode = "0755")
  }

  # prepare files
  de_dds <- file.path(outdir, "DESeq2_dds.rds")
  de_fpkm <- file.path(outdir, "gene_fpkm.csv")
  de_count <- file.path(outdir, "transcripts_deseq2.csv")
  de_xls <- file.path(outdir, "transcripts_deseq2.fix.xls")
  de_plots <- file.path(outdir, c(
    "figure1_MA_plot.png",
    "figure2_MA_plot_LFC.png",
    "figure3_sample_counts.png",
    "figure4_PCA_plot.png",
    "figure5_dispersion.png",
    "figure6_sample_distance.png",
    "figure7_top_genes.png",
    "figure8_volcano.png"
  ))
  # save dds to file
  saveRDS(dds, file = de_dds)
  # dds <- readRDS(de_dds) # restore data

  # Run the DESeq pipeline
  dds <- DESeq2::DESeq(dds)

  # save table
  res <- DESeq2::results(dds)

  # model test
  # Get differential expression results
  resLFC <- DESeq2::lfcShrink(dds, coef = 2, res = res, type = "apeglm")
  rld <- DESeq2::rlogTransformation(dds, blind = FALSE)
  ntd <- DESeq2::normTransform(dds)

  ## save to file
  res <- res[order(res$padj), ] # Order by adjusted p-value
  resSig <- subset(as.data.frame(res), padj < pval_cutoff)
  ncount <- DESeq2::counts(dds, normalized = TRUE) # Normalied counts

  ## --------------------##
  ## Add annotation     ##
  ## --------------------##
  if (isTRUE(readable)) {
    message(paste0("DESeq2 deseq_node() ", organism))
    res$symbol   <- gene_to_symbol(rownames(res), organism, "SYMBOL")
    res$entrezid <- gene_to_symbol(rownames(res), organism, "ENTREZID")
  }

  ## Merge with normalized count data
  resdata <- merge(as.data.frame(ncount),
                   as.data.frame(res),
                   by = "row.names",
                   sort = FALSE)
  resdata <- resdata[order(resdata$padj), ]
  names(resdata)[1] <- "Gene"

  # save data to file
  resdata <- deseq_sig(resdata, fc_cutoff, pval_cutoff)
  write.csv(resdata, de_count, quote = TRUE, row.names = TRUE)

  # save to xls
  df_csv <- read_deseq_csv(de_count)
  readr::write_delim(df_csv, de_xls, delim = "\t", col_names = TRUE)

  # save FPKM to file
  df_fpkm <- DESeq2::fpkm(dds) %>%
    as.data.frame.matrix() %>%
    tibble::rownames_to_column("Gene")
  readr::write_csv(df_fpkm, de_fpkm, quote_escape = "double")

  # MA
  png(de_plots[1], width = 1200, height = 1200, res = 300)
  DESeq2::plotMA(res, ylim = c(-2, 2))
  dev.off()

  # MA for LFC
  png(de_plots[2], width = 1200, height = 1200, res = 300)
  DESeq2::plotMA(resLFC, ylim = c(-2, 2))
  dev.off()

  # Sample counts
  png(de_plots[3], width = 1200, height = 1200, res = 300)
  DESeq2::plotCounts(dds, gene = which.min(res$padj), intgroup = "condition")
  dev.off()

  # PCA
  png(de_plots[4], width = 2000, height = 2000, res = 300)
  print(DESeq2::plotPCA(rld, intgroup = c("condition")))
  dev.off()

  # Dispersion
  png(de_plots[5], width = 1500, height = 1500, res = 300)
  DESeq2::plotDispEsts(dds)
  dev.off()

  # Sample distance
  png(de_plots[6], width = 1000, height = 1000, res = 300)
  sampleDists <- dist(t(SummarizedExperiment::assay(rld)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- rld$condition
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255)
  pheatmap::pheatmap(sampleDistMatrix,
                     clustering_distance_rows = sampleDists,
                     clustering_distance_cols = sampleDists,
                     col = colors
  )
  dev.off()

  # Top genes
  png(de_plots[7], width = 1200, height = 1200, res = 300)
  select <- order(rowMeans(DESeq2::counts(dds, normalized = TRUE)),
                  decreasing = TRUE
  )[1:30]
  ma <- SummarizedExperiment::assay(ntd)[select, ]
  df <- as.data.frame(SummarizedExperiment::colData(dds)[, c("condition")])
  colnames(df) <- "condition"
  rownames(df) <- colnames(ma)
  pheatmap::pheatmap(SummarizedExperiment::assay(ntd)[select, ],
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     show_rownames = FALSE,
                     annotation_col = df
  )
  dev.off()
}


#' import matrix data to DESeq2
#' using DESeqDataSetFromMatrix()
#'
#' @param ma matrix sample vs genes, non-negative integers
#' @param smp_control vector, sample names in ma, as control
#' @param smp_treatment vector, sample names in ma, as treatment
#'
#' output: DESeqDataSet: dds
#'
#' @export
deseq_import_from_matrix <- function(ma, smp_control, smp_treatment) {
  # ma should be matrix
  stopifnot(inherits(ma, "matrix"))

  # convert matrix to non-negative integers
  ma <- round(ma, digits = 0) #

  # remove negative columns
  neg_cols <- apply(ma < 0, 1, any)
  if (sum(neg_cols) > 0) {
    warning(paste0("Remove ", sum(neg_cols), " records, with negative counts"))
    ma <- ma[!neg_cols, ]
  }

  # check samples names in ma
  stopifnot(all(c(smp_control, smp_treatment) %in% colnames(ma)))

  # create design
  colData <- deseq_design(ma, smp_control, smp_treatment)

  # import from matrix
  # DESeq2::DESeqDataSetFromMatrix
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = ma,
    colData = colData,
    design = ~condition
  )

  # output
  return(dds)
}


#' create experiment design for DESeq2 analysis
#'
#' @param ma count matrix for deseq2
#' @param smp_control vector, a list of character, names of control
#' @param smp_treatment vector, a list of character, names of treatment
#'
#' @export
deseq_design <- function(ma, smp_control, smp_treatment) {
  # build design data.frame
  # sample_name/condition
  # rownames (samples)
  # colnames (condition), ctl, exp
  stopifnot(all(c(smp_control, smp_treatment) %in% colnames(ma)))
  # conditions
  df <- data.frame(condition = factor(
    c(
      rep("control", times = length(smp_control)),
      rep("treatment", times = length(smp_treatment))
    ),
    levels = c("control", "treatment")
  ))
  # assign sample names
  rownames(df) <- c(smp_control, smp_treatment)
  return(df)
}


