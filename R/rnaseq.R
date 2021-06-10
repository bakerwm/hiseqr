#' Functions for RNAseq analysis
#'
#' Apply DESeq2 for RNAseq DE analysis
#' require >=2 replicates for each sample
#'
#' input: count.txt, from featureCounts output
#' output: DESeq2 output
#'
#' helper functions:
#'
#' 1. read_fc
#'
#' 2.
#'
#'
#' @import fs
#' @import ggplot2
#' @import ggthemes
#' @import dplyr
#' @import tidyr
#' @import readr
#' @import tibble
#' @import patchwork
#'
#'
#' @name rnaseq


# library(ggplot2)
# library(ggrepel)
# library(patchwork)
# library(dplyr)
# library(ggthemes)


# --Main: Port for pipeline ----------------------------------------------------

#' @describeIn rnaseq_hub A port for RNAseq analysis, DE analysis
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
rnaseq_hub <- function(x, ...) {
  if(! is_hiseq_dir(x, "rnaseq_rx")) {
    msg <- paste0("`x` expect a rnaseq_rx, ", x)
    stop(msg)
  }
  # Prepare data for DESeq2 analysis: count.matrix, ctl_name, tret_name
  pd <- prep_deseq(x)
  # Convert to dds for DESeq2 analysis
  dds <- import_matrix(pd$ma, pd$control_name, pd$treatment_name)
  # Add meta data (Length) to dds, for FPKM/RPKM
  r1_list <- list_hiseq_file(x, "count_sens", "_r1")
  df      <- read_fc(r1_list[1], get_gene_length = TRUE)
  m_data  <- data.frame(basepairs = df$Length)
  S4Vectors::mcols(dds) <- S4Vectors::DataFrame(S4Vectors::mcols(dds), m_data)
  # Run standard DEseq2
  # Add extra arguments: fc_cutoff, pval_cutoff, readable
  deseq2_main(dds,
              outdir = list_hiseq_file(x, "deseq_dir", "_rx"),
              organism = list_hiseq_file(x, "genome", "_rx"))
  # Publish quality figures: scatter, MA, volcano
  make_publish_plots(x)
  # fix_xls <- get_fix_xls(x)
  # if(file.exists(fix_xls)) {
  #   message("Generating publish quality plots")
  #   message(paste0("found DESeq2 ouptut: ", fix_xls))
  #   make_publish_plots(x)
  # }
  # #-- Done in python code, pipeline
  # #-- To-Do: report
  # report_dir <- file.path(x, "report")
  # rnaseq_report(x, report_dir)
}



#' @describeIn rnaseq_enrich A port for RNAseq enrich analysis
#'
#' @param x path to the RNAseqRx directory: a.vs.b
#' @param ... extra arguments for demseq2_main function
#'
#' extra arguments:
#' - overwrite    : FALSE
#' - text_width   : 40
#' - pval_cutoff  : 0.9,
#' - qval_cutoff  : 0.9,
#' - level        : 2
#' - readable     : bool
#' - layout       : nicely, kk
#' - show_category: 12
#'
#' @import dplyr
#' @import DESeq2
#'
#' @example rnaseq_enrich_hub(x)
#'
#' @export
rnaseq_enrich_hub <- function(x, ...) {
  if(! is_hiseq_dir(x, "rnaseq_rx")) {
    msg <- paste0("`x` expect a RNAseqRx, eg: a.vs.b directory, failed: ", x)
    stop(msg)
  }
  # contains: sig_files, orgdb, organism, keytype, enrich_dir
  dots  <- rlang::list2(...)
  edata <- prep_rnaseq_enrich(x)
  arg_vars  <- purrr::list_modify(dots, !!!edata)
  if(is(edata, "list")) {
    #-- Run: enrich analysis
    tmp <- lapply(edata$sig_files, function(s) {
      arg_vars$outdir <- dirname(s) #@@@ update outdir
      df <- read_text(s)
      fold_change <- setNames(df$log2FoldChange, nm = df$Gene)
      arg_vars    <- arg_vars %>%
        purrr::list_modify(
          gene_list   = df %>% dplyr::pull(Gene),
          fold_change = structure(df$log2FoldChange, names = df$Gene))
      tmp <- do.call(enrich_hub, arg_vars)
    })
    #-- Report: for RNAseq
    template    <- system.file("rnaseq",
                               "hiseq_report_rx_enrich.Rmd",
                               package = "hiseqr")
    report_dir  <- file.path(x, "enrich", "report")
    rnaseq_report(x, report_dir, template)
  }
}




# --Main: DESeq2 standard pipeline ---------------------------------------------

#' @describeIn rnaseq The main function do the analysis
#'
#' @param dds count data in matrix
#' @param outdir Directory to save the results
#' @param fc_cutoff float Cutoff for foldchange, default: 2
#' @param pval_cutoff float Cutoff for pvalue, default: 0.05
#' @param organism character, the name of the organism, eg: dm6, hg38, human,
#' @param readable bool, Add symbol and entrezid to table
#'
#' @import DESeq2
#' @import apeglm
#' @import ggplot2
#' @import pheatmap
#' @import RColorBrewer
#' @import SummarizedExperiment
#'
#' @export
deseq2_main <- function(dds,
                        outdir,
                        fc_cutoff = 2,
                        pval_cutoff = 0.05,
                        organism = NULL,
                        readable = TRUE) {
  # parallel
  BiocParallel::register(BiocParallel::MulticoreParam(4))
  # Check dds, DESeq2 object
  stopifnot(inherits(dds, "DESeqDataSet"))
  # Prepare output directory
  if (! dir.exists(outdir)) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE, mode = "0755")
  }
  # Default output files
  dds_rds   <- file.path(outdir, "DESeq2_dds.rds")
  de_fpkm   <- file.path(outdir, "gene_fpkm.csv")
  de_count  <- file.path(outdir, "transcripts_deseq2.csv")
  de_xls    <- file.path(outdir, "transcripts_deseq2.fix.xls")
  plot_name <- setNames(c(
    "figure1_MA_plot.png",
    "figure2_MA_plot_LFC.png",
    "figure3_sample_counts.png",
    "figure4_PCA_plot.png",
    "figure5_dispersion.png",
    "figure6_sample_distance.png",
    "figure7_top_genes.png",
    "figure8_volcano.png"
  ),
  nm = c(
    "ma",
    "ma_lfc",
    "counts",
    "pca",
    "dispersion",
    "sample_distance",
    "top_gene",
    "volcano"
  )
  )
  de_plots <- sapply(plot_name, function(i) {
    file.path(outdir, i)
  }, USE.NAMES = TRUE, simplify = FALSE)
  # Run the DESeq pipeline
  if(file.exists(dds_rds)) {
    dds <- readRDS(dds_rds)
  } else {
    dds <- DESeq2::DESeq(dds)
    # Save dds to file
    saveRDS(dds, file = dds_rds)
  }
  # Get table
  res <- DESeq2::results(dds)
  # Test models
  # Get differential expression results
  res_lfc <- DESeq2::lfcShrink(dds, coef = 2, res = res, type = "apeglm")
  ntd <- DESeq2::normTransform(dds)
  vsd <- DESeq2::vst(dds, blind=FALSE)
  rld <- DESeq2::rlog(dds, blind=FALSE)
  # rld <- DESeq2::rlogTransformation(dds, blind = FALSE)
  res <- res[order(res$padj), ] # Order by adjusted p-value
  res_sig <- subset(as.data.frame(res), padj < pval_cutoff)
  ncount <- DESeq2::counts(dds, normalized = TRUE) # Normalied counts
  # Add annotation, gene_symbol, entrezid
  if (isTRUE(readable) & is.character(organism)) {
    keytype <- guess_keytype(rownames(res), organism = organism) # !!!! error mm10
    df_gene <- convert_id(rownames(res),
                          from_keytype = keytype,
                          to_keytype   = c("ENTREZID", "SYMBOL"),
                          organism     = organism,
                          na_rm        = FALSE)
    ## convert ids
    ## to entrezid
    to_entrez  <- setNames(df_gene$ENTREZID, nm = df_gene[, keytype])
    res$entrez <- dplyr::recode(rownames(res), !!!to_entrez)
    ## to symbol
    to_symbol  <- setNames(df_gene$SYMBOL, nm = df_gene[, keytype])
    res$symbol <- dplyr::recode(rownames(res), !!!to_symbol)
  }
  # Merge with normalized count data
  resdata <- merge(as.data.frame(ncount),
                   as.data.frame(res),
                   by   = "row.names",
                   sort = FALSE
  )
  resdata <- resdata[order(resdata$padj), ]
  names(resdata)[1] <- "Gene"
  # Save results to file, csv
  # resdata2 <- read_deseq_csv(resdata, fc_cutoff, pval_cutoff)
  write.csv(resdata, de_count, quote = TRUE, row.names = TRUE)
  # Save results to file, xls
  df_csv <- read_deseq_csv(de_count, fc_cutoff, pval_cutoff, readable = TRUE,
                           organism = organism)
  readr::write_delim(df_csv, de_xls, delim = "\t", col_names = TRUE)
  # Save FPKM to file, xls
  df_fpkm <- DESeq2::fpkm(dds) %>%
    as.data.frame.matrix() %>%
    tibble::rownames_to_column("Gene")
  readr::write_csv(df_fpkm, de_fpkm, quote_escape = "double")
  # MA plot
  png(de_plots$ma, width = 1200, height = 1200, res = 300)
  DESeq2::plotMA(res, ylim = c(-2, 2))
  dev.off()
  # MA for LFC
  png(de_plots$ma_lfc, width = 1200, height = 1200, res = 300)
  DESeq2::plotMA(res_lfc, ylim = c(-2, 2))
  dev.off()
  # Sample counts
  png(de_plots$counts, width = 1200, height = 1200, res = 300)
  DESeq2::plotCounts(dds, gene = which.min(res$padj), intgroup = "condition")
  dev.off()
  # PCA
  png(de_plots$pca, width = 2000, height = 2000, res = 300)
  print(DESeq2::plotPCA(rld, intgroup = c("condition")))
  dev.off()
  # Dispersion
  png(de_plots$dispersion, width = 1500, height = 1500, res = 300)
  DESeq2::plotDispEsts(dds)
  dev.off()
  # Sample distance
  png(de_plots$sample_distance, width = 1000, height = 1000, res = 300)
  sample_dist <- dist(t(SummarizedExperiment::assay(rld)))
  sample_dist_matrix <- as.matrix(sample_dist)
  rownames(sample_dist_matrix) <- rld$condition
  colnames(sample_dist_matrix) <- NULL
  colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255)
  pheatmap::pheatmap(sample_dist_matrix,
                     clustering_distance_rows = sample_dist,
                     clustering_distance_cols = sample_dist,
                     col = colors
  )
  dev.off()
  # Top genes
  png(de_plots$top_gene, width = 1200, height = 1200, res = 300)
  s <- order(rowMeans(DESeq2::counts(dds, normalized = TRUE)),
             decreasing = TRUE
  )[1:30]
  ma <- SummarizedExperiment::assay(ntd)[s, ]
  df <- as.data.frame(SummarizedExperiment::colData(dds)[, c("condition")])
  colnames(df) <- "condition"
  rownames(df) <- colnames(ma)
  pheatmap::pheatmap(SummarizedExperiment::assay(ntd)[s, ],
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     show_rownames = FALSE,
                     annotation_col = df
  )
  dev.off()
  # output
  de_xls
}







# --Utils: Prepare data --------------------------------------------------------

#' @describeIn import_matrix Construct dds for DESeq2 analysis using matrix
#'
#' using DESeqDataSetFromMatrix()
#'
#' @param ma matrix sample vs genes, non-negative integers
#' @param control_name vector, sample names in ma, as control
#' @param treatment_name vector, sample names in ma, as treatment
#'
#' output: DESeqDataSet: dds
#'
#' @export
import_matrix <- function(ma, control_name, treatment_name) {
  # Check arguments
  stopifnot(inherits(ma, "matrix"))
  # Check sample names in ma
  stopifnot(all(c(control_name, treatment_name) %in% colnames(ma)))
  # Convert to integers, if float exists
  ma <- round(ma, digits = 0)
  # Save only positive values
  neg_cols <- apply(ma < 0, 1, any)
  if (sum(neg_cols) > 0) {
    warning(paste0("Remove ", sum(neg_cols), " records, with negative counts"))
    ma <- ma[! neg_cols, ]
  }
  # create design
  col_data <- data.frame(condition = factor(c(
    rep("control", times = length(control_name)),
    rep("treatment", times = length(treatment_name))
  ),
  levels = c("control", "treatment")
  ))
  rownames(col_data) <- c(control_name, treatment_name)
  # return dds
  DESeq2::DESeqDataSetFromMatrix(
    countData = ma,
    colData   = col_data,
    design    = ~condition
  )
}




#' @describeIn prep_deseq Prepare data for DESeq analysis
#'
#' @param x path to the directory of RNAseqRx
#'
#' @import readr
#' @import configr
#' @import dplyr
#'
#' @export
prep_deseq <- function(x) {
  if(! is_hiseq_dir(x, "rnaseq_rx")) {
    msg <- paste0("`x` expect a RNAseqRx, eg: a.vs.b directory, failed: ", x)
    stop(msg)
  }
  # Load control/treatment, wildtype/mutant files
  pd <- read_hiseq(x)
  # Wildtype dir, smp_name
  wt_dir    <- list_hiseq_file(x, "wt_dir", "_rx")
  wt_r1     <- list_hiseq_dir(wt_dir, "_r1")
  wt_names  <- list_hiseq_file(wt_dir, "smp_name", "_r1")
  wt_count_files <- list_hiseq_file(wt_dir, "count_sens", "r1")
  # Mutant dir, smp_name
  mut_dir   <- list_hiseq_file(x, "mut_dir", "_rx")
  mut_r1    <- list_hiseq_dir(mut_dir, "_r1")
  mut_names <- list_hiseq_file(mut_dir, "smp_name", "r1")
  mut_count_files <- list_hiseq_file(mut_dir, "count_sens", "r1")
  # Parsing the genome
  genome <- list_hiseq_file(x, "genome", "_rx")
  # Parsing count.txt files
  c_files <- c(wt_count_files, mut_count_files)
  df <- hiseqr::read_fc2(c_files)
  colnames(df) <- c("id", wt_names, mut_names) # shorter names?!!!!
  # Convert data.frame to matrix
  ma <- df %>%
    tibble::column_to_rownames("id") %>%
    as.matrix()
  # Output
  list(
    ma              = ma,
    control_name    = wt_names,
    treatment_name  = mut_names,
    genome          = genome,
    control_count   = wt_count_files,
    treatment_count = mut_count_files
  )
}




#' @describeIn prep_rnaseq_enrich Prepare data for Enrich analysis
#'
#' @description Extract significant genes, create dirs
#'
#' @param x path to the directory of RNAseqRx
#'
#' @import readr
#' @import configr
#' @import dplyr
#'
#' @export
prep_rnaseq_enrich <- function(x) {
  if(! is_hiseq_dir(x, "rnaseq_rx")) {
    stop(glue::glue("rnaseq_rx dir expected, faild {x}"))
  }
  #--Parsing config.toml
  px       <- read_hiseq(x[1])
  organism <- list_hiseq_file(x, "genome", "rnaseq_rx")
  if(is.null(organism)) {
    msg <- paste0("`organism` not supported, [", as.character(organism), "], ",
                  "see `Organism.dplyr::supportedOrganisms()` ",
                  "for full list of organisms")
    warning(msg)
    return(NULL)
  }
  orgdb <- get_orgdb(organism)
  # enrich_dir <- px$args$enrich_dir
  enrich_dir <- file.path(x, "enrich")
  report_dir <- file.path(enrich_dir, "report")
  fix_xls    <- get_fix_xls(x)
  if(is(fix_xls, "character") & file.exists(fix_xls)) {
    df <- read_text(fix_xls)
  } else {
    msg <- paste0("`transcripts_deseq2.fix.xls` not found in `x`, ", x)
    warning(msg)
    return(NULL)
  }
  #--Output: required columns
  req_cols <- c("Gene", "sig")
  if(! all(req_cols %in% names(df))) {
    warning("`x` required columns missing: ",
            paste(req_cols, collapse = ", "))
    return(NULL)
  }
  gene_list <- df$Gene
  # guess keytype
  keytype   <- guess_keytype(gene_list, organism)
  # sig genes: up, down, not, sig
  # files: gene_list.xls, require `Gene`
  sig_files  <- sapply(c("sig", "up", "down"), function(s) {
    subdir   <- file.path(enrich_dir, s)
    sig_file <- file.path(subdir, "gene_list.xls")
    if(! dir.exists(subdir)) {
      dir.create(subdir, recursive = TRUE, mode = "0755")
    }
    if(s == "sig") {
      s = c("up", "down")
    }
    sig_df <- df %>%
      dplyr::filter(sig %in% s)
    write.table(sig_df, file = sig_file, sep = "\t", row.names = FALSE,
                col.names = TRUE)
    # output
    normalizePath(sig_file)
    # sig_df$Gene
  })
  list(
    sig_files  = sig_files,
    # sig_genes  = sig_genes,
    # sig_types  = c("sig", "up", "down"),
    organism   = organism,
    orgdb      = orgdb,
    keytype    = keytype,
    enrich_dir = enrich_dir
  )
}



# --Utils: Functions -----------------------------------------------------------

#' @describeIn read_deseq_csv Default csv file
#'
#' cal mean of replicates
#' default, 2 replicates for each condition
#'
#' @param x string Path to the .csv file
#' @param fc_cutoff numeric foldChange cutoff, default: 2
#' @param pval_cutoff numeric pvalue cutoff, default: 0.05
#'
#' @export
read_deseq_csv <- function(x, fc_cutoff = 2, pval_cutoff = 0.05,
                           readable = TRUE, organism = NULL) {
  if(! is(x, "character")) {
    warning("'x' required data.frame")
    return(NULL)
  }
  data <- readr::read_csv(x) %>%
    deseq_csv_mean %>%
    get_sig_gene(fc_cutoff = fc_cutoff, pval_cutoff = pval_cutoff)
  # Add symbol, entrezid
  if(is.character(organism) && isTRUE(readable)) {
    keytype <- guess_keytype(data$Gene, organism = organism)
    df_gene <- convert_id(data$Gene,
                          from_keytype = keytype,
                          to_keytype   = c("ENTREZID", "SYMBOL"),
                          organism     = organism)
    merge(data, df_gene, by.x = "Gene", by.y = keytype, all.x = TRUE)
  } else {
    data
  }
}






#' @describeIn get_sig_gene Extract the significantly changed genes
#'
#' add sig label
#'
#' @param data data.frame, csv file
#' @param fc_cutoff numeric, cutoff for foldchange, default: 2
#' @param pval_cutoff numeric, cutoff for pvalue, default: 0.05
#' @param type string, all, both, up, down, not, default: all
#'
#' @import dplyr
#'
#' @export
get_sig_gene <- function(data, fc_cutoff = 2, pval_cutoff = 0.05,
                         type = "all") {
  stopifnot(is(data, "data.frame"))
  if(! is(fc_cutoff, "integer")) {
    fc_cutoff <- 2
  }
  if(! is(pval_cutoff, "float")) {
    pval_cutoff <- 0.05
  }
  if(! "log2FoldChange" %in% names(data)) {
    stop("'log2FoldChange' not found in data")
  }
  if(type %in% c("up", "not", "down", "sig", "all")) {
    sig_groups <- switch (
      type,
      "sig"  = c("up", "down"),
      "all"  = c("up", "not", "down"),
      "up"   = "up",
      "down" = "down",
      "not"  = "not"
    )
  } else {
    sig_groups <- c("up", "down", "not")
  }
  #--sig: exists
  if("sig" %in% names(data)) {
    if(all(data$sig %in% c("up", "not", "down"))) {
      out <- dplyr::filter(data, sig %in% sig)
      return(out)
    } else {
      data <- data %>%
        dplyr::rename(sig_old = sig)
    }
  }
  #--filter: fc, pval
  log2_fc <- log2(as.numeric(fc_cutoff))
  if("padj" %in% names(data)) {
    data$pc <- data$padj
  } else if("pvalue" %in% names(data)) {
    data$pc <- data$pvalue
  } else {
    stop("padj or pvalue, not found in data")
  }
  # replace pval: NA->1; log2FoldChange: NA->0
  data <- data %>%
    tidyr::replace_na(list(pc = 1, log2FoldChange = 0))
  # add sig
  data <- data %>%
    dplyr::mutate(
      sig = ifelse(pc > pval_cutoff, "not",
                   ifelse(log2FoldChange >= log2_fc, "up",
                          ifelse(log2FoldChange <= -log2_fc, "down", "not"))))
  #--output: sig_groups
  d2 = data %>%
    dplyr::select(-pc) %>%
    dplyr::filter(sig %in% sig_groups)
}





#' @describeIn deseq_mean mean value of replicates
#'
#' @param data data.frame From res(dds) output
#'
#' @export
deseq_csv_mean <- function(data) {
  stopifnot(is(data, "data.frame"))
  col_required <- c("Gene", "baseMean", "padj")
  if(! all(col_required %in% names(data))) {
    stop("Missing required columns: ", paste(col_required, collapse = ", "))
  }
  # sample names
  rep_name <- data %>%
    dplyr::select(Gene:baseMean) %>%
    dplyr::select(-Gene, -baseMean) %>%
    names()
  # groups
  group_name <- fq_name(rep_name, fix_rep = TRUE) %>% unique
  if(! length(group_name) == 2) {
    stop("Only two groups supported")
  }
  g1 <- group_name[1]
  g2 <- group_name[2]
  data %>%
    dplyr::mutate(!! g1 := dplyr::select(., starts_with(g1)) %>% rowMeans(),
                  !! g2 := dplyr::select(., starts_with(g2)) %>% rowMeans()) %>%
    dplyr::select(Gene, all_of(group_name), all_of(rep_name), baseMean:padj)
}









#' @describeIn make_publish_plots Create plots with pubilsh quality
#'
#' @param x deseq directory, file of DEseq2 output
#' @param outdir character Path to the directory, saving the plots
#' @param save2pdf bool Save the plots in PDF file,
#'
#' @import ggplot2
#' @import ggrepel
#'
#' @export
make_publish_plots <- function(x, outdir = NULL, to_pdf = TRUE) {
  fix_xls <- get_fix_xls(x)
  if(is.null(fix_xls)) {
    on.exit(paste0("not a a.vs.b directory: ", x))
  }
  # prepare data
  df <- readr::read_delim(fix_xls, "\t", col_types = readr::cols())
  if("symbol" %in% names(df)) {
    df$label <- df$symbol
  } else if("SYMBOL" %in% names(df)) {
    df$label <- df$SYMBOL
  } else {
    df$label <- df$Gene
  }
  #--args: output dir
  deseq_dir <- list_hiseq_file(x, "deseq_dir", "rnaseq_rx")
  if(is.null(outdir)) {
    outdir <- deseq_dir
  }
  if(! dir.exists(outdir)) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE, mode = "0740")
  }
  #--args: deseq_dir
  xname <- list_hiseq_file(x, "wt_name", "rnaseq_rx")
  yname <- list_hiseq_file(x, "mut_name", "rnaseq_rx")
  #--args: in case xname,yname is NULL (not from deseq_dir)
  if(is.null(xname) | is.null(yname)) {
    xname <- names(df)[2]
    yname <- names(df)[3]
  }
  #--check: required columns
  req_cols <- c("Gene", "baseMean", "log2FoldChange", "padj", "sig", xname,
                yname)
  if(! all(req_cols %in% names(df))) {
    on.exit(message("required columns not found"))
  }
  #--labels: sig genes, label
  labels <- df %>%
    dplyr::filter(sig %in% c("up", "down")) %>%
    dplyr::arrange(padj) %>%
    dplyr::pull(label) %>%
    head(5)
  labels <- purrr::discard(labels, is.na)
  # png files
  p1_png <- file.path(outdir, "publish.1.scatter.png")
  p2_png <- file.path(outdir, "publish.2.ma.png")
  p3_png <- file.path(outdir, "publish.3.volcano.png")
  p4_png <- file.path(outdir, "publish.4.heatmap.png")
  # plot-1: scatter
  p1 <- scatter_plot2(
    df, xname, yname,
    highlight_column = "sig",
    highlight_values = c("up", "down"),
    label_column     = "SYMBOL",
    label_values     = labels,
    point_color      = "grey60")
  png(p1_png, width = 5, height = 4.5, units = "in", res = 150)
  print(p1)
  dev.off()
  # plot-2: ma
  p2 <- ma_plot(df, labels = labels)
  png(p2_png, width = 5, height = 4.5, units = "in", res = 150)
  print(p2)
  dev.off()
  # plot-3: volcano
  df$log10pval <- sapply(df$padj, function(i) {
    ifelse(is.null(i), 0, -log10(i))
  })
  p3 <- volcano_plot(df, x = "log2FoldChange", y = "log10pval", labels = labels)
  png(p3_png, width = 5, height = 4.5, units = "in", res = 150)
  print(p3)
  dev.off()
  # plot-4: heatmap
  # fpkm
  # library(ComplexHeatmap)
  df_sig <- df %>%
    dplyr::filter(sig %in% c("up", "down"))
  f_fpkm <- file.path(outdir, "gene_fpkm.csv")
  if(file.exists(f_fpkm)) {
    if(nrow(df_sig) > 0) {
      df3 <- readr::read_csv(f_fpkm, col_types = readr::cols())
      # only sig genes
      df4 <- df3 %>%
        dplyr::filter(Gene %in% df_sig$Gene) %>%
        tibble::column_to_rownames("Gene") %>%
        as.matrix()
      # check sig
      sig_list <- dplyr::recode(rownames(df4), !!! setNames(df_sig$sig, df_sig$Gene))
      df5 <- log10(df4+1)
      png(p4_png, width = 4, height = 10, units = "in", res = 150)
      p4 <- ComplexHeatmap::Heatmap(
        df5,
        show_row_names = nrow(df5) < 30,
        cluster_columns = FALSE,
        row_split = sig_list,
        heatmap_legend_param = list(
          title = "log10(FPKM+1)",
          legend_height = unit(4, "cm"),
          title_position = "lefttop-rot"))
      print(p4)
      dev.off()
    }
  }
  # Add legend
  legend <- glue::glue(
    "Figure. Differentially expression analysis. ",
    "(A-B) Comparasion of ene expression is shown as rpm from two conditions.",
    "Dashed lines indicate two fold change. A. Scatter plot, B. MA plot.",
    "(C) Volcano plot showing enrichment values and corresponding significance levels.", .sep = "\n")
  p <- patchwork::wrap_plots(list(p1, p2, p3), ncol = 1) +
    patchwork::plot_annotation(tag_levels = "A",
                               title = "RNAseq plot",
                               caption = legend)
  # save plots to RDS
  px <- list(
    p1     = p1,
    p2     = p2,
    p3     = p3,
    p4     = p4,
    data   = df,
    labels = labels,
    xname  = xname,
    yname  = yname)
  px_file <- file.path(outdir, "publish_plot_data.rds")
  saveRDS(px, file = px_file)
  # save plots to PDF
  if(isTRUE(to_pdf)) {
    pdf_fname <- glue::glue("DESeq2.{xname}.vs.{yname}.publish.pdf")
    pdf_file  <- file.path(outdir, pdf_fname)
    pdf(pdf_file, width = 4.5, height = 12, paper = "a4")
    print(p)
    dev.off()
  }
}












#' #' @describeIn get_fix_xls Parsing the fix.xls file from DESeq2
#' #'
#' #'
#' #' @export
#' get_fix_xls <- function(x) {
#'   if(is.character(x)) {
#'     x <- x[1]
#'   } else {
#'     on.exit(message("'x' require character"), add = TRUE)
#'     return(NULL)
#'   }
#'   # RNAseqRx directory
#'   # x/deseq/transcripts_deseq2.fix.xls
#'   hiseq_type <- get_hiseq_type(x)
#'   if(is.null(hiseq_type)) {
#'     msg = paste0("`x` is not hiseq directory, ", x)
#'     on.exit(message(msg), add = TRUE)
#'     return(NULL)
#'   } else if (! hiseq_type == "rnaseq_rx") {
#'     msg = paste0("`x` is not hiseq_rx directory: ", x)
#'     on.exit(message(msg), add = TRUE)
#'     return(NULL)
#'   }
#'   pd        <- read_hiseq(x)
#'   organism  <- pd$args$genome
#'   deseq_dir <- pd$args$deseq_dir
#'   fix_xls   <- file.path(deseq_dir, "transcripts_deseq2.fix.xls")
#'   de_csv    <- file.path(deseq_dir, "transcripts_deseq2.csv")
#'   # convert csv to xls
#'   if(! file.exists(fix_xls)) {
#'     message("Fix deseq csv: add mean, symbol, sig")
#'     fix_df <- read_deseq_csv(de_csv, fc_cutoff = 2, pval_cutoff = 0.05,
#'                              readable = TRUE, organism = organism)
#'     readr::write_delim(fix_df, fix_xls, delim = "\t", col_names = TRUE)
#'   }
#'   # output
#'   if(file.exists(fix_xls)) {
#'     fix_xls
#'   } else {
#'     stop(paste0("fix_xls file not found: ", x))
#'   }
#' }



#' @describeIn get_fix_xls Parsing the fix.xls file from DESeq2
#'
#' @param x character, path to the a.vs.b, rnaseq_rx directory
#'
#' @export
get_fix_xls <- function(x) {
  if(! is_hiseq_dir(x, "rnaseq_rx")) {
    stop(glue::glue("rnaseq_rx dir expected, failed: {x}"))
  }
  # x/deseq/transcripts_deseq2.fix.xls
  if(is_hiseq_dir(x[1], "rnaseq_rx")) {
    organism  <- list_hiseq_file(x, "genome", "rnaseq_rx")
    deseq_dir <- list_hiseq_file(x, "deseq_dir", "rnaseq_rx")
    fix_xls   <- file.path(deseq_dir, "transcripts_deseq2.fix.xls")
    de_csv    <- file.path(deseq_dir, "transcripts_deseq2.csv")
    # convert csv to fix_xls
    if(! file.exists(fix_xls)) {
      if(file.exists(de_csv)) {
        message("Fix deseq csv: add mean, symbol, sig")
        fix_df <- read_deseq_csv(de_csv, fc_cutoff = 2, pval_cutoff = 0.05,
                                 readable = TRUE, organism = organism)
        readr::write_delim(fix_df, fix_xls, delim = "\t", col_names = TRUE)
      }
    }
    if(file.exists(fix_xls)) {
      fix_xls
    } else {
      warning("`x` transcript_deseq2.fix.xls file not exists")
    }
  }
}








#-- Report: Functions ----------------------------------------------------------

#' @describeIn rnaseq_report
#'
#' @param input directory to the sample
#' @param output directory to the html file
#' @param template the template, default from hiseqr
#'
#' @export
rnaseq_report <- function(input, output, template_rmd = NULL) {
  # input
  input <- normalizePath(input)
  if(! dir.exists(output)) dir.create(output, recursive = TRUE)
  output <- normalizePath(output)
  # output
  outhtml <- file.path(output, "HiSeq_report.html")
  # subtype
  if(is_hiseq_dir(input, "rnaseq_r1")) {
    hiseq_subtype = "hiseq_report_r1.Rmd"
  } else if(is_hiseq_dir(input, "rnaseq_rn")) {
    hiseq_subtype = "hiseq_report_rn.Rmd"
  } else if(is_hiseq_dir(input, "rnaseq_rx")) {
    hiseq_subtype = "hiseq_report_rx.Rmd"
  } else {
    hiseq_subtype = "tmp"
  }
  # template
  hiseq_type <- list_hiseq_file(input, "hiseq_type")
  hiseq_type <- gsub("_\\w+$", "", hiseq_type[1])
  if(is.null(template_rmd)) {
    template <- system.file(hiseq_type, hiseq_subtype, package = "hiseqr")
  } else {
    template <- template_rmd
  }
  stopifnot(file.exists(template))
  ## copy template to output
  template_to <- file.path(output, basename(template))
  file.copy(template, template_to, overwrite = TRUE)
  rmarkdown::render(input       = template_to,
                    output_file = outhtml,
                    params      = list(input_dir = input))
}




#' #' @describeIn list_rnaseq_single_dirs
#' #'
#' #' @param x path to the input directories
#' #' @param log boolen, whether print the log status, default: FALSE
#' #'
#' #' @export
#' list_rnaseq_single_dirs <- function(x){
#'   if(is_hiseq_dir(x)) {
#'     if(is_hiseq_single_dir(x)) {
#'       x
#'     } else if(is_hiseq_merge_dir(x)) {
#'       px   <- read_hiseq(x)
#'       dirs <- px$args$rep_list
#'       dirs <- purrr::discard(dirs, is.null)
#'       if(length(dirs) > 0) {
#'         sapply(dirs, list_rnaseq_single_dirs, simplify = TRUE)
#'       }
#'     } else if(is_hiseq_multiple_dir(x)) {
#'       px   <- read_hiseq(x)
#'       dirs <- c(px$args$wildtype_dir, px$args$mutant_dir)
#'       dirs <- purrr::discard(dirs, is.null)
#'       if(length(dirs) > 0) {
#'         sapply(dirs, list_rnaseq_single_dirs, simplify = TRUE)
#'       }
#'     }
#'   }
#' }
#'
#'
#'
#' #' @describeIn list_rnaseq_merge_dirs
#' #' rnaseq_merge_dirs
#' #'
#' #' @param x path to the input directories
#' #'
#' #' @export
#' list_rnaseq_merge_dirs <- function(x){
#'   if(is_hiseq_dir(x)) {
#'     if(is_hiseq_single_dir(x)) {
#'       NULL
#'     } else if(is_hiseq_merge_dir(x)) {
#'       x
#'     } else if(is_hiseq_multiple_dir(x)) {
#'       px   <- read_hiseq(x)
#'       dirs <- c(px$args$wildtype_dir, px$args$mutant_dir)
#'       dirs <- purrr::discard(dirs, is.null)
#'       if(length(dirs) > 0) {
#'         sapply(dirs, list_rnaseq_merge_dirs, simplify = TRUE)
#'       }
#'     }
#'   }
#' }
#'
#'
#'
#'
#' #' @describeIn list_rnaseq_multiple_dirs
#' #' rnaseq_merge_dirs
#' #'
#' #' @param x path to the input directories
#' #'
#' #' @export
#' list_rnaseq_multiple_dirs <- function(x){
#'   if(is_hiseq_dir(x)) {
#'     if(is_hiseq_multiple_dir(x)) {
#'       x
#'     }
#'   }
#' }






#' #' @describeIn get_rnaseq_trim_stat
#' #'
#' #' parsing the trimming status
#' #'
#' #' @param x path to hiseq, single
#' #' @import dplyr
#' #' @import readr
#' #'
#' #' @export
#' get_rnaseq_trim_stat <- function(x) {
#'   f <- list_hiseq_file(x, "trim_summary_json", "r1")
#'   df <- lapply(f, function(i) {
#'     jsonlite::read_json(i) %>%
#'       as.data.frame()
#'   }) %>%
#'     dplyr::bind_rows()
#'   if(nrow(df) > 0) {
#'     df %>%
#'       dplyr::select(name, input, output, out_pct, rm_pct) %>%
#'       dplyr::mutate(across(-name, as.numeric)) %>%
#'       dplyr::mutate(name = forcats::fct_rev(name)) %>%
#'       dplyr::rename(raw       = input,
#'                     clean     = output,
#'                     clean_pct = out_pct,
#'                     short_pct = rm_pct)
#'   }
#' }
#'
#'
#'
#'
#' #' @describeIn get_rnaseq_align_stat, read TOML files
#' #'
#' #' read align from rnaseq seq
#' #' @export
#' get_rnaseq_align_stat <- function(x) {
#'   f <- list_hiseq_file(x, "align_summary_json", "r1")
#'   df <- lapply(f, function(i) {
#'     jsonlite::read_json(i) %>%
#'       as.data.frame()
#'   }) %>%
#'     dplyr::bind_rows()
#'   if(nrow(df) > 0) {
#'     df %>%
#'       dplyr::mutate(
#'         map_pct = round(map / total * 100, 2),
#'         unique_pct = round(unique / total * 100, 2),
#'         unmap_pct = round(unmap / total * 100, 2)
#'       ) %>%
#'       dplyr::select(
#'         name, index, total, map, unique, multi, spikein, rRNA, unmap,
#'         map_pct, unique_pct, unmap_pct
#'       )
#'   }
#' }



#' #' @describeIn get_rnaseq_strandness
#' #'
#' #' extract strandness info
#' #' 1. strandness.json
#' #'
#' #' sens: 1, anti: 2, ++ -- / +- -+
#' #' sens: 2, anti: 1, +- -+ / ++ --
#' #'
#' #' @export
#' get_rnaseq_strandness <- function(x) {
#'   # sense strand
#'   f1 <- list_hiseq_file(x, "count_sens", "r1")
#'   if(is.character(f1)) {
#'     f1_json <- paste0(f1, ".summary.json")
#'     df1 <- lapply(f1_json, function(i) {
#'       m = jsonlite::read_json(i)
#'       as.data.frame(m[[1]])
#'     }) %>%
#'       dplyr::bind_rows() %>%
#'       dplyr::mutate(strand = "sense")
#'   } else {
#'     df1 <- NULL
#'   }
#'   # antisense strand
#'   f2 <- list_hiseq_file(x, "count_anti", "r1")
#'   if(is.character(f2)) {
#'     f2_json <- paste0(f2, ".summary.json")
#'     df2 <- lapply(f2_json, function(i) {
#'       m = jsonlite::read_json(i)
#'       as.data.frame(m[[1]])
#'     }) %>%
#'       dplyr::bind_rows() %>%
#'       dplyr::mutate(strand = "antisense")
#'   } else {
#'     df2 <- NULL
#'   }
#'   # output
#'   dplyr::bind_rows(list(df1, df2))
#' }





#' #' @describeIn get_rnaseq_count_txt
#' #'
#' #' @export
#' get_rnaseq_count_txt <- function(x, strand = "sens") {
#'   dirs <- list_rnaseq_single_dirs(x)
#'   dirs <- purrr::discard(dirs, is.null)
#'   if(length(dirs) > 0) {
#'     sapply(dirs, function(i) {
#'       px <- read_hiseq(i)
#'       if(strand == "anti") {
#'         px$args$count_anti
#'       } else {
#'         px$args$count_sens
#'       }
#'     })
#'   }
#' }



#' #' @export
#' get_rnaseq_report <- function(x) {
#'   # search for single dir
#'   rep_list      <- list_rnaseq_single_dirs(x)
#'   merge_list    <- list_rnaseq_merge_dirs(x)
#'   multiple_list <- list_rnaseq_multiple_dirs(x)
#'   dirs          <- sort(c(rep_list, merge_list, multiple_list))
#'   dirs          <- purrr::discard(dirs, is.null)
#'   report_list <- sapply(dirs, function(i){
#'     px <- read_hiseq(i)
#'     report_dir <- px$args$report_dir
#'     f_html  <- list.files(report_dir, "*.html", full.names = TRUE)
#'     if(length(f_html) > 0) {
#'       f_html[1]
#'     }
#'   }, simplify = TRUE, USE.NAMES = FALSE) %>%
#'     unlist()
#'
#'   purrr::discard(report_list, is.null)
#' }





