## functions for RNAseq


##----------------------------------------------------------------------------##
## For hiseq package, parsing directory
##

#' list_arguments_file
#'
#' @param x path to the directory.
#' expect x/config/{feature}/config/arguments.pickle,
#' x/{feature}/config/arguments.pickle
#'
#' @export
list_arguments_file <- function(x){
  # arguments file
  # RNAseq multiple: x/config/{feature}/config/arguments.txt
  # RNAseq single: x/{feature}/config/arguments.txt
  # RNAseq deseq: x/{feature}/config/arguments.txt
  chk0 <- dir.exists(x)

  # version-1: x/{feature}/...
  dir1 <- list.dirs(x, full.names = TRUE, recursive = FALSE)
  ps1 <- mapply(function(i) file.path(i, "config", "arguments.pickle"), dir1)
  if (length(ps1) > 0) {
    ps1 <- ps1[file.exists(ps1)]
  }
  chk1 <- length(ps1) > 0

  # version-2: x/config/{feature}/...
  dir2 <- list.dirs(file.path(x, "config"), full.names = TRUE, recursive = FALSE)
  ps2 <- mapply(function(i) file.path(i, "config", "arguments.pickle"), dir2)
  if (length(ps2) > 0) {
    ps2 <- ps2[file.exists(ps2)]
  }
  chk2 <- length(ps2) > 0

  ## mission
  ## RNAseq multiple: ps1=0, ps2=1
  ## RNAseq single:   ps1=1, ps2=0
  if(isTRUE(chk1)) {
    plist <- as.list(ps1)
  } else if(isTRUE(chk2)) {
    plist <- as.list(ps2)
  } else {
    plist <- NULL
  }

  ## change name
  if(! is.null(plist)) {
    names(plist) <- basename(dirname(dirname(unlist(plist))))
  }

  return(plist)
}


#' load pickle file
#'
#' @param x path to pickle,
#'
#' @import reticulate
#' @export
#'
load_pickle <- function(x){
  # python version issue !!
  if(file.exists(x) & endsWith(x, ".pickle")) {
    pd <- reticulate::import("pandas")
    tag <- pd$read_pickle(x)
  } else {
    tag <- NULL
  }
}


#' hiseq_type
#'
#' @param x path, the directory of RNAseq output
#' @param log boolen, whether print the log status, default: FALSE
#'
#' @export
#'
is_hiseq_dir <- function(x) {
  chk0 <- dir.exists(x)

  # check arguments.pickle file
  pk_files <- list_arguments_file(x)
  # choose one

  if(is.null(pk_files) | is.na(pk_files)) {
    return(NULL)
  } else {
    # pk_files <- pk_files[1] # 1st item
    pk <- load_pickle(pk_files[[1]])

    # get values
    if("rnaseq_type" %in% names(pk)) {
      tag <- pk$rnaseq_type
    } else if("atacseq_type" %in% names(pk)) {
      tag <- pk$atacseq_type
    } else {
      tag <- NULL # not known
    }
    return(tag) # return
  }
}


#' rnaseq_single_dir
#'
#' reading data from RNAseq Single
#'
#' @param x path to the directory of RNAseq single
#'
#' @export
#'
rnaseq_single_dir <- function(x, feature="gene"){
  # check
  chk0 <- hiseqr::is_hiseq_dir(x) == "rnaseq_single"

  if(! isTRUE(chk0)){
    return(NULL)
  }

  # pick pickle file
  pk_files <- list_arguments_file(x)
  args     <- hiseqr::load_pickle(pk_files[[feature]])

  ## check required args
  required_names <- c("align_stat", "smp_name", "genome", "gtf", "gtf",
                      "count_sens", "count_anti", "strandness_status")

  for(name in required_names){
    # args[[name]] = ifelse(name %in% names(args), c(args[[name]]), NULL)
    if(! name %in% names(args)){
      args[[name]] <- NULL
    }
  }

  ## mapping data
  df <- hiseqr::readAlign1(args$align_stat)
  args$reads_total <- sum(df$count)

  ## mito.u, mito.m, genome.u genome.m
  for(i in seq_len(nrow(df))) {
    name <- as.character(df$group[i])
    args[[name]] <- df$count[[i]]
  }

  ## save df
  args$align_stat_df <- df

  return(args)
}


#' rnaseq_multiple_dir
#'
#' return args for RNAseq_multiple/single/
#'
#' @param x path to the directory of RNAseq multiple
#'
#' @export
#'
#' @export
rnaseq_multiple_dir <- function(x, feature = "gene"){
  # check
  chk0 <- hiseqr::is_hiseq_dir(x) == "rnaseq_multiple"

  if(! isTRUE(chk0)){
    return(NULL)
  }

  # pick pickle file
  pk_files <- list_arguments_file(x)
  args     <- hiseqr::load_pickle(pk_files[[feature]])

  ## check required args
  required_names <- c("smp_name", "outdir", "genome", "gtf")

  for(name in required_names){
    # args[[name]] <- ifelse(name %in% names(args), c(args[[name]]), NULL)
    if(! name %in% required_names) {
      args[[name]] <- NULL
    }
  }

  ## parsing RNAseq single, saved as vector
  args$single_dirs <- mapply(function(i) file.path(args$outdir, i), args$smp_name)
  args$single_args <- lapply(args$single_dirs, rnaseq_single_dir)

  return(args)
}


##----------------------------------------------------------------------------##
## Run DESeq analysis using DESeq2
##


#' deseqHub
#' design for hiseq package, directory structure
#'
#' @param x path to the directory of project directory
#' @param feature character, gene|te|..., default, gene
#'
#' deseq_dir: prj_dir/feature/deseq
#' count_dir: prj_dir/feature/count
#'
deseqHub <- function(x, feature = "gene"){
  x = "results/hiseq_v2/RNAseq_read1/RNAseq_DaGal4XshPiwi-2_6h.vs.RNAseq_DaGal4XshWhite_6h"
  hiseq_type <- is_hiseq_dir(x)
  if(! hiseq_type == "deseq_single"){
    stop(glue::glue("DESeq dir expected, {hseq_type} found, {x}"))
  }

  # call for count files
  pk_files <- list_arguments_file(x)
  args     <- hiseqr::load_pickle(pk_files[[feature]])

  # reading data
  pdata <- countToMatrix(count_ctl = args$count_ctl,
                         count_exp = args$count_exp)

  # for DESeqDataSet
  dds <- deseqImportFromMatrix(ma            = pdata$ma,
                               smp_control   = pdata$smp_ctl,
                               smp_treatment = pdata$smp_exp)
  # for DEseq2
  deseqNode(dds,
            outdir = args$deseqdir,
            pvalue_cutoff = 0.1)
}


#' deseqHub2
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
deseqHub2 <- function(count_ctl, count_exp, outdir){
  # check
  stopifnot(all(file.exists(count_ctl)))
  stopifnot(all(file.exists(count_exp)))

  # import data
  pdata <- countToMatrix(count_ctl, count_exp)

  # for DESeqDataSet
  dds <- deseqImportFromMatrix(ma            = pdata$ma,
                               smp_control   = pdata$smp_ctl,
                               smp_treatment = pdata$smp_exp)
  # for DEseq2
  deseqNode(dds,
            outdir = outdir,
            pvalue_cutoff = 0.1)
}


#' run DESeq2 analysis
#' import: dds
#' output: basic analysis
#' gene-level
#'
#' @param dds count data in matrix
#' @param outdir Directory to save the results
#' @param pvalue_cutoff Cutoff to filt the records, padj for DESeq2 output,
#'   default: 0.1
#'
#' @import DESeq2
#' @import ggplot2
#' @import pheatmap
#' @import RColorBrewer
#' @import SummarizedExperiment
#'
#' @export
#'
deseqNode <- function(dds, outdir, pvalue_cutoff = 0.05) {
  # check input param
  stopifnot(inherits(dds, "DESeqDataSet"))
  # dir
  if(! dir.exists(outdir)) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE, mode = "0755")
  }

  # prepare files
  de_dds   <- file.path(outdir, "DESeq2_dds.rds")
  de_count <- file.path(outdir, "transcripts_deseq2.csv")
  de_xls   <- file.path(outdir, "transcripts_deseq2.fix.xls")
  de_plots <- file.path(outdir, c("figure1_MA_plot.png",
                                  "figure2_MA_plot_LFC.png",
                                  "figure3_sample_counts.png",
                                  "figure4_PCA_plot.png",
                                  "figure5_dispersion.png",
                                  "figure6_sample_distance.png",
                                  "figure7_top_genes.png",
                                  "figure8_volcano.png"))
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
  rld    <- DESeq2::rlogTransformation(dds, blind = FALSE)
  ntd    <- DESeq2::normTransform(dds)

  ## save to file
  res    <- res[order(res$padj), ]   # Order by adjusted p-value
  resSig <- subset(as.data.frame(res), padj < pvalue_cutoff)
  ncount <- DESeq2::counts(dds, normalized = TRUE)   # Normalied counts

  ## Merge with normalized count data
  resdata <- merge(as.data.frame(ncount),
                   as.data.frame(res),
                   by = "row.names",
                   sort = FALSE)
  resdata <- resdata[order(resdata$padj), ]
  names(resdata)[1] <- "Gene"

  # save data to file
  write.csv(resdata, de_count, quote = TRUE, row.names = TRUE)
  df_csv <- deseqCsvReader(de_count)
  readr::write_delim(df_csv, de_xls, delim = "\t", col_names = TRUE)

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
  png(de_plots[5],  width = 1500, height = 1500, res = 300)
  DESeq2::plotDispEsts(dds)
  dev.off()

  # Sample distance
  png(de_plots[6], width = 1000, height = 1000, res = 300)
  sampleDists <- dist(t(SummarizedExperiment::assay(rld)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- rld$condition
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(RColorBrewer::brewer.pal(9, "Blues")) )(255)
  pheatmap::pheatmap(sampleDistMatrix,
                     clustering_distance_rows = sampleDists,
                     clustering_distance_cols = sampleDists,
                     col = colors)
  dev.off()

  # Top genes
  png(de_plots[7], width = 1200, height = 1200, res = 300)
  select <- order(rowMeans(DESeq2::counts(dds, normalized = TRUE)),
                  decreasing = TRUE)[1:30]
  ma <- SummarizedExperiment::assay(ntd)[select, ]
  df <- as.data.frame(SummarizedExperiment::colData(dds)[, c("condition")])
  colnames(df) <- "condition"
  rownames(df) <- colnames(ma)
  pheatmap::pheatmap(SummarizedExperiment::assay(ntd)[select,],
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     show_rownames = FALSE,
                     annotation_col = df)
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
deseqImportFromMatrix <- function(ma, smp_control, smp_treatment){
  # ma should be matrix
  stopifnot(inherits(ma, "matrix"))

  # convert matrix to non-negative integers
  ma <- round(ma, digits = 0) #

  # remove negative columns
  neg_cols <- apply(ma < 0, 1, any)
  if(sum(neg_cols) > 0){
    warning(paste0("Remove ", sum(neg_cols), " records, with negative counts"))
    ma <- ma[! neg_cols, ]
  }

  # check samples names in ma
  stopifnot(all(c(smp_control, smp_treatment) %in% colnames(ma)))

  # create design
  colData <- deseqDesign(ma, smp_control, smp_treatment)

  # import from matrix
  # DESeq2::DESeqDataSetFromMatrix
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = ma,
    colData = colData,
    design = ~condition)

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
deseqDesign <- function(ma, smp_control, smp_treatment) {
  # build design data.frame
  # sample_name/condition
  # rownames (samples)
  # colnames (condition), ctl, exp
  stopifnot(all(c(smp_control, smp_treatment) %in% colnames(ma)))
  # conditions
  df <- data.frame(condition = factor(
    c(rep("control", times = length(smp_control)),
      rep("treatment", times = length(smp_treatment))),
    levels = c("control", "treatment")))
  # assign sample names
  rownames(df) <- c(smp_control, smp_treatment)
  return(df)
}


#' read output of DESeq2_run function
#' csv format
#'
#' @param file Path to csv file
#'
#' @import readr
#' @import dplyr
#'
#' @export
deseqCsvReader <- function(x) {
  # read DESeq2 csv file
  # df <- readr::read_csv(x, col_types = readr::cols())
  df <- read.csv(x) %>%
    dplyr::select(-1) %>%
    dplyr::rename(id = Gene) %>%
    dplyr::mutate(id = as.character(id))

  # add mean columns
  # df2 <- DESeq2_add_mean(df)
  df2 <- deseqCsvMean(df)

  return(df2)
}


#' cal mean of replicates
#' default, 2 replicates for each condition
#'
#' @param data A data fraem of DESeq2 output, first 5 columns:
#'   <id> <a-rep1> <a-rep2> <b-rep1> <b-rep2> ...
#'
#' @export
deseqCsvMean <- function(data){
  # only for DESeq2 output csv
  # insert to: col-2, col-3
  # conserved columns
  stopifnot(is.data.frame(data))
  col_defaults <- c("id", "baseMean", "log2FoldChange", "lfcSE",
                    "stat", "pvalue", "padj")
  stopifnot(all(col_defaults %in% colnames(data)))

  # sample names
  smp_names  <- colnames(data)[! colnames(data) %in% col_defaults]
  smp_groups <- unique(fqName(smp_names, rm_rep = TRUE))
  stopifnot(length(smp_groups) == 2)

  # calculate mean
  ctl_mean <- data %>%
    dplyr::select(contains(smp_groups[1])) %>%
    rowMeans()
  exp_mean <- data %>%
    dplyr::select(contains(smp_groups[2])) %>%
    rowMeans()
  # assemble
  df1 <- data.frame(
    id = data$id,
    ctl = ctl_mean,
    exp = exp_mean,
    stringsAsFactors = FALSE)

  # change names
  colnames(df1) <- c("id", smp_groups)

  # output
  df2 = merge(df1, data, by = "id")

  return(df2)
}


#' import count.txt as matrix
#'
#' read featureCounts, report as matrix
#'
#' @param count_ctl vector, list of count.txt files for control
#' @param count_exp vector, list of count.txt files for treatment
#'
#' @export
countToMatrix <- function(count_ctl, count_exp) {
  if(inherits(count_ctl, "data.frame") & inherits(count_exp, "data.frame")) {
    stopifnot("id" %in% colnames(count_ctl))
    stopifnot("id" %in% colnames(count_exp))
    df  <- merge(count_ctl, count_exp, by = "id")
    ## names
    smp_ctl <- colnames(count_ctl[, -id])
    smp_exp <- colnames(count_exp[, -id])
  } else if(inherits(count_ctl, "character") & inherits(count_exp, "character")) {
    df1 <- fcReader(count_ctl)
    df2 <- fcReader(count_exp)
    df  <- merge(df1, df2, by = "id")
    ## names
    smp_ctl <- colnames(df1)[!colnames(df1) == "id"]
    smp_exp <- colnames(df2)[!colnames(df2) == "id"]
  }
  ## convert df to matrix
  ma <- df %>%
    tibble::column_to_rownames("id") %>%
    as.matrix.data.frame()

  # check names
  return(list(ma = ma, smp_ctl = smp_ctl, smp_exp = smp_exp))
}











#' add sig label to data.frame
#'
#' @param data data.frame
#' @param type all, both, up, down
#' @param fcCutoff numeric, cutoff for foldchange
#' @param pvalCutoff numeric, cutoff for pvalue
#' @param type string, all, both, up, down, not, default: all
#'
#' @import dplyr
#'
#' @export
#'
deseqSig <- function(data, fcCutoff = 2, pvalCutoff = 0.05, type = "all") {
  stopifnot(inherits(data, "data.frame"))
  stopifnot(all(c("id", "log2FoldChange", "pvalue") %in% colnames(data)))
  if ("padj" %in% colnames(data)) {
    data$pvalCheck <- data$padj
  } else {
    data$pvalCheck <- data$pvalue
  }

  type   <- tolower(type) # lower case
  log2FC <- log2(as.numeric(fcCutoff))

  # columns required
  columns1 <- c("id", "log2FoldChange", "pvalue")

  # split group
  df1 <- dplyr::filter(data, log2FoldChange >= log2FC & pvalCheck <= pvalCutoff)
  df2 <- dplyr::filter(data, log2FoldChange <= -log2FC & pvalCheck <= pvalCutoff)
  df3 <- dplyr::filter(data, ! id %in% c(df1$id, df2$id))

  if (nrow(df1) > 0 ) {
    df1$sig <- "up"
  }

  if (nrow(df2) > 0) {
    df2$sig <- "down"
  }

  if (nrow(df3) > 0) {
    df3$sig <- "not"
  }

  # up/down
  if (type == "up") {
    df <- df1
    # df <- df %>% dplyr::filter(sig == "up")
  } else if (type == "down") {
    df <- df2
    # df <- df %>% dplyr::filter(sig == "down")
  } else if (type == "both") {
    df <- dplyr::bind_cols(df1, df2)
    # df <- df %>% dplyr::filter(sig %in% c("up", "down"))
  } else if (type == "not") {
    df <- df3
    # df <- df %>% dplyr::filter(sig %in% c("not"))
  } else {
    df <- dplyr::bind_rows(df1, df2, df3)
  }

  # remove pvalCheck column
  df <- dplyr::select(df, - pvalCheck)

  # return
  return(df)
}






