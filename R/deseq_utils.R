#' Functions for DESeq2 analysis
#'
#' Reading files, directories, pipeline, ...
#' Processing data
#' Prepare for plot
#' Plotting
#'
#' @name deseq_utils


#' @describeIn hiseq_prep_deseq Prepare data for DESeq analysis
#'
#' @param x path to the directory of rnaseq_rx
#' @param strandness character could be "sens", "anti", default "sens"
#' @param fix_batch bool fix batch effect, default: TRUE
#'
#' @import readr
#' @import configr
#' @import dplyr
#'
#' @export
hiseq_prep_deseq <- function(x, strandness = "sens", fix_batch = TRUE) {
  if(!is_hiseq_dir(x, "rnaseq_rx")) {
    fs <- ifelse(is_hiseq_dir(x), read_hiseq(x)$hiseq_type, "NULL")
    warning(glue::glue("x, expect `rnaseq_rx`, got {fs}"))
    return(NULL)
  }
  #-- Check: strand
  count_tag = ifelse(strandness == "sens", "count_sens",
                     ifelse(strandness == "anti", "count_anti", NULL))
  if(is.null(count_tag)) {
    warning(glue::glue("illegal strandness, expect: [\"sens\", \"anti\"], ",
                       "got {strandness}"))
    return(NULL)
  }
  #-- Check: required data : wt
  wt_dir  <- list_hiseq_file(x, "wt_dir", "_rx")
  wt_name <- list_hiseq_file(x, "wt_name", "_rx")
  wt_data <- data.frame(
    files = list_hiseq_file(wt_dir, count_tag, "r1"),
    names = list_hiseq_file(wt_dir, "smp_name", "r1")
  )
  #-- run: sanitize, suffix
  wt_suffix <- deseq_sanitize_str(wt_data$names, 20)
  #-- Check: required data : mut
  mut_dir  <- list_hiseq_file(x, "mut_dir", "_rx")
  mut_name <- list_hiseq_file(x, "mut_name", "_rx")
  mut_data <- data.frame(
    files = list_hiseq_file(mut_dir, count_tag, "r1"),
    names = list_hiseq_file(mut_dir, "smp_name", "r1")
  )
  #-- run: sanitize, suffix
  condition <- deseq_sanitize_str(c(wt_name, mut_name), 20) # wt, mut
  mut_suffix <- deseq_sanitize_str(mut_data$names, 20)
  #-- Check: condition, batch
  fcdata <- rbind(wt_data, mut_data)
  fcdata$condition <- c(rep(condition[1], length(wt_suffix)),
                        rep(condition[2], length(mut_suffix)))
  fcdata$names <- paste0(fcdata$condition, ".", c(wt_suffix, mut_suffix))
  fcdata$condition <- factor(fcdata$condition, levels = condition)
  if(isTRUE(fix_batch)) {
    fcdata$batch <- as.factor(
      c(LETTERS[seq_len(nrow(wt_data))], LETTERS[seq_len(nrow(mut_data))])
      )
  }
  #-- run: load to dds
  dds <- import_featurecounts(fcdata)
  DESeq2::DESeq(dds)
}


# --Utils: Prepare data --------------------------------------------------------

#' @describeIn import_featurecounts Construct dds for DESeq2 analysis using matrix
#'
#' using DESeqDataSetFromMatrix()
#'
#' @param x data.frame require "files", "names", "condition" columns
#'
#' @return `DESeqDataSet`
#'
#' @export
import_featurecounts <- function(x) {
  x <- valid_featurecounts_input(x)
  if(is.null(x)) {
    return(NULL)
  }
  #-- run: load count.txt data
  df <- lapply(x$files, function(i) {
    d <- tryCatch(
      {
        read_fc(i) %>%
          dplyr::select(1:2) # only the first bam record
      },
      error = function(cond) {
        message(glue::glue("failed to read featureCounts file: {i}"))
        return(NULL)
      }
    )
  }) %>%
    bind_cols2()
  #-- run: rename data.frame, conver to matrix
  colnames(df) <- c("id", x$names) # rename
  ma <- tibble::column_to_rownames(df, "id") %>%
    as.matrix() # convert to matrix
  #-- Check: positive values
  rm_rows <- apply(ma < 0, 1, any)
  ma <- ma[!rm_rows, ] # remove rows, with negative values
  if(sum(rm_rows) > 0) {
    message(glue::glue("removing {sum(rm_rows)} rows, with negative values"))
  }
  ma <- round(ma) # to int
  #-- run: coldata
  coldata <- data.frame(
    condition = x$condition,
    row.names = x$names
  )
  if(rlang::has_name(x, "batch")) {
    coldata$batch <- x$batch
  }
  #-- run: design (formula)
  if(rlang::has_name(coldata, "batch")) {
    message("run DESeq2 design: `~ condition + batch`")
    fo <- formula(~ batch + condition)
  } else {
    message("run DESeq2 design: `~ condition`")
    fo <- formula(~ condition)
  }
  tryCatch(
    {
      DESeq2::DESeqDataSetFromMatrix(
        countData = ma,
        colData   = coldata,
        design    = fo
      )
    },
    error = function(cond) {
      warning("failed to import featureCounts")
      return(NULL)
    }
  )
}


#' @describeIn valid_featurecounts_input
#' for function `import_featurecounts()`
#'
#' @param x data.frame
#'
#' @example
#'
#' @export
valid_featurecounts_input <- function(x) {
  #-- Check: x data.frame, required columns: files, names
  if(!inherits(x, "data.frame")) {
    warning(glue::glue("x=, expect `data.frame`, got {class(x)}"))
    return(NULL)
  }
  #-- Check: x column names
  rc <- c("files", "names", "condition")
  if(!all(rc %in% colnames(x))) {
    fb <- rc %in% colnames(x)
    fs <- paste(rc[!fb], collapse = ", ")
    warning(glue::glue("missing columns: []"))
    return(NULL)
  }
  #-- Check: file exists
  if(!all(file.exists(x$files))) {
    fb <- file.exists(x$files)
    fs <- paste(x$files[!fb], collapse = ", ")
    warning(glue::glue("file not exists, check `x$files`: [{fs}]"))
    return(NULL)
  }
  #-- Check: names, unique
  if(any(duplicated(x$names))) {
    fb <- x$names[duplicated(x$names)]
    fs <- paste(fb, collapse = ", ")
    warning(glue::glue("duplicate `names` not allowed, [{fs}]"))
    return(NULL)
  }
  #-- Check: condition, factor
  if(!inherits(x$condition, "factor")) {
    warning("`condition` are {class(x$condition)}, converting to factors")
    x$condition <- as.factor(x$condition)
  }
  #-- Check: at least 2 replicates for each condition
  f <- as.data.frame(table(x$condition))
  f <- f[f$Freq < 2, ]
  if(nrow(f) > 0) {
    fs <- paste(x$condition, collapse = ", ")
    warning(glue::glue("at least 2 rep required, check `x$condition`: [{fs}]"))
    return(NULL)
  }
  #-- Check: batch
  if(rlang::has_name(x, "batch")) {
    f <- as.data.frame(table(x$batch))
    f <- f[f$Freq < 2, ]
    if(nrow(f) > 0) {
      fs <- paste(x$batch, collapse = ", ")
      warning(glue::glue("at least 2 rep required, check `x$batch`: [{fs}]"))
      return(NULL)
    }
  }
  x
}


#' @describeIn filt_sig_gene
#'
#' filt sig genes by `sig` column
#'
#' @param x data.frame, DESeqResults, matrix
#' @param type character could be combination of ["up", "down", "not", "sig"],
#' default: "sig"
#' @param fc numeric, cutoff for foldchange, default: 2
#' @param pvalue numeric, cutoff for pvalue, default: 0.05
#' @param p_adjust bool use p-adjust value instead
#'
#' @import dplyr
#'
#' @export
filt_sig_gene <- function(x, type = "sig", fc = 2, pvalue = 0.05,
                          p_adjust = TRUE) {
  # check input `x`, data.frame
  df <- get_sig_name(x, fc, pvalue, p_adjust, return_dataframe = TRUE)
  if(inherits(df, "data.frame")) {
    sig_type <- switch (
      type,
      "sig"  = c("up", "down"),
      "all"  = c("up", "not", "down"),
      "up"   = "up",
      "down" = "down",
      "not"  = "not"
    )
    out <- df[df$sig %in% sig_type, ]
  } else {
    warning("unknown data, expect `data.frame`, got {class(x)}")
    out <- NULL
  }
  out
}


#' @describeIn get_sig_name
#'
#' add sig name
#'
#' @param data data.frame, csv, xls output of `DESeq2::results(dds)`
#' @param fc numeric, cutoff for foldchange, default: 2
#' @param pvalue numeric, cutoff for pvalue, default: 0.05
#' @param p_adjust bool use p-adjust value instead
#'
#' @import dplyr
#'
#' @return vector, sig names
#'
#' @export
get_sig_name <- function(x, fc = 2, pvalue = 0.05, p_adjust = TRUE,
                         return_dataframe = FALSE) {
  if(inherits(x, "character")) {
    if(file.exists(x)) {
      if(endsWith(x, ".csv")) {
        df <- read.csv(x)
      } else if(endsWith(x, ".xls")) {
        df <- readr::read_delim(x, "\t", col_names = TRUE, col_types = readr::cols)
      } else {
        warning(glue::glue("unknown x, expect `.csv`, `.xls` file, get: {x}"))
        return(NULL)
      }
    } else {
      warning(glue::glue("x, file not exists: {x}"))
    }
  } else if(inherits(x, "data.frame")) {
    df <- x
  } else if(inherits(x, "DESeqResults")) {
    df <- as.data.frame(x)
  } else if(inherits(x, "matrix")) {
    df <- as.data.frame(x)
  } else {
    warning(glue::glue("unknown data, expect `data.frame`, get: {class(x)}"))
    return(NULL)
  }
  df <- tibble::tibble(df)
  #-- Check: required columns
  rc <- c("log2FoldChange", "pvalue", "padj")
  if(!all(rc %in% names(df))) {
    warning(glue::glue(
      "missing required columns, [{paste(rc, collapse = ', ')}]")
    )
    return(NULL)
  }
  #-- run: add sig mark
  # fix NA: log2fc NA -> 0, pval NA -> 1
  log2fc <- dplyr::pull(df, "log2FoldChange")
  if(isTRUE(p_adjust)) {
    pval <- dplyr::pull(df, "padj")
  } else {
    pval <- dplyr::pull(df, "pvalue")
  }
  log2fc[is.na(log2fc)] <- 0 # 1,2,4
  pval[is.na(pval)]     <- 1 # 1,2,4
  #-- run: add sig; up=1, down=-1, not=0
  up   <- as.numeric(pval < pvalue & log2fc >= log2(fc))   # +1
  down <- -as.numeric(pval < pvalue & log2fc <= -log2(fc)) # -1
  sig  <- up + down
  #-- run: convert to marks
  sig[sig == 1]  = "up"
  sig[sig == -1] = "down"
  sig[sig == 0]  = "not"
  #-- return
  if(isTRUE(return_dataframe)) {
    df$sig <- sig
    out <- df
  } else {
    out <- sig
  }
  out
}


#' @describeIn deseq_mean
#' calculate the mean values for each group
#'
#' get the following data from `colData(dds)`
#' condition
#' rownames
#'
#' @param x data.frame norm table, see `deseq()`
#' @param dds DESeqDataSet, parsing design from `colData(dds)`
#'
#' @return data.frame
#'
#' @export
deseq_mean <- function(dds, outdir = NULL) {
  #-- Check: dds
  if(!inherits(dds, "DESeqDataSet")) {
    warning(glue::glue(
      "illegal input, ",
      "dds is {class(dds)} (expect `DESeqDataSet`)"))
    return(NULL)
  }
  #-- Check: norm_table.csv, `run_deseq()`
  df <- NULL #
  if(inherits(outdir, "character")) {
    norm_table <- file.path(outdir, "norm_table.csv")
    if(file.exists(norm_table)) {
      message(glue::glue("loading data from: {norm_table}"))
      df <- read.csv(norm_table)
    }
  }
  #-- run: `run_deseq()`
  if(!inherits(df, "data.frame")) {
    dd  <- run_deseq(dds, outdir)
    df1a <- DESeq2::counts(dds, normalized = TRUE) # normalized counts
    df  <- merge(as.data.frame(df1a), as.data.frame(dd$res), by = "row.names")
    colnames(df)[1] <- "gene_id"
  }
  #-- run: choose columns
  coldata <- SummarizedExperiment::colData(dds)
  wt  <- levels(coldata$condition)[1]
  mut <- levels(coldata$condition)[2]
  wt_names  <- rownames(coldata[coldata$condition == wt,])
  mut_names <- rownames(coldata[coldata$condition == mut,])
  df %>%
    dplyr::mutate(
      !! wt := dplyr::select(., all_of(wt_names)) %>% rowMeans(),
      !! mut := dplyr::select(., all_of(mut_names)) %>% rowMeans()) %>%
    dplyr::select(gene_id, all_of(c(wt, mut)), all_of(names(df)[-1]))
}




#' @describeIn set_readable
#' add symbol to data.frame, by `gene_id` column
#' or use `row.names` values
#'
#' @param x data.frame contains `gene_id`
#' alternative,`Gene`, `id`, or `rown.names()`
#' @param genome character name of the genome,
#' eg: dm6, hg38, mm10
#'
#' @return data.frame
#'
#' @export
set_readable <- function(x, genome, keytype = "auto") {
  #-- Check: arguments
  if(inherits(x, "data.frame")) {
    if(all(c("entrez", "symbol") %in% names(x))) {
      message(glue::glue("column `symbol` and `entrez` already exists"))
      return(x)
    }
  } else {
    message(glue::glue("x is {class(x)}, expect `data.frame`"))
    return(x)
  }
  if(!is_valid_organism(genome)) {
    message(glue::glue("unknown gnome: {genome}"))
    return(x)
  }
  if(inherits(keytype, "character")) {
    if(!is_valid_keytype(keytype, organism = genome)) {
      keytype = NULL
    }
  } else {
    keytype = NULL
  }
  #-- Check: gene_id column
  rc <- c("gene_id", "Gene", "id", "gene_name")
  rc_str <- paste(rc, collapse = ", ")
  rc <- rc[rc %in% names(x)]
  if(length(rc) == 0) {
    if(.row_names_info(x) > 0) {
      x <- cbind(gene_id = rownames(x), x) # row.names to "gene_id"
      gid <- "gene_id"
    } else {
      message(glue::glue("missing gene id column: {rc_str}"))
      return(x)
    }
  } else {
    gid <- rc[1]
  }
  g <- as.character(x[, gid]) # gene names
  g_str <- paste(head(g, 5), collapse = ", ") # example
  message(glue::glue(
    "convert `{gid}` to `ENTREZID` and `SYMBOL`: {g_str} ..."))
  #-- run: guess keytype
  if(!inherits(keytype, "character")) {
    keytype <- tryCatch(
      {
        guess_keytype(g, organism = genome)
      },
      error = function(cond) {
        warning(glue::glue("unknown genes for [{genome}]: {g_str} ..."))
        return(NULL)
      }
    )
  }
  #-- run: convert table
  if(inherits(keytype, "character")) {
    gdf <- convert_id(g, from_keytype = keytype,
                      to_keytype = c("ENTREZID", "SYMBOL"),
                      organism   = genome, na_rm = FALSE)
    #-- run: update genes
    out <- dplyr::left_join(
      x, gdf, by = setNames(keytype, nm = gid)
    )
  } else {
    message(glue::glue("unknown genes for [{genome}]: {g_str} ..."))
  }
  out
}




#' @describeIn sanitize_str
#' for coef, only allow
#' letters, numbers, '_' and '.'
#' convert other characters to '.'
#'
#' fix the sample names by length
#' trim to <= 20 characters
#'
#' @param x character
#'
#' @return character
#'
#' @export
deseq_sanitize_str <- function(x, n_max = 0) {
  x <- as.character(x)
  if(inherits(x, "character")) {
    out <- gsub("[^\\w\\.]", ".", x, perl = TRUE)
    if(nchar(x[1]) > n_max & n_max > 0) {
      # longest prefix, suffix
      lcp <- lcPrefix(out, ignore.case = FALSE)
      lcs <- lcSuffix(out, ignore.case = FALSE)
      if(nchar(lcp) > 0) {
        out <- gsub(lcp, "", out)
      }
      if(nchar(lcs) > 0) {
        out <- gsub(lcs, "", out)
      }
    }
    out
  } else {
    x
  }
}


