#' Functions for DESeq2 analysis
#'
#' Reading files, directories, pipeline, ...
#' Processing data
#' Prepare for plot
#' Plotting
#'
#' @name deseq_utils


#' @describeIn hiseq_prep_deseq4fc
#' Prepare data for DESeq analysis, for featureCounts output
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
hiseq_prep_deseq4fc <- function(x, strandness = "sens", fix_batch = TRUE) {
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
  wt_suffix <- deseq_sanitize_str(wt_data$names, 20, fix_prefix = FALSE)
  wt_suffix <- paste0("rep", wt_suffix) # rep1, rep2
  #-- Check: required data : mut
  mut_dir  <- list_hiseq_file(x, "mut_dir", "_rx")
  mut_name <- list_hiseq_file(x, "mut_name", "_rx")
  mut_data <- data.frame(
    files = list_hiseq_file(mut_dir, count_tag, "r1"),
    names = list_hiseq_file(mut_dir, "smp_name", "r1")
  )
  #-- run: sanitize, suffix
  condition  <- deseq_sanitize_str(c(wt_name, mut_name), 20, fix_prefix = TRUE)
  mut_suffix <- deseq_sanitize_str(mut_data$names, 20, fix_prefix = FALSE)
  mut_suffix <- paste0("rep", mut_suffix) # rep1, rep2
  #-- Check: condition, batch
  fcdata <- rbind(wt_data, mut_data)
  fcdata$smp_name  <- c(wt_data$names, mut_data$names)
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


#' @describeIn hiseq_prep_deseq4salmon
#' Prepare data for DESeq analysis, from salmon results; rnaseq_rx
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
hiseq_prep_deseq4salmon <- function(x, fix_batch = TRUE) {
  if(!is_hiseq_dir(x, "rnaseq_rx")) {
    fs <- ifelse(is_hiseq_dir(x), read_hiseq(x)$hiseq_type, "NULL")
    warning(glue::glue("x, expect `rnaseq_rx`, got {fs}"))
    return(NULL)
  }
  #-- Check: required data : wt
  wt_name  <- list_hiseq_file(x, "wt_name", "_rx")
  wt_dir   <- list_hiseq_file(x, "wt_dirs", "_rx")
  wt_data  <- data.frame(
    files = list_hiseq_file(x, "wt_quant", "_rx"),
    names = sapply(wt_dir, function(i) {
      list_hiseq_file(i, "smp_name",  TRUE)
    }, USE.NAMES = FALSE)
  )
  #-- run: sanitize, suffix
  wt_suffix <- deseq_sanitize_str(wt_data$names, 10, fix_prefix = FALSE)
  wt_suffix <- paste0("rep", wt_suffix) # rep1, rep2
  #-- Check: required data : mut
  mut_name  <- list_hiseq_file(x, "mut_name", "_rx")
  mut_dir   <- list_hiseq_file(x, "mut_dirs", "_rx")
  mut_data  <- data.frame(
    files = list_hiseq_file(x, "mut_quant", "_rx"),
    names = sapply(mut_dir, function(i) {
      list_hiseq_file(i, "smp_name", TRUE)
    }, USE.NAMES = FALSE)
  )
  mut_suffix <- deseq_sanitize_str(mut_data$names, 10, fix_prefix = FALSE)
  mut_suffix <- paste0("rep", mut_suffix) # rep1, rep2
  #-- prepare
  condition  <- deseq_sanitize_str(c(wt_name, mut_name), 10, fix_prefix = TRUE)
  #-- Check: condition, batch
  fcdata <- rbind(wt_data, mut_data)
  fcdata$smp_name  <- c(wt_data$names, mut_data$names)
  fcdata$condition <- c(rep(condition[1], length(wt_suffix)),
                        rep(condition[2], length(mut_suffix)))
  fcdata$names <- paste0(fcdata$condition, ".", c(wt_suffix, mut_suffix))
  fcdata$condition <- factor(fcdata$condition, levels = condition)
  if(isTRUE(fix_batch)) {
    fcdata$batch <- as.factor(
      c(LETTERS[seq_len(nrow(wt_data))], LETTERS[seq_len(nrow(mut_data))])
    )
  }
  #----------------------------------------------------------------------------#
  # tx2gene table, in index
  salmon_index <- list_hiseq_file(x, "salmon_index", "rx")
  if(is(salmon_index, "character")) {
    tx2gene_csv <- file.path(salmon_index, "tx2gene.csv")
    if(file.exists(tx2gene_csv)) {
      tx2gene <- read.csv(tx2gene_csv)
    } else {
      warning(glue::glue(
        "tx2gene.csv is {tx2gene_csv}, file not exists"
      ))
      return(NULL)
    }
  } else {
    warning(glue::glue(
      "salmon_index is {class(salmon_index)}, expect 'character'"
    ))
    return(NULL)
  }
  #-- run: load to dds
  dds <- import_salmon(fcdata, tx2gene)
  DESeq2::DESeq(dds)
}


# --Utils: Prepare data -------------------------------------------------------#

#' @describeIn import_featurecounts Construct dds for DESeq2 analysis using matrix
#'
#' using DESeqDataSetFromMatrix()
#'
#' @param x data.frame require "files", "smp_name", "condition" columns
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
    smp_name  = x$smp_name,
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
  #-- run: import data
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

#' @describeIn import_samlmon
#' Construct dds for DESeq2 analysis using salmon output
#'
#' using DESeqDataSetFromMatrix()
#'
#' @param x data.frame require "files", "names", "condition" columns
#'
#' @return `DESeqDataSet`
#'
#' @export
import_salmon <- function(x, tx2gene) {
  #----------------------------------------------------------------------------#
  #-- Check: arguments
  x <- valid_featurecounts_input(x)
  if(inherits(tx2gene, "data.frame")) {
    t2g <- all(c("TXNAME", "GENEID") %in% names(tx2gene))
  } else {
    t2g <- FALSE
  }
  if(!inherits(x, "data.frame") | !t2g) {
    warning(glue::glue(
      "x is {class(x)}, expect 'data.frame', ",
      "with 'files', 'names', 'smp_name', 'condition', ",
      "tx2gene is {class(tx2gene)}, expect 'data.frame', ",
      "with 'TXNAME', 'GENEID'"
    ))
    return(NULL)
  }
  #----------------------------------------------------------------------------#
  #-- run: coldata
  coldata <- data.frame(
    condition = x$condition,
    smp_name  = x$smp_name,
    row.names = x$names
  )
  if('batch' %in% names(x)) {
    coldata$batch <- x$batch
    message("run DESeq2 design: `~ condition + batch`")
    fo <- formula(~ batch + condition)
  } else {
    message("run DESeq2 design: `~ condition`")
    fo <- formula(~ condition)
  }
  #----------------------------------------------------------------------------#
  #-- run: import salmon data
  salmon_files <- setNames(x[["files"]], nm = row.names(coldata))
  txi <- tximport::tximport(salmon_files, type = "salmon", tx2gene = tx2gene)
  #-- run: to DESeq2
  tryCatch(
    {
      DESeqDataSetFromTximport(txi, colData = coldata, design  = fo)
    },
    error = function(cond) {
      warning("failed to import salmon")
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
#' @param type character could be combination of
#' ["all", "sig", "up", "down", "not"],
#' default: "sig"
#' @param fc numeric, cutoff for foldchange, default: 2
#' @param pvalue numeric, cutoff for pvalue, default: 0.05
#' @param p_adjust bool use p-adjust value instead
#' @param force logical force calculate sig, default: FALSE
#'
#' @import dplyr
#'
#' @export
filt_sig_gene <- function(x, type = "sig", ...) {
  #----------------------------------------------------------------------------#
  #-- Check: default values
  dots <- rlang::list2(...)
  args <- rlang::list2(
    fc         = 2,
    pvalue     = 0.05,
    p_adjust   = TRUE,
    force      = FALSE,
    overwrite  = FALSE,
    return_dataframe = TRUE
  )
  #-- update dots, for child functions
  dots_args <- lapply(names(args), function(i) {
    if(!i %in% names(dots)) {
      args[i]
    }
  })
  dots <- c(dots, unlist(dots_args, recursive = FALSE, use.names = TRUE))
  #-- update global
  for(name in names(dots)) {
    assign(name, dots[[name]])
  }
  #----------------------------------------------------------------------------#
  dots[["return_dataframe"]] <- TRUE # force
  df <- get_sig_name(x, !!!dots)
  # df <- get_sig_name(x, fc, pvalue, p_adjust, return_dataframe = TRUE)
  if(inherits(df, "data.frame")) {
    sig_type <- switch (
      type,
      "sig"  = c("up", "down"),
      "all"  = c("up", "not", "down"),
      "up"   = "up",
      "down" = "down",
      "not"  = "not"
    )
    if(inherits(sig_type, "character")) {
      #-- return: subset
      df[df$sig %in% sig_type, ]
    }
  } else {
    warning("unknown data, expect `data.frame`, got {class(x)}")
  }
}


#' @describeIn get_sig_name
#'
#' add sig name, based on fc:foldchange (not log2), pvalue
#'
#' x could be `data.frame`, `.csv`
#' required columns: `log2FoldChange`, `pvlaue`, `padj`
#'
#'
#' @param x data.frame, csv, xls output of `DESeq2::results(dds)`
#' @param fc numeric, cutoff for foldchange, default: 2
#' @param pvalue numeric, cutoff for pvalue, default: 0.05
#' @param p_adjust logical use p-adjust value instead
#' @param return_dataframe logical return a new data.frame, contain `sig` col
#' @param force logical force calculate sig, default: FALSE
#' @param .col_sig character name of the column for `sig`, default: "sig"
#' @param .col_log2fc character name of the log2fc column, default: "log2FoldChange"
#' @param .col_pvalue character name of the pvalue column, default: "pvalue";
#' @param .col_padj character name of the pvalue column, default: "padj";
#' @param .sig_up character assign name for up-regulated genes, default: "up"
#' @param .sig_down character assign name for down-regulated genes, default: "down"
#' @param .sig_not character assign name for not changed genes, default: "not"
#'
#' @import dplyr
#'
#' @return vector, sig names
#'
#' @export
get_sig_name <- function(x, ...) {
  #----------------------------------------------------------------------------#
  #-- Check: default values
  dots <- rlang::list2(...)
  args <- rlang::list2(
    fc         = 2,
    pvalue     = 0.05,
    p_adjust   = TRUE,
    force      = FALSE,
    overwrite  = FALSE,
    return_dataframe = TRUE,
    .col_sig   = "sig",
    .col_log2fc = "log2FoldChange",
    .col_pvalue = "pvalue",
    .col_padj   = "padj",
    .sig_up     = "up",
    .sig_down   = "down",
    .sig_not    = "not"
  )
  #-- update dots, for child functions
  dots_args <- lapply(names(args), function(i) {
    if(!i %in% names(dots)) {
      args[i]
    }
  })
  dots <- c(dots, unlist(dots_args, recursive = FALSE, use.names = TRUE))
  #-- update global
  for(name in names(dots)) {
    assign(name, dots[[name]])
  }
  #----------------------------------------------------------------------------#
  #-- Check: arguments
  if(inherits(fc, "numeric") & inherits(pvalue, "numeric")) {
    if(fc <= 0 | pvalue <= 0 | pvalue > 1) {
      warning(glue::glue("`fc={fc}`, expect (0-Inf), greater than 0; ",
                         "`pvalue={pvalue}`, expect (0-1)"))
      return(NULL)
    }
  } else {
    warning(glue::glue("`fc` is {class(fc)}, expect `numeric`, ",
                       "`pvalue` is {class(pvalue)}, expect `numeric`"))
    return(NULL)
  }
  if(!inherits(p_adjust, "logical")) {
    warning(glue::glue("`p_value` is {class(pvalue)}, expect logical"))
    return(NULL)
  }
  if(!inherits(return_dataframe, "logical")) {
    warning(glue::glue(
      "`return_dataframe` is {class(return_dataframe)}, expect logical"))
    return(NULL)
  }
  if(!inherits(force, "logical")) {
    warning(glue::glue(
      "`force` is {class(force)}, expect logical"))
    return(NULL)
  }
  #-- Check: arguments
  df <- NULL
  if(inherits(x, "character")) {
    if(file.exists(x)) {
      if(endsWith(x, ".csv")) {
        df <- read.csv(x)
      } else if(endsWith(x, ".xls")) {
        df <- readr::read_delim(x, "\t", col_names = TRUE, col_types = readr::cols)
      } else {
        warning(glue::glue("x is {class(x)}, expect `.csv` file"))
        # return(NULL)
      }
    } else {
      warning(glue::glue("file not exists, `x=`: {x}"))
    }
  } else if(inherits(x, "data.frame")) {
    df <- x
  } else if(inherits(x, "matrix")) {
    df <- as.data.frame(x)
  } else if(inherits(x, "DESeqResults")) {
    df <- as.data.frame(x)
  } else {
    warning(glue::glue("x is {class(x)}, expect `.csv` or `data.frame` file"))
    # return(NULL)
  }
  #-- Check: data.frame
  if(!inherits(df, "data.frame")) {
    return(NULL)
  }
  #-- Check: required columns
  df <- tibble::tibble(df)
  if(.col_sig %in% names(df) & !force) {
    message(glue::glue(
      "column `{.col_sig}` exists, set `force=TRUE` to re-calculate sig"))
    if(return_dataframe) {
      out <- df
    } else {
      out <- df[[.col_sig]]
    }
    return(out)
  }
  #-- run: re-run, sig values
  # .col_pvalue = ifelse(p_adjust, "padj", "pvalue")
  .pvalue_name <- ifelse(p_adjust, .col_padj, .col_pvalue)
  rc <- c(.col_log2fc, .pvalue_name)
  # if(isTRUE(p_adjust)) {
  #   rc <- c(.col_log2fc, .col_padj)
  # } else {
  #   rc <- c(.col_log2fc, .col_pvalue)
  # }
  if(!all(rc %in% names(df))) {
    rc_str <- paste(rc, collapse = ", ")
    warning(glue::glue("missing required columns, [{rc_str}]"))
    return(NULL)
  }
  df_sub <- df[, c(.col_log2fc, .col_pvalue)]
  #-- Check: numeric columns
  if(!all(apply(df_sub, 2, is.numeric))) {
    warning("expect numeric columns, but get:")
    str(df_sub)
    return(NULL)
  }
  #-- run: assign marks
  up   <- df[[.pvalue_name]] < pvalue & df[[.col_log2fc]] > log2(fc)
  down <- df[[.pvalue_name]] < pvalue & df[[.col_log2fc]] < -log2(fc)
  up[is.na(up)] <- FALSE
  down[is.na(down)] <- FALSE
  #-- add
  df[[.col_sig]]       <- .sig_not # init
  df[[.col_sig]][up]   <- .sig_up
  df[[.col_sig]][down] <- .sig_down
  #-- return
  if(return_dataframe) {
    df
  } else {
    df[[.col_sig]]
  }
}


#' @describeIn deseq_mean
#' calculate the mean values for each group
#'
#' get the following data from `colData(dds)`
#' condition
#' rownames
#'
#' @param x DESeqDataSet, or '.csv' parsing design from `colData(dds)`
#' @param outdir for run_deseq_res()
#' @param transform_method character for dds count transformed,
#' 'standard', 'vst', 'rlog'; default: 'standard'
#'
#'
#' @return data.frame
#'
#' @export
deseq_mean <- function(x, ...) {
  #----------------------------------------------------------------------------#
  #-- Check: default values
  dots <- rlang::list2(...)
  args <- rlang::list2(
    outdir     = NULL,
    transform_method = "standard" # vst, rlog
  )
  #-- update dots, for child functions
  dots_args <- lapply(names(args), function(i) {
    if(!i %in% names(dots)) {
      args[i]
    }
  })
  dots <- c(dots, unlist(dots_args, recursive = FALSE, use.names = TRUE))
  #-- update global
  for(name in names(dots)) {
    assign(name, dots[[name]])
  }
  #----------------------------------------------------------------------------#
  # #-- arguments
  # outdir <- NULL
  # transform_method = "standard" # for dds
  # dots <- rlang::list2(...)
  # for(name in names(dots)) {
  #   assign(name, dots[[name]])
  # }
  #-- transform
  if(inherits(transform_method, "character")) {
    if(!transform_method %in% c("standard", "vst", "rlog")) {
      warning(glue::glue(
        "transform_method is {transform_method}, expect ",
        "'standard', 'vst', 'rlog'; set to 'standard' "
      ))
      transform_method <- "standard"
    }
  } else {
    warning(glue::glue(
      "transform_method is {class(transform_method)}, expect 'character'",
      "'standard', 'vst', 'rlog'; set to 'standard' "
    ))
    transform_method <- "standard"
  }
  #-- run: support; csv
  df <- NULL
  if(inherits(x, "data.frame") | inherits(x, "character")) {
    df <- deseq_csv_mean(x) # another option
  } else if(inherits(x, "DESeqDataSet")) {
    dd <- run_deseq_res(x, outdir = outdir, shrink = FALSE, transform = TRUE) #!!!
    if(inherits(dd, "list")) {
      df1a <- dd[["dds_trans"]][[transform_method]] # log10 table
      df1  <- as.data.frame(2^assay(df1a)) # log10 -> counts
      df2  <- as.data.frame(dd$res)
      df   <- merge(df1, df2, by = "row.names")
      colnames(df)[1] <- "gene_id"
      #-- run: choose columns
      coldata <- SummarizedExperiment::colData(x)
      wt  <- levels(coldata$condition)[1]
      mut <- levels(coldata$condition)[2]
      wt_names  <- rownames(coldata[coldata$condition == wt,])
      mut_names <- rownames(coldata[coldata$condition == mut,])
      df <- df %>%
        dplyr::mutate(
          !!wt := dplyr::select(., all_of(wt_names)) %>% rowMeans(),
          !!mut := dplyr::select(., all_of(mut_names)) %>% rowMeans()) %>%
        dplyr::select(gene_id, all_of(c(wt, mut)), all_of(names(df)[-1]))
    }
  } else {
    warning(glue::glue(
      "illegal input, ",
      "x is {class(x)}, expect `DESeqDataSet`, `.csv`"))
  }
  #-- run: `run_deseq()`
  if(inherits(df, "data.frame")) {
    df
  }
}


#' @describeIn deseq_csv_mean
#' calculate the mean values
#' support old version: transcripts_deseq2.csv
#' @param x data.frame
deseq_csv_mean <- function(x) {
  #-- Check: arguments
  df <- NULL
  if(inherits(x, "data.frame")) {
    df <- x
  } else if(inherits(x, "character")) {
    if(file.exists(x)) {
      if(endsWith(x, ".csv")) {
        df <- read.csv(x)
      }
    }
  }
  if(!inherits(df, "data.frame")) {
    warning(glue::glue("'x' is {class(x)}, expect 'data.frame'"))
    return(NULL)
  }
  #-- Check: columns
  # columns before 'baseMean', remove 'gene_id', 'Gene'
  # i: baseMean
  # g: first of sample
  rc <- "baseMean"
  if(rc %in% names(df)) {
    i <- grep(rc, names(df), fixed = TRUE) # baseMean
    g <- ifelse(names(df)[1] == "X", 3, 2) # !!! support for old version, transcripts_deseq2.csv
    ix <- names(df)[g:(i - 1)]
  } else {
    warning(glue::glue("unknown x, column '{rc}' not found"))
    return(NULL)
  }
  col_gene <- names(df)[1:(g-1)] #
  col_smp  <- names(df)[-c(1:(g-1))]
  # column Gene, gene_id
  #-- unique names
  iu <- unique(fq_name(ix, fix_rep = TRUE))
  if(length(iu) != 2) {
    iu_str <- paste(iu, collapse = ", ")
    warning(glue::glue("{length(iu)} samples found, expect 2; get: {iu_str}"))
    return(NULL)
  }
  wt  <- iu[1]
  mut <- iu[2]
  wt_names  <- ix[startsWith(ix, wt)]
  mut_names <- ix[startsWith(ix, mut)]
  #-- Check: wt, mut exists or not
  if(any(c(wt, mut) %in% names(df))) {
    warning(glue::glue(
      "'x' might contains merged columns: {wt}, {mut}, \n",
      "check 'x' again"
    ))
    return(NULL)
  }
  #-- run:
  df %>%
    dplyr::mutate(
      !!wt := dplyr::select(., all_of(wt_names)) %>% rowMeans(),
      !!mut := dplyr::select(., all_of(mut_names)) %>% rowMeans()) %>%
    dplyr::select(all_of(c(col_gene, wt, mut, col_smp)))
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
set_readable <- function(x, ...) {
  #----------------------------------------------------------------------------#
  #-- Check: default values
  dots <- rlang::list2(...)
  args <- rlang::list2(
    genome       = NULL,
    keytype      = "auto",
    gene_table   = NULL,
    overwrite    = FALSE,
    .col_gene_id = 1,
    .cutoff      = 0.8
  )
  #-- update dots, for child functions
  dots_args <- lapply(names(args), function(i) {
    if(!i %in% names(dots)) {
      args[i]
    }
  })
  dots <- c(dots, unlist(dots_args, recursive = FALSE, use.names = TRUE))
  #-- update global
  for(name in names(dots)) {
    assign(name, dots[[name]])
  }
  #----------------------------------------------------------------------------#
  #-- Check: arguments
  if(!inherits(overwrite, "logical")) {
    overwrite = FALSE
  }
  if(inherits(x, "data.frame")) {
    if(all(c("ENTREZID", "SYMBOL") %in% names(x)) & !overwrite) {
      message(glue::glue("column `SYMBOL` and `ENTREZ` already exists"))
      return(x)
    }
  } else {
    message(glue::glue("x is {class(x)}, expect `data.frame`"))
    return(x)
  }
  #-- Check: x, gene_id column
  rc <- c("gene_id", "Gene", "id", "gene_name")
  rc <- rc[rc %in% names(x)]
  if(length(rc) == 0) {
    rc_str <- paste(rc, collapse = ", ")
    if(.row_names_info(x) > 0) {
      x <- cbind(gene_id = rownames(x), x) # row.names to "gene_id"
      gid <- "gene_id"
    } else {
      message(glue::glue("missing 'gene_id' column: {rc_str}"))
      return(x)
    }
  } else {
    gid <- rc[1]
  }
  g <- as.character(x[[gid]]) # gene names
  g_str <- paste(g[1:3], collapse = ", ") # example
  message(glue::glue(
    "convert `{gid}` to `ENTREZID` and `SYMBOL`: {g_str} ..."
  ))
  #----------------------------------------------------------------------------#
  #-- run: load gene table
  #-- level-1: from table
  gdf <- NULL
  if(inherits(gene_table, "character")) {
    if(file.exists(gene_table) & endsWith(gene_table, ".csv")) {
      gdf <- read.csv(gene_table)
      keytype <- names(gdf)[.col_gene_id] # 1-st column
      # check genes available
      pct <- sum(g %in% gdf[[keytype]]) / length(g)
      if(pct < .cutoff) {
        g_str <- paste(g[1:3], collapse = ", ")
        gt_str <- paste(gdf[[keytype]][1:3], collapse = ", ")
        warning(glue::glue(
          "'gene_table' not valid, no more than {.cutoff*100}% genes found; \n",
          "genes in 'x' are: {g_str} ... \n",
          "genes in 'gene_table' are: {gt_str} ... "
        ))
        return(x)
      }
    }
  }
  #-- level-2: from genome
  if(!inherits(gdf, "data.frame")) {
    if(inherits(genome, "character")) {
      if(is_valid_organism(genome)) {
        if(!is_valid_keytype(keytype, organism = genome)) {
          # guess keytype
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
        # load gene_table from org.*.eg.db
        if(inherits(keytype, "character")) {
          if(is_valid_keytype(keytype, organism = genome)) {
            gdf <- convert_id(g, from_keytype = keytype,
                              to_keytype = c("ENTREZID", "SYMBOL"),
                              organism   = genome,
                              rm_na      = FALSE)
          }
        }
      }
    }
  }
  #----------------------------------------------------------------------------#
  #-- run: convert
  #-- force: gdf to character (ENTREZID, as.number)
  out <- x
  if(inherits(gdf, "data.frame")) {
    gdf <- dplyr::mutate_if(gdf, is.numeric, as.character)
    if(keytype %in% names(gdf)) {
      out <- dplyr::left_join(x, gdf, by = setNames(keytype, nm = gid))
    } else {
      message(glue::glue(
        "'set_readable()' skipped, 'genome' or 'gene_table' not valid; \n",
        "'genome' is {genome}, expect 'character', \n",
        "'gene_table' is {gene_table}, expect '.csv', \n",
        "either 'genome' or 'genome_table' should be valid"
      ))
    }
  } else {
    warning(glue::glue(
      "'set_readable()' skipped, ",
      "either 'genome' or 'genome_table' should be valid; \n",
      "genome = '{genome}', genome_table = '{gene_table}'"
    ))
  }
  #----------------------------------------------------------------------------#
  # return
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
#' 1. remove not supported characters
#' - support: [\\w.] letters, numbers, _, .
#'
#' 2. remove common:
#' - lcPrefix(), longest common prefix
#' - lcSuffix(), longest common suffix
#'
#' 3. fix prefix
#' - prefix start with "letters"
#'
#' 4. ignore, if all x are the same
#'
#' @param x character
#'
#' @return character
#'
#' @export
deseq_sanitize_str <- function(x, n_max = 0, fix_prefix = FALSE) {
  x <- as.character(x)
  if(inherits(x, "character")) {
    # 1. supported characters
    out <- gsub("[^\\w\\.]", ".", x, perl = TRUE)
    if(length(unique(out)) == 1) {
      return(out)
    }
    # 2. prefix, suffix
    if(nchar(x[1]) > n_max & n_max > 0) {
      # longest prefix, suffix
      lcp <- Biobase::lcPrefix(out, ignore.case = FALSE)
      lcs <- Biobase::lcSuffix(out, ignore.case = FALSE)
      if(nchar(lcp) > 0) {
        out <- gsub(lcp, "", out)
      }
      if(nchar(lcs) > 0) {
        out <- gsub(lcs, "", out)
      }
    }
    # 3. fix prefix
    if(isTRUE(fix_prefix)) {
      if(any(grepl("^[^A-Za-z]", out, perl = TRUE))) {
        out <- paste0("X", out)
      }
    }
    out
  } else {
    x
  }
}


