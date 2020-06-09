# functions for DESeq2 downstream analysis
#
#
##-----------------------------##


#' read featureCounts output
#'
#' @param x str or list of feqtureCounts files
#'
#' @export
read_fc2 <- function(x){
  # gene
  if(is(x, "character")) {
    b <- lapply(x, function(i) {
      readr::read_delim(i, "\t", comment = "#", col_types = readr::cols())
    })
    names(b) <- basename(dirname(dirname(dirname(x))))
    b
  } else if(is(x, "list")) {
    lapply(x, function(k) {
      read_fc2(k)
    })
  } else {
    warning("unknown input")
  }
}


##------------------------------##
## single deseq dirs
## merge deseq_dir features
## ...
#' read files in DESeq2 output
#'
#' @param deseq_dir path to a.vs.b, contains deseq/enrich/report/... directories
#'
#' @export
merge_deseq <- function(x) {
  # load data
  pd <- read_rnaseq(x)

  ## write to files
  f_suffix  <- ".count.sens.txt"
  count_dir <- file.path(x, "merge", "count")
  dir.create(count_dir, recursive = TRUE, showWarnings = FALSE)

  ##-------------------------------------------##
  # control: ctl
  f1 <- sapply(names(pd$args), function(i) {
    pd$args[[i]]$count_ctl
  }, simplify = FALSE, USE.NAMES = TRUE)

  # read data: features
  p1 <- read_fc2(f1)

  # features, samples
  smps_ctl  <- names(p1[[1]])
  features  <- names(p1)
  count_ctl <- sapply(smps_ctl, function(i) {
    # prepare data
    tmp1 <- sapply(features, function(k) {
      df1 <- p1[[k]][[i]]
      names(df1)[7] <- i # rename bam
      df1 %>%
        dplyr::mutate_if(is.numeric, as.character)
    }, simplify = FALSE, USE.NAMES = TRUE)

    # data.frame
    df2 <- dplyr::bind_rows(tmp1)

    # save to file
    ff <- file.path(count_dir, paste0(i, f_suffix))
    readr::write_delim(df2, ff, "\t", col_names = TRUE)

    # output
    ff
  }, simplify = TRUE, USE.NAMES = TRUE)

  ##-------------------------------------------##
  # control: exp
  f1 <- sapply(names(pd$args), function(i) {
    pd$args[[i]]$count_exp
  }, simplify = FALSE, USE.NAMES = TRUE)

  # read data: features
  p1 <- read_fc2(f1)

  # features, samples
  smp_list  <- names(p1[[1]])
  features  <- names(p1)
  count_exp <- sapply(smp_list, function(i) {
    # prepare data
    tmp1 <- sapply(features, function(k) {
      df1 <- p1[[k]][[i]]
      names(df1)[7] <- i # rename bam
      df1 %>%
        dplyr::mutate_if(is.numeric, as.character)
    }, simplify = FALSE, USE.NAMES = TRUE)

    # data.frame
    df2 <- dplyr::bind_rows(tmp1)

    # save to file
    ff <- file.path(count_dir, paste0(i, f_suffix))
    readr::write_delim(df2, ff, "\t", col_names = TRUE)

    # output
    ff
  }, simplify = TRUE, USE.NAMES = TRUE)

  ##-------------------------------------------##
  # ctl vs exp
  pp <- read_deseq(x, "gene")
  deseq_dir <- file.path(x, "merge", "deseq")
  deseq_hub2(count_ctl, count_exp, deseq_dir, pp$args$genome)

  list(count_ctl = count_ctl,
       count_exp = count_exp,
       deseq_dir = x,
       feature = "merge",
       args = pd)
}


#' read atac directory
#'
#' input dir
#' output files, config
#'
#' @param x string path
#'
#' @export
read_rnaseq <- function(x) {
  x_type <- is_hiseq_dir(x)
  if(isTRUE(is_rnaseq_dir(x))) {
    message(paste0("Reading: ", x_type))
  } else {
    warning(paste0("Not a RNAseq dir: ", x))
    return(NULL)
  }

  # read args
  pk1 <- list_arguments_file(x)
  pk1 <- unlist(pk1) # all list
  if(! is.null(pk1)) {
    args <- sapply(pk1, load_pickle, USE.NAMES = TRUE)
  } else {
    args <- list()
  }

  list(atacseq_type = x_type,
       file_list = file_to_list(x, recursive = TRUE),
       args = args)
}


#' @param x string, path to dir
#'
#' @export
is_rnaseq_dir <- function(x) {
  x_type <- is_hiseq_dir(x)

  if(! is.null(x_type)) {
    grepl("^rnaseq|^deseq", x_type)
  }
}


#' @param x string, path to dir
#'
#' @export
is_rnaseq_single_dir <- function(x) {
  x_type <- is_hiseq_dir(x)
  x_type == "atacseq_single" # sample:1, replicate:n
}


#' @param x string, path to dir
#'
#' @export
is_rnaseq_multiple_dir <- function(x) {
  x_type <- is_hiseq_dir(x)
  x_type == "rnaseq_multiple" # sample:1, replicate:n
}


#' @param x string, path to dir
#'
#' @export
is_deseq_single_dir <- function(x) {
  x_type <- is_hiseq_dir(x)
  x_type == "deseq_single" # sample:1, replicate:n
}


#' get fix xls
#'
#' @param path string, path to the directory of deseq
#' @param feature string, gene|te|...
#'
#' @export
get_deseq_xls <- function(path, feature = "gene") {
  # check dir
  chk0 <- is_hiseq_dir(path) == "deseq_single"
  if(! isTRUE(chk0)) {
    warning("not a deseq_single dir")
    return(NULL)
  }

  deseq_xls <- file.path(path, feature, "deseq",
                         "transcripts_deseq2.fix.xls")

  if(file.exists(deseq_xls)) {
    return(deseq_xls)
  } else {
    warning(paste("file not exists: ", deseq_xls))
    return(NULL)
  }
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
read_deseq_csv <- function(x) {
  # read DESeq2 csv file
  # df <- readr::read_csv(x, col_types = readr::cols())
  df <- read.csv(x) %>%
    dplyr::select(-1) %>%
    dplyr::mutate(Gene = as.character(Gene))
  # dplyr::rename(id = Gene) %>%
  # dplyr::mutate(id = as.character(id))

  # add mean columns
  # df2 <- DESeq2_add_mean(df)
  df2 <- deseq_csv_mean(df)

  return(df2)
}


#' cal mean of replicates
#' default, 2 replicates for each condition
#'
#' @param data A data fraem of DESeq2 output, first 5 columns:
#'   <id> <a-rep1> <a-rep2> <b-rep1> <b-rep2> ...
#'
#' @export
deseq_csv_mean <- function(data) {
  # only for DESeq2 output csv
  # insert to: col-2, col-3
  # conserved columns
  stopifnot(is.data.frame(data))
  col_required <- c(
    "Gene", "baseMean", "log2FoldChange", "lfcSE",
    "stat", "pvalue", "padj"
  )
  col_meta <- c("sig", "symbol", "entrezid")
  stopifnot(all(col_required %in% colnames(data)))

  # sample names
  smp_names <- colnames(data)[!colnames(data) %in% c(col_required, col_meta)]
  smp_groups <- unique(fq_name(smp_names, rm_rep = TRUE))
  # stopifnot(length(smp_groups) == 2)

  # calculate mean
  ctl_mean <- data %>%
    dplyr::select(contains(smp_groups[1])) %>%
    rowMeans()
  exp_mean <- data %>%
    dplyr::select(contains(smp_groups[2])) %>%
    rowMeans()

  # assemble
  df1 <- data.frame(
    id = data$Gene,
    ctl = ctl_mean,
    exp = exp_mean,
    stringsAsFactors = FALSE
  )

  # change names
  colnames(df1) <- c("Gene", smp_groups)

  # output
  df2 <- merge(df1, data, by = "Gene")

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
count_to_matrix <- function(count_ctl, count_exp) {
  if (inherits(count_ctl, "data.frame") & inherits(count_exp, "data.frame")) {
    stopifnot("id" %in% colnames(count_ctl))
    stopifnot("id" %in% colnames(count_exp))
    df <- merge(count_ctl, count_exp, by = "id")
    ## names
    smp_ctl <- colnames(count_ctl[, -id])
    smp_exp <- colnames(count_exp[, -id])
  } else if (inherits(count_ctl, "character") & inherits(count_exp, "character")) {
    df1 <- read_fc(count_ctl)
    df2 <- read_fc(count_exp)
    df <- merge(df1, df2, by = "id")
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



#' list_arguments_file
#'
#' @param x path to the directory.
#' expect x/config/{feature}/config/arguments.pickle,
#' x/{feature}/config/arguments.pickle
#'
#' @export
list_arguments_file <- function(x) {
  # arguments file
  # RNAseq multiple: x/config/{feature}/config/arguments.txt
  # RNAseq single: x/{feature}/config/arguments.txt
  # RNAseq deseq: x/{feature}/config/arguments.txt
  chk0 <- dir.exists(x)

  # version-1: x/{feature}/...
  dir1 <- list.dirs(x, full.names = TRUE, recursive = FALSE)
  ps1 <- sapply(c(x, dir1), function(i) {
    file.path(i, "config", "arguments.pickle")}, simplify = TRUE)
  # ps1 <- mapply(function(i) file.path(i, "config", "arguments.pickle"), dir1)
  if (length(ps1) > 0) {
    ps1 <- ps1[file.exists(ps1)]
  }
  chk1 <- length(ps1) > 0

  # version-2: x/config/{feature}/...
  dir2 <- list.dirs(file.path(x, "config"), full.names = TRUE, recursive = FALSE)
  ps2 <- sapply(dir2, function(i) {
    file.path(i, "config", "arguments.pickle")}, simplify = TRUE)
  # ps2 <- mapply(function(i) file.path(i, "config", "arguments.pickle"), dir2)
  if (length(ps2) > 0) {
    ps2 <- ps2[file.exists(ps2)]
  }
  chk2 <- length(ps2) > 0

  ## mission
  ## RNAseq multiple: ps1=0, ps2=1
  ## RNAseq single:   ps1=1, ps2=0
  if (isTRUE(chk1)) {
    plist <- as.list(ps1)
  } else if (isTRUE(chk2)) {
    plist <- as.list(ps2)
  } else {
    plist <- NULL
  }

  ## change name
  if (!is.null(plist)) {
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
load_pickle <- function(x) {
  # python version issue !!
  if (file.exists(x) & endsWith(x, ".pickle")) {
    pd <- reticulate::import("pandas")
    tag <- pd$read_pickle(x)
  } else {
    tag <- NULL
  }
}


##------------------------------##
## single deseq dirs
## including dataset
## including plots
## including *.rds
## ...
#' read files in DESeq2 output
#'
#' @param deseq_dir path to a.vs.b, contains deseq/enrich/report/... directories
#' @param feature string, gene|te|...
#'
#' @export
read_deseq <- function(deseq_dir, feature = "gene") {
  # check if dirA/dirB is "deseq_single"
  chk0 <- hiseqr::is_hiseq_dir(deseq_dir) == "deseq_single"
  if(! chk0) {
    stop(paste0("not a deseq_single dir, ", deseq_dir))
  }

  # args
  args <- deseq_single_dir(deseq_dir, feature)

  # deseq genes
  sig_gene <- get_sig_gene(deseq_dir, feature)

  # to df
  count_df <- sapply(sig_gene, nrow) %>%
    as.data.frame %>%
    dplyr::rename(count = ".") %>%
    tibble::rownames_to_column("type") %>%
    dplyr::filter(! type %in% "not")


  ##------------------------##
  ## colors
  cc <- scales::hue_pal()(3) # red, green, blue
  cc <- c(cc[1], cc[3], cc[2]) # red, blue, green

  # plot
  count_plot <- count_df %>%
    dplyr::mutate(type = factor(type, levels = c("up", "not", "down"))) %>%
    ggplot(aes(x = count, y = "1", fill = type)) +
    geom_col(color = "grey20") +
    geom_text(aes(label = count), position = position_stack(vjust = 0.5)) +
    scale_fill_manual(values = cc) +
    theme_bw() +
    theme(legend.position = "right",
          legend.title    = element_blank(),
          panel.grid      = element_blank(),
          plot.title      = element_text(size = 10))

  ##-----------------------##
  ## enrich analysis
  # go_data.rds
  # go_plots.rds
  # kegg_data.rds
  # kegg_plots.rds
  # file_path: ...
  enrich_dir <- args$enrichdir
  ff <- file_to_list(enrich_dir, recursive = TRUE)

  ## output
  list(args       = args,
       count_df   = count_df,
       count_plot = count_plot,
       file_list  = ff)
}


#' read enrich dir
#' @param x path to enrich dir, eg: up/down/..., contains files go_data.rds, ...
#'
#' @export
read_enrich <- function(x) {
  # dir of up/down/not/up_and_down ...
  if(! is_enrich_dir(x)) {
    warning("not a enrich directory")
    return(NULL)
  }

  # all files
  ff <- is_enrich_dir(x, return_file = TRUE)

  # show list
  tmp <- sapply(sort(names(ff)), function(i){
    d <- file.path(x, i)
    message(paste("Found enrich:", i, "...",
                  ifelse(dir.exists(d), "directory", "file")))
  })

  ff
}


#' fetch files/directories from directory
#' @param x path to enrich dir, eg: up/down/..., contains files go_data.rds, ...
#'
#' @export
is_enrich_dir <- function(x, return_file = FALSE) {
  # input
  if(! dir.exists(x)) {
    warning("Directory expected")
    return(FALSE)
  }

  # files/folders
  req_files <- c("go_data.rds",
                 "go_plots.rds",
                 # "kegg_data.rds",
                 # "kegg_plots.rds",
                 "go_group",
                 "go_enrich",
                 "go_gsea") # ,
  # "kegg_enrich",
  # "kegg_gsea")

  # level-1: files in x
  f1 <- list.files(x)
  if(! all(req_files %in% f1)) {
    msg <- paste(paste(f1, file.exists(file.path(x, req_files)), sep = "..."),
                 collapse = "\n")
    # warning("not a enrich dir, go_data.rds, kegg_data.rds, ..., expected")
    warning(msg)
    return(FALSE)
  }

  # loop over files
  ff <- file_to_list(x, recursive = TRUE)

  if(isTRUE(return_file)) {
    return(ff)
  } else {
    return(TRUE)
  }
}


#' hiseq_type
#'
#' @param x path, the directory of RNAseq output
#'
#' @export
is_hiseq_dir <- function(x) {
  chk0 <- dir.exists(x)

  # check arguments.pickle file
  pk_files <- list_arguments_file(x)
  # choose one

  if (is.null(pk_files)) { # | is.na(pk_files)) {
    return(NULL)
  } else {
    # pk_files <- pk_files[1] # 1st item
    pk <- load_pickle(pk_files[[1]])

    # get values
    if ("rnaseq_type" %in% names(pk)) {
      tag <- pk$rnaseq_type
    } else if ("atacseq_type" %in% names(pk)) {
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
rnaseq_single_dir <- function(x, feature = "gene") {
  # check
  chk0 <- is_hiseq_dir(x) == "rnaseq_single"

  if (!isTRUE(chk0)) {
    return(NULL)
  }

  # pick pickle file
  pk_files <- list_arguments_file(x)
  args <- load_pickle(pk_files[[feature]])

  ## check required args
  required_names <- c(
    "align_stat", "smp_name", "genome", "gtf", "gtf",
    "count_sens", "count_anti", "strandness_status"
  )

  for (name in required_names) {
    # args[[name]] = ifelse(name %in% names(args), c(args[[name]]), NULL)
    if (!name %in% names(args)) {
      args[[name]] <- NULL
    }
  }

  ## mapping data
  df <- read_align1(args$align_stat)
  args$reads_total <- sum(df$count)

  ## mito.u, mito.m, genome.u genome.m
  for (i in seq_len(nrow(df))) {
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
rnaseq_multiple_dir <- function(x, feature = "gene") {
  # check
  chk0 <- is_hiseq_dir(x) == "rnaseq_multiple"

  if (!isTRUE(chk0)) {
    return(NULL)
  }

  # pick pickle file
  pk_files <- list_arguments_file(x)
  args <- load_pickle(pk_files[[feature]])

  ## check required args
  required_names <- c("smp_name", "outdir", "genome", "gtf")

  for (name in required_names) {
    # args[[name]] <- ifelse(name %in% names(args), c(args[[name]]), NULL)
    if (!name %in% required_names) {
      args[[name]] <- NULL
    }
  }

  ## parsing RNAseq single, saved as vector
  args$single_dirs <- mapply(function(i) file.path(args$outdir, i), args$smp_name)
  args$single_args <- lapply(args$single_dirs, rnaseq_single_dir)

  return(args)
}



#' deseq_single_dir
#'
#' reading data from RNAseq Single
#'
#' @param x path to the directory of RNAseq single
#'
#' @export
deseq_single_dir <- function(x, feature = "gene") {
  # check
  chk0 <- is_hiseq_dir(x) == "deseq_single"

  if (!isTRUE(chk0)) {
    return(NULL)
  }

  # pick pickle file
  pk_files <- list_arguments_file(x)

  if(! feature %in% names(pk_files)) {
    warning(paste0("feature not exists: ", feature))
    return(NULL)
  }

  args <- load_pickle(pk_files[[feature]])

  ## check required args
  required_names <- c("deseqdir", "dirs_ctl", "dirs_exp",
                      "count_ctl", "count_exp")

  for (name in required_names) {
    # args[[name]] = ifelse(name %in% names(args), c(args[[name]]), NULL)
    if (!name %in% names(args)) {
      args[[name]] <- NULL
    }
  }

  return(args)
}


#' get fix xls
#'
#' @param path string, path to the directory of deseq
#' @param feature string, gene|te|...
#' @param gene_name_only boolean, whether get gene name only
#'
#' @export
#'
get_sig_gene <- function(path, feature = "gene",
                         gene_name_only = FALSE) {
  # xls file
  deseq_xls <- get_deseq_xls(path, feature)

  # convert to df
  if(is.null(deseq_xls)) {
    return(NULL)
  }

  # read
  df <- readr::read_delim(deseq_xls, "\t", col_types = readr::cols())

  # prepare list
  p <- lapply(c("up", "down", "not"), function(i) {
    get_sig_gene2(df, i, gene_name_only)
  })

  names(p) <- c("up", "down", "not")
  p
}



#' rnaseq_report
#'
#' @param df data.frame
#' @param type string up,down,not
#' @param gene_name_only boolean, whether get gene name only
#'
#' @export
get_sig_gene2 <- function(df, type = "up", gene_name_only = FALSE) {
  if(inherits(df, "data.frame")) {
    # add sig
    if(! "sig" %in% colnames(df)) {
      df <- deseq_sig(df)
    }

    # check
    if("sig" %in% colnames(df)) {
      df_sub <- subset(df, sig == type)
      # name
      if(isTRUE(gene_name_only)) {
        if("id" %in% colnames(df)) {
          df_sub <- df_sub$id
        } else if("Gene" %in% colnames(df)) {
          df_sub <- df_sub$Gene
        } else {
          warning("Gene|id not found")
        }
      }
    } else {
      warning("column:sig not found")
      df_sub <- NULL
    }
  } else {
    warning("expect data.frame")
    df_sub <- NULL
  }

  df_sub
}


#' add sig label to data.frame
#'
#' @param data data.frame
#' @param fc_cutoff numeric, cutoff for foldchange
#' @param pval_cutoff numeric, cutoff for pvalue
#' @param type string, all, both, up, down, not, default: all
#'
#' @import dplyr
#'
#' @export
#'
deseq_sig <- function(data, fc_cutoff = 2, pval_cutoff = 0.05, type = "all") {
  stopifnot(inherits(data, "data.frame"))

  # required: genes
  if ("id" %in% colnames(data)) {
    data$gene_id <- data$id
  } else if ("Gene" %in% colnames(data)) {
    data$gene_id <- data$Gene
  } else {
    warning("Gene_name column not found, id|Gene expected")
    return(data)
  }

  # required columns
  stopifnot(all(c("log2FoldChange", "pvalue") %in% colnames(data)))
  if ("padj" %in% colnames(data)) {
    data$pval_check <- data$padj
  } else {
    data$pval_check <- data$pvalue
  }

  type <- tolower(type) # lower case
  log2_fc <- log2(as.numeric(fc_cutoff))

  # split group
  df1 <- dplyr::filter(data, log2FoldChange >= log2_fc & pval_check <= pval_cutoff)
  df2 <- dplyr::filter(data, log2FoldChange <= -log2_fc & pval_check <= pval_cutoff)
  df3 <- dplyr::filter(data, !gene_id %in% c(df1$gene_id, df2$gene_id))

  if (nrow(df1) > 0) {
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
  } else if (type == "not") {
    df <- df3
    # df <- df %>% dplyr::filter(sig %in% c("not"))
  } else if (type == "both") {
    df <- dplyr::bind_cols(df1, df2)
    # df <- df %>% dplyr::filter(sig %in% c("up", "down"))
  } else {
    df <- dplyr::bind_rows(df1, df2, df3)
  }

  # remove pval_check column
  df <- dplyr::select(df, -c(gene_id, pval_check))

  # return
  return(df)
}




#' get fix xls
#'
#' @param path string, path to the directory of deseq
#' @param feature string, gene|te|...
#' @param gene_name_only boolean, whether get gene name only
#'
#' @export
#'
get_sig_gene_count <- function(...) {
  sig_gene <- get_sig_gene(...)
  sapply(sig_gene, nrow) %>%
    as.data.frame() %>%
    dplyr::rename(count = ".") %>%
    t()
}





