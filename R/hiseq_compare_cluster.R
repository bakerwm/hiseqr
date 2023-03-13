#' functions for compareCluster



#' hiseq_prep_compare_go
#' prepare data for `compareCluster()` function
#'
#' @param x character path to the rnaseq_rx, parsing the mut_name
#' @param outdir string path to the directory, saving results
#' @param fc numeric fold-change cutoff, default: 2
#' @param pvalue float cutoff for pvalue, default: 0.05
#' @param p_adjust logical, whether use adjusted pvalue, default: TRUE
#' @param shrink_method string how to shrink the fold-change results
#'
#' @return data.frame including columns: gene_id, group, sig, genome
#'
#' arguments
#' @export
hiseq_prep_compare_go <- function(x, ...) {
  #-- run
  dots = rlang::list2(...)
  args <- list(
    outdir   = "./",
    fc       = 2,
    pvalue   = 0.05,
    p_adjust = TRUE,
    shrink_method = "apeglm"
  )
  args <- purrr::list_modify(args, !!!dots)
  #-- make sure, x is rnaseq_rx
  if(all(is_hiseq_dir(x, "rnaseq_rx"))) {
    message("preparing compare_go() data ...")
    #-- update names
    smp <- names(x)
    if(is.null(smp)) {
      smp <- sapply(x, function(i) list_hiseq_file(i, "mut_name", "rx"))
    }
    # simplify
    smp <- deseq_sanitize_str(basename(smp), 10)
    smp <- gsub("\\.", "_", smp) # fix, for formula
    names(x) <- smp
  } else {
    message("compare_go() failed, illegal x")
    return(NULL)
  }
  #-- update genome
  g_list <- sapply(x, function(i) list_hiseq_file(i, "genome", "rx"))
  genome <- g_list[[1]]
  if(inherits(genome, "character")) {
    if(!(length(unique(g_list)) == 1 & is_valid_organism(genome))) {
      g_str <- paste(unique(g_list), collapse = ",")
      message(glue::glue("x, genome not consistent: {g_str}"))
      return(NULL)
    }
  } else {
    message(glue::glue("genome is {class(genome)}, expect character"))
    return(NULL)
  }
  #-- load data
  check_path(args$outdir)
  #-- gene_id, sig, group, genome
  args$genome <- genome
  df_go <- lapply(names(x), function(i) {
    deseq_dir <- list_hiseq_file(x[i], "deseq_dir", "rx")
    deseq_qc_res(deseq_dir, !!!args) %>%
      dplyr::filter(!sig == "not") %>%
      dplyr::select(gene_id, sig) %>%
      dplyr::mutate(group = i, genome = genome)
  }) %>%
    dplyr::bind_rows()
  #-- prep kegg data
  organism <- get_organism_name(genome)
  #-- keytypes for kegg
  kegg_keytype <- ifelse(
    organism == organism,
    "ncbi-geneid", "kegg"
  )
  #-- keytypes for kegg
  keytype = guess_keytype(as.character(df_go$gene_id), organism = genome)
  df_go %>%
    dplyr::mutate(
      kegg_gene = to_kegg_gene_id(
        as.character(df_go$gene_id), organism, keytype,
        simplified = TRUE, rm_na = FALSE, return_table = FALSE
      )
    )
}


#' hiseq_prep_compare_go
#' prepare data for `compareCluster()` function
#'
#' @param x character path to the rnaseq_rx, parsing the mut_name
#' @param outdir string path to the directory, saving results
#' @param fc numeric fold-change cutoff, default: 2
#' @param pvalue float cutoff for pvalue, default: 0.05
#' @param p_adjust logical, whether use adjusted pvalue, default: TRUE
#' @param shrink_method string how to shrink the fold-change results
#'
#' @return data.frame including columns: gene_id, group, sig, genome
#'
#' arguments
#' @export
hiseq_prep_pairwise_go <- function(x, ...) {
  # only the first two of x choosen for down-stream analysis
  if(inherits(x, "character")) {
    if(length(x) >= 2) {
      x <- x[1:2]
    } else {
      message(glue::glue("x length is {length(x)}, expect >= 2"))
      return(NULL)
    }
  } else {
    message(glue::glue("x is {class(x)}, expect character"))
    return(NULL)
  }
  # loading data
  dots <- rlang::list2(...)
  df <- hiseq_prep_compare_go(x, !!!dots)
  # fix, overlap between two samples
  if(inherits(df, "data.frame")) {
    df %>%
      dplyr::group_by(gene_id, kegg_gene, sig) %>%
      dplyr::summarise(smp = paste(group, collapse = "_and_"), .groups = "drop") %>%
      dplyr::mutate(group = ifelse(grepl("_and_", smp), "both", smp),
                    genome = df$genome[1])
  }
}


#' arguments
hiseq_compare_go_enrich <- function(x, ...) {
  #-- default args
  dots <- rlang::list2(...)
  args <- purrr::list_modify(list(outdir = "./"), !!!dots)
  df <- hiseq_prep_compare_go(x, !!!args)
  #-- run
  tmp <- compare_go_enrich(df, !!!args)
}


#' arguments
hiseq_compare_kegg_enrich <- function(x, ...) {
  #-- default args
  dots <- rlang::list2(...)
  args <- purrr::list_modify(list(outdir = "./"), !!!dots)
  #-- loading data
  df <- hiseq_prep_compare_go(x, !!!dots)
  #-- run
  tmp <- compare_kegg_enrich(df, !!!args)
}


#' arguments
hiseq_pairwise_go_enrich <- function(x, ...) {
  #-- default args
  dots <- rlang::list2(...)
  args <- purrr::list_modify(list(outdir = "./"), !!!dots)
  #-- loading data
  df <- hiseq_prep_pairwise_go(x, !!!dots)
  #-- run
  tmp <- compare_go_enrich(df, !!!args)
}


#' arguments
hiseq_pairwise_kegg_enrich <- function(x, ...) {
  #-- default args
  dots <- rlang::list2(...)
  args <- purrr::list_modify(list(outdir = "./"), !!!dots)
  #-- loading data
  df <- hiseq_prep_pairwise_go(x, !!!dots)
  #-- run
  tmp <- compare_kegg_enrich(df, !!!args)
}



#' @param x data.frame list of genes, required columns:
#' gene_id, sig, group, genome
#' @param outdir string path to the directory, saving results
#' @param suffix string name of the group samples, default: NULL
#'
#' @return
compare_go_enrich <- function(x, ...) {
  #-- default args
  dots <- rlang::list2(...)
  args <- list(
    outdir = "./",
    suffix = NULL
  )
  args <- purrr::list_modify(args, !!!dots)
  # #-- update; outdir
  check_path(args$outdir)
  args$outdir <- normalizePath(args$outdir)
  #-- loading data
  #-- check data
  if(!inherits(x, "data.frame")) {
    message(glue::glue("compare_go() failed, x is {class(x)}, expect data.frame"))
    return(NULL)
  }
  c_req <- c("gene_id", "sig", "group", "genome")
  if(!all(c_req %in% names(x))) {
    c_str <- paste(c_req, collapse = ",")
    message(glue::glue("compare_go() failed, missing clolumns: {c_str}"))
    return(NULL)
  }
  #-- output files
  genome   <- x$genome[1]
  if(!inherits(args$suffix, "character")) {
    suffix <- paste(unique(x$group), collapse = "_")
  } else {
    suffix <- args$suffix
  }
  out_name <- glue::glue("compare_GO.{suffix}")
  go_data_rds   <- file.path(args$outdir, paste0(out_name, ".data.rds"))
  #-- guess keytype
  if(!file.exists(go_data_rds)) {
    keytype <- guess_keytype(as.character(x$gene_id), organism = genome)
    onts <- c("ALL", "BP", "CC", "MF")
    tmp <- lapply(onts, function(ont) {
      message(glue::glue("compare_go() for {ont}"))
      args_local <- list(
        geneClusters = gene_id~group+sig,
        fun          = "enrichGO",
        data         = x,
        OrgDb        = get_orgdb(genome),
        keyType      = keytype,
        ont          = ont
      )
      cc <- tryCatch(
        {
          rlang::exec(clusterProfiler::compareCluster, !!!args_local)
        },
        error = function(cond) {
          message(glue::glue("compare_go() failed for {ont}"))
          return(NULL)
        }
      )
    })
    names(tmp) <- onts
    saveRDS(tmp, file = go_data_rds) # saved to file
  }
  #-- load data
  if(file.exists(go_data_rds)) {
    go_data <- readRDS(go_data_rds)
  } else {
    go_data <- NULL
  }
  #-- save plots
  message(glue::glue("go_data is {class(go_data)}"))
  if(inherits(go_data, "list")) {
    message("plot-1")
    # tmp3 <- lapply(names(go_data), function(ont) {
    for(ont in names(go_data)) {
      message(glue::glue("plot-2 {ont}"))
      cc <- go_data[[ont]]
      if(is.null(cc)) {
        next
      }
      out_png_name <- paste0(out_name, ".", ont, ".dotplot.png")
      out_png <- file.path(args$outdir, out_png_name)
      title <- glue::glue("GO: {ont}")
      p <- dotplot(go_data[[ont]], x="group", showCategory = 10, label_format = 40) +
        facet_grid(~sig) +
        ggtitle(title) +
        theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))
      # output
      ggsave(out_png, p, width = 10, height = 16)
    }
  }
}



#' @param x data.frame list of genes, required columns:
#' gene_id, sig, group, genome
#' @param outdir string path to the directory, saving results
#' @param suffix string name of the group samples, default: NULL
#'
#'
#' @return
compare_kegg_enrich <- function(x, ...) {
  #-- default args
  dots <- rlang::list2(...)
  args <- list(
    outdir = "./",
    suffix = NULL
  )
  args <- purrr::list_modify(args, !!!dots)
  # #-- update; outdir
  check_path(args$outdir)
  args$outdir <- normalizePath(args$outdir)
  #-- loading data
  #-- check data
  if(!inherits(x, "data.frame")) {
    message(glue::glue("compare_go() failed, x is {class(x)}, expect data.frame"))
    return(NULL)
  }
  c_req <- c("gene_id", "sig", "group", "genome", "kegg_gene")
  if(!all(c_req %in% names(x))) {
    c_str <- paste(c_req, collapse = ",")
    message(glue::glue("compare_go() failed, missing clolumns: {c_str}"))
    return(NULL)
  }
  #-- kegg args
  genome   <- x$genome[1]
  organism <- get_organism_name(genome)
  kegg_keytype <- get_kegg_keytype(organism) # for OrgDb
  kegg_keyType <- get_kegg_keyType(organism) # for clusterProfiler
  kegg_code    <- get_kegg_code(organism) # organism
  #-- output
  if(!inherits(args$suffix, "character")) {
    suffix <- paste(unique(x$group), collapse = "_")
  } else {
    suffix <- args$suffix
  }
  # suffix   <- paste(unique(x$group), collapse = "_")
  out_name <- glue::glue("compare_KEGG.{suffix}")
  kegg_data_rds   <- file.path(args$outdir, paste0(out_name, ".data.rds"))
  #-- run
  if(!file.exists(kegg_data_rds)) {
    args_local <- list(
      geneClusters = kegg_gene~group+sig,
      fun  = "enrichKEGG",
      data = x,
      keyType       = kegg_keyType,
      organism      = kegg_code,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05
    )
    kegg_data <- tryCatch(
      {
        rlang::exec(clusterProfiler::compareCluster, !!!args_local)
      },
      error = function(cond) {
        message("hiseq_compare_kegg() failed")
        return(NULL)
      }
    )
    saveRDS(kegg_data, file = kegg_data_rds)
  }
  #-- load data
  if(file.exists(kegg_data_rds)) {
    kegg_data <- readRDS(kegg_data_rds)
  } else {
    kegg_data <- NULL
  }
  #-- save plots
  if(inherits(kegg_data, "compareClusterResult")) {
    out_png <- file.path(args$outdir, paste0(out_name, ".dotplot.png"))
    title <- glue::glue("KEGG")
    p <- dotplot(kegg_data, x="group", showCategory = 10) +
      facet_grid(~sig) +
      ggtitle(title) +
      theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))
    # output
    ggsave(out_png, p, width = 9, height = 12)
  } else {
    message(glue::glue("kegg_data is {class(kegg_data)}, expect compareClusterResult"))
  }
}


