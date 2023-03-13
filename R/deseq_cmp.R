#' Functions for RNAseq compare analysis
#'
#' Apply DESeq2 for RNAseq DE analysis
#' require >=2 replicates for each sample
#'
#' input: count.txt, from featureCounts output
#' output: DESeq2 output
#'
#' helper functions:
#'


#' deseq_cmp
#' 
#' @description
#'    A tmp function, help documentation 
deseq_cmp <- function(x) {
  x
}


#' deseq_pair_report
#'
#' @param x string to GroupA "A.vs.B"
#' @param y string to GroupB "A.vs.B"
#' @param feature string default: gene
#' @param outdir string directory to save files
#'
#' @export
deseq_pair_report <- function(x, y, outdir, shrink_method = "apeglm") {
  message("# Run RNAseq pair ...")
  #----------------------------------------------------------------------------#
  subname <- paste0(basename(x), ".compare.", basename(y))
  subdir <- file.path(outdir, subname)
  if(! dir.exists(subdir)) {
    message("Create dir .......")
    dir.create(subdir, mode = "0755", recursive = TRUE, showWarnings = FALSE)
  }
  subdir <- normalizePath(subdir)
  template <- system.file("rnaseq", "deseq_pair_report.Rmd",
                          package = "hiseqr")
  template_to <- file.path(subdir, basename(template))
  # copy Rmd
  file.copy(template, template_to)
  out_html <- file.path(subdir, "HiSeq_report.html")
  ## run Rmd
  if(file.exists(out_html)) {
    message("output exists, skipping ...")
  } else {
    message("render html")
    rmarkdown::render(input       = template_to,
                      output_file = out_html,
                      params      = list(x = x,
                                         y = y,
                                         shrink_method = shrink_method))
  }
}



#' RNAseq cmp, stat
#'
#' @param x string to GroupA "A.vs.B"
#' @param y string to GroupB "A.vs.B"
#' @param ... extra parameters for deseq_qc_res()
#'
#'
#' @export
read_deseq_pair <- function(x, y, ...) {
  if(all(is_hiseq_dir(c(x, y), "rnaseq_rx"))) {
    mut_name <- list_hiseq_file(c(x, y), "mut_name",  "rx") %>% unlist
    mut_name <- deseq_sanitize_str(mut_name, 10)
    title    <- paste(mut_name, collapse = ".vs.")
    message(glue::glue("Compare RNAseq: {title}"))
  } else {
    flag  <- is_hiseq_dir(c(x, y), "rnaseq_rx")
    title <- paste(c("x", "y"), as.character(flag), sep = ": ", collapse = ", ")
    message(glue::glue("read_deseq_pair failed: {title}"))
    return(NULL)
  }
  ## check reference genome
  g <- list_hiseq_file(c(x, y), "genome") %>% unlist()
  if(g[1] != g[2]) {
    flag <- paste(as.character(g), collapse = ", ")
    warning("read_deseq_pair failed: ref genome not identical: {flag}")
    return(NULL)
  }
  ## convert to absolute path
  x <- normalizePath(x)
  y <- normalizePath(y)
  deseq_dir <- list_hiseq_file(c(x, y), "deseq_dir") %>% unlist()
  x_smpname <- list_hiseq_file(x, "smp_name")
  y_smpname <- list_hiseq_file(y, "smp_name")
  ## in case, mut_name identical
  if(mut_name[1] == mut_name[2]) {
    mut_name[1] <- paste0(mut_name[1], ".x")
    mut_name[2] <- paste0(mut_name[2], ".y")
  }
  ## in case, x-name y-name identical
  if(x_smpname == y_smpname) {
    x_smpname <- paste0(x_smpname, ".x")
    y_smpname <- paste0(x_smpname, ".y")
  }
  ## degenes
  # output
  list(
    x = x,
    y = y,
    title = title,
    x_res = deseq_qc_res(deseq_dir[1]) %>% mutate(gene_id = as.character(gene_id)),
    y_res = deseq_qc_res(deseq_dir[2]) %>% mutate(gene_id = as.character(gene_id)),
    x_name = mut_name[1],
    y_name = mut_name[2],
    x_smpname = x_smpname,
    y_smpname = y_smpname,
    genome    = g[1]
  )
}



#' stat_deseq_pair
#'
#' @param x string to GroupA "A.vs.B"
#' @param y string to GroupB "A.vs.B"
#' @param .pd_rds string path to the file, saving read_deseq_pair() output
#' @param kepp_not_sig logical whether keep not sig changed genes, default: TRUE
#' @param group character "all", "gene", "te", "piRC"
#'
#'
#' @export
stat_deseq_pair <- function(x = NULL, y = NULL, .pd_rds = NULL,
                            keep_not_sig = TRUE, group = "all", ...) {
  if(inherits(.pd_rds, "character")) {
    pd = readRDS(.pd_rds)
  } else {
    pd <- read_deseq_pair(x, y, ...)
  }
  if(! inherits(pd, "list")) {
    warning("stat_deseq_pair failed")
    return(NULL)
  }
  #----------------------------------------------------------------------------#
  ## gene_prefix
  ## only for drosophila: dm6
  g_list <- pd$x_res$gene_id
  if(inherits(pd[["genome"]], "character")) {
    if(pd[["genome"]] == "dm6") {
      a    <- table(stringr::str_sub(pd$x_res$gene_id, 1, 3))
      gp   <- names(a[a == max(a)]) # gene prefix
      gene <- pd$x_res$gene_id[startsWith(pd$x_res$gene_id, gp)]
      a2   <- a[a<max(a)] # exlcue gene
      clup <- names(a2[a2 == max(a2)]) # cluster prefix
      clu  <- pd$x_res$gene_id[startsWith(pd$x_res$gene_id, clup)]
      clu  <- c(clu, "42AB", "flam") # for drosophila
      te   <- pd$x_res$gene_id[!pd$x_res$gene_id %in% c(gene, clu)]
      ## filt
      if(group == "gene") {
        g_list <- gene
      } else if(group == "te") {
        g_list <- te
      } else if(group == "cluster") {
        g_list <- clu
      } else {
        g_list <- pd$x_res$gene_id
      }
    }
  }
  pd$x_res <- dplyr::filter(pd$x_res, gene_id %in% g_list)
  pd$y_res <- dplyr::filter(pd$y_res, gene_id %in% g_list)
  #----------------------------------------------------------------------------#
  ## return sig table, sig list
  df1x <- as.data.frame(table(pd$x_res$sig))
  df1y <- as.data.frame(table(pd$y_res$sig))
  ## for 0-row conditions
  if(nrow(df1x) == 0) {
    df1x <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("Var1", "Freq"))
  }
  if(nrow(df1y) == 0) {
    df1y <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("Var1", "Freq"))
  }
  df1  <- merge(df1x, df1y, by = "Var1", all = TRUE)
  df1[is.na(df1)] <- 0 #
  names(df1) <- c("sig", pd$x_name, pd$y_name)
  ## sig list
  df2 <- dplyr::left_join(pd$x_res[, c("gene_id", "sig")],
                          pd$y_res[, c("gene_id", "sig")], by = "gene_id")
  names(df2) <- c("gene_id", pd$x_name, pd$y_name)
  df2 <- dplyr::mutate(df2, across(everything(), as.character))
  df2[is.na(df2)] <- "not"
  ## add trend
  df3 <- df2 %>%
    tidyr::unite("trend", -gene_id, sep = "-", remove = FALSE)
  ## filt by not-nog
  if(! isTRUE(keep_not_sig)) {
    df3 <- df3 %>%
      dplyr::filter(! trend %in% c("not-not"))
  }
  ## output
  list(
    pd = pd,
    df1 = df1,
    df2 = df2,
    df3 = df3
  )
}


#' plot_overlap_deseq_pair
#'
#' @param x string to GroupA "A.vs.B"
#' @param y string to GroupB "A.vs.B"
#' @param .pd_rds string path to the file, saving read_deseq_pair() output
#' @param ... parameters for deseq_qc_res()
#'
#'
#' @export
plot_overlap_deseq_pair <- function(x = NULL, y = NULL, .pd_rds = NULL, group = "all", ...) {
  # load data
  px <- stat_deseq_pair(x, y, .pd_rds, group = group, ...) # px
  if(! inherits(px, "list")) {
    warning("stat_deseq_pair failed")
    return(NULL)
  }
  # plots - overlap2
  p_list <- lapply(c("up", "down", "not"), function(s) {
    df1 <- px$df2 %>%
      dplyr::filter(!!sym(px$pd$x_name) == s | !!sym(px$pd$y_name) == s) %>%
      dplyr::mutate(across(-gene_id, ~ .x == s))
    # make plot
    ggplot(df1, aes(A = !!sym(px$pd$x_name), B = !!sym(px$pd$y_name))) +
      geom_venn(
        show_percentage = FALSE,
        fill_color    = scales::hue_pal()(2),
        fill_alpha    = .8,
        stroke_size   = .5,
        set_name_size = 4
      ) +
      ggtitle(s) +
      theme_void() +
      coord_fixed() +
      theme(plot.title = element_text(hjust = .5))
  })
  names(p_list) <- c("up", "down", "not")
  ## combine, up/down
  p1 <- patchwork::wrap_plots(p_list, nrow = 1)
  c(p_list, list(merge = p1))
}




#' .cluster_ma
#' @param method string see dist(method = ),
#'
#' @export
.cluster_ma <- function(ma, method = "euclidean") {
  if(is.matrix(ma) & is.numeric(ma)) {
    d <- dist(ma, method = "euclidean", upper = F)
    h <- hclust(d, "complete")
    # sort matrix
    a <- ma[h$order, ]
  } else {
    warning("matrix required, skipped...")
    a <- ma
  }
  ## output
  a
}


#' plot_heatmap_deseq_pair
#'
#' @param x string to GroupA "A.vs.B"
#' @param y string to GroupB "A.vs.B"
#' @param feature string default: gene
#' @param .pd_rds string path to the file, saving read_deseq_pair() output
#'
#'
#' @export
plot_heatmap_deseq_pair <- function(x = NULL, y = NULL, .pd_rds = NULL, group = "all", ...) {
  # load data
  px <- stat_deseq_pair(x, y, .pd_rds, group = group, ...) # px
  if(! inherits(px, "list")) {
    warning("stat_deseq_pair failed")
    return(NULL)
  }
  ##------------------------------------------------##
  ## colors
  cc <- scales::hue_pal()(3) # red, green, blue
  cc <- c(cc[1], cc[3], cc[2]) # red, blue, green
  ## data
  df1 <- px$df3 %>%
    dplyr::filter(!trend %in% c("not-not")) %>%
    dplyr::select(-trend) %>%
    tidyr::pivot_longer(-gene_id, names_to = "sample", values_to = "sig") %>%
    dplyr::mutate(sig = factor(sig, levels = c("up", "not", "down")),
                  gene_id = as.character(gene_id))
  ## cluster rows
  ma <- px$df3 %>%
    dplyr::filter(!trend %in% c("not-not")) %>%
    dplyr::select(-trend) %>%
    dplyr::mutate(across(-gene_id, ~ ifelse(.x == "up", 1,
                                            ifelse(.x == "down", -1, 0)))) %>%
    tibble::column_to_rownames("gene_id") %>%
    as.matrix() %>%
    .cluster_ma()
  ## heatmap
  if(nrow(df1) > 0) {
    p1 <- ggplot(df1, aes(sample, gene_id, fill = sig)) +
      geom_tile() +
      scale_fill_manual(values = cc) +
      scale_x_discrete(position = "top") +
      scale_y_discrete(limits = rownames(ma)) +
      ggtitle(px$pd$title) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.title = element_blank()
      )
    ## remove y-text
    if(nrow(df1) > 100) {
      # remove y text
      p1 <- p1 +
        theme(axis.text.y = element_blank())
    }
  } else {
    p1 <- NULL
  }
  #output
  p1
}


#' plot_alluvial_deseq_pair
#' @param x string to GroupA "A.vs.B"
#' @param y string to GroupB "A.vs.B"
#' @param feature string default: gene
#' @param .pd_rds string path to the file, saving read_deseq_pair() output
#'
#'
#' @export
plot_alluvial_deseq_pair <- function(x = NULL, y = NULL, .pd_rds = NULL, group = "all", ...) {
  # load data
  px <- stat_deseq_pair(x, y, .pd_rds, group = group, ...) # px
  if(! inherits(px, "list")) {
    warning("stat_deseq_pair failed")
    return(NULL)
  }
  ##------------------------------------------------##
  ## colors
  cc <- scales::hue_pal()(3) # red, green, blue
  cc <- c(cc[1], cc[3], cc[2]) # red, blue, green
  ## data
  df1 <- px$df3 %>%
    dplyr::filter(!trend %in% c("not-not")) %>%
    dplyr::select(-trend) %>%
    tidyr::pivot_longer(-gene_id, names_to = "sample", values_to = "sig") %>%
    dplyr::mutate(sig = factor(sig, levels = c("up", "not", "down")),
                  gene_id = as.character(gene_id))
  ## plot
  if(nrow(df1) > 0) {
    ## plots - Alluvial plot
    p1 <- df1 %>%
      ggplot(aes(x = sample, y = 1, stratum = sig, fill = sig,
                 alluvium = gene_id, label = sig)) +
      geom_flow() +
      geom_stratum(alpha = .8) +
      scale_fill_manual(values = cc) +
      ggtitle(px$pd$title) +
      geom_text(stat = "stratum") +
      scale_x_discrete(position = "top") +
      scale_y_continuous(expand = c(0, 0)) +
      ylab("Number of Genes") +
      theme_base() +
      theme(
        plot.title = element_text(size = 10),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y  = element_line(color = "black", size = .5),
        axis.ticks.y = element_line(color = "black", size = .5))
  } else {
    p1 <- NULL
  }
  ## output
  p1
}


#' deseq_pair_stat_plot
#' @param x string to GroupA "A.vs.B"
#' @param y string to GroupB "A.vs.B"
#' @param feature string default: gene
#' @param .pd_rds string path to the file, saving read_deseq_pair() output
#'
#'
#' @export
deseq_pair_stat_plot <- function(x    = NULL,
                                 y    = NULL,
                                 .pd_rds = NULL) {
  ## load data
  if(is.null(.pd_rds)) {
    pd <- read_deseq_pair(x, y, feature) # basic data
  } else {
    pd <- readRDS(.pd_rds)
  }

  # plot dat
  pa <- deseq_pair_stat(x, y, feature, .pd_rds)

  # checkpoint
  if(is.null(pa)) {
    warning("x, y required")
    return(NULL)
  }

  ##------------------------##
  ## colors
  cc <- scales::hue_pal()(3) # red, green, blue
  cc <- c(cc[1], cc[3], cc[2]) # red, blue, green
  pa$sig_plot_table %>%
    dplyr::mutate(sig = factor(sig, c("up", "not", "down"))) %>%
    ggplot(aes(count, sample, fill = sig)) +
    geom_col(color = "grey20") +
    geom_text(aes(label = count), position = position_stack(vjust = 0.5)) +
    scale_y_discrete(limits = c("GroupB", "GroupA")) +
    scale_fill_manual(values = cc) +
    ggtitle(pd$title) +
    theme_bw() +
    theme(legend.position = "right",
          legend.title    = element_blank(),
          panel.grid      = element_blank(),
          plot.title      = element_text(size = 10))

}


#' deseq_pair_heatmap_te
#' @param x string to GroupA "A.vs.B"
#' @param y string to GroupB "A.vs.B"
#' @param feature string default: gene
#' @param .pd_rds string path to the file, saving read_deseq_pair() output
#'
#'
#' @export
deseq_pair_heatmap_te <- function(x    = NULL,
                                  y    = NULL,
                                  feature = "te",
                                  topN    = 70,
                                  .pd_rds = NULL) {
  ## load data
  pd <- deseq_pair_stat(x, y, feature, .pd_rds)

  ## check feature: TE
  if(! pd$feature == "te") {
    warning("deseq_pair_heatmap_te() only for TEs, skipped.")
    return(NULL)
  }

  # checkpoint
  if(is.null(pd)) {
    warning("x, y required")
    return(NULL)
  }

  ##------------------------------------------------##
  ## prepare data
  fa <- bind_rows(pd$sig_genes$x) %>%
    dplyr::select(1:3)

  fb <- bind_rows(pd$sig_genes$y) %>%
    dplyr::select(1:3)

  ## set names
  colnames(fa) <- c("Gene", "A_ctl", "A_exp")
  colnames(fb) <- c("Gene", "B_ctl", "B_exp")

  df <- merge(fa, fb, by = "Gene") %>%
    dplyr::mutate(B_exp_norm = B_exp * (A_ctl / (B_ctl + .5))) %>%
    tidyr::separate("Gene", c("fb", "Gene"), sep = "_") %>%
    dplyr::select(Gene, A_ctl, A_exp, B_exp_norm)

  # pick top70 TEs
  topHits <- df %>%
    tidyr::pivot_longer(-Gene, names_to = "sample", values_to = "count") %>%
    dplyr::group_by(sample) %>%
    dplyr::arrange(desc(count)) %>%
    dplyr::top_n(topN) %>%
    dplyr::pull(Gene) %>%
    unique()

  # for plot
  df2 <- df %>%
    dplyr::filter(Gene %in% topHits) %>%
    dplyr::arrange(desc(B_exp_norm)) %>%
    dplyr::mutate_if(is.numeric, log10) %>%
    dplyr::mutate(y = row_number())

  if(nrow(df2) > 0) {

    # plot
    p1 <- df2 %>%
      tidyr::pivot_longer(-c(Gene, y), names_to = "sample", values_to = "RPM") %>%
      ggplot(aes(sample, Gene, fill = RPM)) +
      geom_tile(colour = "grey90") +
      scale_fill_gradient(low = "white", high = "blue2") +
      scale_y_discrete(position = "left", limits = rev(df2$Gene)) +
      scale_x_discrete(position = "top") +
      ggtitle(pd$title) +
      xlab(NULL) + ylab(NULL) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 8, color = "black"),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 8),
        axis.title = element_blank()
      )
  } else {
    p1 <- NULL
  }
  ## output
  p1
}






#' deseq_pair_alluvial_plot
#'
#' @param x string to GroupA "A.vs.B"
#' @param y string to GroupB "A.vs.B"
#' @param feature string default: gene
#' @param .pd_rds string path to the file, saving read_deseq_pair() output
#'
#'
#' @export
deseq_pair_alluvial_plot <- function(x    = NULL,
                                     y    = NULL,
                                     feature = "gene",
                                     .pd_rds = NULL) {
  ## load data
  pd <- deseq_pair_stat(x, y, feature, .pd_rds)

  # checkpoint
  if(is.null(pd)) {
    warning("x, y required")
    return(NULL)
  }

  ##------------------------------------------------##
  ## colors
  cc <- scales::hue_pal()(3) # red, green, blue
  cc <- c(cc[1], cc[3], cc[2]) # red, blue, green

  ## data
  df <- pd$sig_list %>%
    dplyr::filter(!trend %in% c("not-not")) %>%
    dplyr::select(-trend) %>%
    tidyr::pivot_longer(-Gene, names_to = "sample", values_to = "sig") %>%
    dplyr::mutate(sig = factor(sig, levels = c("up", "not", "down")))

  if(nrow(df) > 0) {
    ## plots - Alluvial plot
    p1 <- df %>%
      ggplot(aes(x = sample, y = 1, stratum = sig, fill = sig,
                 alluvium = Gene, label = sig)) +
      geom_flow() +
      geom_stratum(alpha = .8) +
      scale_fill_manual(values = cc) +
      ggtitle(pd$title) +
      geom_text(stat = "stratum") +
      scale_x_discrete(position = "top") +
      scale_y_continuous(expand = c(0, 0)) +
      ylab("Number of Genes") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 10),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y  = element_line(color = "black", size = .5),
        axis.ticks.y = element_line(color = "black", size = .5))
  } else {
    p1 <- NULL
  }

  ## output
  p1
}




#' compare_deseq
#'
#' @param x string to GroupA "A.vs.B"
#' @param y string to GroupB "A.vs.B"
#' @param feature string default: gene
#'
#'
#' @export
compare_deseq <- function(x, y, feature = "gene") {
  # assign group
  deseq_dirs <- setNames(c(x, y), c("GroupA", "GroupB"))
  # config
  pa <- read_hiseq(x)
  pb <- read_hiseq(y)
  # global: title
  ptitle <- paste(c(glue::glue("GroupA: {basename(x)}"),
                    glue::glue("GroupB: {basename(y)}"),
                    "Critera: foldChange > 2 & pvalue < 0.05"),
                  collapse = "\n")
  ##------------------------##
  ## gene list
  df1 <- get_fix_xls(x) %>%
    readr::read_delim("\t", col_types = readr::cols())
  df2 <- get_fix_xls(y) %>%
    readr::read_delim("\t", col_types = readr::cols())

  # gene count
  df_count <- lapply(c(x, y), function(i) {
    get_fix_xls(i) %>%
      readr::read_delim("\t", col_types = readr::cols()) %>%
      dplyr::pull(sig) %>%
      table()
  }) %>%
    dplyr::bind_rows()
  df_count$sample <- c("GroupA", "GroupB")

  # for plot
  df_count_plot <- df_count %>%
    tidyr::pivot_longer(names_to = "group", values_to = "count", down:up) %>%
    dplyr::filter(group %in% c("up", "down"))


  # for table (with gene names)
  df_table <- lapply(names(deseq_dirs), function(i) {
    deseq_dirs[i] %>%
      get_fix_xls() %>%
      readr::read_delim("\t", col_types = readr::cols()) %>%
      dplyr::select(Gene, SYMBOL, sig, log2FoldChange) %>%
      dplyr::mutate(log2FoldChange = round(log2FoldChange, 1)) %>%
      tidyr::unite("fc", c("sig", "log2FoldChange"), sep = ":") %>%
      dplyr::mutate(group = i)
  }) %>%
    dplyr::bind_rows()
  #!!!! to-to
  df_table2 <- df_table %>%
    tidyr::pivot_wider(names_from = "group", values_from = "fc") %>%
    dplyr::filter(! (grepl("not", GroupA) & grepl("not", GroupB))) %>%
    dplyr::mutate(Trend = paste(
      gsub(":.*", "", GroupA),
      gsub(":.*", "", GroupB),
      sep = "-"
    ))

  # for alluvial plot
  df_alluvial <- df_table %>%
    dplyr::mutate(sig = gsub(":.*", "", fc)) %>%
    dplyr::select(-fc) %>%
    dplyr::mutate(sig = factor(sig, levels = c("up", "not", "down")))

  ## for heatmap plot
  df_heatmap <- df_table %>%
    tidyr::separate("fc", c("sig", "log2fc"), sep = ":") %>%
    tidyr::pivot_wider(names_from = "group", values_from = "log2fc")


  dplyr::mutate(score = plyr::mapvalues(sig,
                                        c("up", "down", "not"),
                                        c(1, -1, 0))) %>%
    dplyr::mutate(score = as.numeric(as.character(score))) %>%
    dplyr::select(-sig) %>%
    tidyr::spread("sample", "score") %>%
    tibble::column_to_rownames("Gene")

  ##------------------------##
  ## colors
  cc <- scales::hue_pal()(3) # red, green, blue
  cc <- c(cc[1], cc[3], cc[2]) # red, blue, green

  ## plots - count
  plot_count <- df_count_plot %>%
    dplyr::mutate(sig = factor(sig, c("up", "not", "down"))) %>%
    ggplot(aes(count, sample, fill = sig)) +
    geom_col(color = "grey20") +
    geom_text(aes(label = count), position = position_stack(vjust = 0.5)) +
    scale_y_discrete(limits = c("GroupB", "GroupA")) +
    scale_fill_manual(values = cc) +
    ggtitle(ptitle) +
    theme_bw() +
    theme(legend.position = "right",
          legend.title    = element_blank(),
          panel.grid      = element_blank(),
          plot.title      = element_text(size = 10))

  ## plots - overlap2
  plot_overlap_list <- sapply(c("up", "down", "not"), function(i){
    # # gene list
    list(GroupA = g1[[i]]$Gene,
         GroupB = g2[[i]]$Gene) %>%
      ggvenn::ggvenn(fill_color    = scales::hue_pal()(2),
                     fill_alpha    = .8,
                     stroke_size   = .5,
                     set_name_size = 4) +
      ggtitle(i)
  }, simplify = FALSE, USE.NAMES = TRUE)

  ## combine, up/down
  plot_overlap2 <- cowplot::plot_grid(plot_overlap_list$up,
                                      plot_overlap_list$down,
                                      ncol = 2)
  plot_overlap3 <- cowplot::ggdraw(cowplot::add_sub(plot_overlap2,
                                                    ptitle,
                                                    x     = 0,
                                                    hjust = 0,
                                                    size  = 10))

  ## plots - overlap4
  plot_overlap4 <- list(GroupA_up   = g1[["up"]]$Gene,
                        GroupA_down = g1[["down"]]$Gene,
                        GroupB_up   = g2[["up"]]$Gene,
                        GroupB_down = g2[["down"]]$Gene) %>%
    ggvenn::ggvenn(fill_color      = scales::hue_pal()(4),
                   fill_alpha      = .8,
                   set_name_size   = 3.5,
                   show_percentage = FALSE,
                   stroke_size     = .5)

  ## plots - Alluvial plot
  plot_alluvial <- df_alluvial %>%
    ggplot(aes(x = sample, y = 1, stratum = sig, fill = sig,
               alluvium = Gene, label = sig)) +
    geom_flow() +
    geom_stratum(alpha = .8) +
    scale_fill_manual(values = cc) +
    ggtitle(ptitle) +
    geom_text(stat = "stratum") +
    theme_minimal() +
    theme(plot.title = element_text(size = 10))

  ## plots - heatmap
  plot_heatmap <- pheatmap::pheatmap(as.matrix(df_heatmap),
                                     color          = rev(cc),
                                     silent         = T,
                                     border_color   = "grey60",
                                     show_rownames  = FALSE,
                                     cluster_cols   = FALSE,
                                     legend_labels  = c("down", "not", "up"),
                                     legend_breaks  = c(-1, 0, 1))
  ## output obj
  list(x          = x,
       y          = y,
       ptitle        = ptitle,
       listA         = g1,
       listB         = g2,
       df_count      = df_count,
       df_table      = df_table,
       plot_count    = plot_count,
       plot_overlap2 = plot_overlap2,
       plot_overlap3 = plot_overlap3,
       plot_overlap4 = plot_overlap4,
       plot_heatmap  = plot_heatmap,
       plot_alluvial = plot_alluvial)
}

