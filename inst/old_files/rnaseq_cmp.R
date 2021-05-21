

#' functions for pair deseq
#'
#' @param dirA string to GroupA "A.vs.B"
#' @param dirB string to GroupB "A.vs.B"
#' @param feature string default: gene
#' @param outdir string directory to save files
#'
#' @export
# pair_deseq_report <- function(dirA, dirB, feature, outdir) {
deseq_pair_report <- function(dirA, dirB, outdir) {
  message("# Run RNAseq pair ...")
  ##---------------------##
  # report
  subdir   <- paste0(basename(dirA), ".compare.", basename(dirB))
  outdir   <- normalizePath(outdir)
  # outdir   <- file.path(outdir, subdir)

  if(! dir.exists(outdir)) {
    dir.create(outdir, mode = "0755", recursive = TRUE)
  }

  template <- system.file("rnaseq", "deseq_pair_report.Rmd",
                          package = "hiseqr")
  template_to <- file.path(outdir, basename(template))
  ## copy Rmd
  file.copy(template, template_to)
  out_html <- file.path(outdir, "deseq_pair_report.html")

  ## run Rmd
  if(file.exists(out_html)) {
    message("output exists, skipping ...")
  } else {
    rmarkdown::render(input       = template,
                      output_file = out_html,
                      params      = list(dirA = dirA,
                                         dirB = dirB))
  }
}




#' RNAseq cmp, stat
#' number of sig genes
#'
#' @param dirA string to GroupA "A.vs.B"
#' @param dirB string to GroupB "A.vs.B"
#' @param feature string default: gene
#'
#' @import dplyr
#'
#' @export
read_deseq_pair <- function(dirA, dirB, feature = "gene") {
  # check if dirA/dirB is "deseq_single"
  if(is.null(dirA) | is.null(dirB)) {
    warning("dirA, dirB not detected")
    return(NULL)
  } else if(is_deseq_single_dir(dirA) & is_deseq_single_dir(dirB)) {
    message("input are deseq_single dirs")
  } else {
    warning("input are not deseq_single, skipped")
    return(NULL)
  }

  ## convert to absolute path
  dirA <- normalizePath(dirA)
  dirB <- normalizePath(dirB)

  ##------------------------##
  ## gene list
  sig_genes <- list(dirA = get_sig_gene(dirA, feature, gene_name_only = FALSE),
                    dirB = get_sig_gene(dirB, feature, gene_name_only = FALSE))

  ##------------------------##
  ## labels
  argsA <- deseq_single_dir(dirA)
  argsB <- deseq_single_dir(dirB)

  if(argsA$genome != argsB$genome) {
    warning(paste("reference genome not inconsistent:", argsA$genome, argsB$genome, sep = " "))
    return(NULL)
  }

  # global: title
  title <- paste(c(glue::glue("GroupA: {basename(dirA)}"),
                   glue::glue("GroupB: {basename(dirB)}"),
                   "Critera: foldChange > 2 & pvalue < 0.05"),
                 collapse = "\n")

  # output
  list(sig_genes = sig_genes,
       title     = title,
       dirA      = dirA,
       dirB      = dirB,
       feature   = feature)
}



#' RNAseq cmp, stat
#' number of sig genes
#'
#' @param dirA string to GroupA "A.vs.B"
#' @param dirB string to GroupB "A.vs.B"
#' @param feature string default: gene
#' @param kepp_not_sig logical whether keep not sig changed genes, default: TRUE
#' @param .pd_rds string path to the file, saving read_deseq_pair() output
#'
#' @import dplyr
#'
#' @export
deseq_pair_stat <- function(dirA         = NULL,
                            dirB         = NULL,
                            feature      = "gene",
                            .pd_rds      = NULL,
                            keep_not_sig = TRUE) {
  # loading data from rds file
  if(is.null(.pd_rds)) {
    pd <- read_deseq_pair(dirA, dirB, feature)
  } else {
    pd <- readRDS(.pd_rds) # dirA, dirB
    dirA    <- pd$dirA
    dirB    <- pd$dirB
    feature <- pd$feature
  }

  # checkpoint
  if(is.null(pd)) {
    warning("dirA, dirB required")
    return(NULL)
  }

  # not_sig: off, feature=gene
  if(feature == "gene") {
    keep_not_sig = FALSE
  }

  # sig genes, count
  df1 <- merge(sapply(pd$sig_genes$dirA, nrow) %>%
                 as.data.frame() %>%
                 tibble::rownames_to_column("sig"),
               sapply(pd$sig_genes$dirB, nrow) %>%
                 as.data.frame() %>%
                 tibble::rownames_to_column("sig"),
               by = "sig") %>%
    dplyr::rename(GroupA = "..x", GroupB = "..y")

  # sig genes, plot
  df2 <- df1 %>%
    tidyr::gather("sample", "count", -1) %>%
    dplyr::filter(sig %in% c("up", "down"))

  # for table (with gene names)
  gl <- lapply(list(pd$sig_genes$dirA,
                    pd$sig_genes$dirB), function(i){
                      dplyr::bind_rows(i) %>%
                        dplyr::select(Gene, sig)
                    })

  # # id conversion
  # df3 <- dplyr::bind_cols(gl) %>%
  #   dplyr::select(- Gene...3) %>%
  #   dplyr::rename(GroupA = sig...2,
  #                 GroupB = sig...4,
  #                 Gene   = Gene...1) %>%
  #   dplyr::mutate(trend = paste(GroupA, GroupB, sep = "-"))

  df3 <- merge(gl[[1]], gl[[2]], by = "Gene") %>%
    dplyr::rename(GroupA = sig.x,
                  GroupB = sig.y) %>%
    dplyr::mutate(trend  = paste(GroupA, GroupB, sep = "-"))

  # remove not_sig
  if(! isTRUE(keep_not_sig)) {
    df3 <- df3 %>%
      dplyr::filter(! trend %in% c("not-not"))
      # dplyr::filter(! grepl("not", trend))
  }

  # output
  c(pd,
    list(sig_count_table = df1,
       sig_plot_table  = df2,
       sig_list        = df3)
  )
}


#' RNAseq cmp,
#' number of sig genes
#' .pd_rds is read_deseq_pair() output
#'
#' @param dirA string to GroupA "A.vs.B"
#' @param dirB string to GroupB "A.vs.B"
#' @param feature string default: gene
#' @param .pd_rds string path to the file, saving read_deseq_pair() output
#'
#' @import dplyr
#'
#' @export
deseq_pair_stat_plot <- function(dirA    = NULL,
                                 dirB    = NULL,
                                 feature = "gene",
                                 .pd_rds = NULL) {
  ## load data
  if(is.null(.pd_rds)) {
    pd <- read_deseq_pair(dirA, dirB, feature) # basic data
  } else {
    pd <- readRDS(.pd_rds)
  }

  # plot dat
  pa <- deseq_pair_stat(dirA, dirB, feature, .pd_rds)

  # checkpoint
  if(is.null(pa)) {
    warning("dirA, dirB required")
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



#' RNAseq cmp,
#' number of sig genes
#' .pd_rds is read_deseq_pair() output
#'
#' @param dirA string to GroupA "A.vs.B"
#' @param dirB string to GroupB "A.vs.B"
#' @param feature string default: gene
#' @param .pd_rds string path to the file, saving read_deseq_pair() output
#'
#' @import dplyr
#'
#' @export
deseq_pair_overlap_plot <- function(dirA    = NULL,
                                    dirB    = NULL,
                                    feature = "gene",
                                    .pd_rds = NULL) {
  ## load data
  pd <- deseq_pair_stat(dirA, dirB, feature, .pd_rds)

  # checkpoint
  if(is.null(pd)) {
    warning("dirA, dirB required")
    return(NULL)
  }

  ## plots - overlap2
  p_list1 <- sapply(c("up", "down", "not"), function(i) {
    list(GroupA = pd$sig_genes$dirA[[i]]$Gene,
         GroupB = pd$sig_genes$dirB[[i]]$Gene) %>%
      ggvenn::ggvenn(fill_color    = scales::hue_pal()(2),
                     fill_alpha    = .8,
                     stroke_size   = .5,
                     set_name_size = 4) +
      ggtitle(i)
  }, simplify = FALSE, USE.NAMES = TRUE)

  ## combine, up/down
  p1 <- patchwork::wrap_plots(p_list1, nrow = 1) +
    plot_annotation(
      title = pd$title,
      caption = 'made with patchwork'
    )

  ## plots - overlap4
  ## up/down, pair
  p2 <- list(A_up = pd$sig_genes$dirA[["up"]]$Gene,
             B_up = pd$sig_genes$dirB[["up"]]$Gene,
             A_down = pd$sig_genes$dirA[["down"]]$Gene,
             B_down = pd$sig_genes$dirB[["down"]]$Gene) %>%
    ggvenn::ggvenn(fill_color    = scales::hue_pal()(4),
                   show_percentage = FALSE,
                   fill_alpha    = .8,
                   stroke_size   = .5,
                   set_name_size = 4) +
    ggtitle(pd$title)
  ## output
  c(p_list1, list(merge = p1, up_down = p2))
}


#' cluster
#' data.frame
#'
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


#' to-do: hclust rows !!!!
#' RNAseq cmp,
#' number of sig genes
#'
#' @param dirA string to GroupA "A.vs.B"
#' @param dirB string to GroupB "A.vs.B"
#' @param feature string default: gene
#' @param .pd_rds string path to the file, saving read_deseq_pair() output
#'
#' @import dplyr
#'
#' @export
deseq_pair_heatmap <- function(dirA    = NULL,
                               dirB    = NULL,
                               feature = "gene",
                               .pd_rds = NULL) {
  ## load data
  pd <- deseq_pair_stat(dirA, dirB, feature, .pd_rds)

  # checkpoint
  if(is.null(pd)) {
    warning("dirA, dirB required")
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

  ## for cluster
  ma <- df %>%
    dplyr::mutate(sig = plyr::mapvalues(sig, c("up", "not", "down"),
                                        c(1, 0, -1),
                                        warn_missing = FALSE)) %>%
    dplyr::mutate(sig = as.numeric(as.character(sig))) %>%
    pivot_wider(names_from = "sample", values_from = "sig") %>%
    tibble::column_to_rownames("Gene") %>%
    as.matrix() %>%
    .cluster_ma()

  ## plot
  if(nrow(df) > 0) {
    p1 <- df %>%
      ggplot(aes(sample, Gene, fill = sig)) +
      geom_tile() +
      scale_fill_manual(values = cc) +
      scale_x_discrete(position = "top") +
      scale_y_discrete(limits = rownames(ma)) +
      ggtitle(pd$title) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.title = element_blank()
      )

    ## remove y-text
    if(nrow(df) > 100) {
      # remove y text
      p1 <- p1 +
        theme(axis.text.y = element_blank())
    }
  } else {
    p1 <- NULL
  }

  ## output
  p1
}



#' to-do: hclust rows !!!!
#' RNAseq cmp,
#' number of sig genes
#'
#' @param dirA string to GroupA "A.vs.B"
#' @param dirB string to GroupB "A.vs.B"
#' @param feature string default: gene
#' @param .pd_rds string path to the file, saving read_deseq_pair() output
#'
#' @import dplyr
#'
#' @export
deseq_pair_heatmap_te <- function(dirA    = NULL,
                                  dirB    = NULL,
                                  feature = "te",
                                  topN    = 70,
                                  .pd_rds = NULL) {
  ## load data
  pd <- deseq_pair_stat(dirA, dirB, feature, .pd_rds)

  ## check feature: TE
  if(! pd$feature == "te") {
    warning("deseq_pair_heatmap_te() only for TEs, skipped.")
    return(NULL)
  }

  # checkpoint
  if(is.null(pd)) {
    warning("dirA, dirB required")
    return(NULL)
  }

  ##------------------------------------------------##
  ## prepare data
  fa <- bind_rows(pd$sig_genes$dirA) %>%
    dplyr::select(1:3)

  fb <- bind_rows(pd$sig_genes$dirB) %>%
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






#' RNAseq cmp,
#' number of sig genes
#'
#' @param dirA string to GroupA "A.vs.B"
#' @param dirB string to GroupB "A.vs.B"
#' @param feature string default: gene
#' @param .pd_rds string path to the file, saving read_deseq_pair() output
#'
#' @import dplyr
#'
#' @export
deseq_pair_alluvial_plot <- function(dirA    = NULL,
                               dirB    = NULL,
                               feature = "gene",
                               .pd_rds = NULL) {
  ## load data
  pd <- deseq_pair_stat(dirA, dirB, feature, .pd_rds)

  # checkpoint
  if(is.null(pd)) {
    warning("dirA, dirB required")
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



##-----------------------------##
## functions for two group pair
## 1. sig genes
## 2. overlap genes
## 3. plot: sig genes
## 4. plot: overlap genes
## 5. TE scatter
## 6. TE heatmap

#' functions for compare two deseq
#'
#' @param dirA string to GroupA "A.vs.B"
#' @param dirB string to GroupB "A.vs.B"
#' @param feature string default: gene
#'
#' @import ggplot2
#' @import dplyr
#' @import pheatmap
#' @import ggalluvial
#' @import RColorBrewer
#' @import clusterProfiler
#'
#' @export
compare_deseq <- function(dirA, dirB, feature = "gene") {
  # check if dirA/dirB is "deseq_single"
  chk0 <- is_hiseq_dir(dirA) == "deseq_single"
  chk1 <- is_hiseq_dir(dirB) == "deseq_single"

  if(! chk0) {
    stop(paste0("not a deseq_single dir, ", dirA))
  }

  if(! chk1) {
    stop(paste0("not a deseq_single dir, ", dirB))
  }

  # args
  argsA <- deseq_single_dir(dirA)
  argsB <- deseq_single_dir(dirB)

  if(argsA$genome != argsB$genome) {
    stop(paste("reference genome not inconsistent:", argsA$genome, argsB$genome, sep = " "))
  }

  # global: title
  ptitle <- paste(c(glue::glue("GroupA: {basename(dirA)}"),
                    glue::glue("GroupB: {basename(dirB)}"),
                    "Critera: foldChange > 2 & pvalue < 0.05"),
                  collapse = "\n")

  ##------------------------##
  ## gene list
  g1 <- get_sig_gene(dirA, feature, gene_name_only = FALSE) # up/not/down
  g2 <- get_sig_gene(dirB, feature, gene_name_only = FALSE) # up/not/down

  # for count
  df_count <- merge(sapply(g1, nrow) %>%
                      as.data.frame() %>%
                      tibble::rownames_to_column("sig"),
                    sapply(g2, nrow) %>%
                      as.data.frame() %>%
                      tibble::rownames_to_column("sig"),
                    by = "sig") %>%
    dplyr::rename(GroupA = "..x", GroupB = "..y")

  # for plot
  df_count_plot <- df_count %>%
    tidyr::gather("sample", "count", -1) %>%
    dplyr::filter(sig %in% c("up", "down"))

  # for table (with gene names)
  a = lapply(list(g1, g2), function(i){
    dplyr::bind_rows(i) %>%
      dplyr::select(Gene, sig)
  })

  # id conversion
  df_table <- merge(a[[1]], a[[2]], by = "Gene") %>%
    dplyr::filter(! (sig.x == "not" & sig.y == "not")) %>%
    dplyr::rename(GroupA = sig.x,
                  GroupB = sig.y) %>%
    dplyr::mutate(trend = paste(GroupA, GroupB, sep = "-"))

  # convert to symbol
  df_table$symbol <- gene_to_symbol(df_table$Gene, argsA$genome)

  # for alluvial plot
  df_alluvial <- df_table %>%
    dplyr::select(Gene, GroupA, GroupB) %>%
    tidyr::gather("sample", "sig", -1) %>%
    dplyr::mutate(sig = factor(sig, levels = c("up", "not", "down")))

  ## for heatmap plot
  df_heatmap <- df_alluvial %>%
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
                                     # treeheight_row = 2,
                                     legend_labels  = c("down", "not", "up"),
                                     legend_breaks  = c(-1, 0, 1))
  ## output obj
  list(dirA          = dirA,
       dirB          = dirB,
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



#'
#' ##-----------------------------##
#' ## functions for two group compare
#' ## 1. sig genes
#' ## 2. overlap genes
#' ## 3. plot: sig genes
#' ## 4. plot: overlap genes
#'
#' #' functions for compare two deseq
#' #'
#' #' @param dirA string to GroupA "A.vs.B"
#' #' @param dirB string to GroupB "A.vs.B"
#' #' @param feature string default: gene
#' #'
#' #' @import ggplot2
#' #' @import dplyr
#' #' @import pheatmap
#' #' @import ggalluvial
#' #' @import RColorBrewer
#' #' @import clusterProfiler
#' #'
#' #' @export
#' compare_deseq <- function(dirA, dirB, feature = "gene") {
#'   # check if dirA/dirB is "deseq_single"
#'   chk0 <- hiseqr::is_hiseq_dir(dirA) == "deseq_single"
#'   if(! chk0) {
#'     stop(paste0("not a deseq_single dir, ", dirA))
#'   }
#'
#'   chk1 <- hiseqr::is_hiseq_dir(dirB) == "deseq_single"
#'   if(! chk1) {
#'     stop(paste0("not a deseq_single dir, ", dirB))
#'   }
#'
#'   # global: title
#'   ptitle <- paste(c(glue::glue("GroupA: {basename(dirA)}"),
#'                     glue::glue("GroupB: {basename(dirB)}"),
#'                     "Critera: foldChange > 2 & pvalue < 0.05"),
#'                   collapse = "\n")
#'
#'   ##------------------------##
#'   ## gene list
#'   ## get sig gene list
#'   g1 <- get_sig_gene(dirA, feature, gene_name_only = TRUE) # up/not/down
#'   g2 <- get_sig_gene(dirB, feature, gene_name_only = TRUE) # up/not/down
#'
#'   # list to data.frame
#'   to_df <- function(x, summary = FALSE) {
#'     df <- data.frame(id = unlist(x)) %>%
#'       tibble::rownames_to_column("sig") %>%
#'       dplyr::mutate(sig = gsub("[0-9]", "", sig))
#'
#'     if(isTRUE(summary)) {
#'       df <- df %>%
#'         dplyr::group_by(sig) %>%
#'         dplyr::summarise(count = n())
#'     }
#'
#'     return(df)
#'   }
#'
#'   # for summary
#'   df_sig_count <- merge(to_df(g1, TRUE), to_df(g2, TRUE), by = "sig") %>%
#'     dplyr::rename(GroupA = count.x,
#'                   GroupB = count.y) %>%
#'     tidyr::gather("sample", "count", -sig) %>%
#'     tidyr::spread("sig", "count")
#'
#'   # for plot
#'   df2 <- merge(to_df(g1), to_df(g2), by = "id") %>%
#'     dplyr::rename(GroupA = sig.x,
#'                   GroupB = sig.y) %>%
#'     dplyr::filter(! (GroupA == "not" & GroupB == "not"))
#'
#'   # for table
#'   df_sig_both <- df2 %>%
#'     dplyr::mutate(group = paste(GroupA, GroupB, sep = "-")) %>%
#'     dplyr::mutate(symbol = AnnotationDbi::mapIds(org.Dm.eg.db::org.Dm.eg.db,
#'                                                  keys = id,
#'                                                  column = "SYMBOL",
#'                                                  keytype = "ENSEMBL",
#'                                                  multiVals = "first"))
#'
#'   ##------------------------##
#'   ## colors
#'   cc <- scales::hue_pal()(3) # red, green, blue
#'   cc <- c(cc[1], cc[3], cc[2]) # red, blue, green
#'
#'   ## plots - count
#'   plot_count <- df_sig_count %>%
#'     tidyr::gather("sig", "count", -sample) %>%
#'     dplyr::mutate(sig = factor(sig, c("up", "not", "down"))) %>%
#'     ggplot(aes(count, sample, fill = sig)) +
#'     geom_col(color = "grey20") +
#'     geom_text(aes(label = count), position = position_stack(vjust = 0.5)) +
#'     scale_y_discrete(limits = c("GroupB", "GroupA")) +
#'     scale_fill_manual(values = cc) +
#'     ggtitle(ptitle) +
#'     theme_bw() +
#'     theme(legend.position = "right",
#'           legend.title    = element_blank(),
#'           panel.grid      = element_blank(),
#'           plot.title      = element_text(size = 10))
#'
#'   ## plots - overlap2
#'   plot_overlap_list <- lapply(c("up", "down", "not"), function(i){
#'     # gene list
#'     a <- list(GroupA = g1[[i]],
#'               GroupB = g2[[i]])
#'
#'     # plot
#'     p <- ggvenn::ggvenn(a,
#'                         fill_color = scales::hue_pal()(2),
#'                         fill_alpha = .8,
#'                         set_name_size = 4,
#'                         stroke_size = .5) +
#'       ggtitle(i)
#'   })
#'   names(plot_overlap_list) <- c("up", "down", "not")
#'   ## combine, up/down
#'   plot_overlap2 <- cowplot::plot_grid(plot_overlap_list$up,
#'                                       plot_overlap_list$down,
#'                                       ncol = 2)
#'   plot_overlap3 <- cowplot::ggdraw(cowplot::add_sub(plot_overlap2,
#'                                                     ptitle,
#'                                                     x     = 0,
#'                                                     hjust = 0,
#'                                                     size  = 10))
#'
#'   ## plots - overlap4
#'   g_list <- list(GroupA_up = g1$up,
#'                  GroupA_down = g1$down,
#'                  GroupB_up = g2$up,
#'                  GroupB_down = g2$down)
#'   plot_overlap4 <- ggvenn::ggvenn(g_list,
#'                                   fill_color      = scales::hue_pal()(4),
#'                                   fill_alpha      = .8,
#'                                   set_name_size   = 3.5,
#'                                   show_percentage = FALSE,
#'                                   stroke_size     = .5)
#'
#'   ## plots - Alluvial plot
#'   plot_alluvial <- df2 %>%
#'     tidyr::gather("sample", "sig", -id) %>%
#'     dplyr::mutate(sig = factor(sig, levels = c("up", "not", "down"))) %>%
#'     ggplot(aes(x = sample, y = 1, stratum = sig, fill = sig,
#'                alluvium = id, label = sig)) +
#'     geom_flow() +
#'     geom_stratum(alpha = .8) +
#'     scale_fill_manual(values = cc) +
#'     ggtitle(ptitle) +
#'     geom_text(stat = "stratum") +
#'     theme_minimal() +
#'     theme(plot.title = element_text(size = 10))
#'
#'   ## plots - heatmap
#'   ma <- df2 %>%
#'     dplyr::mutate(GroupA = plyr::mapvalues(GroupA, c("up", "down", "not"), c(1, -1, 0)),
#'                   GroupB = plyr::mapvalues(GroupB, c("up", "down", "not"), c(1, -1, 0))) %>%
#'     dplyr::mutate(GroupA = as.numeric(GroupA),
#'                   GroupB = as.numeric(GroupB)) %>%
#'     tibble::column_to_rownames("id") %>%
#'     as.matrix()
#'
#'   plot_heatmap <- pheatmap::pheatmap(ma,
#'                                      color          = rev(cc),
#'                                      silent         = T,
#'                                      border_color   = "grey60",
#'                                      show_rownames  = FALSE,
#'                                      cluster_cols   = FALSE,
#'                                      # treeheight_row = 2,
#'                                      legend_labels  = c("down", "not", "up"),
#'                                      legend_breaks  = c(-1, 0, 1))
#'   ## output obj
#'   list(dirA          = dirA,
#'        dirB          = dirB,
#'        ptitle        = ptitle,
#'        listA         = g1,
#'        listB         = g2,
#'        df_sig_count  = df_sig_count,
#'        df_sig_both   = df_sig_both,
#'        plot_count    = plot_count,
#'        plot_overlap2 = plot_overlap2,
#'        plot_overlap3 = plot_overlap3,
#'        plot_overlap4 = plot_overlap4,
#'        plot_heatmap  = plot_heatmap,
#'        plot_alluvial = plot_alluvial)
#' }
#'
