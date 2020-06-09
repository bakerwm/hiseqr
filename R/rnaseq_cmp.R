
##-----------------------------##
## functions for two group compare
## 1. sig genes
## 2. overlap genes
## 3. plot: sig genes
## 4. plot: overlap genes

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
    # gene list
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




#' functions for compare two deseq
#'
#' @param dirA string to GroupA "A.vs.B"
#' @param dirB string to GroupB "A.vs.B"
#' @param feature string default: gene
#' @param outdir string directory to save files
#'
#' @export
compare_deseq_report <- function(dirA, dirB, feature, outdir) {
  message("# Run RNAseq compare ...")
  if(! dir.exists(outdir)) {
    dir.create(outdir, mode = "0755")
  }

  ##---------------------##
  # report
  outdir <- normalizePath(outdir)
  template <- system.file("rnaseq", "RNAseq_cmp_report.Rmd",
                          package = "hiseqr")
  template_to <- file.path(outdir, basename(template))
  ## copy Rmd
  file.copy(template, template_to)
  out_html <- file.path(outdir, "RNAseq_cmp_report.html")

  ## run Rmd
  if(file.exists(out_html)) {
    message("output exists, skipping ...")
  } else {
    rmarkdown::render(input       = template,
                      output_file = out_html,
                      params      = list(dirA = dirA,
                                         dirB = dirB,
                                         feature    = feature))
  }
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
