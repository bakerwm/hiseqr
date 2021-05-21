#' Functions for check P7 structure
#'
#' Check the structure of P7, TruSeq, Nextera
#' Example:
# --------------------------------------------------------------------------------
#   fastq: data/raw_data/smRNAseq_shWhite_2_4h_rep4_1.fq.gz
# total:      101,513
# group   total   count   pct     sample
# #P7a_i7 101513  97494   96.0    smRNAseq_shWhite_2_4h_rep4_1
# #P7a_7b 101513  95178   93.8    smRNAseq_shWhite_2_4h_rep4_1
# #i7_P7b 101513  98869   97.4    smRNAseq_shWhite_2_4h_rep4_1
# --------------------------------------------------------------------------------
#   demx    demx_rc demx_name       demx_n  demx_pct        group   total   i7      i7_pct  sample
# >CACGAT CACGAT  TruSeq_Index31  94135    96.6   i7      101513  97494   96.0    smRNAseq_shWhite_2_4h_rep4_1
# >CAAAAG CAAAAG  TruSeq_Index28  2237      2.3   i7      101513  97494   96.0    smRNAseq_shWhite_2_4h_rep4_1
# >CACGAA CACGAA  CACGAA  203       0.2   i7      101513  97494   96.0    smRNAseq_shWhite_2_4h_rep4_1
# >CGATGT ACATCG  P7_3B   97166    99.7   bc      101513  97494   96.0    smRNAseq_shWhite_2_4h_rep4_1
# >CGAAGT ACTTCG  ACTTCG  167       0.2   bc      101513  97494   96.0    smRNAseq_shWhite_2_4h_rep4_1
# >CGATGA TCATCG  TCATCG  27        0.0   bc      101513  97494   96.0    smRNAseq_shWhite_2_4h_rep4_1
#'
#'
#' @name hiseq_p7


#' @describeIn processing the P7 structure
#' @import patchwork
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#'
#' @param group character options ["i7", "bc"]
hiseq_p7 <- function(x) {
  lapply(x, function(f) {
    ##--------------------##
    a <- readLines(f)
    ## load i7+bc
    ## line startswith ">"
    a <- a[startsWith(a, ">")] # only lines startswith ">"
    a <- gsub(" |>", "", a)    # remove ">" and "white spaces"
    df1 <- stringr::str_split(a, "\t", simplify = TRUE) %>%
      as.data.frame() %>%
      dplyr::mutate_at(c(4:5, 7:9), as.numeric)
    colnames(df1) <- c("demx", "demx_rc", "demx_name", "demx_n", "demx_pct",
                       "group", "total", "i7", "i7_pct", "sample")
    df1 <- df1 %>%
      dplyr::group_by(group) %>%
      dplyr::mutate(rank = row_number()) %>%
      dplyr::mutate(rank = as.character(rank))
    ##--------------------##
    # plots
    df2 <- df1 %>%
      dplyr::filter(group == "i7")
    p1 <- df2 %>%
      hiseqr::bar_plot(x = "demx_pct", y = "rank", label = "demx_name",
                       group = "rank", direction = "horizontal") +
      geom_vline(xintercept = c(50, 100), color = "blue", linetype = 2) +
      scale_x_continuous(limits = c(0, 150), position = "top") +
      ggtitle(df2$sample[1]) +
      theme(axis.title.x = element_blank())
    # bc
    df3 <- df1 %>%
      dplyr::filter(group == "bc")
    p2 <- df3 %>%
      hiseqr::bar_plot(x = "demx_pct", y = "rank", label = "demx_name",
                       group = "rank", direction = "horizontal") +
      geom_vline(xintercept = c(50, 100), color = "blue", linetype = 2) +
      scale_x_continuous(limits = c(0, 150), position = "top") +
      ggtitle(df3$sample[1]) +
      theme(axis.title.x = element_blank())
    ##--------------------##
    ## load p7
    ## line starswith "#"
    b <- readLines(f)
    b <- b[startsWith(b, "#")]
    b <- gsub(" |#", "", b) # remove "#" and "white spaces"
    df4 <- stringr::str_split(b, "\t", simplify = TRUE) %>%
      as.data.frame() %>%
      dplyr::mutate_at(2:4, as.numeric)
    colnames(df4) <- c("group", "total", "count", "pct", "sample")
    df4$group <- factor(
      df4$group,
      levels = c("no_ad", "p7b", "p7a_p7b", "p7a"))
    # ## structure
    # df5 <- df4 %>%
    #   dplyr::select(-pct) %>%
    #   tidyr::pivot_wider(names_from = "group", values_from = "count") %>%
    #   dplyr::mutate(P7a = P7a_i7 - P7a_7b,
    #                 P7b = i7_P7b - P7a_7b,
    #                 no_ad = total - P7a - P7b - P7a_7b) %>%
    #   dplyr::select(total, sample, P7a, P7a_7b, P7b, no_ad) %>%
    #   tidyr::pivot_longer(names_to = "group", values_to = "count", P7a:no_ad) %>%
    #   dplyr::mutate(pct = round(count / total * 100, 1))
    ## plot
    p3 <- df4 %>%
      hiseqr::bar_plot(x = "pct", y = "sample", label = NA,
                       group = "group", direction = "horizontal",
                       fill = "group") +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank())
    # ## output
    list(i7 = p1, bc = p2, p7 = p3,
         i7_table = df2, bc_table = df3, p7_table = df4)
  })
}






#' @describeIn Generate the report for hiseq_p7
hiseq_p7_report <- function(input, output) {
  input  <- normalizePath(input)
  if(! dir.exists(output)) {
    dir.create(output, recursive = TRUE)
  }
  output  <- normalizePath(output)
  outhtml <- file.path(output, "HiSeq_P7_report.html")
  template <- system.file("qc", "hiseq_qc_p7.Rmd", package = "hiseqr")
  # check files, template
  if(! file.exists(template)) {
    stop(paste0("template not exists: ", template))
  }
  # check stat.log files
  f_list <- list.files(input, "*stat.log$", full.names = TRUE)
  if(length(f_list) == 0) {
    stop("No *stat.log files detected")
  }
  ## copy template to output
  template_to <- file.path(output, basename(template))
  file.copy(template, template_to, overwrite = TRUE)
  rmarkdown::render(input       = template_to,
                    output_file = outhtml,
                    params      = list(input_dir = input))
}




























