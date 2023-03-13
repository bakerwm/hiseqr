#' Functions for check P7 structure
#'
#' Check the structure of P7, TruSeq, Nextera



#' hiseq_p7_report
#' @export
hiseq_p7_report <- function(input, output=NULL) {
  input  <- normalizePath(input)
  if(is.null(output)) output <- input
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
  f_list <- list.files(input, "*.config.yaml$", full.names = TRUE)
  if(length(f_list) == 0) {
    stop("No *.config.yaml files detected")
  }
  ## copy template to output
  template_to <- file.path(output, basename(template))
  file.copy(template, template_to, overwrite = TRUE)
  rmarkdown::render(input       = template_to,
                    output_file = outhtml,
                    params      = list(input_dir = input))
}





#' hiseq_p7
#' @param x character path to the directory of hiseq p7 output
#'
#' @export
hiseq_p7 <- function(x) {
  if(! dir.exists(x)) {
    on.exit("Not a directory")
  }
  #-- 1. list config.yaml
  config_list <- list.files(x, "*.config.yaml", full.names = TRUE)
  if(length(config_list) == 0) {
    on.exit(glue::glue("No *.config.yaml found in {x}"))
  }
  #-- 2. list. i7, barcode, p7 .json
  lapply(config_list, function(i) {
    #-- 2.1 load config
    pd <- configr::read.config(i)
    #-- 2.2 hiseq_type
    st <- ifelse(pd$is_nextera, "Nextera",
                 ifelse(pd$is_truseq, "TruSeq",
                        ifelse(pd$is_smallrna, "smallRNA", "unknown")))
    if(pd$is_nsr) {
      st <- paste0(st, "-NSR")
    }
    #-- 2.3 json data
    p1 <- hiseq_lib_p7(pd$p7_json)[[1]] +
      ggtitle(paste0(pd$fname, ": ", st))
    p2 <- hiseq_lib_i7(pd$i7_json)[[1]]
    if(is(p2, "ggplot")) {
      p2 <- p2 + ggtitle("i7 index")
    }
    p3 <- hiseq_lib_barcode(pd$barcode_json)[[1]] +
      ggtitle("barcode")
    #-- 2.4 output
    list(
      smp_name   = pd$fname,
      hiseq_type = st,
      p7         = p1,
      i7         = p2,
      barcode    = p3
    )
  })
}






#' hiseq_lib_i7
#'
hiseq_lib_i7 <- function(x) {
  lapply(x, function(f) {
    fname <- gsub(".i7.json", "", basename(f))
    #-- 1. load data
    l <- jsonlite::read_json(f)
    if(length(l) == 0) {
      return(NULL)
    }
    df <- lapply(l, as.data.frame.list) %>%
      dplyr::bind_rows() %>%
      dplyr::mutate(sample = fname) %>%
      dplyr::arrange(desc(pct)) %>%
      dplyr::mutate(rank = dplyr::row_number())
    #-- 2. add others
    df_other <- data.frame(
      count = round(mean(df$count / df$pct * 100), 0) - sum(df$count),
      name  = "other",
      pct   = 100 - sum(df$pct),
      seq   = "-",
      seq_revcomp = "-",
      sample = fname,
      rank   = nrow(df) + 1
    )
    #-- 3. combine
    df2 <- dplyr::bind_rows(df, df_other) %>%
      dplyr::mutate(rank = as.character(rank))
    #-- 4. barplot
    df2 %>%
      hiseqr::bar_plot(
        x = "pct", y = "rank", label = "name",
        group = "name", direction = "horizontal") +
      geom_vline(xintercept = c(50, 100), color = "blue", linetype = 2) +
      scale_x_continuous(limits = c(0, 150),
                         breaks = c(0, 50, 100),
                         position = "top") +
      ggtitle(fname) +
      theme(axis.title.x = element_blank(),
            legend.position = "None")
  })
}



#' hiseq_lib_barcode
#'
hiseq_lib_barcode <- function(x) {
  lapply(x, function(f) {
    fname <- gsub(".barcode.json", "", basename(f))
    #-- 1. load data
    l <- jsonlite::read_json(f)
    df <- lapply(l, as.data.frame.list) %>%
      dplyr::bind_rows() %>%
      dplyr::mutate(sample = fname) %>%
      dplyr::arrange(desc(pct)) %>%
      dplyr::mutate(rank = dplyr::row_number())
    #-- 2. add others
    df_other <- data.frame(
      count = round(mean(df$count / df$pct * 100), 0) - sum(df$count),
      name  = "other",
      pct   = 100 - sum(df$pct),
      seq   = "-",
      seq_revcomp = "-",
      sample = fname,
      rank   = nrow(df) + 1
    )
    #-- 3. combine
    df2 <- dplyr::bind_rows(df, df_other) %>%
      dplyr::mutate(rank = as.character(rank))
    #-- 4. barplot
    df2 %>%
      hiseqr::bar_plot(
        x = "pct", y = "rank", label = "name",
        group = "name", direction = "horizontal") +
      geom_vline(xintercept = c(50, 100), color = "blue", linetype = 2) +
      scale_x_continuous(limits = c(0, 150),
                         breaks = c(0, 50, 100),
                         position = "top") +
      ggtitle(fname) +
      theme(axis.title.x = element_blank(),
            legend.position = "None")
  })
}




#' hiseq_lib_p7
#'
hiseq_lib_p7 <- function(x) {
  lapply(x, function(f) {
    fname <- gsub(".p7.json", "", basename(f))
    #-- 1. load data
    df <- jsonlite::read_json(f) %>%
      as.data.frame.list() %>%
      dplyr::mutate(sample = fname,
                    p7     = "p7") %>%
      tidyr::pivot_longer(no_ad:p7b,
                          names_to = "group",
                          values_to = "count") %>%
      dplyr::group_by(sample) %>%
      dplyr::mutate(pct = round(count / sum(count) * 100, 1)) %>%
      dplyr::mutate(group = forcats::fct_relevel(group, "p7a", "p7a_p7b", "p7b", "no_ad"))
    #-- 2. barplot
    df %>%
      hiseqr::bar_plot(
        x = "pct", y = "p7", label = NA,
        group = "group", direction = "horizontal") +
      geom_vline(xintercept = c(50, 100), color = "blue", linetype = 2) +
      scale_x_continuous(limits = c(0, 150),
                         breaks = c(0, 50, 100),
                         position = "top") +
      ggtitle(fname) +
      theme(axis.title.x = element_blank(),
            legend.position = "top")
  })
}
