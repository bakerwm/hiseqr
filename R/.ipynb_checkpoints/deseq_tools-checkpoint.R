
#' @param x string deseq_dir
get_align_stat3 <- function(x) {
  x_type <- hiseqr::is_hiseq_dir(x)

  if(! x_type == "deseq_single") {
    warning(paste0("deseq_single expected, ", x_type, " got"))
    return(NULL)
  }

  tmp <- lapply(c("gene", "te", "piRNA_cluster"), function(f){
    args <- hiseqr::deseq_single_dir(x, f)
    if(is.null(args)) {
      return(NULL)
    }

    tmp2 <- lapply(c(args$dirs_ctl, args$dirs_exp), function(i){
      get_align_stat2(i, f)
    })
    dplyr::bind_rows(tmp2)
  })

  ## total mapping
  df <- dplyr::bind_rows(tmp)

  ## total reads
  df2 <- df %>%
    dplyr::filter(index_name == "1.rRNA") %>%
    dplyr::select(fqname, total)

  ## mapped
  df_tmp <- df %>%
    dplyr::select(fqname, index_name, map) %>%
    tidyr::spread("index_name", "map")

  df_table <- merge(df2, df_tmp, by = "fqname") %>%
    dplyr::mutate(unmap = total * 2 - select_if(., is.numeric) %>% rowSums())
  df_table
}




#' @param x string deseq_dir
get_align_stat <- function(x, feature = "gene") {
  x_type <- hiseqr::is_hiseq_dir(x)

  if(! x_type == "deseq_single") {
    warning(paste0("deseq_single expected, ", x_type, " got"))
    return(NULL)
  }

  args <- hiseqr::deseq_single_dir(x)

  tmp <- lapply(c(args$dirs_ctl, args$dirs_exp), function(i){
    get_align_stat2(i, feature)
  })

  dplyr::bind_rows(tmp)
}



#' @param x string, rnaseq_single dir
get_align_stat2 <- function(x, feature = "gene") {
  x_type <- hiseqr::is_hiseq_dir(x)

  if(! x_type == "rnaseq_single") {
    warning(paste0("rnaseq_single expected, ", x_type, " got"))
    return(NULL)
  }
  # x is rnaseq_single
  align_dir <- file.path(x, feature, "align")
  stat_list <- list.files(align_dir, "*.json", T, T, T)

  tmp <- lapply(stat_list, function(f){
    d <- jsonlite::read_json(f)
    as_tibble(d)
  })

  dplyr::bind_rows(tmp)
}

