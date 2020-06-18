




##----------------------------##
## DT table
get_DT_table <- function(df, mode = 1) {

  if(mode == 1) {
    DT::datatable(df,
                  extensions = 'Buttons',
                  options = list(
                    pageLength = 10,
                    scrollX = TRUE,
                    dom = 'Bfrtip',
                    buttons =
                      list('copy', #'print',
                           list(extend = 'collection',
                                buttons = c('excel', 'csv'),
                                text = 'Download')))
    )
  } else if(mode == 2) {
    DT::datatable(df, rownames = TRUE, filter = "top",
                  options = list(pageLength = 10, scrollX = TRUE))
  } else if(mode == 3){
    DT::datatable(df, rownames = TRUE, filter = "top", escape = FALSE,
                  options = list(pageLength = 10, scrollX = TRUE))
  } else {
    DT::datatable(df)
  }
}


#' constructure list for file/directory
#' @param x path to a directory/file
#' @param recursive boolean walk in directories
#'
#' @export
file_to_list <- function(x, recursive = FALSE) {
  sapply(x, function(i){
    file_to_list_single(i, recursive)
  }, USE.NAMES = FALSE)
}


#' constructure list for file/directory
#' @param x path to a directory, single object
#' @param recursive boolean walk in directories
#'
#' @export
file_to_list_single <- function(x, recursive = FALSE) {
  # input
  if(length(x) > 1) {
    message("Multiple elements found, only the first one selected")
    x <- x[1]
  }

  if(dir.exists(x) & recursive) {
    # a1 <- split(unname(x), basename(x))
    x_list <- list.files(x, full.names = TRUE)
    a2 <- sapply(x_list, function(i) {
      # split(unname(i), basename(i))
      file_to_list_single(i, recursive)
    }, USE.NAMES = FALSE)
    # set name
    setNames(list(a2), nm = basename(x))
  } else {
    split(unname(x), basename(x))
  }
}


#' fq names
#'
#' remove .fa, .fa.gz, .fq, .fq.gz
#' remove, "_1" or "_2" suffix
#' remove, _rep1, _rep2, ...
#'
#' @param x character, name of fastq file
#' @param fix_pe bool, whether trim _1, _2 suffix, default: FALSE
#'
#' @export
fq_name <- function(x, fix_pe = TRUE, rm_rep = FALSE){
  name <- gsub("[.](fast|f)[aq](.gz)?$",
               "",
               basename(x),
               ignore.case = TRUE,
               perl = TRUE)

  ## for paired names: _1, _2
  if(isTRUE(fix_pe)) {
    name <- gsub("[._][12]$", "", name, perl = TRUE)
  }

  ## for replicates
  if(isTRUE(rm_rep)) {
    name <- gsub("[._](rep|r)[0-9]+$",
                 "",
                 name,
                 ignore.case = TRUE,
                 perl = TRUE)
  }

  return(name)
}




#' read_fc
#' read featureCount output, .txt
#'
#' @param x path to count.txt file, featureCounts output
#' @param digits integes, Number of digits keep, default: 0
#'
#' @import dplyr
#'
#' @export
read_fc <- function(x, digits = 0, gene_length = FALSE){
  if(isTRUE(gene_length)) {
    # return the length of genes, column
    df <- readr::read_delim(x[1], "\t", col_types = readr::cols(), comment = "#") %>%
      dplyr::select(Geneid, Length) %>%
      tibble::column_to_rownames("Geneid")

  } else {
    tmp <- lapply(x, function(f){
      dfx <- readr::read_delim(f, "\t", col_types = readr::cols(), comment = "#") %>%
        dplyr::select(-c(2:6)) %>%
        dplyr::rename(id = Geneid) %>%
        tidyr::gather("sample", "count", -1) %>%
        dplyr::mutate(count = round(count, digits = digits)) %>%
        tidyr::spread("sample", "count")

      ## basename
      colnames(dfx) <- gsub(".bam$", "", basename(colnames(dfx)))

      return(dfx)
    })

    # combine col
    df <- merge_list_of_dataframes(tmp, by = "id")
  }
  return(df)
}



#' fc_summary
#'
#' @param x path to the summary of featureCounts file
#'
#' @export
#'
fc_summary <- function(x, fix_names = TRUE) {
  df <- readr::read_delim(x, "\t", col_types = readr::cols()) %>%
    tidyr::gather(key = "sample", value = "count", -1)
  # fix the names of bam
  if(isTRUE(fix_names)){
    df$sample = gsub(".bam$", "", basename(df$sample))
  }
  return(df)
}



#' bed_venn
#'
#' @param blist list, a list of BED files
#' @param names character, the names for each group
#'
#' @export
bed_venn <- function(blist, names = NULL){
  if(length(blist) == 2){
    # x <- bedIntersect2(blist)
    x <- intersect_bed2(blist)
  } else if(length(blist) == 3) {
    # x <- bedIntersect3(blist)
    x <- intersect_bed3(blist)
  } else if(length(blist) == 4) {
    # x <- bedIntersect4(blist)
    x <- intersect_bed4(blist)
  } else {
    stop("only accept narrowpeaks: 2-4 files")
  }

  p <- vennplot(x, names)

  return(p)
}



#' read_align1
#'
#' @param x path to stat file
#' @import readr
#' @import dplyr
#'
#' @export
read_align1 <- function(x){
  if(! file.exists(x)) {
    warning(paste0("file not exists: ", x))
    # return(data.frame(id = character(0),
    #                   group = character(0),
    #                   count = numeric(0)))
    df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)),
                   c("id", "group", "count"))
    return(df)
  }
  # hd <- c("total", "unmap", "unique", "multi", "map", "id", "index")
  df <- readr::read_delim(x, "\t", col_names = TRUE,
                          col_types = readr::cols()) %>%
    dplyr::filter(grepl("^\\d+$", map)) %>%
    dplyr::rename(total = `#total`,
                  id    = fqname,
                  index = index_name)

  ## rRNA/chrM
  df_chrM <- df %>%
    dplyr::filter(grepl("^1\\.", index)) %>%
    dplyr::select(id, total, unique, multiple)
  df_genome <- df %>%
    dplyr::filter(grepl("^2\\.", index)) %>%
    dplyr::select(id, unique, multiple, unmap)
  df_out <- merge(df_chrM, df_genome, by = "id")
  names(df_out) <- c("id", "total", "mito.u",
                     "mito.m", "genome.u", "genome.m",
                     "unmap")
  # plot
  groups <- c("mito.u", "mito.m", "genome.u", "genome.m", "unmap")

  # format
  df2 <- df_out %>%
    dplyr::mutate(id = factor(id, levels = rev(id))) %>%
    dplyr::select(-total) %>%
    tidyr::gather("group", "count", -1) %>%
    dplyr::mutate(group = factor(group, levels = rev(groups)),
                  id    = as.character(id),
                  count = as.numeric(count))

  return(df2)
}


#' read_align2
#'
#' @param x path to stat file
#' @import readr
#' @import dplyr
#'
#' @export
read_align2 <- function(x){
  df <- readr::read_delim(x, ",", col_names = TRUE,
                          col_types = readr::cols()) %>%
    dplyr::select(1, 2, 5, 6, 9:11)
  colnames(df) <- c("id", "total", "mito.u", "mito.m",
                    "genome.u", "genome.m", "unmap")
  groups <- c("mito.u", "mito.m", "genome.u", "genome.m", "unmap")

  # convert
  df2 <- df %>%
    dplyr::mutate(id = factor(id, levels = rev(id))) %>%
    dplyr::select(-total) %>%
    tidyr::gather("group", "count", -1) %>%
    mutate(group = factor(group, levels = rev(groups)))
  return(df)
}




#' readStat
#'
#' @param flist list, multiple files of alignment stats
#' @import readr
#' @import dplyr
#'
#' @export
stat_align <- function(flist){
  # check format
  if(all(grepl("*.align.txt", flist))) {
    tmp <- lapply(flist, read_align1)
  } else if(all(grepl("*mapping_stat.csv", flist))) {
    tmp <- lapply(flist, read_align2)
  } else {
    stop("File type unknown: *.align.txt, or *mapping_stat.csv")
  }
  # df  <- dplyr::bind_rows(tmp) %>%
  #   mutate(mito.pct = (mito.u + mito.m) / total)
  df <- dplyr::bind_rows(tmp)
  return(df)
}




#' read_frag
#'
#' @param x path to fragment length
#' @import readr
#' @import dplyr
#' @import ggplot2
#'
#' @export
read_frag <- function(x){
  # tmp <- lapply(x, function(f){
  #   df <- read.delim(f, header = F, sep = "\t",
  #                    col.names = c("length", "count"))
  #   return(df)
  # })

  # # merge data.frame
  # df <- bind_rows(tmp)
  df <- read.delim(x[1], header = F, sep = "\t",
                   col.names = c("length", "count"))
  return(df)
}



#' corPlot
#'
#' @param x path to file, count matrix, bedtools computeMatrix output
#' @import readr
#' @import dplyr
#' @import RColorBrewer
#'
#' @export
cor_matrix <- function(x, id = "all"){
  df1  <- readr::read_delim(x, "\t", col_names = TRUE,
                            col_types = readr::cols()) %>%
    dplyr::select(-c(1:3))
  # all samples
  all_samples <- colnames(df1)

  # choose sample
  if(id == "all") {
    df2 <- df1
  } else {
    df2 <- df1 %>%
      dplyr::select(contains(id))
  }

  # check
  if(ncol(df2) == 0){
    stop("samples not found")
  }

  # names
  colnames(df2) <- gsub("ATACseq_DaGal4X|'|.not_MT_trRNA|.map_dm6|_1|\\.1", "", colnames(df2))

  return(df2)
}


#' merge a list of data.frames
#'
#' @param list a list of data.frames
#' @param by column name for function merge()
#'
#' @export
merge_list_of_dataframes <- function(list, by) {
  stopifnot(is.list(list))
  if (length(list) == 1) {
    return(list[[1]])
  } else if (length(list) == 2) {
    df <- merge(list[[1]], list[[2]], by = by, all = TRUE)
    df[is.na(df)] <- 0
    return(df)
  } else {
    df <- merge(list[[1]], list[[2]], by = by, all = TRUE)
    df[is.na(df)] <- 0
    # drop 1, 2
    list <- rlist::list.remove(list, c(1, 2))
    list_new <- rlist::list.insert(list, 1, df)
    return(merge_list_of_dataframes(list_new, by = by))
  }
}




## default colors in ggplot2
#' ggplot2_colors
#'
#' @param n integer
#'
#' @import scales
#'
#' @export
gg_color <- function(n = 3) {
  # version 1
  # n = 3
  # scales::hue_pal()(3)
  #
  # scales::show_col(scales::hue_pal()(3))

  # version2
  # https://stackoverflow.com/a/8197703/2530783
  #
  # gg_color_hue <- function(n) {
  #   hues = seq(15, 375, length = n + 1)
  #   hcl(h = hues, l = 65, c = 100)[1:n]
  # }
  #
  # scales::show_col(gg_color_hue(3))

  # version3
  # https://stackoverflow.com/a/8197706/2530783
  # ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  #   if ((diff(h) %% 360) < 1) {
  #     h[2] <- h[2] - 360/n
  #   }
  #   hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
  # }
  #
  # scales:::show_col(ggplotColours(n=3))

  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#' choose organism
#'
#'
#' @export
choose_orgdb <- function(organism) {
  switch(
    organism,
    "dm3" = org.Dm.eg.db::org.Dm.eg.db,
    "dm6" = org.Dm.eg.db::org.Dm.eg.db,
    "hs19" = org.Hs.eg.db::org.Hs.eg.db,
    "hg38" = org.Hs.eg.db::org.Hs.eg.db,
    "mm9" = org.Mm.eg.db::org.Mm.eg.db,
    "mm9" = org.Mm.eg.db::org.Mm.eg.db,
    "*" = NULL
  )
}











#'
#' ##------------------##
#' ## deprecated,
#' ## duplicated list levels
#'
#' #' constructure list for file/directory
#' #' @param x path to a directory/file
#' #' @param recursive boolean walk in directories
#' #'
#' #' @export
#' file_to_list <- function(x, recursive = FALSE) {
#'   sapply(x, function(i) {
#'     file_to_list_single(i, recursive)}, USE.NAMES = FALSE)
#' }
#'
#' #' constructure list for file/directory
#' #' @param x path to a directory, single object
#' #' @param recursive boolean walk in directories
#' #'
#' #' @export
#' file_to_list_single <- function(x, recursive = FALSE) {
#'   # input
#'   if(length(x) > 1) {
#'     warning("Multiple elements found, only the first one selected")
#'     x <- x[1]
#'   }
#'
#'   # file/dir to list
#'   if(dir.exists(x) & isTRUE(recursive)) {
#'     # all files in directory
#'     if(isTRUE(recursive)) {
#'       f <- list.files(x)
#'     } else {
#'       f <- x
#'     }
#'     setNames(object = list(
#'       sapply(f, function(i) {
#'         file_to_list(file.path(x, i), recursive)
#'       }, simplify = FALSE, USE.NAMES = TRUE)),
#'       nm = basename(x))
#'   } else if(file.exists(x)) {
#'     # single file, set filename as names
#'     split(unname(x), basename(x))
#'   } else {
#'     warning("file/directory not exists")
#'   }
#' }
