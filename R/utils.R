

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
fqName <- function(x, fix_pe = TRUE, rm_rep = FALSE){
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




#' fcReader
#' read featureCount output, .txt
#'
#' @param x path to count.txt file, featureCounts output
#' @param digits integes, Number of digits keep, default: 0
#'
#' @import dplyr
#'
#' @export
fcReader <- function(x, digits = 0){
  tmp <- lapply(x, function(f){
    df <- readr::read_delim(f, "\t", col_types = readr::cols(), comment = "#") %>%
      dplyr::select(-c(2:6)) %>%
      dplyr::rename(id = Geneid) %>%
      tidyr::gather("sample", "count", -1) %>%
      dplyr::mutate(count = round(count, digits = digits)) %>%
      tidyr::spread("sample", "count")

    ## basename
    colnames(df) <- gsub(".bam$", "", basename(colnames(df)))

    return(df)
  })

  # combine col
  df2 <- merge_list_of_dataframes(tmp, by = "id")
  return(df2)
}



#' fc_summary
#'
#' @param x path to the summary of featureCounts file
#'
#' @export
#'
fcSummary <- function(x, fix_names = TRUE) {
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
bedVenn <- function(blist, names = NULL){
  if(length(blist) == 2){
    x <- bedIntersect2(blist)
  } else if(length(blist) == 3) {
    x <- bedIntersect3(blist)
  } else if(length(blist) == 4) {
    x <- bedIntersect4(blist)
  } else {
    stop("only accept narrowpeaks: 2-4 files")
  }

  p <- vennplot(x, names)

  return(p)
}



#' readAlign1
#'
#' @param x path to stat file
#' @import readr
#' @import dplyr
#'
#' @export
readAlign1 <- function(x){
  if(! file.exists(x)) {
    warning(paste0("file not exists: ", x))
    # return(data.frame(id = character(0),
    #                   group = character(0),
    #                   count = numeric(0)))
    df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("id", "group", "count"))
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
    filter(grepl("^1\\.", index)) %>%
    select(id, total, unique, multiple)
  df_genome <- df %>%
    filter(grepl("^2\\.", index)) %>%
    select(id, unique, multiple, unmap)
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


#' Deprecated: column order dependent.
#'
#' #' readAlign1
#' #'
#' #' @param x path to stat file
#' #' @import readr
#' #' @import dplyr
#' #'
#' #' @export
#' readAlign1 <- function(x){
#'   hd <- c("total", "unmap", "unique", "multi", "map", "id", "index")
#'   df <- readr::read_delim(x, "\t", col_names = hd,
#'                           col_types = readr::cols(), comment = "#")
#'   df_chrM   <- df %>%
#'     # filter(grepl("chrM", index)) %>%
#'     filter(grepl("^1\\.", index)) %>%
#'     select(id, total, unique, multi)
#'   df_genome <- df %>%
#'     # filter(!grepl("chrM", index)) %>%
#'     filter(grepl("^2\\.", index)) %>%
#'     select(id, unique, multi, unmap)
#'   df_out <- merge(df_chrM, df_genome, by = "id")
#'   names(df_out) <- c("id", "total", "mito.u",
#'                      "mito.m", "genome.u", "genome.m",
#'                      "unmap")
#'   # plot
#'   groups <- c("mito.u", "mito.m", "genome.u", "genome.m", "unmap")
#'
#'   # format
#'   df2 <- df_out %>%
#'     mutate(id = factor(id, levels = rev(id))) %>%
#'     select(-total) %>%
#'     tidyr::gather("group", "count", -1) %>%
#'     mutate(group = factor(group, levels = rev(groups)))
#'
#'   return(df2)
#' }




#' readAlign2
#'
#' @param x path to stat file
#' @import readr
#' @import dplyr
#'
#' @export
readAlign2 <- function(x){
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
alignStat <- function(flist){
  # check format
  if(all(grepl("*.align.txt", flist))) {
    tmp <- lapply(flist, readAlign1)
  } else if(all(grepl("*mapping_stat.csv", flist))) {
    tmp <- lapply(flist, readAlign2)
  } else {
    stop("File type unknown: *.align.txt, or *mapping_stat.csv")
  }
  # df  <- dplyr::bind_rows(tmp) %>%
  #   mutate(mito.pct = (mito.u + mito.m) / total)
  df <- dplyr::bind_rows(tmp)
  return(df)
}




#' fragReader
#'
#' @param x path to fragment length
#' @import readr
#' @import dplyr
#' @import ggplot2
#'
#' @export
fragReader <- function(x){
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
corMatrix <- function(x, id = "all"){
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

