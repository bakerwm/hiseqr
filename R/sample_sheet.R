# Example:



#' read sample info from Excel.xlsx
#'
#' @param x file xlsx
#' @n_max int max number rows to read
#' @autofix bool auto fix the samplename, etc, default: TRUE
#'
#' @export
sheet_load <- function(x, n_max=10000, autofix = TRUE, index_list = NULL,
                       sampleid_start = 1) {
  df <- readxl::read_xlsx(x, sheet = "sample_sheet", trim_ws = TRUE, n_max = n_max)
  colnames(df) <- gsub("[^\\w]", "", colnames(df), perl = TRUE)
  colnames(df) <- tolower(colnames(df))

  # add date
  f_date <- stringr::str_split(basename(x), "_", simplify = TRUE)[1, 2] #date
  df$date <- as.Date(f_date, format="%Y%m%d")

  # remove blank rows
  df <- df %>%
    dplyr::filter(! is.na(lib_number))

  # auto_fix
  if(isTRUE(autofix)) {
    df <- sheet_autofix(df)
  }

  # convert index from id to sequence
  index <- sheet_read_index(index_list)

  if(is.null(index)) {
    index <- NULL
  } else {
    p7 <- index[df$p7_index_id]
    bc <- index[df$barcode_id]
    # convert NULL to NA
    p7[is.na(p7)] <- "NULL"
    bc[is.na(bc)] <- "NULL"
    # assign
    df$p7_index    <- p7
    df$barcode_seq <- bc
  }

  # output
  df
}


#' check samplename
#' < 60 chars
#' no duplicates
#'
#' @export
is_sheet_samplename <- function(x) {
  # not longer than: 60 chars
  #
  all(grepl("^[\\w|\\-|\\.|]{1,60}$", x, perl = TRUE, ignore.case = TRUE))
}


#' check duplication
#'
#' @export
is_sheet_duplicate <- function(x, y = NULL) {
  if(is.null(y)) {
    length(x) > length(unique(x))
  } else {
    if(length(x) == length(y)) {
      xy <- paste0(x, y)
      length(xy) > length(unique(xy))
    } else {
      rep(FALSE, length(x)) # all FALSE
    }
  }
}


#' check libname
#' @export
is_sheet_libname <- function(x) {
  # YY00, YS00
  all(grepl("^\\w{2,4}\\d{2,3}$", x, perl = TRUE))
}


#' check index id
#' @export
is_sheet_index_id <- function(x) {
  # TruSeq_index1-48
  # Next_Ad2.1-24
  # Null
  # no duplicate names: p7_index + barcode
  all(grepl("^(truseq_index\\d+)|(next_ad2.\\d+)|(null)$", df$p7_index_id, perl = TRUE, ignore.case = TRUE))

}


#' check barcode id
#' @export
is_sheet_barcode_id <- function(x) {
  all(grepl("^(p7_\\d+A|B)|(iclip\\d+)|(null)$", df$barcode_id, perl = TRUE, ignore.case = TRUE))
}


#' check index name/ barcode name
#' @param x string, id of the index
#' @param index string/vector, named vector,
#'
#' @export
is_sheet_valid_index <- function(x, index = NULL) {
  if(is.null(index)) {
    TRUE
  } else {
    x <- x[!grepl("NULL", x, ignore.case = TRUE)] # remove NULL
    x <- x[! is.na(x)]
    all(x %in% names(index))
  }
}


#' read index sequence
#' @param x string path to the index file, could be *.rds, *.txt, *.csv
#'
#' @export
sheet_read_index <- function(x) {
  # check index
  if(! is.null(x)) {
    if(endsWith(x, ".rds")) {
      l <- readRDS(x) # named sequence
      # format: truseq.TruSeq_Index1
      a <- unlist(l)
      # format: TruSeq_Index1
      b <- stringr::str_split(names(a), "\\.", n = 2, simplify = T) #
      names(a) <- b[, 2]
      a
    } else {
      if(endsWith(x, ".txt")) {
        di <- readr::read_delim(x, "\t", col_types = readr::cols())
      } else if(endsWith(x, "*.csv")) {
        di <- readr::read_csv(x, col_types = readr::cols())
      } else {
        di <- setNames(data.frame(matrix(ncol = 2, nrow = 0)),
                       c("name", "sequence"))
      }
      # name, sequence required
      if(all(c("name", "sequence") %in% colnames(di))) {
        setNames(di$sequence, di$name)
      } else {
        warning("name, sequence, columns not found")
        NULL
      }
    }
  }
}


#' Check the fields in sheet
#' @param data.frame sample data
#' @param index_list list, named index sequence
#'
#' @export
sheet_check <- function(df, index_list = NULL) {
  # required columns
  required_cols <- c("lib_number", "lib_user", "sample_name", "rbp",
                     "cell_line", "species", 'spikein', "p7_index_id",
                     "barcode_id", "seq_type", "lib_type", "readsm")
  chk1 <- all(required_cols %in% colnames(df))

  # load index
  index <- sheet_read_index(index_list)

  # chk
  if(chk1) {
    chk2 <- is_sheet_libname(df$lib_number)
    chk3 <- is_sheet_libname(df$lib_user)
    chk4 <- is_sheet_samplename(df$sample_name)
    chk5 <- is_sheet_index_id(df$p7_index_id)
    chk6 <- is_sheet_barcode_id(df$barcode_id)
    ## check duplication
    chk7 <- is_sheet_duplicate(df$sample_name)
    chk8 <- is_sheet_duplicate(df$p7_index_id, df$barcode_id)
    chk9 <- is_sheet_valid_index(df$p7_index_id, index)
    chk10 <- is_sheet_valid_index(df$barcode_id, index)
  } else {
    chk2 <- chk3 <- chk4 <- chk5 <- chk6 <- chk7 <- chk8 <- TRUE
  }

  ## message
  df_msg <- data.frame(
    "content" = c("required columns",
                  "Lib_number",
                  "Lib_user",
                  "Sample_name",
                  "P7_index_id",
                  "Barcode_id"),
    "status"  = c(chk1 & !chk7,
                  chk2,
                  chk3,
                  chk4,
                  chk5 & !chk8 & chk9,
                  chk6 & !chk8 & chk10),
    "Expected" = c("see examples above",
                   "eg: YY11",
                   "eg: MW03",
                   "no duplication, shorter than 60",
                   "TruSeq_index1-48/Next_Ad2.1-24",
                   "eg: P7_1A, iCLIP2"),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(status = ifelse(status, "ok", "failed"))

  ## output
  list(status = all(df_msg$status == "ok"),
       msg    = df_msg)
}



#' Auto fix sample_sheet
#'
#' replace, "[^\\w\\-\\.] by "_"
#'
#' lib_number: [A-Z]{2}\\d{2}
#' lib_user:   [A-Z]{2,4}\\d{2}
#' sample_name: [\\w\\-\\.]+
#' p7_index_id:
#' barcode_id:
#'
#' @param df data.frame
#'
#' @export
sheet_autofix <- function(df) {
  # auto fix columns

  ## no blanks
  df$lib_number  <- gsub("[^\\w\\-\\.]", "_", df$lib_number, per = TRUE)
  df$lib_user    <- gsub("[^\\w\\-\\.]", "_", df$lib_user, per = TRUE)
  df$sample_name <- gsub("[^\\w\\-\\.]", "_", df$sample_name, per = TRUE)
  df$p7_index_id <- gsub("[^\\w\\-\\.]", "_", df$p7_index_id, per = TRUE)
  df$barcode_id  <- gsub("[^\\w\\-\\.]", "_", df$barcode_id, per = TRUE)

  ## upper case
  df$lib_number  <- toupper(df$lib_number)
  df$lib_user    <- toupper(df$lib_user)
  df$seq_type    <- toupper(df$seq_type)
  df$p7_index_id <- gsub("null", "NULL", df$p7_index_id, per = TRUE, ignore.case = TRUE)
  df$barcode_id  <- gsub("null", "NULL", df$barcode_id, per = TRUE, ignore.case = TRUE)

  ## remove double __
  df$sample_name <- gsub("\\_+", "_", df$sample_name)

  ## seqname: DNAseq, RNAseq, ChIPseq, smallRNAseq, GoldCLIP, GROseq
  df$sample_name <- gsub("(RNA|DNA|ChIP|Gold|GRO)\\_?(seq)", "\\1\\2", df$sample_name, perl = TRUE, ignore.case = TRUE)
  df$sample_name <- gsub("RNAseq", "RNAseq", df$sample_name, ignore.case = TRUE)
  df$sample_name <- gsub("DNAseq", "DNAseq", df$sample_name, ignore.case = TRUE)
  df$sample_name <- gsub("ChIPseq", "ChIPseq", df$sample_name, ignore.case = TRUE)
  df$sample_name <- gsub("smallRNAseq", "smallRNAseq", df$sample_name, ignore.case = TRUE)
  df$sample_name <- gsub("small_RNAseq", "smallRNAseq", df$sample_name, ignore.case = TRUE)
  df$sample_name <- gsub("GROseq", "GROseq", df$sample_name, ignore.case = TRUE)

  ## rep1:
  df$sample_name <- gsub("(r|rep)(\\d+)", "rep\\2", df$sample_name, perl = TRUE)

  df
}



