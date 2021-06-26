#' Example:



#' read sample info from Excel.xlsx
#'
#' @param x file xlsx
#' @n_max int max number rows to read
#' @autofix bool auto fix the samplename, etc, default: TRUE
#'
#' @export
read_sheet <- function(x, n_max=10000, autofix = TRUE,
                       sampleid_start = 1) {
  df <- tryCatch(
    {
      df <- readxl::read_xlsx(x, sheet = "sample_sheet",
                              trim_ws = TRUE, n_max = n_max)
      colnames(df) <- gsub("[^\\w]", "", colnames(df), perl = TRUE)
      colnames(df) <- tolower(colnames(df))
      df <- dplyr::filter(df, ! is.na(lib_number))
    },
    error=function(cond) {
      message(glue::glue("Failed reading sheet: {x}"))
      message(cond)
      # Choose a return value in case of error
      return(NA)
    }
  )
  # retrieve date from filename
  # YY01_20170114_GJQ_...
  f_date <- stringr::str_extract(basename(x), "20[\\d]{6}")
  if(is.na(f_date)) {
    f_date <- Sys.Date()
  } else {
    f_date <- as.Date(f_date, format = "%Y%m%d")
  }
  # add date
  df$date <- f_date
  # fix lane
  df$lib_number  <- stringr::str_extract(basename(x), "^Y[A-Z]\\d+")
  # auto_fix
  if(isTRUE(autofix)) {
    df <- autofix_sheet(df, sampleid_start)
  }
  df
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
autofix_sheet <- function(df, sampleid_start = 1) {
  if(is_valid_sheet_df(df)) {
    message("Auto fix sample sheet")
  } else {
    message("autofix_sheet() failed, expect data.frame")
    return(NA)
  }
  required_cols <- c("sampleid", "lib_number", "lib_sub",  "lib_user",
                     "sample_name", "rbp",
                     "cell_line", "species", 'spikein', "p7_index_id",
                     "barcode_id", "seq_type", "lib_type", "readsm",
                     "fcid", "lane", "reference", "spikein_ref", "p7_index",
                     "barcode_seq", "reminder", "date")
  df1 <- df %>%
    dplyr::mutate(
      lib_number  = gsub("[\\W]", "_", lib_number, perl = TRUE),
      lib_number  = toupper(lib_number),
      lib_user    = .sanitize_str(lib_user),
      lib_user    = toupper(lib_user),
      sample_name = .sanitize_str(sample_name),
      sample_name = autofix_hiseq_name(sample_name),
      p7_index_id = .sanitize_str(p7_index_id),
      p7_index_id = gsub("null", "NULL", p7_index_id, ignore.case = TRUE),
      barcode_id  = .sanitize_str(barcode_id),
      barcode_id  = gsub("null", "NULL", barcode_id, ignore.case = TRUE),
      seq_type    = toupper(seq_type), #SE|PE
      lib_type    = .sanitize_str(lib_type),
      readsm      = ifelse(is.na(readsm), 3, readsm)
    )
  # add lib_sub
  if(! "lib_sub" %in% colnames(df)) {
    df$lib_sub <- "G1"
  }
  # fix sampleid
  df <- df %>%
    dplyr::mutate(
      s = stringr::str_pad(
        as.character(sampleid_start + row_number() - 1),
        side = "left", width = 6, pad = "0"
      ),
      sampleid = paste0("YYs", s)
    )
  # output
  dplyr::select(df, all_of(required_cols))
}





#' @describeIn is_valid_sheet_df
#'
#' Check the fields in sheet
#' @param data.frame sample data
#' @param index_list list, named index sequence
#'
#' @export
is_valid_sheet_df <- function(df, verbose = FALSE) {
  # required columns
  required_cols <- c("lib_number", "lib_user", "sample_name", "rbp",
                     "cell_line", "species", 'spikein', "p7_index_id",
                     "barcode_id", "seq_type", "lib_type", "readsm")
  k1 <- all(required_cols %in% colnames(df))
  # chk
  if(k1) {
    k2 <- is_valid_index(df$p7_index_id, "p7") | grepl("NULL", df$p7_index_id, ignore.case = TRUE)
    k3 <- is_valid_index(df$barcode_id, "barcode") | grepl("NULL", df$barcode_id, ignore.case = TRUE)
    # check duplication
    k4 <- ! any(duplicated(df$sample_name))
    k5 <- ! any(duplicated(paste(df$p7_index_id, df$barcode_id)))
  } else {
    chk2 <- chk3 <- chk4 <- chk5 <- TRUE
  }
  # message
  msg <- data.frame(
    Content = c("required columns", "Sample_name", "P7_index_id", "Barcode_id"),
    Status  = c(k1, k4, all(k2), all(k3)),
    Expected = c("see examples above", "no dup", "TruSeq18", "P7_1A"),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(Status = ifelse(Status, "ok", "failed")) %>%
    dplyr::select(Status, Content, Expected)
  if(isTRUE(verbose)) {
    print(msg)
  }
  all(c(k1, k2, k3, k4, k5))
}


#' replace non-character by "_"
.sanitize_str <- function(x, replace = "_") {
  sapply(x, function(i) {
    # gsub("[^\\w\\.\\-]", replace,  i)
    i <- gsub("[^A-Za-z0-9.-_]", replace,  i)
    gsub("(_)+", "_", i, perl = TRUE)
  })
}


#' lib_number, auto_fix


#' eg: YY01, WM03
#'
#' @param x character lib_user
#'
#' @export
autofix_lib_user <- function(x) {
  i <- toupper(x)
  gsub("([A-Z]{2,4})([\\W])?(\\d+)", "\\1\\3", i)
}



#' mission:
#' 1. replace non-letters, by "_"
#' 2. fix hiseq_type
#' 3. add rep1/rep2 in the tail
#'
#' @param x character
#'
#' @export
autofix_sample_name <- function(x) {
  sapply(x, function(i) {
    j <- .sanitize_str(i)
    j <- autofix_hiseq_name(j)
    if(grepl("_(r|rep)(\\d+)$", j, ignore.case = TRUE)) {
      gsub("_(r|rep)(\\d+)$", "_rep\\2", j)
    } else {
      paste0(j, "_rep1")
    }
  })
}



#' in the filename, HiSeq located in the first part
#'
#'
#' Standard hiseq names in this package
#'
#' RNAseq
#' ChIPseq
#' ATACseq
#' DNAseq
#' CnR (CUT&RUN)
#' CnT (CUT&TAG)
#' smRNAseq (small RNAseq)
#' STACCseq
#' HiC
#' RiboSeq
#' GoldCLIP
#' CLIPseq
#' RIPseq
#' GridSeq
#'
#'
#' @param x character the filename or lib_type
#'
#' @export
autofix_hiseq_name <- function(x) {
  j <- gsub("(ATAC|ChIP|DNA|RNA|GoldCLIP|GRO|STACC)(_seq)?", "\\1seq", i,
            perl = TRUE, ignore.case = TRUE)
  j <- gsub("^mRNAseq", "RNAseq", j, ignore.case = TRUE)
  j <- gsub("DNAseq", "DNAseq", j, ignore.case = TRUE)
  j <- gsub("ATACseq", "ATACseq", j, ignore.case = TRUE)
  j <- gsub("ChIPseq", "ChIPseq", j, ignore.case = TRUE)
  j <- gsub("CLIP(_)?(seq)?(_)?", "ChIPseq_", j, ignore.case = TRUE)
  j <- gsub("(small|sm)(_)?RNAseq", "smRNAseq", j, ignore.case = TRUE)
  j <- gsub("(CUTRUN|CUT_RUN|CUT_and_RUN)", "CnR", j, ignore.case = TRUE)
  j <- gsub("(CUTTAG|CUT_TAG|CUT_and_TAG)", "CnT", j, ignore.case = TRUE)
  j <- gsub("STACCseq", "STACCseq", j, ignore.case = TRUE)
  j <- gsub("GROseq", "GROseq", j, ignore.case = TRUE)
  j
}




#' load index,
#'
#' @export
load_hiseq_index <- function(hiseq_type = TRUE) {
  f <- system.file("data", "hiseq_index.rds", package = "hiseqr")
  l <- readRDS(f)
  if(isTRUE(hiseq_type)) {
    hiseq_type <- names(l)
  }
  l <- l[hiseq_type]
  ## output
  # lapply(hiseq_type, function(i) {
  #   l[[i]] %>%
  #     as.data.frame() %>%
  #     dplyr::rename(index = ".") %>%
  #     tibble::rownames_to_column("name")
  # }) %>%
  #   dplyr::bind_rows()
  ## remove names for 1st level
  names(l) <- NULL
  unlist(l)
}



#' check index id
#' @param x character
#' @param hiseq_type character, p7, barcode
#'
#'
#' @export
is_valid_index <- function(x, hiseq_type = "p7", skip_null = TRUE) {
  if(hiseq_type %in% c("p7", "truseq", "nextera", "barcode")) {
    if(hiseq_type == "p7") {
      hiseq_type <- c("truseq", "nextera")
    }
  } else {
    warning(glue::glue("unknown hiseq type: {hiseq_type}, expect: p7, truseq, nextera, barcode"))
    return(FALSE)
  }
  # p7 index: TruSeq, Nextera
  p7 <- load_hiseq_index(hiseq_type)
  sapply(x, function(i) {
    if(isTRUE(skip_null) & tolower(i) == "null") {
      out <- TRUE
    } else {
      t <- ifelse(i %in% p7, names(p7)[p7 == i],
                  ifelse(i %in% names(p7), p7[i], NA))
      if(is.na(t)) {
        message(glue::glue("unknown index: {i}"))
      }
      out <- ! is.na(t)
    }
    out
  })
}







#' #' read index sequence
#' #' @param x string path to the index file, could be *.rds, *.txt, *.csv
#' #'
#' #' @export
#' sheet_read_index <- function(x) {
#'   # check index
#'   if(is.null(x)) {
#'     return(NULL)
#'   }
#'
#'   if(is_sheet_valid_index(names(x))) {
#'     return(x)
#'   } else {
#'     x <- x[1] # the first one
#'     if(endsWith(x, ".rds")) {
#'       l <- readRDS(x) # named sequence
#'       # format: truseq.TruSeq_Index1
#'       a <- unlist(l)
#'       # format: TruSeq_Index1
#'       b <- stringr::str_split(names(a), "\\.", n = 2, simplify = T) #
#'       names(a) <- b[, 2]
#'       a
#'     } else {
#'       if(endsWith(x, ".txt")) {
#'         di <- readr::read_delim(x, "\t", col_types = readr::cols())
#'       } else if(endsWith(x, "*.csv")) {
#'         di <- readr::read_csv(x, col_types = readr::cols())
#'       } else {
#'         di <- setNames(data.frame(matrix(ncol = 2, nrow = 0)),
#'                        c("name", "sequence"))
#'       }
#'       # name, sequence required
#'       if(all(c("name", "sequence") %in% colnames(di))) {
#'         setNames(di$sequence, di$name)
#'       } else {
#'         warning("name, sequence, columns not found")
#'         NULL
#'       }
#'     }
#'   }
#' }



#'
#' #' check duplication
#' #'
#' #' @export
#' is_sheet_duplicate <- function(x, y = NULL) {
#'   if(is.null(y)) {
#'     length(x) > length(unique(x))
#'   } else {
#'     if(length(x) == length(y)) {
#'       xy <- paste0(x, y)
#'       length(xy) > length(unique(xy))
#'     } else {
#'       rep(FALSE, length(x)) # all FALSE
#'     }
#'   }
#' }


#'
#' #' check index id
#' #' @param x string
#' #'
#' #' @export
#' is_sheet_index_id <- function(x) {
#'   # TruSeq_index1-48
#'   # Next_Ad2.1-24
#'   # Null
#'   # no duplicate names: p7_index + barcode
#'   all(grepl("^(truseq_index\\d+)|(next_ad2.\\d+)|(null)$", x, perl = TRUE, ignore.case = TRUE))
#' }


#'
#' #' check barcode id
#' #' @export
#' is_sheet_barcode_id <- function(x) {
#'   all(grepl("^(p7_\\d+A|B)|(iclip\\d+)|(null)$", x, perl = TRUE, ignore.case = TRUE))
#' }

#'
#' #' check index name/ barcode name
#' #' @param x string, id of the index
#' #' @param index string/vector, named vector,
#' #'
#' #' @export
#' is_sheet_valid_index <- function(x, index = NULL) {
#'   if(is.null(index)) {
#'     TRUE
#'   } else {
#'     x <- x[!grepl("NULL", x, ignore.case = TRUE)] # remove NULL
#'     x <- x[! is.na(x)]
#'     all(x %in% names(index))
#'   }
#' }






