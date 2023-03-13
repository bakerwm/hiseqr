#' Example:

# f  <- list.files("/data/biodata/mydb/hiseq_sample_sheet/archived-2021/", "*.xlsx", full.names = TRUE)
# t  <- lapply(f, read_hiseq_sheet)
# df <- dplyr::bind_rows(t)
# df$sampleid <- fix_hiseq_sampleid(df$sampleid)
#
# df1 <- df %>%
#   dplyr::mutate(tmp_user = gsub("\\d+", "", lib_user)) %>%
#   tidyr::unite(sp, c(lib_number, tmp_user), sep = "-", remove = FALSE)
#
# # specific number
# df2 <- dplyr::select(df, lib_number, lib_user) %>%
#   dplyr::mutate(lib_user = gsub("\\d+", "", lib_user)) %>%
#   unique() %>%
#   dplyr::group_by(lib_user) %>%
#   dplyr::mutate(num = stringr::str_pad(row_number(), 2, pad = 0)) %>%
#   tidyr::unite(sp, c(lib_number, lib_user), sep = "-")
#
# # combine
# df3 <- merge(df1, df2, by = "sp") %>%
#   dplyr::mutate(lib_user = gsub("\\d+", "", lib_user),
#                 lib_user = paste0(lib_user, num)) %>%
#   dplyr::select(-c(sp, num, tmp_user))
#
# s <- "/data/biodata/mydb/hiseq_sample_sheet/summary_2021.xlsx"
# writexl::write_xlsx(df3, s, col_names = TRUE)


#' load hiseq sample sheet files
#'
#'
load_hiseq_sheet_dir <- function(x) {
  if(inherits(x, "character")) {
    message(glue::glue("loading hiseq sample_sheet from dir: {x}"))
    f_list <- list.files(x[1], "*.xlsx", full.names = TRUE)
    if(length(f_list) == 0) {
      message(paste0("No *.xlsx files found in : ", x[1]))
      return(NULL)
    }
    df <- lapply(f_list, read_hiseq_sheet) %>%
      dplyr::bind_rows() %>%
      dplyr::mutate(lib_number = fix_hiseq_lib_number(lib_number)) %>%
      dplyr::arrange(lib_number) %>%
      dplyr::mutate(sampleid = fix_hiseq_sampleid(sampleid))
    # fix user
    df1 <- df %>%
      dplyr::mutate(tmp_user = gsub("\\d+", "", lib_user)) %>%
      tidyr::unite(sp, c(lib_number, tmp_user), sep = "-", remove = FALSE)
    # fix lib_user
    df2 <-  dplyr::select(df, lib_number, lib_user) %>%
        dplyr::mutate(lib_user = gsub("\\d+", "", lib_user)) %>%
        unique() %>%
        dplyr::group_by(lib_user) %>%
        dplyr::mutate(num = stringr::str_pad(row_number(), 2, pad = 0)) %>%
        tidyr::unite(sp, c(lib_number, lib_user), sep = "-")
    # combine
    df3 <- merge(df1, df2, by = "sp") %>%
      dplyr::mutate(lib_user = gsub("\\d+", "", lib_user),
                    lib_user = paste0(lib_user, num)) %>%
      dplyr::select(-c(sp, num, tmp_user))
    # show
    message(glue::glue(
      "{length(f_list)} files found, including {nrow(df3)} samples."
    ))
    df3
  }
}


#' read sample info from Excel.xlsx
#'
#' @param x file xlsx
#' @n_max int max number rows to read
#' @fix bool auto fix the samplename, etc, default: TRUE
#'
#' @export
read_hiseq_sheet <- function(x, n_max = 10000, fix = TRUE,
                             sampleid_start = 1) {
  df <- tryCatch(
    {
      df1 <- readxl::read_xlsx(x, sheet = "sample_sheet",
                               trim_ws = TRUE, n_max = n_max)
      colnames(df1) <- gsub("[^\\w]", "", colnames(df1), perl = TRUE)
      colnames(df1) <- tolower(colnames(df1))
      df1 <- dplyr::filter(df1, ! is.na(sample_name))
    },
    error=function(cond) {
      message(glue::glue("Failed reading sheet: {x}"))
      message(cond)
      return(NA)
    }
  )
  # auto_fix
  if(isTRUE(fix)) {
    df <- fix_hiseq_sheet(df, sampleid_start)
  }
  # fix: date
  f_date <- stringr::str_extract(basename(x), "20[\\d]{6}")
  if(is.na(f_date)) {
    f_date <- Sys.Date()
  } else {
    f_date <- as.Date(f_date, format = "%Y%m%d")
  }
  # add date
  df$date <- f_date
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
fix_hiseq_sheet <- function(df, sampleid_start = 1) {
  required_cols <- c("lib_number", "lib_user", "sample_name", "cell_line",
                     "species", 'spikein', "p7_index_id", "barcode_id",
                     "seq_type", "lib_type")
  if(is(df, "data.frame")) {
    m <- required_cols %in% colnames(df)
    m_cols <- required_cols[! m]
    if(length(m_cols) > 0) {
      m_err <- paste(m_cols, collapse = ", ")
      warning(glue::glue("missing columns: {m_err}"))
      return(NULL)
    }
  } else {
    warning(glue::glue("unknown `df`, expect data.frame, got {class(df)}"))
    return(NULL)
  }
  # hiseq index
  hidx <- load_hiseq_index()
  # fix columns
  if(! "lib_sub" %in% colnames(df)) {
    df$lib_sub <- "G1"
  }
  if(! "rbp" %in% colnames(df)) {
    df$rbp <- "NULL"
  }
  df1 <- df %>%
    dplyr::mutate(
      sampleid    = fix_hiseq_sampleid(sampleid, sampleid_start),
      lib_number  = fix_hiseq_lib_user(lib_number),
      lib_user    = fix_hiseq_lib_user(lib_user),
      lib_sub     = fix_hiseq_lib_sub(lib_sub),
      reminder    = sample_name, # save old names
      sample_name = fix_hiseq_sample_name(sample_name),
      species     = fix_hiseq_species(species),
      spikein     = fix_hiseq_species(spikein),
      p7_index_id = fix_hiseq_index(p7_index_id),
      barcode_id  = fix_hiseq_index(barcode_id),
      seq_type    = toupper(seq_type), #SE|PE
      lib_type    = fix_hiseq_lib_type(lib_type), # .sanitize_str(lib_type),
      readsm      = fix_hiseq_reads(readsm),
      p7_index    = ifelse(is.na(hidx[p7_index_id]), "NULL", hidx[p7_index_id]),
      barcode_seq = ifelse(is.na(hidx[barcode_id]), "NULL", hidx[barcode_id]),
      spikein_ref = ifelse(is.na(spikein_ref), "NULL", spikein_ref),
      fcid        = "",
      lane        = "",
      status      = "",
      date        = ""
    )
  # output
  out_cols <- c("sampleid", "lib_number", "lib_sub",  "lib_user",
                "sample_name", "rbp", "cell_line", "species", 'spikein',
                "p7_index_id", "barcode_id", "seq_type", "lib_type", "readsm",
                "fcid", "lane", "reference", "spikein_ref", "p7_index",
                "barcode_seq", "reminder", "status")
  dplyr::select(df1, dplyr::all_of(out_cols))
}


#' is_valid_sheet_df
#'
#' Check the fields in sheet
#' @param data.frame sample data
#' @param index_list list, named index sequence
#'
#' @export
is_valid_sheet_df <- function(df, verbose = FALSE) {
  required_cols <- c("lib_number", "lib_user", "sample_name",
                     "species", 'spikein', "p7_index_id", "barcode_id",
                     "seq_type", "lib_type", "readsm")
  k1 <- k2 <- k3 <- k4 <- k5 <- FALSE
  if(is(df, "data.frame")) {
    k1 <- all(required_cols %in% colnames(df))
    if(k1) {
      k2 <- ! any(duplicated(df$sample_name))
      k3 <- is_valid_index(df$p7_index_id)
      k4 <- is_valid_index(df$barcode_id)
      k5 <- ! any(duplicated(paste(df$p7_index_id, df$barcode_id)))
    }
  } else {
    warning(glue::glue("Expect data.frame, got {class(df)}"))
  }
  # show msg
  df_msg <- data.frame(
    Fileds = c(
      "header",
      "sample_name",
      "P7_index",
      "barcode",
      "P7+barcode"
    ),
    Status = c(
      k1, k2, k3, k4, k5
    ),
    Expect = c(
      "see above examples",
      "not duplicated",
      "TruSeq_Index1,Next_Ad2.1",
      "P7_1A,P7_1B",
      "not duplicated"
    )
  )
  if(isTRUE(verbose)) {
    print(df_msg)
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



#' fix_hiseq_sampleid
#'
#' @export
fix_hiseq_sampleid <- function(x, start = 1) {
  s <- stringr::str_pad(
    as.character(start + seq_len(length(x)) - 1),
    side = "left", width = 6, pad = "0"
  )
  paste0("YYs", s)
}


#' lib_number, auto_fix
#' eg: YY01, WM03
#'
#' @param x character lib_user
#'
#' @export
fix_hiseq_lib_number <- function(x) {
  j <- toupper(x)
  num <- stringr::str_extract(j, "\\d+$")
  num <- stringr::str_pad(num, 3, pad = "0")
  prefix <- stringr::str_extract(j, "^[A-Z]{2,4}")
  paste0(prefix, num)
  # gsub("([A-Z]{2,4})([\\W])?(\\d+)", "\\1\\3", j)
}


#' lib_number, auto_fix
#' eg: YY01, WM03
#'
#' @param x character lib_user
#'
#' @export
fix_hiseq_lib_user <- function(x) {
  j <- toupper(x)
  ## for specific usage ##
  j <- gsub("^JG(\\d+)", "GJQ\\1", j)
  j <- gsub("^LW(\\d+)", "LWW\\1", j)
  j <- gsub("^LW(\\d+)", "LWW\\1", j)
  j <- gsub("^WL(\\d+)", "WLT\\1", j)
  ## END ##
  gsub("([A-Z]{2,4})([\\W])?(\\d+)", "\\1\\3", j)
}


#' lib_sub
#' eg: G1, G2
#'
#' @param x character lib_sub
#'
#' @export
fix_hiseq_lib_sub <- function(x) {
  i <- toupper(x)
  # gsub("(G)([\\W])?([1-9])$", "\\1\\3", i)
  gsub("G0+", "G", i)
}



#' mission:
#' 1. replace non-letters, by "_"
#' 2. fix hiseq_type
#' 3. add rep1/rep2 in the tail
#'
#' @param x character
#'
#' @export
fix_hiseq_sample_name <- function(x) {
  j <- .sanitize_str(x)
  j <- fix_hiseq_lib_type(j)
  sapply(j, function(k) {
    if(grepl("_(r|rep)(\\d+)$", k, ignore.case = TRUE)) {
      gsub("_(r|rep)(\\d+)$", "_rep\\2", k)
    } else {
      paste0(k, "_rep1")
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
#' @param x character the filename or lib_type
#'
#' @export
fix_hiseq_lib_type <- function(x) {
  j <- gsub("(ATAC|ChIP|DNA|RNA|GoldCLIP|GRO|STACC)([-_]?seq)?", "\\1seq", x,
            perl = TRUE, ignore.case = TRUE)
  j <- gsub("^mRNAseq", "RNAseq", j, ignore.case = TRUE)
  j <- gsub("DNAseq", "DNAseq", j, ignore.case = TRUE)
  j <- gsub("ATAC[_-]?seq", "ATACseq", j, ignore.case = TRUE)
  j <- gsub("ChIPseq", "ChIPseq", j, ignore.case = TRUE)
  # j <- gsub("CLIP(_)?(seq)?(_)?", "CLIPseq_", j, ignore.case = TRUE)
  j <- gsub("(small|sm)(_)?RNAseq", "smRNAseq", j, ignore.case = TRUE)
  j <- gsub("(CUTRUN|CUT_RUN|CUT_and_RUN)", "CnR", j, ignore.case = TRUE)
  j <- gsub("(CUTTAG|CUT_TAG|CUT_and_TAG)", "CnT", j, ignore.case = TRUE)
  j <- gsub("STACCseq", "STACCseq", j, ignore.case = TRUE)
  j <- gsub("GROseq", "GROseq", j, ignore.case = TRUE)
  # move hiseq to prefix
  sapply(j, function(k) {
    if(grepl("[A-Za-z0-9]+seq_", k, ignore.case = TRUE)) {
      s <- stringr::str_extract(k, "([A-Za-z0-9]+seq_)")
      k2 <- gsub("([A-Za-z0-9]+seq_)", "", k)
      paste0(s, k2)
    } else {
      k
    }
  })
}


#' to-do
#'
#' @export
fix_hiseq_species <- function(x) {
  x
}


#' non character
#' @export
fix_hiseq_index <- function(x) {
  j <- .sanitize_str(x)
  gsub("null", "NULL", j, ignore.case = TRUE)
}


#' make sure the number of reads, numeric
#'
#' @export
fix_hiseq_reads <- function(x) {
  j <- gsub("[^0-9]", "", x)
  j[is.na(j) | length(j) < 1] <- 1 # default
  as.numeric(j)
}


#' load index
#' @param hiseq_type character, [truseq, nextera, barcode]
#'
#' @export
load_hiseq_index <- function(hiseq_type = TRUE) {
  f <- system.file("data", "hiseq_index.rds", package = "hiseqr")
  l <- readRDS(f)
  if(isTRUE(hiseq_type)) {
    hiseq_type <- names(l)
  }
  l <- l[hiseq_type]
  names(l) <- NULL
  unlist(l)
}


#' check index id
#' @param x character
#' @param hiseq_type character, p7, barcode
#'
#' @export
is_valid_index <- function(x, hiseq_type = "p7", skip_null = TRUE) {
  # index: TruSeq, Nextera, barcode
  p7 <- load_hiseq_index()
  f1 <- x %in% p7 | x %in% names(p7)
  f2 <- isTRUE(skip_null) & tolower(x) == "null"
  # show msg
  x_invalid <- x[! (f1|f2)]
  if(length(x_invalid) > 0) {
    x_err <- paste(x_invalid, collapse = ", ")
    warning(glue::glue("Unknown index: {x_err}"))
  }
  all(f1 | f2)
}


#' #' read index sequence
#' #' @param x string path to the index file, could be *.rds, *.txt, *.csv
#' #'
#' #' @export
#' read_sheet_index <- function(x) {
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






