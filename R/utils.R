#' Functions for HiSeq pipeline, utils for general data processing
#'
#'
#'
#'
#' @name utils




#' @describeIn read_toml Parse TOML files as data.frame
#'
#'
#' @param x path to the toml file
#'
#' @export
read_toml <- function(x) {
  if(! is.character(x)) {
    stop("`x`: Non character detected")
  }
  x <- x[file.exists(x)]
  if(length(x) > 0) {
    lapply(x, function(i) {
      if(configr::is.toml.file(i)) {
        configr::read.config(i) %>%
          as.data.frame.list()
      }
    }) %>%
      dplyr::bind_rows()
  }
}




#' @describeIn read_align Parse Alignment stat, toml
#' map, unique, multiple, unmap
#'
#' @param x path to stat file
#' @import readr
#' @import dplyr
#'
#' @export
read_align <- function(x){
  if(! is.character(x)) {
    stop("'x': non character detected")
  }
  is_file <- fs::file_exists(x)
  x <- x[is_file]
  lapply(x, function(i) {
    configr::read.config(i) %>%
      as.data.frame.list()
  }) %>%
    bind_rows() %>%
    tidyr::pivot_longer(cols      = where(is.numeric),
                        names_to  = "group",
                        values_to = "count")
}





#' @describeIn read_text parse the file, eg: fragment size
#'
#' @param x path to fragment length
#'
#' @import dplyr
#'
#' @export
read_text <- function(x, ...){
  if(! is.character(x)) {
    stop("'x': non character detected")
  }
  x <- x[file.exists(x)]
  lapply(x, function(i) {
    s <- guess_separator(i)
    read.delim(x, sep = s, ...)
  }) %>%
    dplyr::bind_rows()
}








#' @describeIn guess_separater Guess the separater of the plain text file
#' @description Guess the format of the plain text file
#' ,    csv
#' tab  table
#' :
#' white space [default]
#'
#' Check the first 100 lines of the file, if the separater exists
#'
#' @param x character path to files, csv, txt, ...
#'
#' @export
guess_separator <- function(x, nmax = 100) {
  if(is.character(x)) {
    if(file.exists(x)) {
      x = x[1] # first file
      lines <- readLines(x, n = nmax)
      #
      if(all(grepl("\t", lines))) {
        sep = "\t"
      } else if(all(grepl(",", lines))) {
        sep = ","
      } else if(all(grepl(":", lines))) {
        sep = ":"
      } else if(all(grepl(" ", lines))) {
        sep = " "
      } else {
        sep = " "
      }
    } else {
      warning("file not exists")
    }
  } else {
    warning("Expect a filename, failed")
  }
}






#' @describeIn path_to_list Convert directory structure into nested list
#'
#' @param path character Path to the file or directory
#' @param recursive bool Constructure the nested list recursive
#'
#'
#' @export
path_to_list <- function(path, recursive = FALSE) {
  if(is.null(path)) {
    stop("'path': null detected")
  }
  if(is.character(path)) {
    # only 1 item
    if(length(path) > 1) {
      warning("Multiple items detected, only support 1 item")
    }

    # file exists
    if(! file.exists(path)) {
      stop("'path': not exists")
    }

    # convert to absolute path
    path <- normalizePath(path)
    path <- gsub("\\\\/$", "", path) # trim the ends

    # construct the list-of-list
    if(dir.exists(path) && isTRUE(recursive)) {
      files <- list.files(path, full.names = FALSE, recursive = FALSE)
      sapply(files, function(i) {
        f <- file.path(path, i)
        if(dir.exists(f)) {
          path_to_list(f, recursive = recursive)
        } else {
          f
        }
      }, USE.NAMES = TRUE, simplify = FALSE)
    } else if(file.exists(path)) {
      setNames(list(path), nm = basename(path))
      # path
    } else {
      stop("`path` expected path to a file or directory, failed")
    }
  } else {
    stop("'path': non character detected")
  }
}






#' @describeIn fq_name Retrive the names of fastq files
#'
#' remove .fa, .fa.gz, .fq, .fq.gz
#' remove, "_1" or "_2" suffix
#' remove, _rep1, _rep2, ...
#'
#' @param x character, name of fastq files
#' @param fix_pe bool, trim _1, _2 suffix, default: FALSE
#' @param fix_rep bool, trim _rep1, _r1 in tha file name
#'
#' @export
fq_name <- function(x, fix_pe = TRUE, fix_rep = FALSE){

  # Retrive the file name
  name <- gsub("[.](fast|f)[aq](.gz)?$",
               "",
               basename(x),
               ignore.case = TRUE,
               perl = TRUE)

  # Fix the read1/2 in the tail
  if(isTRUE(fix_pe)) {
    name <- gsub("[._][12]$", "", name, perl = TRUE)
  }

  # Fix the rep name: r1 rep1
  if(isTRUE(fix_rep)) {
    name <- gsub("[._](rep|r)[0-9]+$",
                 "",
                 name,
                 ignore.case = TRUE,
                 perl = TRUE)
  }
  name
}





#' @describeIn read_fc2 Parsing multiple featureCounts files
#'
#' @param x character featureCount files, or data.frames
#'
#' @export
read_fc2 <- function(x) {
  # if character
  chk1 <- sapply(x, is.character)
  chk2 <- sapply(x, file.exists)

  if(all(chk1, chk2)) {
    lapply(x, read_fc) %>%
      bind_cols2()
  }
}



#' @describeIn read_fc Parsing single featureCounts file
#'
#' @param x path to count.txt file, featureCounts output
#' @param digits integers, Number of digits keep, default: 0
#' @param get_gene_length bool, save gene_lengthin column 5
#'
#' @import readr
#' @import dplyr
#'
#' @export
read_fc <- function(x, digits = 0, get_gene_length = FALSE){
  if(is.character(x)) {
    x <- x[1] # the first item
  } else {
    warning(paste0("Reading file failed: ", x))
    return(NULL)
  }

  # Read file
  if(isTRUE(get_gene_length)) {
    # readr::read_delim(x, "\t", col_types = readr::cols(), comment = "#") %>%
    read.table(x, header = TRUE, sep = "\t", comment.char = "#") %>%
      dplyr::select(Geneid, Length) %>%
      dplyr::rename(id = Geneid) %>%
      tibble::as_tibble()
  } else {
    # readr::read_delim(x, "\t", col_types = readr::cols(), comment = "#") %>%
    read.table(x, header = TRUE, sep = "\t", comment.char = "#") %>%
      dplyr::select(-c(2:6)) %>%
      dplyr::rename(id = Geneid) %>%
      dplyr::mutate_if(is.numeric, round, digits = digits) %>%
      tibble::as_tibble()
  }
}



#' @describeIn read_fc_summary Parsing the summary file of featureCounts output
#'
#' @param x path to the count.txt.summary file
#' @param fix_name bool extract file name,
#'
#' @export
read_fc_summary <- function(x, fix_name=TRUE) {
  # function to fix filename
  fname <- function(x, fix_name = TRUE) {
    sapply(x, function(i) {
      if(isTRUE(fix_name)) {
        gsub("\\.\\w+$", "", basename(i))
      } else {
        i
      }
    }, simplify = TRUE)
  }

  if(is.character(x)) {
    lapply(x, function(i){
      readr::read_delim(x, "\t", col_types = readr::cols()) %>%
        tidyr::pivot_longer(names_to = "sample", values_to = "count", -1)
    }) %>%
      dplyr::bind_rows() %>%
      dplyr::mutate(sample = fname(sample, fix_name))
  }
}






#' @describeIn bind_cols2 Merge multiple data.frame by column
#'
#' @param x a list of data.frames
#' @param by column name for function merge()
#'
#' @description see example on https://www.r-bloggers.com/2018/10/how-to-perform-merges-joins-on-two-or-more-data-frames-with-base-r-tidyverse-and-data-table/
#' see dplyr::bind_cols, for how to process dots
#'
#' @export
bind_cols2 <- function(..., by = NULL) {
  dots <- rlang::list2(...)
  dots <- rlang::squash_if(dots, vctrs::vec_is_list)
  dots <- purrr::discard(dots, is.null)
  is_data_frame <- purrr::map_lgl(dots, is.data.frame)
  names(dots)[is_data_frame] <- ""
  if(! all(is_data_frame)) {
    stop("`...`, Input are not all data.frames")
  }
  if(length(dots) > 1) {
    if(is.null(by)) {
      xnames <- lapply(dots, names)
      by <- Reduce(function(a, b, ...) intersect(a, b), xnames)
    }
    out <- Reduce(function(a, b, ...) merge(a, b, by = by, ...), dots)
  } else if(length(dots) == 1) {
    out <- dots
  } else {
    stop("`...`, No data.frame input")
  }
  out
}



#' @describeIn str_similar Extract similar string by Levenshtein-distance,
#'
#' function: utils::adist
#' alternative: stringdist package
#'
#' @param x
#' @param y
#'
#'
#' @import stringdist
#'
#' @export
str_similar <- function(x, y, ignore_case = FALSE) {
  if(is.character(x) && is.character(y)) {
    x <- x[1]
  } else {
    msg <- "require character"
    stop("'x' and 'y': ", msg)
  }
  d <- sapply(y, function(i) {
    adist(x, i, ignore.case = ignore_case)
  }, simplify = TRUE)
  di <- which.min(d)
  y[di]
}




# --Functions for html report --------------------------------------------------

#' @describeIn to_DT convert data.frames to DT table
#'
#' @param df data.frame
#' @param mode integer default 1
#' @param page_length integer default 10
#'
#' @description see https://rstudio.github.io/DT/ for details
#'
#' @export
to_DT <- function(df, mode = 1, pageLength = 10) {
  if(! is.data.frame(df)) {
    stop("`df`, not a data.frame")
  }
  if(mode == 1) {
    DT::datatable(
      df,
      extensions   = 'Buttons',
      options      = list(
        pageLength = pageLength,
        scrollX    = TRUE,
        dom        = 'Bfrtip',
        buttons    =
          list('copy', #'print',
               list(extend  = 'collection',
                    buttons = c('excel', 'csv'),
                    text    = 'Download')))
    )
  } else if(mode == 2) {
    DT::datatable(
      df,
      rownames = TRUE,
      filter   = "top",
      options  = list(
        pageLength = pageLength,
        scrollX    = TRUE
      )
    )
  } else if(mode == 3){
    DT::datatable(
      df,
      rownames = TRUE,
      escape   = FALSE, # show html code
      options  = list(
        dom        = 'Bftp',
        pageLength = pageLength,
        scrollX    = TRUE,
        extensions = 'Buttons',
        buttons    = list('copy',
                          list(extend  = "collection",
                               buttons = "excel",
                               text    = "Download"))
      )
    )
  } else {
    DT::datatable(df)
  }
}





