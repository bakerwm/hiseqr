
## Alignment


#' read_hiseq_align
#' 
#' @param x character or json
#'
#' @export
read_hiseq_align <- function(x) {
  if(is(x, "character")) {
    if(endsWith(x, ".json")) {
      df <- tryCatch(
        error = function(cnd) NULL,
        jsonlite::read_json(x) %>%
          as.data.frame
      )
      # basic
      required_cols <- c("name", "index", "total", "map",
                         "unique", "multi", "unmap")
      if(all(rlang::has_name(df, required_cols))) {
        dplyr::select(df, all_of(required_cols))
      } else {
        warning(glue::glue("unknown .json file:"))
      }
    } else if(is_hiseq_dir(x, TRUE)) {
      pd <- read_hiseq(x)
      # if(pd$hiseq_type == "alignment") {
      # if(is_hiseq_dir(x, "alignment")) {
      j <- list_hiseq_file(x, "align_json", TRUE) %>%
        unique()
      lapply(j, read_hiseq_align) %>%
        dplyr::bind_rows()
      # } else {
      #   warning(glue::glue("unknown hiseq dir: {pd$hiseq_type}"))
      # }
    } else {
      warning(glue::glue("not a hiseq dir: {x}"))
    }
  } else {
    warning(glue::glue("unknown x: {x}"))
  }
}


##----------------------------------------------------------------------------##
## For hiseq package, parsing directory
##
read_align_prj <- function(x) {
  x_type <- is_hiseq_dir(x)
  if(is_align_dir(x)) {
    message(paste0("Reading: ", x_type))
  } else {
    warning(paste0("Not a align dir: ", x))
    return(NULL)
  }

  # read args
  pk <- list_align_arguments_file(x)
  pk <- unlist(pk)

  # read pickle
  load_pickle(pk[[1]]) # the first one
}


#' load all pickles
read_align_prj2 <- function(x) {
  if(is_align_single_dir(x)) {
    read_align_prj(x)
  } else if(is_align_multiple_dir(x)) {
    m <- read_align_prj(x)
    # sub directories
    sub <- file.path(m$outdir, m$smp_name)
    p   <- lapply(sub, read_align_prj)
    names(p) <- m$smp_name
    p # return
  } else {
    warning(paste0("Not a align dir: ", x))
  }
}


#' list_arguments_file
#'
#' @param x path to the directory.
#'   expect x/config/{feature}/config/arguments.pickle,
#'   x/{feature}/config/arguments.pickle
#'
#' @export
list_align_arguments_file <- function(x){
  # arguments file
  # Alignment multiple: x/config/arguments.txt
  # Alignment single: x/config/arguments.txt
  chk0 <- dir.exists(x)

  # version-1: x/{feature}/...
  f1 <- file.path(x, "config", "arguments.pickle")
  chk1 <- file.exists(f1)

  if(isTRUE(chk1)) {
    return(f1)
  } else {
    return(NULL)
  }
}

#' Check if align dir
#' 
#' @param x string, path to dir
#'
#' @export
is_align_dir <- function(x) {
  x_type <- is_hiseq_dir(x)

  if(! is.null(x_type)) {
    grepl("^align", x_type)
  }
}


#' Check if align_r1
#' 
#' @param x string, path to dir
#'
#' @export
is_align_single_dir <- function(x) {
  if(is_align_dir(x)){
    p <- read_align_prj(x) # 2 2 (fq, index)
    p$align_type[1] == 1
  }
}


#' Check if align_rn
#' 
#' @param x string, path to dir
#'
#' @export
is_align_multiple_dir <- function(x) {
  if(is_align_dir(x)){
    p <- read_align_prj(x) # 2 2 (fq, index)
    p$align_type[1] > 1
  }
}


#' align stat
#'
#' @param x string path to the dir
#' @export
get_align_prj_stat <- function(x) {
  if(is_align_single_dir(x)) {
    p <- read_align_prj2(x)
    read_align1(p$align_stat)
  } else if(is_align_multiple_dir(x)) {
    p <- read_align_prj2(x)
    tmp <- lapply(p, function(i) {read_align1(i$align_stat)} )
    dplyr::bind_rows(tmp)
  } else {
    warning(paste0("Not a align dir: ", x))
  }
}



#' align stat
#'
#' @param x string path to the dir
#' @export
get_align_prj_stat2 <- function(x) {
  if(is_align_single_dir(x)) {
    p <- read_align_prj2(x)
    read_align3(p$align_stat)
  } else if(is_align_multiple_dir(x)) {
    p <- read_align_prj2(x)
    tmp <- lapply(p, function(i) {read_align3(i$align_stat)} )
    dplyr::bind_rows(tmp)
  } else {
    warning(paste0("Not a align dir: ", x))
  }
}

