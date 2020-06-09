
## Alignment

##----------------------------------------------------------------------------##
## For hiseq package, parsing directory
##

#' list_arguments_file
#'
#' @param x path to the directory.
#' expect x/config/{feature}/config/arguments.pickle,
#' x/{feature}/config/arguments.pickle
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


#' align
#'
#' @param x path, the directory of RNAseq output
#' @param log boolen, whether print the log status, default: FALSE
#'
#' @export
#'
is_align_dir <- function(x) {
  chk0 <- dir.exists(x)

  # check arguments.pickle file
  pk_files <- list_align_arguments_file(x)

  if(is.null(pk_files) | is.na(pk_files)) {
    return(NULL)
  } else {
    # pk_files <- pk_files[1] # 1st item
    pk <- load_pickle(pk_files[[1]])

    # get values
    if("align_type" %in% names(pk)) {
      tag <- pk$align_type #'alignment'
    } else if("rnaseq_type" %in% names(pk)) {
      tag <- pk$rnaseq_type
    } else if("atacseq_type" %in% names(pk)) {
      tag <- pk$atacseq_type
    } else {
      tag <- NULL # not known
    }
    return(tag) # return
  }
}



#' align_single_dir
#'
#' reading data from RNAseq Single
#'
#' @param x path to the directory of RNAseq single
#'
#' @export
#'
align_single_dir <- function(x, feature="gene"){
  # check
  chk0 <- is_align_dir(x)[1] == 1

  if(! isTRUE(chk0)){
    return(NULL)
  }

  # pick pickle file
  pk_files <- list_align_arguments_file(x)
  args     <- load_pickle(pk_files)

  ## check required args
  required_names <- c("align_stat", "smp_name", "genome")

  for(name in required_names){
    # args[[name]] = ifelse(name %in% names(args), c(args[[name]]), NULL)
    if(! name %in% names(args)){
      args[[name]] <- NULL
    }
  }

  ## mapping data
  df <- hiseqr::readAlign1(args$align_stat)
  args$reads_total <- sum(df$count)

  ## mito.u, mito.m, genome.u genome.m
  for(i in seq_len(nrow(df))) {
    name <- as.character(df$group[i])
    args[[name]] <- df$count[[i]]
  }

  ## save df
  args$align_stat_df <- df

  return(args)
}





#' rnaseq_multiple_dir
#'
#' return args for RNAseq_multiple/single/
#'
#' @param x path to the directory of RNAseq multiple
#'
#' @export
#'
#' @export
align_multiple_dir <- function(x, feature = "gene"){
  # check
  chk0 <- is_align_dir(x)[1] > 1

  if(! isTRUE(chk0)){
    return(NULL)
  }

  # pick pickle file
  pk_files <- list_align_arguments_file(x)
  args     <- hiseqr::load_pickle(pk_files)

  ## check required args
  required_names <- c("smp_name", "outdir", "genome")

  for(name in required_names){
    if(! name %in% required_names) {
      args[[name]] <- NULL
    }
  }

  ## parsing RNAseq single, saved as vector
  args$single_dirs <- mapply(function(i) file.path(args$outdir, i), args$smp_name)
  args$single_args <- lapply(args$single_dirs, align_single_dir)

  return(args)
}


