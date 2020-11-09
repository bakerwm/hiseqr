#' Functions for HiSeq pipeline
#'
#'
#'



#' hiseq_type
#'
#' @param x path, the directory of RNAseq output
#'
#' @export
is_hiseq_dir <- function(x) {
  chk0 <- dir.exists(x)

  # check arguments.pickle file
  pk_files <- list_arguments_file(x)
  # choose one

  if (is.null(pk_files)) { # | is.na(pk_files)) {
    return(NULL)
  } else {
    # pk_files <- pk_files[1] # 1st item
    pk <- load_pickle(pk_files[[1]])

    # get values
    if ("rnaseq_type" %in% names(pk)) {
      tag <- pk$rnaseq_type
    } else if ("atacseq_type" %in% names(pk)) {
      tag <- pk$atacseq_type
    } else if("align_type" %in% names(pk)) {
      tag <- "align"
    } else if("chipseq_type" %in% names(pk)) {
      tag <- pk$chipseq_type
    } else {
      tag <- NULL # not known
    }
    # return(tag) # return
    tag
  }
}


