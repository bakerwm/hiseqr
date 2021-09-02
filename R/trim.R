#' functions for hiseq trim
#'
#' read data from hiseq.trim output
#'
#' @name trim



#' @describeIn read_hiseq_trim
#'
#' @param x character or data.frame
#'
#' @export
read_hiseq_trim <- function(x) {
  if(is(x, "character")) {
    if(endsWith(x, ".json")) {
      df <- tryCatch(
        error = function(cnd) NULL,
        jsonlite::read_json(x) %>%
          as.data.frame
      )
      # basic
      df_cols       <- colnames(df)
      required_cols <- c("name", "total", "clean", "too_short", "dup")
      extra_cols    <- df_cols[! df_cols %in% required_cols]
      if(all(rlang::has_name(df, required_cols))) {
        dplyr::select(df, all_of(c(required_cols, extra_cols))) %>%
          dplyr::rename(
            input = total,
            output = clean
          )
      } else {
        warning(glue::glue("unknown .json file:"))
      }
    } else if(is_hiseq_dir(x, "trim")) {
      pd <- read_hiseq(x)
      if(pd$hiseq_type == "trim_r1") {
        j_list <- list_hiseq_file(x, "trim_json")
      } else if(pd$hiseq_type == "trim_rn") {
        j_list <- list_hiseq_file(x, "trim_json", "r1")
        # project_dir <- list_hiseq_file(x, "project_dir", TRUE)
        # j_list <- sapply(list_hiseq_file(x, "smp_name", TRUE), function(i) {
        #   r1 <- file.path(project_dir, i)
        #   list_hiseq_file(r1, "trim_json")
        # })
      }
      # load data
      lapply(j_list, read_hiseq_trim) %>%
        dplyr::bind_rows()
    } else {
      warning(glue::glue("unknown x: {x}"))
    }
  } else {
    warning(glue::glue("expect character, got {x}"))
  }
}













