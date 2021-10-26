#' Functions for HiSeq report
#'
#' Summary data in project_dir.
#' Processing data
#' Prepare for plot
#' Plotting
#'
#' @name hiseq_report




#' hiseq_report
#'
#' @param input directory to the sample
#' @param output directory to the html file
#' @param template the template, default from hiseqr
#'
#' @export
hiseq_report <- function(input, output = NULL, template_rmd = NULL) {
  if(!is_hiseq_dir(input)) {
    warning(paste0("input is not hiseq dir: ", input))
    return(NULL)
  }
  input <- normalizePath(input)
  if(!inherits(output, "character")) {
    output <- file.path(input, "report")
  }
  output <- file.path(normalizePath(dirname(output)), basename(output))
  check_path(output)
  outhtml <- file.path(output, "HiSeq_report.html")
  # check input
  pd <- read_hiseq(input)
  ps <- unlist(strsplit(pd$hiseq_type, "_")) # cnr r1
  if(length(ps) != 2) {
    warning(paste0("unknown hiseq_type: ", pd$hiseq_type))
    return(NULL)
  }
  if(inherits(template_rmd, "character")) {
    template <- template_rmd[1]
  } else {
    template_name <- paste0("hiseq_report_", ps[2], ".Rmd")
    template <- system.file(ps[1], template_name, package = "hiseqr")
  }
  if(!file.exists(template)) {
    warning(paste0("template rmarkdown file not exists: ", template))
    return(NULL)
  }
  # run
  template_to <- file.path(output, basename(template))
  file.copy(template, template_to, overwrite = TRUE)
  rmarkdown::render(input       = template_to,
                    output_file = outhtml,
                    params      = list(input_dir = input))
}



#' csv_to_html
#'
#' @param input csv file
#' @param output html file
#' @param template the template, default from hiseqr
#'
#' @export
csv_to_html_report <- function(input, output) {
  input <- normalizePath(input)
  output <- file.path(normalizePath(dirname(output)), basename(output))
  # output <- normalizePath(output)
  template <- system.file('utils', 'csv_to_html.Rmd', package = "hiseqr")
  stopifnot(file.exists(template))
  rmarkdown::render(input       = template,
                    output_file = output,
                    params      = list(input_csv = input))
}



#' @describeIn fix samples names
#'
#' fix hiseq sample names, remove the most common string in names
#'
#' @param x names
#'
#' @export
fix_hiseq_names <- function(x, ...) {
  x <- "~/work/devel_pipeline/hiseq/atac/output/fruitfly/"
  s <- list_hiseq_file(x, "smp_name", "r1")
  Biobase::lcPrefix(s)
}




#' @describeIn fig_to_panel
#'
#' Generate panel in Xaringan slides using XaringanExtra package
#'
#' one figure in each page
#'
#' Example:
#' .panelset[
#' .panel[.panel-name[R Code]
#'
#'  ```{r panel-chunk, fig.show='hide'}
#'  # ... r code ...
#'  ```
#'  ]
#' .panel[.panel-name[Plot]
#'  ![](README_files/figure-gfm/panel-chunk-1.png)
#'  ]
#'  ]
#'
#' @param x character path to the figures
#' @param nm name of the figures
#'
#' @export
fig_to_panel <- function(x, nm = NULL, ...) {
  if(! is(x, "character")) {
    warning("Only characters supported")
  }
  # for name
  if(is.null(nm)) {
    nm <- gsub("\\.\\w+$", "", basename(x))
  } else if(is(x, "numeric")) {
    nm <- seq_len(length(x))
  } else if(is(x, "character")) {
    if(! length(x) == length(nm)) {
      nm <- seq_len(length(x))
      warning("invalid nm=, length not consistent with files")
    } else {
      a <- 1
    }
  } else {
    warning("unknown nm, use int instead")
    nm <- seq_len(length(x))
  }
  # panel-body
  p_body <- lapply(seq_len(length(x)), function(i) {
    .single_panel(x[i], i, ...)
  }) %>%
    unlist %>%
    paste(collapse = "\n")
  # panel-frame
  glue::glue(
    ".panelset[",
    {p_body},
    "]",
    sep = "\n"
  )
}


#' Convert url to link in html
#' resize the image by scale
#' <img src="img.jpg" alt="a img" style="width:500px;height:600px;">
.url_to_link <- function(x, ...) {
  args <- rlang::list2(...)
  style <- ifelse("style" %in% names(args), args$style, "height:100%")
  alt   <- ifelse("alt" %in% names(args), args$alt, "figure")
  glue::glue('<img src="{x}" alt="{alt}" style="{style}">')
}




# add image to panel
.single_panel <- function(x, n=NULL, ...) {
  if(is.null(n)) {
    n <- gsub("\\.\\w+$", "", basename(x))
  }
  # format
  # image format: markdown, html
  #
  # "![]({x})",
  glue::glue(
    ".panel[",
    ".panel-name[{n}]",
    "## {basename(x)}",
    .url_to_link(x, ...),
    "]",
    sep = "\n"
  )
}
#--subfunctions--#












