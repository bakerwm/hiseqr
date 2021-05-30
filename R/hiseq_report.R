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
hiseq_report <- function(input, output, template_rmd = NULL) {
  ## input
  input  <- normalizePath(input)
  if(! dir.exists(output)) {
    dir.create(output, recursive = TRUE)
  }
  output <- normalizePath(output)
  ## output
  outhtml <- file.path(output, "HiSeq_report.html")

  ## determine the template
  hiseq_type = get_hiseq_type(input)

  # hiseq (hiseq, atacseq, chipseq, rnaseq, ...)
  if(is_hiseq_dir(input)) {
    if(startsWith(hiseq_type, "hiseq_")) {
      hiseq_dir = "hiseq"
    } else if(startsWith(hiseq_type, "cnr_")) {
      hiseq_dir = "cnr"
    } else if(startsWith(hiseq_type, "atacseq_")) {
      hiseq_dir = "atacseq"
    } else if(startsWith(hiseq_type, "chipseq_")) {
      hiseq_dir = "chipseq"
    } else if(startsWith(hiseq_type, "rnaseq_")) {
      hiseq_dir = "rnaseq"
    } else {
      hiseq_dir = "hiseq"
    }
  }

  # subtype
  if(is_hiseq_single_dir(input)) {
    hiseq_subtype = "hiseq_report_single.Rmd"
  } else if(is_hiseq_merge_dir(input)) {
    hiseq_subtype = "hiseq_report_merge.Rmd"
  } else if(is_hiseq_multiple_dir(input)) {
    hiseq_subtype = "hiseq_report_multiple.Rmd"
  } else {
    hiseq_subtype = "tmp"
  }

  # template
  if(is.null(template_rmd)) {
    template <- system.file(hiseq_dir, hiseq_subtype, package = "hiseqr")
  } else {
    template <- template_rmd
  }
  stopifnot(file.exists(template))

  ## copy template to output
  template_to <- file.path(output, basename(template))
  file.copy(template, template_to, overwrite = TRUE)
  rmarkdown::render(input       = template_to,
                    output_file = outhtml,
                    params      = list(input_dir = input))
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












