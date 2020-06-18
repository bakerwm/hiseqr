
# make report

#' atac_report_single
#'
#' @param input directory to the sample
#' @param output directory to the html file
#' @param template the template, default from hiseqr
#'
#' @export
atac_report <- function(input, output) {
  ## input
  input <- normalizePath(input)
  if(! dir.exists(output)) dir.create(output, recursive = TRUE)
  output <- normalizePath(output)
  ## output
  outhtml <- file.path(output, "atac_report.html")

  if(is_atac_single_dir(input)){
    template <- system.file('atacseq', 'atac_report_single.Rmd', package = "hiseqr")
  } else if(is_atac_merge_dir(input)) {
    template <- system.file('atacseq', 'atac_report_merge.Rmd', package = "hiseqr")
  } else if(is_atac_multiple_dir(input)) {
    template <- system.file("atacseq", "atac_report_multiple.Rmd", package = "hiseqr")
  } else {
    warning("unknown input:")
    stop("Stop report")
  }

  ## copy template to output
  template_to <- file.path(output, basename(template))
  file.copy(template, template_to)

  rmarkdown::render(input       = template_to,
                    output_file = outhtml,
                    params      = list(input_dir = input))
}




