
# make report

#' atac_report_single
#'
#' @param input directory to the sample
#' @param output directory to the html file
#' @param template the template, default from hiseqr
#'
#' @export
atac_report <- function(input, output, template_rmd = NA) {
  ## input
  input <- normalizePath(input)
  if(! dir.exists(output)) dir.create(output, recursive = TRUE)
  output <- normalizePath(output)
  ## output
  outhtml <- file.path(output, "atac_report.html")

  ## check input is atac_single, or atac_multi
  if(is_atac_dir(input)){
    template <- system.file("atac_report_single.Rmd", package = "hiseqr")
    input_type = "single"
  } else if(is_atac_merge_dir(input)) {
    template <- system.file("atac_report_merge.Rmd", package = "hiseqr")
  } else {
    dirs <- atac_dirs(input)
    if(length(dirs) >= 1){
      template <- system.file("atac_report_multiple.Rmd", package = "hiseqr")
      input_type = "multi"
    } else {
      input_type = NULL
      stop(paste0("ATACseq files not found in: ", input))
    }
  }

  ##
  if(! is.na(template_rmd)){
    template <- template_rmd # custome
  }

  ## copy template to output
  template_to <- file.path(output, basename(template))
  file.copy(template, template_to)

  rmarkdown::render(input       = template_to,
                    output_file = outhtml,
                    params      = list(input_dir = input))
}


#' rnaseq_report
#'
#' @param input directory to the sample
#' @param output directory to the html file
#' @param template the template, default from hiseqr
#'
#' @export
rnaseq_report <- function(input, output, feature = "gene", template_rmd = NA) {
  ## input
  input <- normalizePath(input)
  if(! dir.exists(output)) dir.create(output, recursive = TRUE)
  output <- normalizePath(output)
  ## output
  outhtml <- file.path(output, "rnaseq_report.html")

  ## check input is atac_single, or atac_multi
  hiseq_type <- is_hiseq_dir(input)

  if(hiseq_type == "rnaseq_single"){
    template <- system.file("rnaseq_single_report.Rmd", package = "hiseqr")
  } else if(hiseq_type == "rnaseq_multiple") {
    template <- system.file("rnaseq_multiple_report.Rmd", package = "hiseqr")
  } else if(hiseq_type == "deseq_single") {
    template <- system.file("rnaseq_deseq_single_report.Rmd", package = "hiseqr")
  } else if(hiseq_type == "deseq_single") {
    template <- system.file("rnaseq_deseq_multiple_report.Rmd", package = "hiseqr")
  } else {
    dirs <- atac_dirs(input)
    if(length(dirs) >= 1){
      template <- system.file("atac_report_multiple.Rmd", package = "hiseqr")
      input_type = "multi"
    } else {
      input_type = NULL
      stop(paste0("ATACseq files not found in: ", input))
    }
  }

  ##
  if(! is.na(template_rmd)){
    template <- template_rmd # custome
  }

  ## copy template to output
  template_to <- file.path(output, basename(template))
  file.copy(template, template_to)

  rmarkdown::render(input       = template_to,
                    output_file = outhtml,
                    params      = list(input_dir = input,
                                       feature   = feature))
}

