
#' rnaseq pipe
#' for python package: hiseq rnaseq
#'
#' @param deseq_dir string path to the a.vs.b directory
#' @param feature string gene|te|...
#'
#' @export
#'
rnaseq_pipe <- function(deseq_dir, feature = "gene") {
  # check input (deseq_single)
  dir_type <- is_hiseq_dir(deseq_dir) #

  if(is.null(dir_type)) {
    warning("not a RNAseq directory, wt.vs.mut expected")
    return(NULL)
  } else if(dir_type == "deseq_single") {
    message(paste0("Found DESeq2 output dir: ", deseq_dir))
  } else {
    warning("deseq_single dir, eg: wt.vs.mut expected, ", dir_type, " got")
    return(NULL)
  }

  # read config from config (pickle)
  args   <- deseq_single_dir(deseq_dir, feature)

  # check required objects
  message("Checking required data:")
  required_names <- c("count_ctl", "count_exp", "prefix_ctl",
                      "prefix_exp", "deseqdir", "genome")

  chk0 <- sapply(required_names, function(i){
    tag <- ifelse(i %in% names(args), "ok", "failed")
    message(paste0("check ", i, "... ", tag))
  })
  stopifnot(all(unlist(chk0)))

  ##--------------------------------------------------------------------------##
  ## RNAseq: deseq
  outdir1 <- args$deseqdir
  deseq_hub(deseq_dir, feature, outdir1)

  ##--------------------------------------------------------------------------##
  ## RNAseq: GO analysis, KEGG analysis
  if(feature == "gene") {
    go_pipe(deseq_dir, feature)
  }
}



#' rnaseq_report
#'
#' @param input directory to the sample
#' @param output directory to the html file
#' @param template the template, default from hiseqr
#'
#' @export
rnaseq_report <- function(input, output, feature = "gene", template_rmd = NULL) {
  ## input
  input <- normalizePath(input)
  if(! dir.exists(output)) dir.create(output, recursive = TRUE)
  output <- normalizePath(output)
  ## output
  outhtml <- file.path(output, "rnaseq_report.html")

  ## check input is atac_single, or atac_multi
  hiseq_type <- is_hiseq_dir(input)

  if(hiseq_type == "rnaseq_single"){
    template <- system.file("rnaseq", "rnaseq_single_report.Rmd", package = "hiseqr")
  } else if(hiseq_type == "rnaseq_multiple") {
    template <- system.file("rnaseq", "rnaseq_multiple_report.Rmd", package = "hiseqr")
  } else if(hiseq_type == "deseq_single") {
    template <- system.file("rnaseq", "rnaseq_deseq_single_report.Rmd", package = "hiseqr")
  } else if(hiseq_type == "deseq_single") {
    template <- system.file("rnaseq", "rnaseq_deseq_multiple_report.Rmd", package = "hiseqr")
  } else {
    dirs <- atac_dirs(input)
    if(length(dirs) >= 1){
      template <- system.file("rnaseq", "rnaseq_report_multiple.Rmd", package = "hiseqr")
      input_type = "multi"
    } else {
      input_type = NULL
      stop(paste0("RNAseq files not found in: ", input))
    }
  }

  ##
  if(! is.null(template_rmd)){
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



