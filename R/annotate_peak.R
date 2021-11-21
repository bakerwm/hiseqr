#' annotate peaks using ChIPseeker package
#'
#'
#'
#' @name annotate_peak


#' to-do
#' 1. annotation, plotAnnoBar
#' 2. tss, plotDistToTSS
#' 3. tss, region, plotAvgProf
#' 4. annotation, vennpie, upsetplot(peakAnno, vennpie=TRUE)



#' @describeIn read_peak_annotation
#'
#' @param x peak file or hiseq_dir path
#' @param genome string, optional if x is hiseq_dir
#' @param tss_region numeric, 3000
#'
#' @export
read_peak_annotation <- function(x, genome = NULL, tss_region = 3000,
                                 ...) {
  message("read peak annotation ...")
  # check parameters
  if(!inherits(x, "character")) {
    message(glue::glue(
      "read_peak_annotation() failed, x is {class(x)}, expect character"
    ))
    return(NULL)
  }
  # check peak
  if(all(is_hiseq_dir(x))) {
    peak <- setNames(list_hiseq_file(x, "peak"),
                     nm = list_hiseq_file(x, "smp_name"))
    genome <- list_hiseq_file(x[1], "genome")
  } else if(all(file.exists(x))) {
    peak <- x
  }
  # check genome
  if(!is_valid_organism(genome)) {
    message(glue::glue(
      "read_peak_annotation() failed, genome={genome} is not valid"
    ))
    return(NULL)
  }
  txdb <- get_txdb(genome)
  # check peak have names
  if(is.null(names(peak))) {
    names(peak) <- basename(peak)
  }
  # run
  lapply(peak, annotatePeak, TxDb = txdb,
         tssRegion = c(-tss_region, tss_region), verbose = FALSE)
}



#'
#' make plot
#' @param x csAnno object or list of csAnno
#' @param plot_type string, bar, tss, venn, upset
#'
#' @export
plot_peak_annotation <- function(x, genome = NULL, tss_region = 3000,
                                 plot_type = "bar", fish = "Trimma_lantana",
                                 ...) {
  message("plot peak annotation ...")
  # run annotate
  if(inherits(x, "csAnno") | all(sapply(x, function(i) inherits(i, "csAnno")))) {
    pa <- x
  } else if(inherits(x, "character")) {
    pa <- tryCatch(
      {
        read_peak_annotation(x, genome, tss_region)
      },
      error=function(cond) {
        message("Failed reading peak profile")
        message(cond)
        return(NA)
      }
    )
  } else {
    message(glue::glue("x is {class(x)}, expect str,list"))
    pa <- NULL
  }
  # annotate results
  if(is.null(pa)) {
    message("failed to read peak annotation")
    return(NULL)
  }
  # plot
  if(plot_type == "bar") {
    p <- plotAnnoBar(pa)
  } else if(plot_type == "tss") {
    p <- plotDistToTSS(pa)
  } else if(plot_type == "venn") {
    p <- lapply(pa, vennpie) # multiple plots
  } else {
    message(glue::glue("unknown plot_type: {plot_type}"))
    p <- NULL
  }
  # output
  if(inherits(p, "ggplot")) {
    out <- p +
      fishualize::scale_fill_fish_d(option = fish)
  } else if(inherits(p, "list")) {
    if(all(sapply(p, function(i) inherits(i, "ggplot")))) {
      out <- lapply(p, function(i) {
        i +
          fishualize::scale_fill_fish_d(option = fish)
      })
    } else {
      out <- NULL
    }
  } else {
    out <- NULL
  }
  # return
  return(out)
}


#' @describeIn read_peak_profile
#'
#' on promoter
#'
#' @param x peak file or hiseq_dir path
#' @param genome string, optional if x is hiseq_dir
#' @param tss_region numeric, 3000
#'
#' @export
read_peak_profile <- function(x, genome = NULL, flanking = 2000, bed = NULL,
                              on_summit = FALSE) {
  message("reading peak profile")
  if(!inherits(x, "character")) {
    message(glue::glue(
      "read_peak_profile() faild, x is {class(x)}, expect character"
    ))
    return(NULL)
  }
  if(all(is_hiseq_dir(x))) {
    genome <- list_hiseq_file(x[1], "genome")
    peak   <- list_hiseq_file(x, "peak")
    summit <- sub("_peaks.narrowPeak", "_summits.bed", peak)
  } else {
    peak <- x
    summit <- NULL
  }
  # check peak have names
  if(is.null(names(peak))) {
    names(peak) <- basename(peak)
  }
  # check genome
  if(!is_valid_organism(genome)) {
    message(glue::glue(
      "read_peak_profile() failed, genome={genome} is not valid"
    ))
    return(NULL)
  }
  txdb <- get_txdb(genome)
  # choose region: promoter, extra_bed, summit
  promoter <- getPromoters(TxDb = txdb,
                           upstream = flanking,
                           downstream = flanking)
  if(isTRUE(on_summit) & inherits(summit, "character")) {
    region_bed <- summit
    gr <- read_peak(summit, upstream = flanking, downstream = flanking)
  } else if(inherits(bed, "character")) {
    region_bed <- bed
    gr <- read_peak(bed, upstream = flanking, downstream = flanking)
  } else {
    region_bed <- NULL # skipped
    gr <- promoter
  }
  if(length(gr) > 1 & !length(gr) == length(peak)) {
    message(glue::glue(
      "bed and peak not in the same length, ",
      "peak={length(peak)}, bed={length(region_bed)}, ",
      "use [promoter]"
    ))
    gr <- promoter
  }
  # convert gr to list, match peak
  if(inherits(gr, "GRanges")) {
    tag_list <- lapply(peak, getTagMatrix, windows = promoter)
  } else if(inherits(gr, "list")) {
    tag_list <- lapply(seq_len(length(peak)), function(i) {
      getTagMatrix(peak[[i]], windows = gr[[i]])
    })
  } else {
    message("read_peak_profile() failed, illegal bed file")
    return(NULL)
  }
  # add names to tag list
  names(tag_list) <- names(peak)
  tag_list
}


#' @describeIn plot_peak_profile
#'
#' @export
plot_peak_profile <- function(x, genome = NULL, flanking = 2000, bed = NULL,
                              on_summit = FALSE, facet = "none", conf = NA,
                              fish = "Trimma_lantana", ...) {
  # auto-recognize x, peaks or matrix
  if(inherits(x, "matrix")) {
    tag_matrix <- x
  } else if(inherits(x, "list")) {
    if(all(sapply(x, is.matrix))) {
      tag_matrix <- x
    } else {
      tag_matrix <- NULL
    }
  } else {
    tag_matrix <- tryCatch(
      {
        read_peak_profile(x, genome, flanking, bed, on_summit)
      },
      error=function(cond) {
        message("Failed reading peak profile")
        message(cond)
        return(NA)
      }
    )
  }
  # check profile matrix
  if(is.null(tag_matrix)) {
    message("plot_peak_profile() failed")
    return(NULL)
  }
  # plot
  p <- tryCatch(
    {
      plotAvgProf(tag_matrix, xlim = c(-flanking, flanking), conf = conf,
                  resample = 1000, facet = facet) +
        fishualize::scale_color_fish_d(option = fish)
    },
    error=function(cond) {
      message("plotAvgProf() failed")
      message(cond)
      return(NULL)
    }
  )
  # output
  p
}


read_peak <- function(x, upstream = 2000, downstream = 20000) {
  if(inherits(x, "character")) {
    gr_list <- lapply(x, function(i) {
      gr <- rtracklayer::import(i)
      start(gr) = start(gr) - upstream
      end(gr) = end(gr) + downstream
      gr
    })
    # single
    if(length(x) == 1) {
      gr_list[[1]]
    } else {
      gr_list
    }
  }
}


#' @describeIn get_txdb
#'
#' @param x character The name of organism, or build name, eg: dm3, fruitfly
#' Support human, mouse and fruitfly
#'
#' to-do: fetch from bioconductor
#'   laod orgdb
#'   GOSemSim::load_OrgDb()
#'
#' @export
get_txdb <- function(x) {
  #----------------------------------------------------------------------------#
  if(!inherits(x, "character")) {
    warning(glue::glue("x is '{class(x)}', expect character"))
    return(NULL)
  }
  #-- empty, na
  if(length(x) == 0 | is.na(x)) {
    message(glue::glue("x is {x}, zero in length, or is NA"))
    return(NULL)
  }
  #-- multiple items
  if(length(x) > 1) {
    x_str <- paste(x[1:3], collapse = ", ")
    message(glue::glue("x, multiple items found, choose the first: {x_str}"))
    x <- x[1]
  }
  #----------------------------------------------------------------------------#
  #-- NA, 0 length
  sci_name <- get_organism_name(x) # to scientific name
  sup_org  <- Organism.dplyr::supportedOrganisms() # OrgDb, TxDb
  r_idx    <- sup_org$organism == sci_name
  if(sum(r_idx) > 0) {
    txdb <- unique(sup_org[r_idx, ][["TxDb"]])
    # txdb <- purrr::keep(txdb, function(i) endsWith(i, "knownGene"))
    # latest build
    b <- as.numeric(deseq_sanitize_str(txdb, 3))
    bi <- match(max(b), b)
    txdb <- txdb[bi]
    # load
    if(!require(txdb, character.only = TRUE)) {
      BiocManager::install(txdb)
    }
    eval(parse(text = txdb)) # output
  } else {
    warning(glue::glue(
      "get_txdb() failed, invalid x={x}, ",
      "Using function `Organism.dplyr::supportedOrganisms()` ",
      "to list supported organisms,",
    ))
  }
}

