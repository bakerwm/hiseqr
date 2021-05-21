#' Functions for DiffBind analysis
#'
#' See documentation: http://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf
#'
#'
#' ! to-do
#' peak annotation
#' Support: dm6, only (up-to-date)
#'
#'
#' Input design, outdir
#' Output diff peaks, plots
#'
#' helper functions:
#'
#' @import fs
#' @import ggplot2
#' @import ggthemes
#' @import dplyr
#' @import tidyr
#' @import readr
#' @import tibble
#' @import patchwork
#'
#'
#' @import DiffBind
#' @import ggplot2
#' @import dplyr
#' @import ChIPseeker
#' @import rtracklayer
#' @import org.Dm.eg.db
#' @import TxDb.Dmelanogaster.UCSC.dm6.ensGene
#'
#'
#' @name diffbind


#' get_hiseq_data
#'
#' retrive data from hiseq directory
#'
#' @param dir path to the hiseq
#' @param group attribute of the file, peak, bam, ...
#'
#' @export
get_hiseq_data <- function(dir, group = "peak") {
  sapply(dir, function(i) {
    if(is_hiseq_dir(i)) {
      pd <- read_hiseq(i)
      pd$args[[group]]
    }
  })
}



#' Convert narrowPeak to bed4
#'
#'
#' np2bed4
#'
#' Convert narrowPeak to BED4
#'
#' save column-7: signalValue
#'
#' @param np narrowPeak file
#' @param outdir path to save the bed file, if NULL, save to the same directory
#'   of np
#'
#' @export
np2bed4 <- function(np, outdir=NULL) {
  sapply(np, function(i) {
    if(is.null(outdir)) {
      outdir <- dirname(i)
    }
    bname <- gsub("_peaks.narrowPeak$", ".bed", basename(i))
    bed4 <- file.path(outdir, bname)

    if(file.exists(bed4)) {
      message(paste0("narrowPeak-to-bed4: file exists, ", bed4))
    } else {
      message(paste0("narrowPeak-to-bed4: ", bed4))
      gr <- hiseqr::read_narrowpeak(i)
      data.frame(seqnames=seqnames(gr),
                 starts=start(gr)-1,
                 ends=end(gr),
                 score = gr$signalValue) %>%
        write.table(file = bed4, quote = FALSE, sep = "\t", row.names = FALSE,
                    col.names = FALSE)
    }
    bed4
  }, simplify = TRUE)
}




#' ATACseqDiff
#'
#' DiffBind analysis for ATACseq results
#'
#' @param dirA path to directory of merge sample, B
#' @param dirB path to directory of merge sample, B
#' @param outdir path to the output, dirA.vs.dirB
#'
#' @export
atacseq_diff <- function(dirA, dirB, outdir = NULL) {
  #--Input dirs
  if(! is_hiseq_merge_dir(dirA) | ! is_hiseq_merge_dir(dirB)) {
    stop("Require hiseq_merge directory, exiting...")
  }

  # default: in dirA
  if(is.null(outdir) | outdir == "NULL") {
    outdir <- dirname(dirA)
  }
  # update outdir
  smp_name <- paste0(basename(dirB), ".vs.", basename(dirA))
  outdir   <- file.path(outdir, smp_name)

  if(! dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }

  #--Prepare design
  # required: 9 columns
  req_cols <- c("SampleID", "Tissue", "Factor", "Condition",
                "Treatment", "Replicate", "bamReads",
                "Peaks", "PeakCaller")

  # single replicate
  dirAx <- c(list_hiseq_single_dirs(dirA), dirA)
  dirBx <- c(list_hiseq_single_dirs(dirB), dirB)

  # peak_files
  peak_np_files <- get_hiseq_data(c(dirAx, dirBx), "peak")
  peak_files    <- np2bed4(peak_np_files)

  # bam_files
  bam_files <- get_hiseq_data(c(dirAx, dirBx), "bam_rmdup")

  # sample infor
  df <- data.frame(Peaks = peak_files,
                   bamReads = bam_files) %>%
    dplyr::mutate(
      SampleID = gsub(".bed", "", basename(Peaks)),
      Tissue   = stringr::str_extract(SampleID, "sh[A-Za-z0-9]+(_\\d+h|L\\d)?"),
      Factor   = "ATACseq",
      Condition  = Tissue,
      Treatment  = Tissue, # stringr::str_extract(SampleID, "\\d+h|L\\d"),
      Replicate  = stringr::str_extract(SampleID, "\\d$"),
      PeakCaller = "MACS"
    )
  df$Replicate[is.na(df$Replicate)] <- 0 # replicates

  #--Run DiffBind
  diffbind_hub(df, outdir)
}



#' Run DiffBind for peaks (groups)
#'
#' Standard
#'
#' @import DiffBind
#' @import ChIPseeker
#'
#' @param design data.frame or csv file, for the peak and bam files
#' @param outdir output dirs
#'
#' @export
diffbind_hub <- function(design, outdir) {
  #--Loading data
  pd <- prep_diffbind(design, outdir)

  #--Run DiffBind:
  if(file.exists(pd$out_files$diffbind_res)) {
    message(glue::glue("Loading DiffBind_res from: {pd$out_files$diffbind_res}"))
    rp <- readRDS(pd$out_files$diffbind_res)
  } else {
    message(glue::glue("Running DiffBind, save as: {pd$out_files$diffbind_res}"))
    rp <- DiffBind::dba(sampleSheet = pd$design, minOverlap = 4)
    rp <- DiffBind::dba.count(rp, bParallel = TRUE, summits = 100, minOverlap = 4)
    rp <- DiffBind::dba.normalize(rp)
    rp <- DiffBind::dba.contrast(rp, categories = DBA_CONDITION)
    rp <- DiffBind::dba.analyze(rp, method = DBA_DESEQ2, bParallel = TRUE)
    saveRDS(rp, pd$out_files$diffbind_res)
  }

  #--Visualization
  #--01.box
  message(glue::glue("Generating plot: {pd$out_files$box}"))
  tryCatch(
    {
      png(pd$out_files$box, width = 1024, height = 1024, res = 200)
      DiffBind::dba.plotBox(rp, contrast = 1)
      dev.off()
    },
    error = function(cond) { message("Failed, error ...") },
    warning = function(cond) { message("Failed, warning ...") }
  )

  #--02.PCA
  message(glue::glue("Generating plot: {pd$out_files$pca}"))
  tryCatch(
    {
      png(pd$out_files$pca, width = 1024, height = 1024, res = 200)
      DiffBind::dba.plotPCA(rp, contrast = 1, label = DBA_REPLICATE)
      dev.off()
    },
    error = function(cond) { message("Failed, error ...") },
    warning = function(cond) { message("Failed, warning ...") }
  )

  #--03.volcano
  message(glue::glue("Generating plot: {pd$out_files$volcano}"))
  tryCatch(
    {
      png(pd$out_files$volcano, width = 1024, height = 1024, res = 200)
      DiffBind::dba.plotVolcano(rp, contrast = 1)
      dev.off()
    },
    error = function(cond) { message("Failed, error ...") },
    warning = function(cond) { message("Failed, warning ...") }
  )

  #--04.heatmap
  message(glue::glue("Generating plot: {pd$out_files$heatmap}"))
  png(pd$out_files$heatmap, width = 1024, height = 1024, res = 200)
  DiffBind::dba.plotHeatmap(rp)
  dev.off()

  #--05.diff peaks
  message(glue::glue("Saving DiffBind peaks: {pd$out_files$diffbind_peaks}"))
  rp_db   <- DiffBind::dba.report(rp)
  rp_up   <- rp_db[rp_db$Fold > 0, ]
  rp_down <- rp_db[rp_db$Fold < 0, ]
  # labels
  label_up   <- paste0(basename(dirB), "_(", length(rp_up), ")")
  label_down <- paste0(basename(dirA), "_(", length(rp_down), ")")
  if(length(rp_db) > 0) {
    readr::write_delim(data.frame(rp_db),
                       pd$out_files$diffbind_peaks,
                       "\t",
                       col_names = FALSE)
    readr::write_delim(data.frame(rp_up),
                       pd$out_files$diffbind_peaks_up,
                       "\t",
                       col_names = FALSE)
    readr::write_delim(data.frame(rp_down),
                       pd$out_files$diffbind_peaks_down,
                       "\t",
                       col_names = FALSE)
  }

  #--06.Annotation peaks
  message(glue::glue("Generating annotation plots: {pd$out_files$anno_png}"))

  #--Annotation--
  #--Extra packages required
  message("000")
  pk <- setNames(list(rp_up, rp_down), c(label_up, label_down))
  pk <- purrr::discard(pk, function(i) {length(i) == 0})
  p_anno_1 <- lapply(
    pk, annotatePeak,
    TxDb = TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene,
    verbose = FALSE)

  # Annotation - all
  message("001")
  pa <- setNames(pd$design$Peaks, pd$design$SampleID)
  p_anno_2 <- lapply(
    pa, annotatePeak,
    TxDb = TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene,
    tssRegion=c(-3000, 3000),
    verbose=FALSE)

  # Combine annotations
  message("002")
  png(pd$out_files$anno_png, width = 1600, height = 800, res = 200)
  p1 <- plotAnnoBar(c(p_anno_1, p_anno_2))
  print(p1)
  dev.off()

  #--TSS
  png(pd$out_files$anno_tss_png, width = 1600, height = 800, res = 200)
  p2 <- plotDistToTSS(c(p_anno_1, p_anno_2))
  print(p2)
  dev.off()

  #--Finish
  message(glue::glue("DiffBind finished!"))
}




#' Prepare data for diffbind
#'
#' @param design data.frame or csv file, for the peak and bam files
#' @param outdir output dirs
#'
#' @export
prep_diffbind <- function(design, outdir) {
  #--Check df
  # required columns
  # peaks exists
  # bam exists
  if(is(design, "character")) {
    df <- read.csv(design, comment.char = "#")
  } else if(is(design, "data.frame")) {
    # load
    df <- design
  } else {
    stop(glue::glue("design=, expect data.frame or csv file, failed"))
  }

  # required: 9 columns
  req_cols <- c("SampleID", "Tissue", "Factor", "Condition",
                "Treatment", "Replicate", "bamReads",
                "Peaks", "PeakCaller")
  chk_req_cols <- sapply(req_cols, function(i) {
    ifelse(i %in% names(df), "ok", "missing")
  }, USE.NAMES = TRUE)
  chk_log = paste(paste(chk_req_cols, req_cols, sep = ": "), collapse = "\n")
  if(! all(req_cols %in% names(df))) {
    stop("design= unknown format:", chk_log)
  }

  # check peak file exists
  chk_peak <- sapply(df$Peaks, function(i) {
    ifelse(file.exists(i), "ok", "missing")
  }, USE.NAMES = TRUE)
  chk_peak_log = paste(paste(chk_peak, df$PeakCaller, sep = ": "),
                       collapse = "\n")
  if(! all(file.exists(df$Peaks))) {
    stop("Peaks, 1 or more file missing:", chk_peak_log)
  }

  # check bam file exists
  chk_bam <- sapply(df$bamReads, function(i) {
    ifelse(file.exists(i), "ok", "missing")
  }, USE.NAMES = TRUE)
  chk_bam_log = paste(paste(chk_bam, df$bamReads, sep = ": "),
                      collapse = "\n")
  if(! all(file.exists(df$bamReads))) {
    stop("bamReads, 1 or more file missing:", chk_bam_log)
  }

  # check Treatment (condition, >=2)
  if(length(unique(df$Treatment)) < 2) {
    stop("Treatment, require at least 2 groups")
  }

  # to-do
  # 3 replicates, ...

  #--Output files
  if(is(outdir, "character")) {
    if(! dir.exists(outdir)) {
      dir.create(outdir, showWarnings = TRUE, recursive = TRUE)
    }
  } else {
    stop("outdir, expect path, failed")
  }

  #--Default: output files
  f_names <- setNames(
    c("config.yaml",
      "00.diffbind_analysis_res.rds",
      "01.plot_box.png",
      "02.plot_pca.png",
      "03.plot_volcano.png",
      "04.plot_heatmap.png",
      "05.diffbind_peaks.bed",
      "06.diffbind_peaks_up.bed",
      "07.diffbind_peaks_down.bed",
      "08.diffbind_peaks.anno.png",
      "09.diffbind_peaks.anno.tss.png"),
    c("config",
      "diffbind_res",
      "box",
      "pca",
      "volcano",
      "heatmap",
      "diffbind_peaks",
      "diffbind_peaks_up",
      "diffbind_peaks_down",
      "anno_png",
      "anno_tss_png"))

  f_files <- sapply(f_names, function(i) {
    file.path(outdir, i)
  }, simplify = FALSE, USE.NAMES = TRUE)

  # return
  list(
    design    = df,
    out_files = f_files
  )
}

















