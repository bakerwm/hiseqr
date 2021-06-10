
#' fastqc_report
#'
#' This function parsing *_fastqc.zip files from input path.
#' It assumes that input path is a directory or a single *_fastqc.zip
#' file.
#'
#' @description Create html report for fastqc output
#' @param indir Path to the input directory, Default is the
#'   current working directory.
#' @param outdir path to the result file prefix (e.g., path/to/qc-result).
#'   Don't add the file extension.
#' @param template a character vector specifying the path to an Rmd template.
#'  file.
#' @param preview logical value. If TRUE, shows a preview of the report.
#' @examples
#' \donotrun{
#' # Demo
#' qc.path <- system.file("fastqc_results", package = "fastqcr")
#'
#' fastqc_report(qc.path, result.file = "demo")
#' }
#'
#' @import fastqcr
#' @import rmarkdown
#'
#' @return A list of paths of zip files
#' @export
fastqc_report <- function (indir, outdir, preview = FALSE) {
  ## input
  indir  <- normalizePath(indir)
  outdir <- normalizePath(outdir)
  report_template <- system.file("qc", "fastqc_report.Rmd",
                                 package = "hiseqr")
  if(! dir.exists(outdir)){
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }
  output_html <- file.path(outdir, "fastqc_report.html")
  rmarkdown::render(input = report_template,
                    output_file = output_html,
                    params = list(input_dir = indir))
  if (preview) {
    utils::browseURL(output_html)
  }
}




#' fastqcFiles
#'
#' This function parsing *_fastqc.zip files from input path.
#' It assumes that input path is a directory or a single *_fastqc.zip
#' file.
#'
#' @param x Path to the input directory
#' @return A list of paths of zip files
#' @export
fastqc_files <- function(x) {
  stopifnot(length(x) == 1)
  if (x == ".")
    x <- getwd()

  # input *.zip, dir/
  if(dir.exists(x)) {
    # search *.zip files in x/
    qcFiles <- list.files(path = x,
                          pattern = "*_fastqc.zip$",
                          all.files = TRUE,
                          full.names = TRUE,
                          recursive = FALSE)
  } else if(file.exists(x)) {
    # check x is *.zip
    qcFiles <- x[endsWith(x, "_fastqc.zip")]
  } else {
    # x is the prefix of *.zip file
    qcFiles <- list.files(path = dirname(x),
                          pattern = "*_fastqc.zip",
                          all.files = TRUE,
                          full.names = TRUE,
                          recursive = FALSE)
  }

  # check output
  if(length(qcFiles) == 0){
    return(NULL)
  }

  return(qcFiles)
}




#' fastqcFiles
#'
#' This function parsing *_fastqc.zip files from input path.
#' It assumes that input path is a directory or a single *_fastqc.zip
#' file.
#'
#' @param x Path to the input directory
#' @return A list of paths of zip files
#' @export
fastqc_files <- function(x) {
  stopifnot(length(x) == 1)
  if (x == ".")
    x <- getwd()

  # input *.zip, dir/
  if(dir.exists(x)) {
    # search *.zip files in x/
    qcFiles <- list.files(path = x,
                          pattern = "*_fastqc.zip$",
                          all.files = TRUE,
                          full.names = TRUE,
                          recursive = FALSE)
  } else if(file.exists(x)) {
    # check x is *.zip
    qcFiles <- x[endsWith(x, "_fastqc.zip")]
  } else {
    # x is the prefix of *.zip file
    qcFiles <- list.files(path = dirname(x),
                          pattern = "*_fastqc.zip",
                          all.files = TRUE,
                          full.names = TRUE,
                          recursive = FALSE)
  }

  # check output
  if(length(qcFiles) == 0){
    return(NULL)
  }

  return(qcFiles)
}




#' @param module string, base quality, base content, length distribution,
#' @export
fastqc_plot <- function(r1_path, r2_path = NULL, module = "base quality") {
  # read1
  r1_qc <- fastqcr::qc_read(r1_path)

  # read2
  if(is.null(r2_path)) {
    r2_qc <- NULL
  } else {
    r2_qc <- fastqcr::qc_read(r2_path)
  }

  # determine the function
  p1 = p2 = NULL
  r1_df = r2_df = NULL
  if(module %in% c(1, "base quality")) {
    plot_function <- plot_base_quality
    r1_df <- r1_qc$per_base_sequence_quality
    # p1    <- plot_base_quality(r1_df)
    if(! is.null(r2_qc)) {
      r2_df <- r2_qc$per_base_sequence_quality
      # p2    <- plot_base_quality(r2_df)
    }
  } else if(module %in% c(2, "base content")) {
    plot_function <- plot_base_content
    r1_df <- r1_qc$per_base_sequence_content
    # p1    <- plot_base_content(r1_df)
    if(! is.null(r2_qc)) {
      r2_df <- r2_qc$per_base_sequence_content
      # p2    <- plot_base_content(r2_df)
    }
  } else if(module %in% c(3, "read length")) {
    plot_function <- plot_len_dis
    r1_df <- r1_qc$sequence_length_distribution
    # p1    <- plot_len_dis(r1_df)
    if(! is.null(r2_qc)) {
      r2_df <- r2_qc$sequence_length_distribution
      # p2    <- plot_len_dis(r2_df)
    }
  } else {
    return(NULL)
  }

  # make plot
  fname <- gsub("_fastqc.zip$", "", basename(r1_path))
  p1 <- plot_function(r1_df, fname)
  p2 <- plot_function(r2_df, "read2")

  # # for modules
  # p1 <- plot_function(r1_df, "read1")
  # p2 <- plot_function(r2_df, "read2")

  p <- cowplot::plot_grid(p1, p2, nrow = 1)

  return(p)
}



#' Per base sequence quality plot
#' @import ggplot2
#' @import dplyr
#' @import scales
#'
#' @export
plot_base_quality <- function(df, title = "Per base quality") {
  required_cols <- c("Base", "Mean", "Median", "Lower Quartile",
                     "Upper Quartile", "10th Percentile", "90th Percentile")
  if(! inherits(df, "data.frame")) {
    return(NULL)
  } else if(! all(required_cols %in% colnames(df))) {
    return(NULL)
  }

  # format data.frame
  df2 <- df %>%
    dplyr::rename(
      y50 = Median,
      y25 = `Lower Quartile`,
      y75 = `Upper Quartile`,
      y10 = `10th Percentile`,
      y90 = `90th Percentile`
    ) %>%
    dplyr::mutate(position = row_number()) %>%
    dplyr::mutate(Base = factor(Base, levels = Base))

  # global vars
  df_nrow <- nrow(df2)

  # determin x labels
  # x_breaks <- c(1:8, seq(9, df_nrow, by = 3), df_nrow)
  # x_labels <- df$Base[x_breaks]
  x_breaks <- df2$position
  x_labels <- df2$Base

  # determine background of x axis
  xmin <- seq(0, df_nrow, by = 2) + 0.5
  xmin <- xmin[xmin < df_nrow]
  xmax <- xmin + 1

  # make basic plot
  # background colors
  p <- ggplot(df2)  +
    annotate("rect",
             xmin = 0,
             xmax = df_nrow + 1,
             ymin = c(0, 20, 28),
             ymax = c(20, 28, 41),
             fill = c("red", "orange", "green2"),
             alpha = .2) +
    annotate("rect",
             xmin = xmin,
             xmax = xmax,
             ymin = 0,
             ymax = 41,
             fill = "grey60",
             alpha = .2)

  # add boxplot + line
  p <- p +
    geom_boxplot(
      aes(
        x = position,
        ymin = y10,
        lower = y25,
        middle = y50,
        upper = y75,
        ymax = y90,
        group = position
      ),
      stat = "identity",
      fill = "yellow2",
      color = "grey30",
      outlier.color = "grey50",
      size = .3
    ) +
    geom_line(aes(position, Mean),
              color = "blue", size = 0.5)

  # modify axis
  p <- p +
    scale_x_continuous(
      expand = c(0, 0),
      limits = c(0, df_nrow + 1),
      breaks = x_breaks,
      labels = x_labels,
      guide = guide_axis(check.overlap = TRUE)
    ) +
    scale_y_continuous(
      expand = c(0, 0),
      breaks = seq(0, 42, by = 2),
      labels = seq(0, 42, by = 2),
    ) +
    ggtitle(title) +
    xlab("position in read (bp)") +
    ylab(NULL) +
    theme_classic() +
    theme(
      axis.text.x = element_text(color = "black")
    )

  return(p)
}



#' @import ggplot2
#' @import dplyr
#' @import scales
#'
#' @export
plot_base_content <- function(df, title = "Per base content") {
  required_cols <- c("Base", "A", "C", "G", "T")
  if(! inherits(df, "data.frame")) {
    return(NULL)
  } else if(! all(required_cols %in% colnames(df))) {
    return(NULL)
  }

  # global vars
  df_nrow <- nrow(df)

  # format data.frame
  df2 <- df %>%
    dplyr::mutate(position = row_number()) %>%
    tidyr::gather(key = "base", value = "score", 2:5) %>%
    dplyr::mutate(base = factor(base, levels = c("A", "C", "G", "T")))

  # determin x labels
  # x_breaks <- c(1:8, seq(9, df_nrow, by = 3), df_nrow)
  # x_labels <- df$Base[x_breaks]
  x_breaks <- seq_len(df_nrow)
  x_labels <- df$Base

  # determine background of x axis
  xmin <- seq(0, df_nrow, by = 2) + 0.5
  xmin <- xmin[xmin < df_nrow]
  xmax <- xmin + 1

  # make basic plot
  # background colors
  p <- ggplot(df2, aes(position, score, color = base)) +
    annotate("rect",
             xmin = xmin,
             xmax = xmax,
             ymin = 0,
             ymax = 100,
             fill = "grey60",
             alpha = .2)

  # add lines
  p <- p +
    geom_line(size = .6)

  # modify axis, layers
  p <- p +
    scale_color_manual(values = c("green", "blue", "grey20", "red")) +
    scale_x_continuous(
      expand = c(0, 0),
      limits = c(0, df_nrow + 1),
      breaks = x_breaks,
      labels = x_labels,
      guide = guide_axis(check.overlap = TRUE)
    ) +
    scale_y_continuous(
      expand = c(0, 0),
      n.breaks = 6
    ) +
    ggtitle(title) +
    xlab("position in read (bp)") +
    ylab(NULL) +
    guides(color = guide_legend(title = NULL)) +
    theme_classic() +
    theme(
      panel.grid.major.y = element_line(color = "grey60", size = .3),
      legend.background = element_rect(fill = "white"),
      legend.position = c(.9, .8)
    )

  return(p)
}



# df <- qc$sequence_length_distribution

#' @import ggplot2
#' @import dplyr
#' @import scales
#'
#' @export
plot_len_dis <- function(df, title = "Distribution of Lengths") {
  required_cols <- c("Length", "Count")
  if(! inherits(df, "data.frame")) {
    return(NULL)
  } else if(! all(required_cols %in% colnames(df))) {
    return(NULL)
  }

  # for 0, 1 row
  if(nrow(df) < 2) {
    df <- df %>%
      tibble::add_row(Length = df$Length - 1, Count = 0, .before = 1) %>%
      tibble::add_row(Length = df$Length + 1, Count = 0)
  }

  # add position
  df <- df %>%
    dplyr::mutate(position = row_number())

  df_nrow <- nrow(df)

  # determine background of x axis
  xmin <- seq(0, df_nrow, by = 2) + 0.5
  xmin <- xmin[xmin < df_nrow]
  xmax <- xmin + 1

  # determine x labels
  # x_breaks <- scales::extended_breaks(n = min(5, df_nrow))(seq_len(df_nrow))
  # x_breaks <- unique(sort(c(1, x_breaks, df_nrow)))
  # x_breaks <- x_breaks[x_breaks %in% c(1:df_nrow)]
  # x_labels <- df$Length[x_breaks]
  x_breaks <- seq_len(df_nrow)
  x_labels <- df$Length

  ymax <- max(df$Count)

  # basic plot
  p <- df %>%
    ggplot(aes(position, Count)) +
    annotate("rect",
             xmin = xmin,
             xmax = xmax,
             ymin = 0,
             ymax = ymax,
             fill = "grey60",
             alpha = .2)

  # point + line
  p <- p +
    geom_point(color = "red", size = 1) +
    geom_line(color = "red", size = .5)

  # layers
  p <- p +
    scale_x_continuous(
      expand = c(0, 0),
      limits = c(0, df_nrow + 1),
      breaks = x_breaks,
      labels = x_labels,
      guide = guide_axis(check.overlap = TRUE)
    ) +
    ggtitle(title) +
    xlab("Sequence Length (bp)") +
    ylab(NULL) +
    theme_classic() +
    theme(
      panel.grid.major.y = element_line(color = "grey60", size = .3)
    )

  return(p)
}

