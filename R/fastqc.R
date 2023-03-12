

#' @param r1 List output of `qc_read()` of read1
#' @param r2 List output of `qc_read()` of read2
#' @param module Character name of the modules
#'
#' @export
plot_fastqc <- function(r1, r2 = NULL, module = "Per base sequence content", ...) {
  # choose plot function
  if(.is_valid_fastqc_module(module)) {
    func <- .plot_func(module)
  } else {
    warning(glue::glue("unknown module: {module}"))
    return(NULL)
  }
  # fix title
  args <- rlang::list2(...)
  title <- args[["title"]]
  if(is.null(title)) {
    title = module
  }
  # read1:
  p1 <- func(r1)
  # update title
  p1$labels$title <- paste0(p1$labels$title, " (read1)")
  # read2:
  if(inherits(r2, "list")) {
    p2 <- func(r2)
    # update title
    p2$labels$title <- paste0(p2$labels$title, " (read2)")
  } else {
    p2 <- patchwork::plot_spacer()
  }
  # combine
  p1 + p2 +
    patchwork::plot_annotation(title = title)
}



#' @param modules Character name of the fastqc modules
.plot_func <- function(module) {
  switch(module,
         "Per base sequence quality" = plot_base_seq_quality,
         "Per tile sequence quality" = plot_tile_seq_quality,
         "Per sequence quality scores" = plot_seq_quality_scores,
         "Per base sequence content" = plot_base_content,
         "Per sequence GC content" = NULL,
         "Per base N content" = NULL,
         "Sequence Length Distribution" = plot_seq_length_distribution,
         "Sequence Duplication Levels"  = NULL,
         "Overrepresented sequences"    = NULL,
         "Adapter Content" = NULL
  )
}


#' @param modules Character name of the fastqc modules
.is_valid_fastqc_module <- function(modules) {
  m <- c(
    "Basic Statistics",
    "Per base sequence quality",
    "Per tile sequence quality",
    "Per sequence quality scores",
    "Per base sequence content",
    "Per sequence GC content",
    "Per base N content",
    "Sequence Length Distribution",
    "Sequence Duplication Levels",
    "Overrepresented sequences",
    "Adapter Content"
  )
  if(isTRUE(modules)) {
    modules <- m
  }
  all(sapply(modules, function(i) i %in% m))
}



#' @param x Character path to the json file of fastqc/falco
#' @param module Character name of the fastqc module
#'
#' @export
qc_read <- function(x, module = TRUE) {
  qc <- jsonlite::read_json(x) # modules
  if(!.is_valid_fastqc_module(module)) {
    mx <- module[! sapply(module, .is_valid_fastqc_module)]
    msg <- paste(mx, collapse = ",")
    stop(paste0("unknown modules found:", msg))
  }
  # software, version
  # ... table to data.frame
  df <- lapply(qc, function(i) {
    if(inherits(i, "list")) {
      if("table" %in% names(i)) {
        dt <- lapply(i[["table"]], unlist) %>%
          as.data.frame
        names(dt) <- names(i[["table"]])
        dt
      }
    } else {
      i
    }
  })
  # remove first two elements: software, version
  df <- within(df, rm(list=c("software", "version")))
  # selection
  if(inherits(module, "character")) {
    out <- sapply(module, function(i) {
      df[[i]]
    }, simplify = F, USE.NAMES = TRUE)
    if(length(module) == 1) {
      out[[1]]
    } else {
      out
    }
  } else {
    df
  }
}


#' @param x Character path to the json file of fastqc/falco
#' @param module Character name of the fastqc module
#'
#' @export
qc_status <- function(x, module = TRUE) {
  qc <- jsonlite::read_json(x) # modules
  if(!.is_valid_fastqc_module(module)) {
    mx <- module[! sapply(module, .is_valid_fastqc_module)]
    msg <- paste(mx, collapse = ",")
    stop(paste0("unknown modules found:", msg))
  }
  df <- lapply(qc, function(i) {
    if(inherits(i, "list")) {
      i[["status"]]
    } else {
      i
    }
  })
  # remove first two elements: software, version
  df <- within(df, rm(list=c("software", "version")))
  # selection
  if(inherits(module, "character")) {
    sapply(module, function(i) {
      df[[i]]
    })
  } else {
    df
  }
}


#' @import ggplot2
#' @import dplyr
#' @import scales
#'
#' @param x List output of `qc_read(x, module = TRUE)`
#' @export
plot_base_content <- function(x, ...) {
  module <- "Per base sequence content"
  if(inherits(x, "list")) {
    df <- x[[module]]
  } else {
    warning("unknown data x")
    return(NULL)
  }
  # bar-plot:
  df1 <- df %>%
    dplyr::mutate(position = row_number()) %>%
    tidyr::pivot_longer(2:5, names_to = "base", values_to = "score") %>%
    dplyr::mutate(base = factor(base, levels = c("A", "C", "G", "T")),
                  score = as.numeric(score))
  # # fix breaks
  x_breaks <- seq_len(nrow(df))
  x_labels <- df$Base
  x_min <- seq(0, nrow(df), by = 2) + 0.5
  x_min <- x_min[x_min < nrow(df)]
  x_max <- x_min + 1
  ## add background ##
  ggplot(df1, aes(position, score, color = base)) +
    annotate("rect", xmin = x_min, xmax = x_max, ymin = 0, ymax = 100,
             fill = "grey60", alpha = .2) +
    geom_line(size = .6) +
    scale_color_manual(values = c("green", "blue", "grey20", "red")) +
    scale_x_continuous(
      expand = c(0, 0),
      limits = c(0, nrow(df) + 1),
      breaks = x_breaks,
      labels = x_labels,
      guide = guide_axis(check.overlap = TRUE)
    ) +
    scale_y_continuous(
      expand = c(0, 0),
      n.breaks = 6
    ) +
    ggtitle(module) +
    xlab("Position in read (bp)") +
    ylab("Nucleotide frequency (%)") +
    guides(color = guide_legend(title = NULL)) +
    theme_classic() +
    theme(
      panel.grid.major.y = element_line(color = "grey60", size = .3),
      legend.background = element_rect(fill = "white"),
      legend.position = c(.9, .8),
      axis.text = element_text(size = 12, color = "grey20")
    )
}


#' @import ggplot2
#' @import dplyr
#' @import scales
#'
#' @param x List output of `qc_read(x, module = TRUE)`
#' @export
plot_seq_length_distribution <- function(x, ...) {
  module <- "Sequence Length Distribution"
  if(inherits(x, "list")) {
    df <- x[[module]] %>%
      dplyr::mutate(across(everything(), as.numeric))
  } else {
    warning("unknown data x")
    return(NULL)
  }
  # for 0, 1 row
  if(nrow(df) < 2) {
    df <- df %>%
      tibble::add_row(Length = df$Length - 1, Count = 0, .before = 1) %>%
      tibble::add_row(Length = df$Length + 1, Count = 0)
  }
  # bar-plot:
  df1 <- df %>%
    dplyr::mutate(position = row_number(),
                  Count = as.numeric(Count))
  # fix breaks
  x_breaks <- seq_len(nrow(df))
  x_labels <- df$Base
  x_min <- seq(0, nrow(df), by = 2) + 0.5
  x_min <- x_min[x_min < nrow(df)]
  x_max <- x_min + 1
  y_max <- max(df1$Count)
  # plot
  ggplot(df1, aes(position, Count)) +
    annotate("rect", xmin = x_min, xmax = x_max, ymin = 0, ymax = y_max,
             fill = "grey60", alpha = .2) +
    geom_point(color = "red", size = 1) +
    geom_line(color = "red", size = .5) +
    scale_x_continuous(
      expand = c(0, 0),
      limits = c(0, nrow(df) + 1),
      breaks = x_breaks,
      labels = x_labels,
      guide  = guide_axis(check.overlap = TRUE)
    ) +
    ggtitle(module) +
    xlab("Sequence Length (bp)") +
    ylab("Number of reads") +
    theme_classic() +
    theme(
      panel.grid.major.y = element_line(color = "grey60", size = .3)
    )
}


#' @import ggplot2
#' @import dplyr
#' @import scales
#'
#' @param x List output of `qc_read(x, module = TRUE)`
#' @export
plot_base_seq_quality <- function(x, ...) {
  module <- "Per base sequence quality"
  if(inherits(x, "list")) {
    df <- x[[module]]
  } else {
    warning("unknown data x")
    return(NULL)
  }
  if(is.null(df)) {
    return(NULL)
  }
  # fix data.frame
  df1 <- df %>%
    dplyr::rename(
      y50 = Median,
      y25 = `Lower Quartile`,
      y75 = `Upper Quartile`,
      y10 = `10th Percentile`,
      y90 = `90th Percentile`
    ) %>%
    dplyr::mutate(across(-1, as.numeric)) %>%
    dplyr::mutate(position = row_number()) %>%
    dplyr::mutate(Base = factor(Base, levels = Base))
  # fix breaks
  x_breaks <- df1$position
  x_labels <- df1$Base
  x_min <- seq(0, nrow(df1), by = 2) + 0.5
  x_min <- x_min[x_min < nrow(df1)]
  x_max <- x_min + 1
  # basic plot
  ggplot(df1) +
    annotate("rect", xmin = 0, xmax = nrow(df1) + 1, ymin = c(0, 20, 28),
             ymax = c(20, 28, 41), fill = c("red", "orange", "green2"),
             alpha = .2) +
    annotate("rect", xmin = x_min, xmax = x_max, ymin = 0, ymax = 41,
             fill = "grey60", alpha = .2) +
    geom_boxplot(
      aes(
        x = position, ymin = y10, lower = y25, middle = y50, upper = y75,
        ymax = y90, group = position
      ),
      stat = "identity", fill = "yellow2", color = "grey30",
      outlier.color = "grey50", size = .3
    ) +
    geom_line(aes(position, Mean), color = "blue", size = 0.5) +
    scale_x_continuous(
      expand = c(0, 0),
      limits = c(0, nrow(df1) + 1),
      breaks = x_breaks,
      labels = x_labels,
      guide = guide_axis(check.overlap = TRUE)
    ) +
    scale_y_continuous(
      expand = c(0, 0),
      breaks = seq(0, 42, by = 2),
      labels = seq(0, 42, by = 2),
      guide = guide_axis(check.overlap = TRUE)
    ) +
    ggtitle(module) +
    xlab("position in read (bp)") +
    ylab(NULL) +
    theme_classic() +
    theme(
      axis.text.x = element_text(color = "black")
    )
}


#' @import ggplot2
#' @import dplyr
#' @import scales
#'
#' @param x List output of `qc_read(x, module = TRUE)`
#' @export
plot_tile_seq_quality <- function(x, ...) {
  module <- "Per tile sequence quality"
  if(inherits(x, "list")) {
    df <- x[[module]]
  } else {
    warning("unknown data x")
    return(NULL)
  }
  # data
  df1 <- df %>%
    dplyr::mutate(Mean = as.numeric(Mean))
  # plot
  ggplot(df1, aes(x = Base, y = Tile, fill = Mean))+
    ggplot2::geom_tile() +
    labs(title = "Per tile sequence quality",
         subtitle = "Quality per tile",
         x = "Position in read (bp)",
         y = NULL)+
    theme_minimal()
}


#' @import ggplot2
#' @import dplyr
#' @import scales
#'
#' @param x List output of `qc_read(x, module = TRUE)`
#' @export
plot_seq_quality_scores <- function(x, ...) {
  module <- "Per sequence quality scores"
  if(inherits(x, "list")) {
    df <- x[[module]]
  } else {
    warning("unknown data x")
    return(NULL)
  }
  # fix data
  df1 <- dplyr::mutate(df, across(everything(), as.numeric))
  #
  x_breaks <- df1$Quality
  x_labels <- df1$Quality
  x_min <- seq(min(df1$Quality), max(df1$Quality), by = 2) + 0.5
  x_min <- x_min[x_min < max(df1$Quality)]
  x_max <- x_min + 1
  # plot
  ggplot(df1, aes(Quality, Count)) +
    annotate("rect", xmin = x_min, xmax = x_max, ymin = 0, ymax = y_max,
             fill = "grey60", alpha = .2) +
    geom_point(color = "red", size = 1) +
    geom_line(color = "red", size = .5) +
    scale_x_continuous(
      expand = c(0, 0),
      limits = c(min(df1$Quality), max(df1$Quality) + 1),
      breaks = x_breaks,
      labels = x_labels,
      guide  = guide_axis(check.overlap = TRUE)
    ) +
    ggtitle(module) +
    xlab("Mean sequence quality (Phred Score)") +
    ylab(NULL) +
    theme_classic() +
    theme(
      panel.grid.major.y = element_line(color = "grey60", size = .3)
    )
}


