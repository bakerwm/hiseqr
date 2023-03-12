#' Functions for DESeq2 down-stream analysis
#'
#' Reading files, directories, pipeline, ...
#' Processing data
#' Prepare for plot
#' Plotting
#'
#' @name scatter_plot


#' @describeIn scatter_plot3 Create scatter plot, version3
#'
#' @param data A data frame for the full version of data, required
#' @param x character Column name show on x-axis
#' @param y character Column name chow on y-axis
#' @param labels character Labels to show in plot
#' @param add_label_point bool Add points for labels
#' @param show_abline bool Add ablines
#' @param trans character Function to convert num log10 or log2
#' @param xmin float Min on x-axis
#' @param xmax float Max x on x-axis
#' @param ymin float Min on y-axis
#' @param ymax float Max y on y-axis
#' @param title character Title, default: "Scatter plot"
#'
#' @example
#' df = readr::read_delim("transcripts_deseq2.fix.xls", "\t")
#' x  = names(df)[2]
#' y  = names(df)[3]
#' highlight_column = "sig"
#' highlight_values = c("up", "down")
#' label_column = "SYMBOL"
#' label_values = df %>% pull(sig) %>% head(10)
#'
#' p1 <- scatter_plot2(
#'  df, x, y,
#' highlight_column=highlight_column,
#' highlight_values=highlight_values,
#' label_column=label_column,
#' label_values=label_values,
#' density_point = FALSE,
#' point_color = "grey50")
#'
#'
#' @import ggplot2
#' @import ggrepel
#' @rlang
#'
#' @export
scatter_plot3 <- function(data, x, y, ...) {
  args <- prep_scatter3(data, x, y, ...)
  for(name in names(args)) {
    assign(name, args[[name]])
  }
  # --main: p0
  p0 <- df_point %>%
    ggplot(aes_string(as.name(x), as.name(y))) #, label = as.name(label_column)))  # !!! is_label
  # --density: p1
  if(isTRUE(density_point)) {
    p1 <- p0 +
      stat_density_2d(
        # aes(fill = after_stat(density)),
        # geom = "raster",
        aes(fill = ..density..^0.25),
        geom = "tile",
        contour = FALSE,
        show.legend = FALSE,
        n = 300) +
      scale_fill_gradient(low = "white", high = "#003366") +
      ggpubr::stat_cor()
    # stat_density2d(
    #   aes(fill = ..density..^0.25),
    #   show.legend = F,
    #   geom = "tile", contour = FALSE, n = 300)
  } else {
    p1 <- p0 +
      geom_point(size = point_size, color = point_color) +
      ggpubr::stat_cor()
  }
  # --lines:
  if(isTRUE(add_fc_lines)) {
    if(trans_axis %in% c("log2", "log10")) {
      tf <- ifelse(trans_axis == "log2", log2, log10)
      i <- c(-tf(2), tf(2))
    } else {
      i <- c(-2, 2)
    }
    p1 <- p1 +
      geom_abline(slope = 1, intercept = 0, color = "grey10") +
      geom_abline(
        slope = 1,
        intercept = i,
        color = "grey40", linetype = 2)
  }
  # --highlight points:
  if(isTRUE(hightligh_column %in% names(data))) {
    p1 <- p1 +
      geom_point(
        aes_string(as.name(x), as.name(y), color = highlight_column),
        data = df_highlight, shape = 16,
        size = point_size_highlight) # !!!! tmp
  }
  # --labels: add label in aes()
  if(isTRUE(label_column %in% names(data))) {
    p1 <- p1 +
      geom_point(
        aes_string(as.name(x), as.name(y)),
        data = df_label, shape = 16,
        size = point_size_highlight, # !!!! tmp
        color = "grey30") +
      ggrepel::geom_text_repel(
        aes_string(as.name(x), as.name(y), label = label_column), # !!! is_label
        data = df_label, inherit.aes = FALSE,
        color              = "black",
        # size               = 3,
        force              = .2,
        max.overlaps       = 80,
        direction          = "both",
        point.padding      = .2,
        box.padding        = .2,
        segment.color      = "grey20",
        segment.size       = .4,
        min.segment.length = 0)
  }
  # --themes:
  p1 <- p1 +
    scale_x_continuous(
      limits = c(xmin, xmax),
      breaks = scales::pretty_breaks()(xmin:xmax),
      labels = scales::pretty_breaks()(xmin:xmax),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(ymin, ymax),
      breaks = scales::pretty_breaks()(ymin:ymax),
      labels = scales::pretty_breaks()(ymin:ymax),
      expand = c(0, 0)
    )
  # --theme:
  p1 <- p1 +
    guides(fill = "none") + # remove density
    ggtitle(title) +
    xlab(xtitle) +
    ylab(ytitle) +
    labs(color = "Group") +
    theme_bw() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = .5),
      plot.title   = element_text(color = "black", hjust = .5, size = 14),
      panel.grid   = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.line    = element_line(color = "black", size = .5),
      axis.ticks   = element_line(color = "black", size = .5),
      axis.text    = element_text(color = "black", size = 10),
      axis.title   = element_text(color = "black", size = 12),
      axis.ticks.length = unit(.2, "cm"),
      aspect.ratio = 0.8
    )
  p1
}



#' prep_scatter3
#'
#' for function scatter_plot3()
prep_scatter3 <- function(data, x, y, ...) {
  args <- list(
    "label_column"     = "label",
    "label_values"     = NULL,
    "hightligh_column" = "sig",
    "highlight_values" = c("up", "down"),
    "point_size_highlight" = 2,
    "add_fc_lines"  = TRUE,
    "density_point" = FALSE,
    "point_size"    = 0.5,
    "point_color"   = "grey50",
    "trans_axis"    = "log10",
    "trans_axis_x"  = "log10",
    "trans_axis_y"  = "log10",
    "xmin"  = 0,
    "xmax"  = 6,
    "ymin"  = 0,
    "ymax"  = 6,
    "title" = "scatter"
  )
  dots <- rlang::list2(...)
  args <- purrr::list_modify(args, !!!dots)
  for(name in names(args)) {
    assign(name, args[[name]])
  }
  # --arguments: data
  if(! inherits(data, "data.frame")) {
    on.exit(message("data not data.frame"), add = TRUE)
  }
  if(! inherits(c(x, y), "character")) {
    on.exit(message("x, y not character"), add = TRUE)
  }
  x <- x[1] # update
  y <- y[1] # update
  if(! all(c(x, y) %in% names(data))) {
    on.exit(message("x, y not found in data.frame"), add = TRUE)
  }
  # --arguments: x, y axis, log-trans
  if(trans_axis %in% c("log2", "log10")) {
    trans_func <- ifelse(trans_axis == "log2", log2, log10)
    data <- dplyr::mutate(data, across(all_of(c(x, y)), ~ trans_func(.x + 1)))
    tx   <- gsub("log", "", trans_axis)
    xtitle <- bquote(.(x) ~ "[" * log[(tx)] ~ "(rpm+1)]")
    ytitle <- bquote(.(y) ~ "[" * log[(tx)] ~ "(rpm+1)]")
  } else {
    xtitle <- x
    ytitle <- y
  }
  # --arguments:
  df_point <- data
  # --highlight:
  if(isTRUE(hightligh_column %in% names(data))) {
    h1 <- dplyr::pull(data, hightligh_column)
    highlight_values <- highlight_values[highlight_values %in% h1]
    df_highlight <- data %>%
      dplyr::filter(!!as.symbol(hightligh_column) %in% highlight_values) %>%
      dplyr::mutate(!!as.symbol(hightligh_column) := factor(
        (!! as.symbol(hightligh_column)), highlight_values
      ))
    # remove from data
    df_point <- data %>%
      dplyr::filter(! (!!as.symbol(highlight_column)) %in% highlight_values)
  } else {
    highlight_values <- NA
    hightligh_column <- NA
    df_highlight     <- NA
  }
  # --labels: add label column
  if(isTRUE(label_column %in% names(data))) {
    h2 <- dplyr::pull(data, label_column)
    label_values <- label_values[label_values %in% h2]
    # update label
    data <- data %>%
      dplyr::mutate(
        label    = !!(as.name(label_column)),
        is_label = ifelse(label %in% label_values, label, NA))
    # update label
    df_label <- data %>%
      dplyr::filter(label %in% label_values) %>%
      dplyr::filter(! is.na(!! as.symbol(label_column))) # NA values
    # remove from data
    df_point <- data %>%
      dplyr::filter(! label %in% label_values)
  } else {
    label_column <- NA
    label_values <- NA
    df_label     <- NA
  }
  # --Updated values
  args2 <- list(
    data   = data,
    xtitle = xtitle,
    ytitle = ytitle,
    label_values = label_values,
    df_label     = df_label,
    df_highlight = df_highlight,
    df_point     = df_point
  )
  # --Output:
  purrr::list_modify(args, !!!args2)
}

