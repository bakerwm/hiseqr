#' Functions for HiSeq plots
#'
#' Including bar_plot, scatter_plot, venn_plot, ma_plot,  box_plot, ...
#'
#'
#'
#'
#' @name plot


#--Global: General functions ---------------------------------------------------

#' @describeIn bar_plot Make barplot
#'
#' @param data data.frame require id, count column
#' @param x character Column on x-axis
#' @param y character Column on y-axis
#' @param direction character Options: vertical, horizontal, default: vertical
#' @param ... extra arguments for plots
#'
#' @import readr
#' @import dplyr
#' @import ggplot2
#'
#' @export
bar_plot <- function(data = NULL, x, y, label = y, group = NULL, fill = group,
                     direction = "vertical", position = "stack", ...) {
  # data <- ggplot2::fortify(data)
  if(! is.character(x) && is.character(y)) {
    stop("`x`, `y` not character")
  }
  x <- x[1]
  y <- y[1]
  if(! all(c(x, y) %in% names(data))) {
    stop("`x`, `y` not found in data.frame")
  }
  if(is(group, "character") && is(fill, "character")) {
    out <- ggplot(data, aes_(as.name(x), as.name(y), group = as.name(group),
                             fill = as.name(fill)))
  } else if(is(group, "character")) {
    out <- ggplot(data, aes_(as.name(x), as.name(y), group = as.name(group)))
  } else if(is(fill, "character")) {
    out <- ggplot(data, aes_(as.name(x), as.name(y), fill = as.name(fill)))
  } else {
    out <- ggplot(data, aes_(as.name(x), as.name(y)))
  }
  out <- out +
    geom_col(position = position) +
    theme_bw()
  if(direction == "vertical") {
    out <- out +
      theme(
        axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
        panel.grid  = element_blank(),
        plot.title  = element_text(hjust = 0.5)
      )
    if(is.character(label)) {
      out <- out +
        geom_text(aes_(label = as.name(label)),
                  vjust = -.1,
                  color = "black",
                  angle = 0)
    }
  } else {
    out <- out +
      scale_x_continuous(position = "top") +
      scale_y_discrete(limits = rev) +
      theme(
        axis.text  = element_text(color = "grey20"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
      )
    if(is.character(label)) {
      out <- out +
        geom_text(aes_(label = as.name(label)),
                  hjust = 0,
                  color = "black",
                  angle = 0)
    }
  }
  out
}



#' length distribution
#'
#' @describeIn barplot for length distribution
#'
#'
#' @export
bar_plot_lendist <- function(data, x = "length", y = "count",
                             range = c(18:40)) {
  # columns
  if(all(c("length", "count", "id") %in% names(data))) {
    data <- data %>%
      dplyr::filter(length %in% range)
    if("strand" %in% names(data)) {
      # v2: strand
      out <- data %>%
        dplyr::mutate(count = ifelse(strand == "+", count, -count)) %>%
        bar_plot(x = x, y = y, label = NA, group = "strand",
                 fill = "strand") +
        scale_fill_manual(values = c("red", "blue")) +
        scale_y_continuous(labels = abs)
    } else {
      # v3: no-strand
      out <- data %>%
        bar_plot(x = x, y = y, label = NA)
    }
  } else {
    stop("file format unknown")
  }
  out
}










#' @describeIn fragsize_plot line plot for fragment size
#'
#' @param data data.frame for plotting
#' @param xmin int Min of the fragsize, default: 0
#' @param xmax int Max of the fragsize, defualt: 1000
#' @param log10_y bool Convert y-axis to log10 scale
#' @param position character, fill percentage, others show the count
#'
#' @import readr
#' @import dplyr
#' @import ggplot2
#'
#' @export
fragsize_plot <- function(data, xmin = 0, xmax = 1000,
                          log10_y = FALSE, position = "fill") {
  data <- ggplot2::fortify(data)
  if(! all(c("id", "length", "count") %in% names(data))) {
    on.exit(message("id, length, count columns missing"), add = TRUE)
  }
  data <- data %>%
    group_by(id) %>%
    mutate(frac = count / sum(count))
  y    <- ifelse(position == "fill", "frac", "count")
  xmin <- ifelse(xmin < 0, 0, xmin)
  xmax <- ifelse(xmax > 1000, 1000, xmax)
  out  <- ggplot(data, aes_(as.name("length"), as.name(y), color = as.name("id"))) +
    geom_line(size = .5) +
    scale_x_continuous(limits = c(xmin, xmax),
                       expand = c(0, 0)) +
    theme_bw() +
    theme(
      axis.text  = element_text(color = "grey20"),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  # log10 of y axis
  if(isTRUE(log10_y)) {
    out <- out +
      scale_y_log10()
  }
  out
}




#' @describeIn venn_plot Create venn plot using ggVennDiagram
#'
#' @param ... list, unique ids for groups
#' @param names character, the names for each group
#'
#' @import ggVennDiagram
#' @import ggplot2
#'
#' @export
venn_plot <- function(..., names = NULL){
  dots <- rlang::list2(...)
  dots <- rlang::squash_if(dots, vctrs::vec_is_list)
  dots <- purrr::discard(dots, is.null)
  is_data_frame <- purrr::map_lgl(dots, is_vector)

  # check category names
  if(is.null(names)) {
    names <- names(x)
  }

  # further
  if(is.null(names)) {
    names <- letters[seq_len(length(x))]
  }
  # main
  ggVennDiagram::ggVennDiagram(
    x,
    label          = "count",
    label_alpha    = 0,
    category.names = names
  ) +
    ggplot2::guides(fill = FALSE) +
    ggplot2::scale_fill_gradient(low = "white", high = "firebrick")
}



#' @describeIn volcano_plot Create volcano plot
#'
#' @param data A data frame for the full version of data
#' @param x character Name of the column, show on X-axis, eg: fold-change
#' @param y character Name of the column, show on Y-axis, eg: p-value, padj
#' @param xintercept numeric Dot lines on x-axis default: c(0.5, 2)
#' @param yintercept numeric Dot lines on y-axis, default: 0.05
#' @param labels characters Text labels on the plot
#' @param show_sig bool Show significant points
#' @param title character Title, default: "Volcano plot"
#'
#' @import ggplot2
#' @import ggrepel
#'
#' @export
volcano_plot <- function(data, x, y,
                         labels     = NULL,
                         xintercept = c(0.5, 2),
                         yintercept = 0.05,
                         show_sig   = TRUE,
                         title      = "Volcano plot") {
  #--args: range of x-axis, y-axis
  fc   <- data %>% pull(x) %>% purrr::keep(is.finite) %>% summary()
  xmin <- floor(fc['Min.'])
  xmax <- ceiling(fc['Max.'])
  pval <- data %>% pull(y) %>% purrr::keep(is.finite) %>% summary()
  ymin <- floor(pval['Min.'])
  ymax <- ceiling(pval['Max.'])

  #--format: range of axis
  scatter_plot(data, x = x, y = y, labels = labels, show_sig = show_sig,
               show_abline = FALSE, trans_x = NULL, trans_y = NULL,
               xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
               title = title) +
    geom_vline(xintercept = 0, size = .5, color = "black") +
    geom_vline(
      xintercept = log2(xintercept),
      size       = .5,
      linetype   = 2,
      color      = "grey60"
    ) +
    geom_hline(
      yintercept = -log10(yintercept),
      size       = .5,
      linetype   = 2,
      color      = "grey60"
    )
}






#' @describeIn ma_plot Create MA plot for deseq output
#'
#' @param data A data frame for the full version of data, required
#' @param x character Column name show on x-axis, eg: baseMean
#' @param y character Column name chow on y-axis, eg: log2FoldChange
#' @param labels character Labels to show in plot
#' @param title character Title, default: "MA plot"
#'
#' @export
ma_plot <- function(data, x = "baseMean", y = "log2FoldChange",
                    labels = NULL, title = "MA plot") {
  #--args: range of x-axis, y-axis
  fc_val <- data %>% pull(x)
  fc   <- log10(fc_val + 1) %>% purrr::keep(is.finite) %>% summary()
  xmin <- floor(fc['Min.'])
  xmax <- ceiling(fc['Max.'])
  pval <- data %>% pull(y) %>% purrr::keep(is.finite) %>% summary()
  ymin <- floor(pval['Min.'])
  ymax <- ceiling(pval['Max.'])

  #--format: range of axis
  scatter_plot(data,
               x = x,
               y = y,
               show_abline = FALSE,
               show_sig    = TRUE,
               labels      = labels,
               add_label_point = TRUE,
               trans_x     = "log10",
               trans_y     = NULL,
               ymin        = ymin,
               ymax        = ymax,
               xmin        = xmin,
               xmax        = xmax,
               title       = title)
}






#' @describeIn scatter_plot Create scatter plot
#'
#' @param data A data frame for the full version of data, required
#' @param x character Column name show on x-axis
#' @param y character Column name chow on y-axis
#' @param labels character Labels to show in plot
#' @param add_label_point bool Add points for labels
#' @param show_sig bool Add colors for sig
#' @param show_abline bool Add ablines
#' @param trans character Function to convert num log10 or log2
#' @param xmin float Min on x-axis
#' @param xmax float Max x on x-axis
#' @param ymin float Min on y-axis
#' @param ymax float Max y on y-axis
#' @param title character Title, default: "Scatter plot"
#'
#' @import ggplot2
#' @import ggrepel
#'
#' @export
scatter_plot <- function(data, x, y, labels = NULL, add_label_point = TRUE,
                         show_sig = TRUE, show_abline = TRUE,
                         trans = "log10", trans_x = trans, trans_y = trans_x,
                         xmin = 0, xmax = 6, ymin = 0, ymax = 6,
                         title = "Scatter plot") {
  data <- ggplot2::fortify(data)

  #--check: input
  if(! is.character(x) && is.character(y)) {
    on.exit(message("x, y not character"), add = TRUE)
  }
  x <- x[1]
  y <- y[1]
  if(! all(c(x, y) %in% names(data))) {
    on.exit(message("x, y not found in data.frame"), add = TRUE)
  }

  #--transform: x, y axis, log-trans
  trans_func <- function(x) {
    if(is.character(x)) {
      switch (
        x,
        "log10" = log10,
        "log2"  = log2,
        "*"     = NULL,
      )
    }
  }

  # trans functions
  trans_func2 <- function(x, trans = "log10") {
    sapply(x, function(i) {
      ifelse(is.null(trans_func(trans)), i, trans_func(trans)(i + 1))
    })
  }

  data <- data %>%
    mutate(!! x := trans_func2(!!sym(x), trans_x),
           !! y := trans_func2(!!sym(y), trans_y),
           is_label = ifelse(label %in% labels, label, ""))

  # --format: axis title
  xtitle <- x
  ytitle <- y
  if(! is.null(trans_func(trans_x))) {
    xtitle <- bquote(.(x) ~ "[" * log["10"] ~ "rpm]")
  }
  if(! is.null(trans_func(trans_y))) {
    ytitle <- bquote(.(y) ~ "[" * log["10"] ~ "rpm]")
  }

  # --format: sig data
  if(! "sig" %in% names(data)) {
    data$sig  <- "sig"
  }
  if(! is.factor(data$sig)) {
    data$sig  <- forcats::as_factor(data$sig)
  }
  sig_vars <- levels(data$sig)
  sig_cols <- setNames(gg_color(length(sig_vars)), nm = sig_vars)
  if(isTRUE(show_sig)) {
    # red, green, blue, gray
    sig_v1 <- c("up", "down", "sig", "not")
    col_v1 <- setNames(c(gg_color(3)[1:2], "black", "grey40"), nm = sig_v1)
    if(all(sig_vars %in% sig_v1)) {
      sig_cols <- col_v1[sig_vars]
    }
    out <- ggplot(
      data, aes_(
        as.name(x), as.name(y),
        color = as.name("sig"),
        label = as.name("is_label"))) +
      scale_color_manual(values = sig_cols)
  } else {
    out <- ggplot(
      data, aes_(
        as.name(x), as.name(y),
        label = as.name("is_label")))
  }

  # --format: add labels
  data_label <- data %>%
    dplyr::filter(label %in% labels)

  if(isTRUE(add_label_point) && nrow(data_label) > 0) {
    out <- out +
      geom_point(aes_(as.name(x), as.name(y), color = as.name("sig")),
                 data = data_label, shape = 16, size = 2)
  }

  out <- out +
    ggrepel::geom_text_repel(
      color              = "grey10",
      size               = 3,
      force              = .2,
      max.overlaps       = 30,
      # nudge_x            = .3,
      # nudge_y            = .3,
      direction          = "both",
      point.padding      = .1,
      box.padding        = .1,
      segment.color      = "grey20",
      segment.size       = .4,
      min.segment.length = 0
    )

  # --format: abline
  if(isTRUE(show_abline)) {
    out <- out +
      geom_abline(slope = 1, intercept = 0, color = "grey10") +
      geom_abline(
        slope = 1, intercept = c(-trans_func(trans_x)(2), trans_func(trans_x)(2)),
        color = "grey40", linetype = 2
      )
  }

  # --format: axis
  out <- out +
    scale_x_continuous(
      limits = c(xmin, xmax),
      breaks = scales::pretty_breaks()(xmin:xmax), # seq(xmin, xmax, 1),
      labels = scales::pretty_breaks()(xmin:xmax), # seq(xmin, xmax, 1),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(ymin, ymax),
      breaks = scales::pretty_breaks()(ymin:ymax), # seq(ymin, ymax, 1),
      labels = scales::pretty_breaks()(ymin:ymax),# seq(ymin, ymax, 1),
      expand = c(0, 0)
    ) +
    geom_point(size = .3)

  # --format: add theme
  out +
    ggtitle(title) +
    xlab(xtitle) +
    ylab(ytitle) +
    labs(color = "Group") +
    theme_clean() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = .5),
      plot.title   = element_text(color = "black", hjust = .5, size = 14),
      panel.grid   = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.line    = element_line(color = "black", size = .5),
      axis.ticks   = element_line(color = "black", size = .5),
      axis.text    = element_text(color = "black", size = 10),
      axis.title   = element_text(color = "black", size = 12),
      axis.ticks.length = unit(.2, "cm")
    )
}


# # px = readRDS("publish_plot_data.rds")
# # df = px$p1$data
# df = readr::read_delim("transcripts_deseq2.fix.xls", "\t")
# x  = names(df)[2]
# y  = names(df)[3]
# highlight_column = "sig"
# highlight_values = c("up", "down")
# label_column = "SYMBOL"
# # label_values = c("Rin2", "Cd36")
# label_values = df %>%
#   filter(sig %in% highlight_values) %>%
#   pull(SYMBOL) %>%
#   head(10)
#
# p1 <- scatter_plot2(
#   df, x, y,
#   highlight_column=highlight_column,
#   highlight_values=highlight_values,
#   label_column=label_column,
#   label_values=label_values,
#   density_point = FALSE,
#   point_color = "grey50")


#' @describeIn scatter_plot2 Create scatter plot, version2
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
scatter_plot2 <- function(data, x, y, ...) {
  args <- rlang::list2(...)
  # density_point <- TRUE
  # point_size  <- 0.3
  # point_color <- "grey30"
  # trans_axis  <- "log10"
  # xmin <- 0
  # xmax <- 6
  # ymin <- xmin
  # ymax <- xmax
  # add_fc_lines <- TRUE
  ## --default values:
  add_fc_lines  <- ifelse(rlang::has_name(args, "add_fc_lines"),
                          args$add_fc_lines, TRUE)
  density_point <- ifelse(rlang::has_name(args, "density_point"),
                          args$density_point, FALSE)
  point_size    <- ifelse(rlang::has_name(args, "point_size"),
                          args$point_size, 0.5)
  point_size_highlight <- ifelse(rlang::has_name(args, "point_size_highlight"),
                          args$point_size_highlight, 2)
  point_color   <- ifelse(rlang::has_name(args, "point_color"),
                          args$point_color, "grey50")
  trans_axis    <- ifelse(rlang::has_name(args, "trans_axis"),
                          args$trans_axis, "log10")
  trans_axis_x  <- ifelse(rlang::has_name(args, "trans_axis_x"),
                          args$trans_axis_x, trans_axis)
  trans_axis_y  <- ifelse(rlang::has_name(args, "trans_axis_y"),
                          args$trans_axis_y, trans_axis_x)
  xmin  <- ifelse(rlang::has_name(args, "xmin"), args$xmin, 0)
  xmax  <- ifelse(rlang::has_name(args, "xmax"), args$xmax, 6)
  ymin  <- ifelse(rlang::has_name(args, "ymin"), args$ymin, xmin)
  ymax  <- ifelse(rlang::has_name(args, "ymax"), args$ymax, xmax)
  title <- ifelse(rlang::has_name(args, "title"), args$title, "Scatter")
  # --plot data: pd
  pd <- prep_scatter_data(data, x, y, ...) # conversion
  # --main: p0
  p0 <- pd$df_point %>%
    ggplot(aes_string(as.name(x), as.name(y), label = "is_label"))
  # --points to density
  if(isTRUE(density_point)) {
    out <- p0 +
      stat_density_2d(
        aes(fill = after_stat(density)),
        geom = "raster",
        contour = FALSE) +
      scale_fill_gradient(low = "white", high = "black")
  } else {
    out <- p0 +
      geom_point(size = point_size, color = point_color)
  }
  # --highlight points:
  if(! is.na(pd$highlight_column)) {
    out <- out +
      geom_point(
        aes_string(as.name(pd$x), as.name(pd$y), color = pd$highlight_column),
        data = pd$df_highlight, shape = 16,
        size = point_size_highlight) # !!!! tmp
  }
  # --labels: add label in aes()
  if(! is.na(pd$label_column)) {
    out <- out +
      geom_point(
        aes_string(
          as.name(pd$x), as.name(pd$y)),
        data = pd$df_label, shape = 16, size = point_size_highlight, # !!!! tmp
        color = "grey30") +
      geom_text_repel(
        aes_string(as.name(pd$x), as.name(pd$y), label = "is_label"),
        data = pd$df_label, inherit.aes = FALSE,
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
  # --lines:
  if(isTRUE(add_fc_lines)) {
    out <- out +
      geom_abline(slope = 1, intercept = 0, color = "grey10") +
      geom_abline(
        slope = 1,
        intercept = c(-pd$trans_func(trans_axis)(2),
                      pd$trans_func(trans_axis)(2)),
        color = "grey40", linetype = 2)
  }
  # --format axis:
  out <- out +
    scale_x_continuous(
      limits = c(xmin, xmax),
      breaks = scales::pretty_breaks()(xmin:xmax), # seq(xmin, xmax, 1),
      labels = scales::pretty_breaks()(xmin:xmax), # seq(xmin, xmax, 1),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(ymin, ymax),
      breaks = scales::pretty_breaks()(ymin:ymax), # seq(ymin, ymax, 1),
      labels = scales::pretty_breaks()(ymin:ymax), # seq(ymin, ymax, 1),
      expand = c(0, 0)
    )
  # --theme:
  out <- out +
    guides(fill = FALSE) + # remove density
    ggtitle(title) +
    xlab(pd$xtitle) +
    ylab(pd$ytitle) +
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
      axis.ticks.length = unit(.2, "cm")
    )
  out
}




## helper functions
prep_scatter_data <- function(data, x, y, ...) {
  # --default args:
  args <- rlang::list2(...)
  if("trans_axis" %in% names(args)) {
    trans_axis <- args$trans_axis
  } else {
    trans_axis <- "log10"
  }
  highlight_column <- ifelse(
    "highlight_column" %in% names(args), args$highlight_column, NA)
  label_column <- ifelse(
    "label_column" %in% names(args), args$label_column, NA)
  if("highlight_values" %in% names(args)) {
    highlight_values <- args$highlight_values
  } else {
    highlight_values <- NA
  }
  if("label_values" %in% names(args)) {
    label_values <- args$label_values
  } else {
    label_values <- NA
  }
  # --arguments: data
  if(! is.character(x) && is.character(y)) {
    on.exit(message("x, y not character"), add = TRUE)
  }
  x <- x[1]
  y <- y[1]
  if(! all(c(x, y) %in% names(data))) {
    on.exit(message("x, y not found in data.frame"), add = TRUE)
  }
  # --arguments: x, y axis, log-trans
  trans_func <- function(x) {
    if(is.character(x)) {
      switch (
        x,
        "log10" = log10,
        "log2"  = log2,
        "*"     = NULL,
      )
    }
  }
  # --trans functions
  trans_func2 <- function(x, trans = "log10") {
    sapply(x, function(i) {
      ifelse(is.null(trans_func(trans)), i, trans_func(trans)(i + 1))
    })
  }
  # --convert:
  df <- data %>%
    dplyr::mutate(
      !! x := trans_func2(!!sym(x), trans_axis),
      !! y := trans_func2(!!sym(y), trans_axis))
  # --arguments:
  xtitle <- x
  ytitle <- y
  if(! is.null(trans_func(trans_axis))) {
    xtitle <- bquote(.(x) ~ "[" * log["10"] ~ "rpm]")
  }
  if(! is.null(trans_func(trans_axis))) {
    ytitle <- bquote(.(y) ~ "[" * log["10"] ~ "rpm]")
  }
  # --highlight_points:
  if(highlight_column %in% names(df)) {
    h1 <- df %>%
      dplyr::pull(highlight_column)
    if(! all(highlight_values %in% h1)) {
      highlight_values <- NA
    }
  } else {
    highlight_values <- NA
    highlight_column <- NA
  }
  # --labels: add label column
  if(is.na(label_column)) {
    df$label     <- NULL
    df$is_label  <- NULL
    label_values <- NA
    label_column <- NA
  } else {
    h2 <- df %>%
      dplyr::pull(label_column)
    if(! all(label_values %in% h2)) {
      label_values <- NA
    }
    # add label column: label
    # add label check: is_label
    df <- df %>%
      dplyr::mutate(
        label = !!(as.name(label_column)),
        is_label = ifelse(label %in% label_values, label, NA))
  }
  # --highlight: data
  if(is.na(highlight_column)) {
    df_highlight <- NA
  } else {
    df_highlight <- df %>%
      dplyr::filter(!!as.symbol(highlight_column) %in% highlight_values) %>%
      dplyr::mutate(!!as.symbol(highlight_column) := factor(
        (!! as.symbol(highlight_column)), highlight_values
      ))
  }
  # --label: data
  if(is.na(label_column)) {
    df_label <- NA
  } else {
    df_label <- df %>%
      dplyr::filter(label %in% label_values)
  }
  # --update: highlight, remove highlight data
  if(! is.na(highlight_column)) {
    df1 <- df %>%
      dplyr::filter(! (!!as.symbol(highlight_column)) %in% highlight_values)
  }
  if(! is.na(label_column)) {
    df1 <- df1 %>%
      dplyr::filter(! label %in% label_values)
  }
  # --Output:
  list(
    data = df,
    x = x,
    y = y,
    xtitle = xtitle,
    ytitle = ytitle,
    highlight_column = highlight_column,
    highlight_values = highlight_values,
    label_column = label_column,
    label_values = label_values,
    df_point     = df1,
    df_highlight = df_highlight,
    df_label = df_label,
    trans_func = trans_func,
    trans_func2 = trans_func2
  )
}





#--Custome: for specific mission: RNAseq ---------------------------------------



#' @describeIn rnaseq_trim_stat_plot Create bar_plot for trim stat
#'
#' output from get_rnaseq_trim_stat(),
#' including columns:
#' id, raw, clean, clean_pct, short_pct
#'
#'
#' @param data data.frame From get_trim_stat
#' @param fish string Name of the fish, use `fishualize::fish_palettes()` list
#' all available fish names
#'
#' @export
rnaseq_trim_stat_plot <- function(data, fish = "Trimma_lantana") {
  if(is(data, "data.frame")) {
    col_required <- c("id", "clean_pct", "short_pct")
    if(! all(col_required %in% names(data))) {
      stop("`data` required columns missing: ",
           paste(col_required, collapse = ", "))
    }
  } else {
    stop("`data` expect data.frame, failed")
  }
  data %>%
    dplyr::select(id, clean_pct, short_pct) %>%
    tidyr::pivot_longer(names_to  = "group",
                        values_to = "count",
                        c(clean_pct, short_pct)) %>%
    bar_plot(x = "count", y = "id", fill = "group", label = "count") +
    fishualize::scale_fill_fish(discrete = TRUE, option = fish) +
    ggtitle("Trim reads")

}




#' @describeIn rnaseq_align_stat_plot Create bar_plot for align_stat
#'
#' Output from get_rnaseq_align_stat()
#' including columns:
#' fqname, map, unique, multiple
#'
#' @param data data.frame From get_rnaseq_align_stat()
#' @param mode integer map=1, unique/multiple=2, default: 1
#' @param fish string Name of the fish, use `fishualize::fish_palettes()` list
#'  all availabel fish names
#'
#' @export
rnaseq_align_stat_plot <- function(data,
                                   mode = 0,
                                   columns = NULL,
                                   fish = "Trimma_lantana",
                                   title = "Barplot") {
  if(! is(data, "data.frame")) {
    stop("`data` expect data.frame, failed")
  }
  if(! "fqname" %in% names(data)) {
    stop("`data` missing the column: fqname")
  }
  # columns
  g1 <- switch(mode,
               c("map", "unmap"),
               c("unique", "multiple", "unmap"),
               c("total"))
  if(length(g1)) {
    columns <- g1
  }
  columns <- purrr::discard(columns, is.null)
  if(! all(columns %in% names(data))) {
    d_line <- paste(names(data), collapse = ", ")
    stop("`columns` missing, expect: ", d_line)
  }
  # subset data.frame
  data %>%
    dplyr::select(any_of(c("fqname", columns))) %>%
    tidyr::pivot_longer(names_to = "group", values_to = "count", -fqname) %>%
    dplyr::mutate(group = factor(group, levels = columns)) %>%
    dplyr::group_by(fqname) %>%
    dplyr::mutate(pct = round(count / sum(count) * 100, 2)) %>%
    bar_plot(x = "pct", y = "fqname", direction = "horizontal",
             fill = "group", label = NULL) +
    fishualize::scale_fill_fish(discrete = TRUE, option = fish) +
    ylab(NULL) +
    xlab("Percentage%") +
    # theme(legend.position = "top") +
    ggtitle(title)
}




#' @describeIn rnaseq_sig_stat_plot Create bar_plot for sig count
#'
#' @param data data.frame From deseq_dir, fix.xls
#'
#' @export
rnaseq_sig_stat_plot <- function(data, fish = "Trimma_lantana") {
  if(is(data, "data.frame")) {
    col_required <- c("sig")
    if(! all(col_required %in% names(data))) {
      stop("`data` required columns missing: ",
           paste(col_required, collapse = ", "))
    }
  } else {
    stop("`data` expect data.frame, failed")
  }
  data %>%
    group_by(sig) %>%
    summarize(count = n()) %>%
    bar_plot(x = "sig", y = "count", fill = "sig", label = "count") +
    fishualize::scale_fill_fish(discrete = TRUE, option = fish) +
    ggtitle("DE genes")
}





#' @describeIn rnaseq_count_summary
#'
#' @param
#'
#'
#' @export
rnaseq_count_summary_plot <- function(data, fish = "Trimma_lantana") {
  if(is(data, "data.frame")) {
    col_required <- c("Status", "sample", "count")
    if(! all(col_required %in% names(data))) {
      stop("`data` required columns missing: ",
           paste(col_required, collapse = ", "))
    }
  } else {
    stop("`data` expect data.frame, failed")
  }
  data %>%
    tidyr::pivot_wider(names_from = "Status", values_from = "count") %>%
    dplyr::mutate(total      = dplyr::select_if(., is.numeric) %>% rowSums,
                  Unassigned = total - Assigned) %>%
    dplyr::select(sample, total, Assigned, Unassigned) %>%
    tidyr::pivot_longer(names_to = "assign",
                        values_to = "count",
                        c(Assigned, Unassigned)) %>%
    bar_plot(x = "count", y = "sample", fill = "assign", label = NULL,
             direction = "horizontal") +
    fishualize::scale_fill_fish(discrete = TRUE, option = fish) +
    ggtitle("FeatureCounts Assign Summary")
}













#' corPlot
#'
#' @param x path to file, count matrix
#' @import readr
#' @import dplyr
#' @import RColorBrewer
#'
#' @export
cor_plot <- function(df){
  ggcor::ggcor(df, type = "lower", show.diag = TRUE,
               cor.test.method = "pearson") +
    ggcor::geom_color(data = get_data(type = "lower", show.diag = TRUE)) +
    ggcor::geom_num(aes(num = r), data = get_data(type = "lower", show.diag = TRUE),
                    colour = "grey90", size = 3) +
    ggplot2::ggtitle("Pearson coefficient correlation")
}
























