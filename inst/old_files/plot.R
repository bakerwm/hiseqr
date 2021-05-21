#' Functions for HiSeq plots
#'
#' Including bar_plot, scatter_plot, venn_plot, ma_plot,  box_plot, ...
#'
#'
#'
#'
#' @name plot



#' @describeIn bar_plot Make barplot
#'
#' @param df data.frame require id, count column
#' @import readr
#' @import dplyr
#' @import ggplot2
#'
#' @export
bar_plot <- function(
  data = NULL,
  x, y,
  label = y,
  direction = "vertical", ...
) {
  data <- fortify(data)
  if(! is_character(x) && is_character(y)) {
    on.exit(message("x, y not character"), add = TRUE)
  }
  x <- x[1]
  y <- y[1]
  if(! all(c(x, y) %in% names(data))) {
    on.exit(message("x, y not found in data.frame"), add = TRUE)
  }
  if(direction == "vertical") {
    out <- data %>%
      ggplot(aes_(as.name(x), as.name(y))) +
      geom_col() +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
      )
    if(is_character(label)) {
      out <- out +
        geom_text(aes_(label = as.name(label)),
                  hjust = -.1,
                  color = "black",
                  angle = 90)
    }
  } else {
    out <- data %>%
      ggplot(aes_(as.name(y), as.name(x))) +
      geom_col() +
      scale_x_continuous(position = "top") +
      scale_y_discrete(limits = rev) +
      theme_bw() +
      theme(
        axis.text  = element_text(color = "grey20"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
      )
    if(is_character(label)) {
      out <- out +
        geom_text(aes_(label = as.name(label)), hjust = -.1, color = "black")
    }
  }
  out
}



#' @describeIn fragsize_plot line plot for fragment size
#'
#' @param df data.frame for plotting
#' @param range
#' @param type string, line or bar
#' @param y_axis string, fraction, count
#'
#' @import readr
#' @import dplyr
#' @import ggplot2
#'
#' @export
fragsize_plot <- function(
  data,
  xmin = 0,
  xmax = 1000,
  log10_y = FALSE,
  position = "fill"
) {
  data <- fortify(data)
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
#' @param x list, unique ids for groups
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












#' readAlign1
#'
#' @param df data.frame, mapping stats
#' @import readr
#' @import dplyr
#' @import ggplot2
#' @import fishualize
#'
#' @export
align_plot <- function(df) {
  stopifnot(all(c("id", "count", "group") %in% names(df)))
  group_colors <- c("darkgreen", "green2", "orange4", "orange", "grey50")

  p <- df %>%
    ggplot(aes(id, count, fill = group)) +
    geom_bar(stat = "identity", position = "fill") +
    # scale_fill_manual(values = rev(group_colors)) +
    scale_fill_fish_d(option = "Bodianus_rufus") + # Scarus_quoyi
    scale_y_continuous(position = "right") +
    # geom_hline(yintercept = mito_mean) +
    xlab(NULL) + ylab("Percentage") +
    coord_flip() +
    guides(fill = guide_legend(title = NULL,
                               label.position = "top",
                               reverse = TRUE)) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 12),
      legend.position = "top"
    )

  return(p)
}



#' readAlign1
#'
#' @param df data.frame, mapping stats
#' @import readr
#' @import dplyr
#' @import ggplot2
#'
#' @export
align_plot2 <- function(df) {
  stopifnot(all(c("id", "count", "group") %in% names(df)))

  # group_colors <- c("darkgreen", "green2", "orange4", "orange",  "grey50")
  group_colors <- RColorBrewer::brewer.pal(9, "Set1")

  df %>%
    ggplot(aes(id, count, fill = group)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = rev(group_colors)) +
    scale_y_continuous(position = "right") +
    # geom_hline(yintercept = mito_mean) +
    xlab(NULL) + ylab("Percentage") +
    coord_flip() +
    guides(fill = guide_legend(title = NULL,
                               label.position = "top",
                               reverse = TRUE)) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 12),
      legend.position = "top"
    )
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
  p <- ggcor::ggcor(df, type = "lower", show.diag = TRUE,
                    cor.test.method = "pearson") +
    ggcor::geom_color(data = get_data(type = "lower", show.diag = TRUE)) +
    ggcor::geom_num(aes(num = r), data = get_data(type = "lower", show.diag = TRUE),
                    colour = "grey90", size = 3) +
    ggplot2::ggtitle("Pearson coefficient correlation")

  return(p)
}







