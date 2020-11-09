

#' vennplot
#'
#' @param x list, unique ids for groups
#' @param names character, the names for each group
#'
#' @import ggVennDiagram
#' @import ggplot2
#'
#' @export
vennplot <- function(x, names = NULL){
  # check category names
  if(is.null(names)) {
    names <- names(x)
  }
  # further
  if(is.null(names)) {
    names <- letters[seq_len(length(x))]
  }
  # main
  p <- ggVennDiagram::ggVennDiagram(x, label = "count", label_alpha = 0,
                                    category.names = names) +
    ggplot2::guides(fill = FALSE) +
    ggplot2::scale_fill_gradient(low = "white", high = "firebrick")
  return(p)
}




#' barplot2
#'
#' @param df data.frame
#' @import readr
#' @import dplyr
#' @import ggplot2
#'
#' @export
barplot2 <- function(df, label = FALSE){
  stopifnot(all(c("id", "count") %in% names(df)))

  p <- ggplot(df, aes(id, count)) +
    geom_bar(stat = "identity") +
    scale_y_continuous(position = "right") +
    coord_flip() +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text  = element_text(color = "grey20"),
      plot.title = element_text(hjust = 0.5)
    )

  if(isTRUE(label)) {
    p <- p +
      geom_text(aes(label = count), hjust = 1.2, color = "grey90")
  }

  return(p)
}




#' fragPlot
#'
#' @param df data.frame for plotting
#' @param range
#' @param type string, line or bar
#' @param y_axis string, fraction, count
#' @import readr
#' @import dplyr
#' @import ggplot2
#'
#' @export
frag_plot <- function(df, type = "line", x_axis = "auto",
                      x_min = 0, x_max = 1000, y_axis = "fraction") {
  stopifnot(all(c("id", "length", "count") %in% names(df)))

  # cal frac%
  df1 <- df %>%
    group_by(id) %>%
    mutate(frac = count / sum(count))

  # plot
  if(y_axis %in% c("count")) {
    p <- ggplot(df1, aes(length, count)) +
      ylab("Number of reads")
  } else {
    p <- ggplot(df1, aes(length, frac)) +
      ylab("Normalized read density")
  }

  # check type: line, bar
  if(type == "line") {
    p <- p +
      geom_line(aes(color = id), color = "red3", size = .5)
  } else {
    p <- p +
      geom_bar(aes(fill = id), stat = "identity")
  }

  # range
  if(x_axis == "auto" & x_max == 1000) {
    p <- p +
      scale_x_continuous(breaks = seq(0, 1000, by = 200),
                     labels = seq(0, 1000, by = 200),
                     limits = c(0, 1000),
                     expand = c(0, 0))
  } else {
    p <- p +
      scale_x_continuous(limits = c(x_min, x_max),
                         expand = c(0, 0))
  }

  # theme
  p +
    theme_bw() +
    theme(
      panel.grid = element_blank()
    )
}



#' fragPlot
#'
#' @param df data.frame for plotting
#' @import readr
#' @import dplyr
#' @import ggplot2
#'
#' @export
frag_plot2 <- function(df) {
  stopifnot(all(c("length", "count", "id") %in% names(df)))

  df2 <- df %>%
    group_by(id) %>%
    mutate(freq = count / sum(count))

  df2 %>%
    ggplot(aes(length, freq, color = id)) +
    geom_line(size = .5) +
    xlab("Fragment length (bp)") +
    ylab("Normalized read density") +
    scale_x_continuous(breaks = seq(0, 1000, by = 200),
                       labels = seq(0, 1000, by = 200),
                       limits = c(0, 1000),
                       expand = c(0, 0)) +
    theme_bw() +
    theme(
      panel.grid = element_blank()
    )

}





#' readAlign1
#'
#' @param df data.frame, mapping stats
#' @import readr
#' @import dplyr
#' @import ggplot2
#'
#' @export
align_plot <- function(df) {
  stopifnot(all(c("id", "count", "group") %in% names(df)))

  group_colors <- c("darkgreen", "green2", "orange4", "orange", "grey50")

  p <- df %>%
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







