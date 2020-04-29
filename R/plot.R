

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
#' @import readr
#' @import dplyr
#' @import ggplot2
#'
#' @export
fragPlot <- function(df) {
  stopifnot(all(c("id", "length", "count") %in% names(df)))

  p <- df %>%
    group_by(id) %>%
    mutate(frac = count / sum(count)) %>%
    ggplot(aes(length, frac, color = id)) +
    geom_line(color = "red3", size = .5) +
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
alignPlot <- function(df) {
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




#' corPlot
#'
#' @param x path to file, count matrix
#' @import readr
#' @import dplyr
#' @import RColorBrewer
#'
#' @export
corPlot <- function(df){
  p <- ggcor::ggcor(df, type = "lower", show.diag = TRUE,
             cor.test.method = "pearson") +
    ggcor::geom_color(data = get_data(type = "lower", show.diag = TRUE)) +
    ggcor::geom_num(aes(num = r), data = get_data(type = "lower", show.diag = TRUE),
             colour = "grey90", size = 3) +
    ggplot2::ggtitle("Pearson coefficient correlation")

  return(p)
}







