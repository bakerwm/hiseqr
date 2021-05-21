#' Functions for colors
#'
#' This is an implementation of pick_colors, from ggplot2 default, or
#' RColorBrewer, fishualize
#'
#'
#' @name colors




#' @describeIn gg_color Pick ggplot2 default colors
#' ggplot2_colors
#'
#' @param n integer
#'
#' @import scales
#'
#' @export
gg_color <- function(n = 3) {
  # version 1
  # n = 3
  # scales::hue_pal()(3)
  #
  # scales::show_col(scales::hue_pal()(3))

  # version2
  # https://stackoverflow.com/a/8197703/2530783
  #
  # gg_color_hue <- function(n) {
  #   hues = seq(15, 375, length = n + 1)
  #   hcl(h = hues, l = 65, c = 100)[1:n]
  # }
  #
  # scales::show_col(gg_color_hue(3))

  # version3
  # https://stackoverflow.com/a/8197706/2530783
  # ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  #   if ((diff(h) %% 360) < 1) {
  #     h[2] <- h[2] - 360/n
  #   }
  #   hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
  # }
  #
  # scales:::show_col(ggplotColours(n=3))

  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}














