
## demx


edit_distance <- function(x, y) {
  ma <- adist(x, y, partial = TRUE, ignore.case = TRUE)

  library(plotly)
  plotly::plot_ly(z = ma, x = c("a", "b", "c"),
                  y = c("a", "b", "c"),
                  type = "heatmap")

}




