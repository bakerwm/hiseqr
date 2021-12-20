#! create overlap plot, euler, venn




## make euler plot

## pairwise
# a, b, ab

overlap_p2 <- function(n1, n2, n12, name1 = "A", name2 = "B",
                       png_file = NULL) {
  # to-do: check input arguments
  # n1 = 40940
  # n2 = 38313
  # n12 = 33731
  # name1 = "A"
  # name2 = "B"
  # construct data
  ia = paste0("A", seq_len(n1 - n12))
  ib = paste0("B", seq_len(n2 - n12))
  iab = paste0("AB", seq_len(n12))
  # for matrix
  ma <- rbind(data.frame(group = name1, id = c(ia, iab)),
              data.frame(group = name2, id = c(ib, iab))) %>%
    dplyr::mutate(score = 1) %>%
    tidyr::pivot_wider(names_from = "group", values_from = score) %>%
    mutate(across(where(is.numeric), ~tidyr::replace_na(., 0))) %>%
    tibble::column_to_rownames("id") %>%
    as.matrix()
  # compute
  f <- eulerr::euler(ma) # compute overlap
  p <- plot(
    f,
    fill = "transparent",
    lwd  = 2,
    col  = c("red3", "blue3", "green3", "darkorange3", "darkorchid", "black"),
    quantities = list(type = c("counts")),
    labels = list(font = 1, cex = 0.8))
  # save to file
  if(inherits(png_file, "character")) {
    if(dir.exists(normalizePath(dirname(png_file)))) {
      png(png_file, width = 5, height = 5, res = 200, units = "in")
      print(p)
      dev.off()
    } else {
      warning(glue::glue("could not write to file: {png_file}"))
    }
  }
  # output
  p
}




