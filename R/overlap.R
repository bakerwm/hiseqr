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





## pairwise
# a, b, ab
overlap_p3 <- function(n1, n2, n3, n12, n13, n23, n123,
                       group_name = NULL,
                       png_file = NULL) {
  # to-do: check input arguments
  # out list
  nums <- setNames(c(n1, n2, n3, n12, n13, n23, n123),
                   nm = c("n1", "n2", "n3", "n12", "n13", "n23", "n123"))
  p <- lapply(seq_len(7), function(i){
    paste(letters[i], seq_len(nums[i]), sep = "")
  })
  names(p) <- names(nums)
  # combine
  out <- list(
    A = c(p$n1, p$n12, p$n13, p$n123),
    B = c(p$n2, p$n12, p$n23, p$n123),
    C = c(p$n3, p$n13, p$n23, p$n123)
  )
  if(inherits(group_name, "character")) {
    group_name <- group_name[1:3]
  } else {
    group_name <- LETTERS[1:3]
  }
  # to matrix
  ma <- rbind(
    data.frame(group = names[1], id = out[[1]]),
    data.frame(group = names[2], id = out[[2]]),
    data.frame(group = names[3], id = out[[3]])) %>%
    dplyr::mutate(score = 1) %>%
    tidyr::pivot_wider(names_from = "group", values_from = score) %>%
    mutate(across(where(is.numeric), ~tidyr::replace_na(., 0))) %>%
    tibble::column_to_rownames("id") %>%
    as.matrix()
  # rename
  colnames(ma) <- group_name
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




