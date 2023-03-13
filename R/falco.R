#' falco
#'
#' read falco output
#'

## parse falco output
# f <- "fastqc_data.txt"
# f <- "../fastqc_data.txt"
# s <- readr::read_file(f)
# s <- unlist(strsplit(s, ">>"))
# split into modules
# s <- purrr::discard(s, grepl("^##|^END", s))
# f <- list.files(".", "fastqc_data.txt$")


#' falco_plot_base_content
#'
#' @param x output of falco_read_txt()
#' @param mode character line, bar
#'
#' # module = 4, Per base sequence content
#'
#' @export
falco_plot_base_content <- function(x, mode = "line", ...) {
  df <- falco_fetch_data(x, module = 4, group = "data") %>%
    tidyr::pivot_longer(-1, names_to = "base", values_to = "freq") %>%
    dplyr::rename(pos = `#Base`)
  if(mode == "line") {
    p <- df %>%
      ggplot(aes(pos, freq, color = base)) +
      geom_line(size = 0.5, ...) +
      scale_color_manual(
        values = setNames(c("green3", "blue3", "grey30", "red3"),
                          nm = c("A", "C", "G", "T"))) +
      scale_y_continuous(breaks = seq(0, 100, 20), labels = seq(0, 100, 20)) +
      xlab("Position in read (bp)") +
      ylab("Per base sequence content (%)") +
      theme_bw() +
      theme(legend.position = c(0.5, 0.9), legend.direction = "horizontal")
  } else if(mode == "bar") {
    p <- df %>%
      ggplot(aes(pos, freq, fill = base)) +
      geom_bar(position = "fill", stat = "identity", ...) +
      scale_fill_manual(values = list(
        A = "green3", C = "blue3", G = "grey30", T = "red3")) +
      scale_y_continuous(breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) +
      xlab("Position in read (bp)") +
      ylab("Per base sequence content (%)") +
      theme_bw() +
      theme(legend.position = "top")
  }
  p
}


#' falco_plot_lendist
#' bar plot or line plot
#' module = 7, Sequence Length Distribution
#'
#' @param x output of falco_read_txt()
#' @param mode character line, bar
#'
#' @export
falco_plot_lendist <- function(x, mode = "line", ...) {
  # Length Count
  df <- falco_fetch_data(x, module = 7, group = "data") %>%
    dplyr::rename(Len = `#Length`)
  # draw plots
  if(mode == "line") {
    p <- df %>%
      ggplot(aes(Len, Count)) +
      geom_point(size = 1, color = "red3") +
      geom_line(size = .5, color = "red3") +
      xlab("Length of reads (nt)") +
      ylab("Number of reads") +
      theme_classic()
  } else if(mode == "bar") {
    p <-   df %>%
      ggplot(aes(Len, Count)) +
      geom_col(...) +
      xlab("Length of reads (nt)") +
      ylab("Number of reads") +
      theme_classic()
  }
  p
}


#' falco_plot_lendist
#' bar plot or line plot
#' module = 7, Sequence Length Distribution
#'
#' @param x output of falco_read_txt()
#' @param mode character line, bar
#'
#' @export
falco_plot_base_content <- function(x, mode = "line") {
  # Length Count
  df <- falco_fetch_data(x, module = 7, group = "data") %>%
    dplyr::rename(Len = `#Length`)
  # draw plots
  if(mode == "line") {
    p <- df %>%
      ggplot(aes(Len, Count)) +
      geom_point(size = 1, color = "red3") +
      geom_line(size = .5, color = "red3") +
      xlab("Length of reads (nt)") +
      ylab("Number of reads") +
      theme_classic()
  } else if(mode == "bar") {
    p <-   df %>%
      ggplot(aes(Len, Count)) +
      geom_col() +
      xlab("Length of reads (nt)") +
      ylab("Number of reads") +
      theme_classic()
  }
  p
}


#' get_falco_data
#'
#' @param x full list of falco_read_txt
#'
#' @export
falco_fetch_data <- function(x, module = "ALL", group = "ALL") {
  if(inherits(x, "list")) {
    # modules
    if(module %in% c("ALL", "all", TRUE)) {
      m <- names(x)
    } else {
      m <- falco_guess_module(module)
    }
    # groups
    if(group[1] %in% c("ALL", "all", TRUE)) {
      group <- c("id", "data", "flag", "extra")
    }
    # extract group
    if(length(m) > 0) {
      out <- lapply(m, function(i) {
        if(length(group) > 1) {
          x[[i]][group]
        } else {
          x[[i]][[group]]
        }
      })
      # add names
      names(out) <- m
      if(length(m) == 1) {
        out <- out[[1]]
      }
      out
    }
  }
}


#' read_falco
#'
#' @param x string fastqc_data.txt file
#' @param module int or character the name of the modules
#' @param group character data, id, flag, extra
#'
#' @export
falco_read_txt <- function(x, module = "ALL", group = "ALL") {
  out <- tryCatch(
    {
      s <- readr::read_file(x)
      s <- unlist(strsplit(s, ">>"))
      s <- purrr::discard(s, grepl("^##|^END", s)) # 10 modules
      out <- lapply(s, falco_read_module)
      names(out) <- unlist(sapply(out, function(i){ i["id"] }, simplify = TRUE))
      # only groups
      if(group %in% c("id", "data", "flag", "extra")) {
        out <- lapply(out, function(i) i[[group]])
      }
      out
    },
    error = function(cond) {
      message(paste0("Failed to read file: ", x))
      return(NA)
    }
  )
  if(inherits(out, "list")) {
    if(module %in% c("ALL", "all", TRUE)) {
      m <- names(out)
    } else {
      m <- falco_guess_module(module)
    }
    out[m]
  }
}


#' read_falco_module
#'
#' @param x string
#'
#' @export
falco_read_module <- function(x) {
  s  <- unlist(strsplit(x, "\n"))
  s2 <- unlist(strsplit(s[1], "\t"))
  if(grepl("#Total Deduplicated", x)) {
    extra <- s[2]
    x2    <- paste(s[-c(1, 2)], collapse = "\n")
  } else {
    extra <- NULL
    x2 <- paste(s[-1], collapse = "\n")
  }
  if(any(grepl("\n", x2))) {
    df <- readr::read_delim(x2, "\t", col_types = readr::cols(),
                            skip_empty_rows = TRUE)
  } else {
    df <- NULL
  }
  list(
    id    = s2[1],
    flag  = s2[2],
    extra = extra,
    data  = df
  )
}


#' falco_guess_module
#'
#' @export
falco_guess_module <- function(x) {
  # default modules in falco v0.2.4
  m <- c("Basic Statistics",
         "Per base sequence quality",
         "Per sequence quality scores",
         "Per base sequence content",
         "Per sequence GC content",
         "Per base N content",
         "Sequence Length Distribution",
         "Sequence Duplication Levels",
         "Overrepresented sequences",
         "Adapter Content"
  )
  # check modules
  if(inherits(x, "numeric")) {
    i <- x
  } else if(inherits(x, "character")) {
    i <- grep(x, m, ignore.case = TRUE)
  } else {
    message(glue::glue("Expect str or int, got {class(x)}"))
    return(NULL)
  }
  # check exceptions
  purrr::discard(m[i], is.na)
}



