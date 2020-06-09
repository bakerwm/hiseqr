
#' check x, go anslysis output
#'
#' support:
#' groupGOResult,
#' enrichResult,
#' ...
#'
#' @param x object
#' @param recursive boolean, whether loop over list
#'
#' @export
is_go <- function(x, recursive = FALSE) {
  # required
  go_class <- c("groupGOResult", "enrichResult", "gseaResult")

  # check
  if(recursive) {
    if(is.list(x)) {
      b_list <- sapply(x, function(i){
        is_go(i, TRUE)
      })
      return(all(b_list))
    } else {
      return(class(x) %in% go_class)
    }

  } else {
    return(class(x) %in% go_class)
  }
}




## input is foldChange
## names, genename
is_gene_fc <- function(x) {
  # numeric
  # names = gene
  return(is.numeric(x) & is.character(names(x)))
}



## input is character
is_gene_list <- function(x){
  # names = gene
  return(is.character(x))
}


