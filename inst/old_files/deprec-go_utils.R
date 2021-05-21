#' Deprecated go_utils functions
#'
#'

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



#' get the orgdb for organism
#'
#' @param x name of the organism, Drosophila melanogaster
#'
#' @import AnnotationHub
#'
#' ah    <- AnnotationHub::AnnotationHub()
#' orgdb <- AnnotationHub::query(ah, "OrgDb") # full
#'
#' sp <- orgdb$species
#' h  <- sp[grep("Drosophila melanogaster", sp)]
#' db <- query(orgdb, h)
#' obj_orgdb <- ah[[db$ah_id]]
#'
#' @export
get_orgdb <- function(x) {
  # for regular speciexs:
  # human, mouse, fruit fly
  db_list <- list(human = org.Hs.eg.db::org.Hs.eg.db,
                  mouse = org.Mm.eg.db::org.Mm.eg.db,
                  fruitfly = org.Dm.eg.db::org.Dm.eg.db)
  if(grepl("Drosophila melanogaster|fruitfly|dm6|dm3", x)) {
    orgdb <- db_list$fruitfly
  } else if(grepl("Homo sapiens|human|hg38|hg19", x)) {
    orgdb <- db_list$human
  } else if(grepl("Mus musculus|mouse|mm10|mm9", x)) {
    orgdb <- db_list$mouse
  } else {
    message("Fetching data using AnnotationHub() from online database")
    f <- "~/.cache/AnnotationHub/hub_local.rds"
    if(! dir.exists(dirname(f))) {
      dir.create(dirname(f))
    }

    if(file.exists(f)) {
      obj <- readRDS(f)
    } else {
      # hub
      ah = AnnotationHub::AnnotationHub()
      orgdb = AnnotationHub::query(ah, "OrgDb")
      obj <- list(ah = ah,
                  OrgDb = orgdb)
      saveRDS(obj, file = f)
    }

    # load
    obj <- readRDS(f)

    # search
    orgdb <- NULL #output
    hits  <- obj$OrgDb$species[grep(x, obj$OrgDb$species)]
    if(length(hits) < 1) {
      message(paste0("organism not found: ", x))
    } else if (length(hits) > 1) {
      df <- data.frame(organism = hits)
      message(paste0("multiple organism hit: ", length(hits)))
      print(df)
    } else {
      print(paste0("organism found: ", hits))
      db <- AnnotationHub::query(obj$OrgDb, hits)
      orgdb <- obj$ah[[db$ah_id]]
    }
  }
  orgdb
}
















