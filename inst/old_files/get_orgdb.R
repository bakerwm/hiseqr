
## check organism



## using AnnotationHub, orgdb/kegg
## KEGG list
#' @param group string, orgdb or kegg, default: orgdb
#' @import clusterProfiler
#' @import AnnotationHub
#' @export
get_species <- function(group = "orgdb") {
  if(group == "orgdb") {
    ah = AnnotationHub::AnnotationHub()
    orgdb = AnnotationHub::query(ah, "OrgDb")
    orgdb$species
  } else if(group == "kegg") {
    clusterProfiler::search_kegg_organism("*") %>%
      dplyr::pull(scientific_name)
  } else {
    warning(paste0("support: orgdb|kegg, input = ", group))
  }
}


##------------------------------------------------------------------------------


#' deprecated, re-write
#' search kegg_code
#' from: http://www.genome.jp/kegg/catalog/org_list.html
#' using tools: clusterProfiler::search_kegg_organism()
#'
#' @param x string organism name
#'
#' @import clusterProfiler
#'
#' @export
get_kegg_code <- function(x) {
  # for regular species
  db_list <- list(human = "hsa",
                  mouse = "mmu",
                  fruitfly = "dme")
  if(grepl("Drosophila melanogaster|fruitfly|dm6|dm3", x)) {
    kegg_code <- db_list$fruitfly
  } else if(grepl("Homo sapiens|human|hg38|hg19", x)) {
    kegg_code <- db_list$human
  } else if(grepl("Mus musculus|mouse|mm10|mm9", x)) {
    kegg_code <- db_list$mouse
  } else {
    message("search KEGG list")
    # search by "scientific_name", "common_name"
    h1 <- clusterProfiler::search_kegg_organism(x, by = "scientific_name")
    h2 <- clusterProfiler::search_kegg_organism(x, by = "common_name")
    # combine
    h <- rbind(h1, h2) %>%
      unique()
    # kegg_code
    kegg_code <- setNames(h$kegg_code, h$scientific_name)
    flag <- paste(kegg_code, collapse = " ")
    # check
    if(nrow(h) == 0) {
      warning(paste0("KEGG organism found:", " [", nrow(h), "] ", " for: ", x))
    } else if(nrow(h) > 1) {
      message(paste0("KEGG organism found:", " [", nrow(h), "] for: ", x, " : ", flag))
    } else {
      message(paste0("KEGG organism found:", " [", nrow(h), "] for: ", x, " : ", flag))
    }
  }
  kegg_code
}





#' deprecated: re-write
## this function from: AnnotationDbi, .testForValidKeys
##
## fks is an alternate vector of keys to consult for validity.
## Normally this will be NULL and the test function should consult
## keys for the supplied keytype
#' @export
is_valid_keys <- function(orgdb, keys, keytype, fks=NULL){
  if (!is.character(keys)){
    warning("'keys' must be a character vector")
    return(FALSE)
  }
  if (length(keys) == 0L) {
    return(FALSE)
  }
  if(is.null(fks)){  ## Normally, fks is just NULL and so we will call keys()
    ktKeys <- AnnotationDbi::keys(orgdb, keytype)
  }else{             ## This lets the caller say wait: use these keys instead
    ktKeys <- fks
  }

  # check pct
  pct <- sum(ktKeys %in% keys) / length(keys)
  if(pct == 0) {
    msg <- paste0("None of the keys are valid keys for ", keytype)
    warning(msg)
    return(FALSE)
  } else if(pct < 0.5){
    msg <- paste0("Less than half of the keys are valid (", pct, ") for ", keytype)
    warning(msg) ## later when things are better, demote this to a warning()
    return(FALSE)
  } else {
    return(TRUE)
  }
}




#' deprecated, see: organism_name
## return scientific name
search_species <- function(x, group = "orgdb") {
  # predata
  sp <- list(fruitfly = "Drosophila melanogaster",
             human    = "Homo sapiens",
             mouse    = "Mus musculus")

  # search
  if(x %in% c("dm3", "dm6", "fruitfly", "Drosophila melanogaster")) {
    sp$fruitfly
  } else if(x %in% c("hg38", "hg19", "human", "Homo sapiens")) {
    sp$human
  } else if(x %in% c("mm10", "mm9", "mouse", "Mus musculus")) {
    sp$mouse
  } else {
    species_list <- get_species(group)
    hits  <- species_list[grep(x, species_list)]
    if(length(hits) < 1) {
      message(paste0("organism not found: ", x))
    } else if (length(hits) > 1) {
      message(paste0("multiple organism matched [",
                     length(hits), "] : ",
                     paste0(hits, collapse = ", ")))
    } else {
      message(paste0("organism matched: ", hits))
    }
    hits
  }
}




#' deprecated, see: guess_keytype
# organism = "Drosophila melanogaster"
#' input is gene names
#' check orgdb, keytypes
#'
#' orgdb
#' @param x string, gene name
#' @param organism string, name of the organism, eg: Homo sapiens
#' @import AnnotationDbi
#'
#' @export
get_gene_keytype <- function(x, organism) {
  orgdb <- get_orgdb(organism)

  # keytypes
  k <- setNames(AnnotationDbi::keytypes(orgdb), AnnotationDbi::keytypes(orgdb))

  # check keytypes
  message("Guess keytype of input genes")
  m <- lapply(k, function(i){
    is_valid_keys(orgdb, x, i)
  })
  m <- unlist(m)
  m <- m[m]

  # output
  if(length(m) == 0) {
    msg <- paste("Unkown gene types, expect x:",
                 paste(AnnotationDbi::keys(orgdb)[1:5], collapse = " "),
                 "...", collapse = " ")
    warning(msg)
    return(NULL)
  } else if(length(m) >= 1) {
    message(paste0("matched keytypes: ",
                   paste0(names(m), collapse = ", "),
                   "\n",
                   "selected: [ ", names(m[1]), " ]"))
    return(names(m[1])) # choose 1
  }
}



#' deprecated, see: convert_id
#' @param x gene names
#' @param organism name of the genome, scientific names, eg: Drosophila melanogaster
#'
#' @export
gene_to_symbol <- function(x, organism, to = "SYMBOL") {
  # pick org.db
  # orgdb <- choose_orgdb(organism)
  orgdb <- get_orgdb(organism)
  keytype <- get_gene_keytype(x, organism)

  # convert
  if(is.null(orgdb)) {
    return(NULL)
  } else if(is.null(keytype)) {
    return(NULL)
  }

  hits <- AnnotationDbi::mapIds(orgdb,
                                keys = x,
                                column = to,
                                keytype = keytype,
                                multiVals = "first")
  # output
  unlist(hits)
}





#' deprecated: is_valid_organism
#' @param x string, supported organism
#' @export
is_valid_organism <- function(x) {
  # orgdb organism
  f1 <- system.file("data", "orgdb_organism.rds", package = "hiseqr")
  pd <- readRDS(f1) # orgdb, kegg

  # alias
  g_alias <- c("dm3", "dm6", "hg19", "hg38", "mm10", "mm9", "fruitfly",
               "fruit fly", "human", "mouse")

  x %in% c(pd$orgdb, g_alias)
}



