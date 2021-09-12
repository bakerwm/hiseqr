#' Functions for GO analysis
#'
#' including functions for AnnotationDbi, AnnotationHub, ...
#'
#'
#'
#' @name go_utils

## @deprecated
#' @describeIn is_named_num Check if input is named numbers
#' eg: the input for GSEA analysis
#'
#'
#' @export
is_named_num <- function(x) {
  is.vector(x) & is.numeric(x) & !is.null(names(x)) & !any(is.na(names(x)))
}


#' @describeIn is_go_result check x, go analysis output
#'
#' support:
#' groupGOResult,
#' enrichResult,
#' enrichKEGG,
#' ...
#'
#' @param x object
#' @param recursive boolean, whether loop over list
#'
#' @export
is_go_result <- function(x, type = TRUE) {
  supp <- c("groupGOResult", "enrichResult", "gseaResult",
            "compareClusterResult")
  if(isTRUE(type)) {
    type <- supp
  } else if(inherits(type, "character")) {
    if(!type[1] %in% supp) {
      supp_str <- paste(supp, collapse = ", ")
      warning(glue::glue(
        "illegal type '{type}', expect {supp_str}"
      ))
      return(FALSE)
    }
  } else {
    warning(glue::glue(
      "illegal type '{class(type})', expect {supp_str}"
    ))
    return(FALSE)
  }
  # class(x) %in% type
  # sapply(x, function(i) inherits(i, type))
  inherits(x, type)
}


#' @describeIn get_orgdb Pick the organism db
#'
#' @param x character The name of organism, or build name, eg: dm3, fruitfly
#' Support human, mouse and fruitfly
#'
#' to-do: fetch from bioconductor
#'   laod orgdb
#'   GOSemSim::load_OrgDb()
#'
#' @export
get_orgdb <- function(x) {
  #----------------------------------------------------------------------------#
  if(!inherits(x, "character")) {
    warning(glue::glue("x is '{class(x)}', expect character"))
    return(NULL)
  }
  #-- empty, na
  if(length(x) == 0 | is.na(x)) {
    message(glue::glue("x is {x}, zero in length, or is NA"))
    return(NULL)
  }
  #-- multiple items
  if(length(x) > 1) {
    x_str <- paste(x[1:3], collapse = ", ")
    message(glue::glue("x, multiple items found, choose the first: {x_str}"))
    x <- x[1]
  }
  #----------------------------------------------------------------------------#
  #-- NA, 0 length
  sci_name <- get_organism_name(x) # to scientific name
  sup_org  <- Organism.dplyr::supportedOrganisms() # OrgDb, TxDb
  r_idx    <- sup_org$organism == sci_name
  if(sum(r_idx) > 0) {
    orgdb <- unique(sup_org[r_idx, ][["OrgDb"]])
    if(length(orgdb) == 1) {
      require(orgdb, character.only = TRUE)
      eval(parse(text = orgdb))
    }
  } else {
    warning(glue::glue(
      "invalid x={x}, ",
      "Using function `Organism.dplyr::supportedOrganisms()` ",
      "to list supported organisms,",
    ))
  }
}


#' @describeIn get_organism_name Extract the organism name
#'
#' @param x character the names used in other project
#' @param group character The name return, ["organism", "OrgDb", "TxDb"]
#'
#' @export
get_organism_name <- function(x, group = "organism") {
  #-- common name to scientific name
  # common organisms
  dm_name <- c("dm", "dm3", "dm6", "fruitfly", "drosophila_melanogaster",
               "drosophila melanogaster")
  hs_name <- c("hs", "hg19", "hg38", "GRCh37", "GRCh38", "human",
               "homo_sapiens", "homo sapiens")
  mm_name <- c("mm", "mm9", "mm10", "GRCm38", "mouse", "mus_musculus",
               "mus musculus")
  sup_org1 <- data.frame(
    organism = c(rep("Drosophila melanogaster",
                     times = length(dm_name)),
                 rep("Homo sapiens",
                     times = length(hs_name)),
                 rep("Mus musculus",
                     times = length(mm_name))),
    name     = c(dm_name, hs_name, mm_name))
  # AnnotationDbi supported organism
  sup_org2 <- Organism.dplyr::supportedOrganisms()
  # Convert to scientific name
  if(is.character(x)) {
    x2 <- tolower(x[1])
    # convert, common name to scientific name
    if(x2 %in% sup_org1$name) {
      x2 <- sup_org1 %>%
        dplyr::filter(name == x2) %>%
        dplyr::pull(organism) %>%
        unique %>%
        head(1)
    }
    # check, scientific name
    if(x2 %in% sup_org2$organism) {
      sup_org2 %>%
        dplyr::filter(organism == x2) %>%
        dplyr::pull(group) %>%
        unique()
    } else {
      sup_org_list <- paste(
        dplyr::pull(sup_org1, organism) %>% unique,
        collapse = ", "
      )
      msg <- paste(
        "`x` [", x, "] unknown;",
        "Using function `Organism.dplyr::supportedOrganisms()`",
        "to list supported organisms:",
        sup_org_list,
        sep = " "
      )
      warning(msg)
      NULL
    }
  }else {
    warning("`x` require character, failed")
    NULL
  }
}


#' to-do: grepl() options, ignore.case
#' @describeIn get_orgdb_info
#'
#' @param x OrgDb
#' @param options character name of metadata, eg: ORGANISM
#'
#' @export
get_orgdb_metadata <- function(x, options = "ORGANISM") {
  if(is(x, "OrgDb")) {
    meta <- metadata(x)
    if(all(options %in% meta$name)) {
      meta %>%
        dplyr::filter(name %in% options) %>%
        dplyr::pull(value)
    } else {
      name_line <- paste(meta$name, collapse = ", ")
      msg <- glue::glue("`options` not valid: [{options}],",
                        "choose from: {name_line}")
      stop(msg)
    }
  } else {
    stop("`x` is not OrgDb")
  }
}


#' @describeIn is_valid_organism
#'
#' @param x string, organism name
#'
#' @export
is_valid_organism <- function(x) {
  !is.null(get_organism_name(x))
}


#' @describeIn guess_keytype Guess the keytype of the genes
#'
#' orgdb
#' @param x string, gene name
#' @param organism string, name of the organism, eg: Homo sapiens
#' @import AnnotationDbi
#'
#' @export
guess_keytype <- function(x, organism = NULL, orgdb = NULL) {
  if(! is(x, "character")) {
    stop("`x` not character, guessing keytype failed.")
  }
  if(! is(orgdb, "OrgDb") & is(organism, "character")) {
    orgdb <- get_orgdb(organism)
  }
  if(! is(orgdb, "OrgDb")) {
    stop("`organism` and `orgdb` are not valid, guessing keytype failed.")
  }
  # Supported keytype
  kt <- AnnotationDbi::keytypes(orgdb)
  message("Guessing keytype:")
  # kt2 <- sapply(kt, function(i){
  #   kt_valid <- is_valid_keys(orgdb, x, i, return_pct = TRUE)
  #   if(isTRUE(kt_valid)) {
  #     i
  #   }
  # }, simplify = TRUE)
  # kt2 <- purrr::discard(kt2, is.null)
  # kt2 <- unlist(kt2)
  # if(length(kt2) > 1) {
  #   kt_line <- paste(kt2, collapse = ", ")
  #   msg <- glue::glue("More than 1 keytypes matched: {kt_line}, choose: {kt2[1]}")
  #   message(msg)
  #   kt2[1]
  # } else if(length(kt2) == 1) {
  #   kt2
  # } else {
  #   keys_line <- paste(keys[1:5], collapse = ", ")
  #   keys_eg   <- paste(AnnotationDbi::keys(orgdb)[1:5], collapse = ", ")
  #   msg <- glue::glue("Unknown keys: {keys_line}; Expect gene names: {keys_eg}")
  #   warning(msg)
  # }
  kt2 <- sapply(kt, function(i) {
    is_valid_keys(orgdb, x, i, return_pct = TRUE)
  })
  kt2 <- kt2[kt2 > 0]
  if(length(kt2) > 1) {
    kt_line <- paste(names(kt2), ":", kt2, "%", sep = "") %>%
      paste(collapse = ", ")
    km <- kt2[which.max(kt2)]
    msg <- glue::glue("Multiple keytypes matched: {kt_line}, choose: [{names(km)}]")
    message(msg)
    names(km)
  } else if(length(kt2) == 1) {
    names(kt2)
  } else {
    x_line    <- x %>% head %>% paste(collapse = ", ")
    keys_line <- paste(AnnotationDbi::keys(orgdb)[1:5], collapse = ", ")
    msg <- glue::glue("`x` is not valid: [{x_line}]; expect example: {keys_line}")
    warning(msg)
    NULL
  }
}


#' @describeIn is_valid_keys check keys is correct keytype in OrgDb
#'
#' Function from AnnotationDbi, .testForValidKeys
#'
#'
#' fks is an alternate vector of keys to consult for validity.
#' Normally this will be NULL and the test function should consult
#' keys for the supplied keytype
#'
#' @param orgdb OrgDb
#' @param keys character gene names
#' @param keytype character keytype
#' @param fks vector
#'
#' @export
is_valid_keys <- function(orgdb, keys, keytype, fks = NULL,
                          return_pct = FALSE){
  if (!is.character(keys)){
    stop("'keys' must be a character vector")
  }
  if (length(keys) == 0L) {
    return()
  }
  if(is.null(fks)){  ## Normally, fks is just NULL and so we will call keys()
    ktKeys <- AnnotationDbi::keys(orgdb, keytype)
  }else{             ## This lets the caller say wait: use these keys instead
    ktKeys <- fks
  }
  keys_valid <- keys[keys %in% ktKeys]
  pct        <- round(length(keys_valid) / length(keys) * 100, 1)
  msg1       <- glue::glue(
    "keytype: {keytype}, {length(keys_valid)} of {length(keys)} ({pct}%)"
  )
  if(pct < 0.5) {
    warning("Less than half of the keys are valid")
    if(return_pct) {
      pct
    }
  } else {
    message(msg1)
    if(return_pct) {
      pct
    } else {
      TRUE
    }
  }
}


#' @describeIn is_valid_keytype Check keytype, from OrgDb, select()
#'
#' @export
is_valid_keytype <- function(x, orgdb = NULL, organism = NULL) {
  out <- FALSE
  if(inherits(x, "character")) {
    # confirm: OrgDb
    if(inherits(organism, "character")) {
      orgdb <- get_orgdb(organism)
    }
    if(is(orgdb, "OrgDb")) {
      out <- x %in% AnnotationDbi::keytypes(orgdb)
    }
  }
  out
}


#' @describeIn convert_id Convert gene ids between keytypes, using AnnotationDbi
#'
#' @param x gene names
#' @param organism name of the genome, eg: "dm6"
#'
#' @export
convert_id <- function(x, from_keytype = NULL, to_keytype = "SYMBOL",
                       organism = NULL, orgdb = NULL, rm_na = FALSE, ...) {
  #----------------------------------------------------------------------------#
  #-- arguments
  dots <- rlang::list2(...)
  args <- list(multi_vars = "first", simplify = FALSE) # multiVars
  dots <- purrr::list_modify(args, !!!dots)
  for(name in names(dots)) {
    assign(name, dots[[name]])
  }
  #----------------------------------------------------------------------------#
  #-- check orgdb, required; (or from organism)
  if(!(inherits(x, "character") & length(x) > 0)) {
    x_str <- paste(x[1:3], collapse = ", ")
    warning(glue::glue(
      "x is '{class(x)}', {x_str} ..., expect 'character'"
    ))
    return(NULL)
  }
  if(!inherits(orgdb, "OrgDb") & is_valid_organism(organism)) {
    orgdb <- get_orgdb(organism)
  }
  if(!inherits(orgdb, "OrgDb")) {
    warning(glue::glue(
      "organism is {organism}, orgdb is {class(orgdb)}, not valid"
    ))
    return(NULL)
  }
  #-- check keytypes
  tk <- sapply(to_keytype, function(i) is_valid_keytype(i, orgdb = orgdb))
  if(!all(tk)) {
    tk_str <- paste(to_keytype, collapse = ", ")
    message(glue::glue(
      "to_keytype is '{tk_str}', not valid"
    ))
    return(NULL)
  }
  #-- guesss, force? in case, from_keytype is not correct
  if(!is_valid_keytype(from_keytype, orgdb)) {
    from_keytype <- guess_keytype(x, orgdb = orgdb)
  }
  if(!is_valid_keytype(from_keytype, orgdb)) {
    x_str <- paste(x[1:3], collapse = ", ")
    warning(glue::glue(
      "x is {x_str} ..., could not determine the keytype"
    ))
    return(NULL)
  }
  #-- run
  if(from_keytype %in% to_keytype & length(to_keytype) == 1) {
    message(glue::glue(
      "'from_keytype' and 'to_keytype' are identical: {from_keytype}"
    ))
    #-- return
    if(simplified) {
      out <- setNames(x, nm = x)
    } else {
      out <- setNames(data.frame(a = x, b = x), nm = c(from_keytype, to_keytype))
    }
  } else {
    #-- run
    args <- purrr::list_modify(
      dots,
      x = orgdb,
      keys = x,
      keytype = from_keytype,
      columns = c(from_keytype, to_keytype),
      multiVals = multi_vars
    )
    out <- rlang::exec(AnnotationDbi::select, !!!args)
    # out <- AnnotationDbi::select(
    #   orgdb,
    #   keys      = x,
    #   keytype   = from_keytype,
    #   columns   = c(from_keytype, to_keytype),
    #   multiVals = "first"
    # )
    # failed rows, keys
    tk_na <- rowSums(as.matrix(apply(out[-1], 2, is.na))) == ncol(out[-1])
    pct   <- round(sum(tk_na) / length(tk_na) * 100, 2)
    # tk_na <- which(is.na(out[[2]])) # na
    if(sum(!tk_na) > 0) {
      message(glue::glue(
        "{sum(tk_na)} of {length(tk_na)} ({pct}%) genes not convert to new keytype"
      ))
      if(rm_na) {
        out <- out[!tk_na, ] # remove na rows
      }
      #-- return
      if(simplify & length(to_keytype) == 1) {
        setNames(out[[to_keytype]], nm = out[[keytype]])
      } else {
        out
      }
    }
  }
}


#' @describeIn convert URL to link
#'
#' markdown: [name](url)
#' html: <a href=url target="_blank">name</a>
#'
#' conflict with: hiseq_reprot.R/.url_to_link2
#'
#' @export
.url_to_link2 <- function(url, name, style = "markdown") {
  if(! is(url, "character") | ! is(url, "character")) {
    stop("`url` and `name` only accept character")
  }
  if(length(url) != length(name)) {
    stop("`url` and `name` not in same length")
  }
  sapply(seq_len(length(url)), function(i) {
    if(style %in% c("md", "markdown")) {
      paste0("[", name[i], "](", url[i], ")")
    } else if(style %in% "html") {
      paste0("<a href='", url[i], "' target='_blank'>", name[i], "</a>")
    } else if(style %in% "url") {
      url[i]
    }
  })
}


#' @describeIn gene_to_link Link to the gene on database (ENSEMBL)
#'
#' create emsembl geneID links
#' format [site]/id/[stable_id]
#' site1: https://asia.ensembl.org
#' site2: https://www.ensembl.org
#'
#' @param x string, gene names
#' @param organism string, scientific name of the organism, eg: Drosophila melanogaster,
#' @param style string, url, markdown, html, ...
#'
#' ensembl: https://asia.ensembl.org/id/FBgn0004872
#' flybase: http://flybase.org/reports/FBgn0004872
#'
#' @export
gene_to_link <- function(x, organism, style = "url",
                         site = "https://asia.ensembl.org",
                         readable = FALSE) {
  # organism, database
  organism_name <- get_organism_name(organism)
  if(is(organism_name, "character")) {
    # redirect to flybase, for Drosophila melanogaster
    if(organism_name == "Drosophila melanogaster") {
      host <- "http://flybase.org/reports"
    } else {
      host <- "https://asia.ensembl.org/id"
    }
  } else {
    warning("`organism` unknown")
    if(is(x, "data.frame")) {
      x$link = NA
    }
    return(x)
  }
  # support, character/data.frame
  if(is(x, "character")) {
    gene_list <- x
  } else if(is(x, "data.frame")) {
    gene_list <- x[[1]] # first column
  } else {
    stop("`x` expect character, data.frame, failed")
  }
  # get id/url
  kt <- guess_keytype(gene_list, organism) # ENSEMBL, FLYBASE
  if(kt %in% c("ENSEMBL", "FLYBASE")) {
    id     <- gene_list
    symbol <- gene_list
    if(isTRUE(readable)) {
      x2 <- convert_id(gene_list,
                       from_keytype = kt,
                       to_keytype   = "SYMBOL",
                       organism     = organism,
                       rm_na        = TRUE) %>%
        unique()
      symbol <- plyr::mapvalues(gene_list, from = x2[[1]], to = x2[[2]],
                                warn_missing = FALSE)
    }
  } else {
    if(organism_name == "Drosophila melanogaster") {
      to_keytype <- c("FLYBASE", "SYMBOL")
    } else {
      to_keytype <- c("ENSEMBL", "SYMBOL")
    }
    x2 <- convert_id(gene_list,
                     from_keytype = kt,
                     to_keytype   = to_keytype,
                     organism     = organism,
                     rm_na        = TRUE) %>%
      unique()
    id     <- plyr::mapvalues(gene_list, from = x2[[1]], to = x2[[2]],
                              warn_missing = FALSE)
    symbol <- plyr::mapvalues(gene_list, from = x2[[1]], to = x2[[3]],
                              warn_missing = FALSE)
  }
  url  <- paste(host, id, sep = "/")
  # readable
  if(isTRUE(readable)) {
    name <- symbol
  } else {
    name <- gene_list
  }
  .url_to_link2(url, name, style = style)
}


