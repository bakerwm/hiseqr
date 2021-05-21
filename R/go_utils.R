#' Functions for GO analysis
#'
#' including functions for AnnotationDbi, AnnotationHub, ...
#'
#'
#'
#' @name go_utils




#' @describeIn check_go_input Check the input for GO analysis
#'
#' Required:
#' - gene_list, chr, characters (group, enrich)
#' - outdir, chr, character
#' - organism, chr, character
#' - orgdb, OrgDb, AnnotationDbi, one of organism, orgdb required
#'
#' optional:
#' - for_gsea, logical, gene_list, num, sorted, named (GSEA)
#'
#'
#'
#' @export
is_valid_go_input <- function(gene_list, organism, ...) {
  arg_vars <- list(
    outdir   = getwd(),
    for_gsea = FALSE,
    orgdb    = NULL,
    keytype  = NULL
  )
  arg_vars <- purrr::list_modify(arg_vars, !!!list(...))
  arg_vars <- purrr::list_modify(arg_vars, gene_list = gene_list,
                                 organism = organism)
  #----------------------------------------------------------------------------#
  #--Check organism, orgdb, force: orgdb, organism
  if(! is(arg_vars$orgdb, "OrgDb") & is(arg_vars$organism, "character")) {
    arg_vars$orgdb <- get_orgdb(arg_vars$organism)
  }
  if(! is(arg_vars$orgdb, "OrgDb")) {
    msg <- paste(c("`orgdb` and `organism`, not valid;",
                   "orgdb (OrgDb):", class(arg_vars$orgdb), ";",
                   "organism (chr):", class(arg_vars$organism)), collapse = " ")
    stop(msg)
  }
  chk_organism <- is(arg_vars$organism, "character") & length(organism) == 1
  chk_orgdb    <- is(arg_vars$orgdb, "OrgDb")
  #----------------------------------------------------------------------------#
  #--Check: keytype
  if(is.null(arg_vars$keytype)) {
    arg_vars$keytype  <- guess_keytype(arg_vars$gene_list, arg_vars$organism)
  }
  if(! is_valid_keytype(arg_vars$keytype, arg_vars$orgdb)) {
    kt_list <- paste(AnnotationDbi::keytypes(arg_vars$orgdb), collapse = ", ")
    msg     <- paste("`keytype` not valid: [", arg_vars$keytype, "]",
                     "choose from: ", kt_list, sep = " ")
    stop(msg)
  }
  #----------------------------------------------------------------------------#
  #--Check: status
  chk <- c(ifelse(isTRUE(arg_vars$for_gsea),
                  is_named_num(gene_list),
                  is(gene_list, "character")),
           chk_organism | chk_orgdb,
           is(arg_vars$outdir, "character") & length(arg_vars$outdir) == 1
  )
  msg <- paste(
    paste(c("gene_list (chr|num):",
            "organism (chr) or orgdb (OrgDb):",
            "outdir (chr):"),
          chk, sep = " "),
    collapse = "; ")
  if(all(chk)) {
    TRUE
  } else {
    warning(msg)
    FALSE
  }
}




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
#' ...
#'
#' @param x object
#' @param recursive boolean, whether loop over list
#'
#' @export
is_go_result <- function(x) {
  go_class <- c("groupGOResult", "enrichResult", "gseaResult")
  class(x) %in% go_class
}



# is_go <- function(x, recursive = FALSE) {
#   # required
#   go_class <- c("groupGOResult", "enrichResult", "gseaResult")
#
#   # check
#   if(recursive) {
#     if(is.list(x)) {
#       b_list <- sapply(x, function(i){
#         is_go(i, TRUE)
#       })
#       return(all(b_list))
#     } else {
#       return(class(x) %in% go_class)
#     }
#
#   } else {
#     return(class(x) %in% go_class)
#   }
# }



#' @describeIn get_orgdb Pick the organism db
#'
#' @param organism character The name of organism, or build name, eg: dm3, fruitfly
#' Support human, mouse and fruitfly
#'
#' to-do: fetch from bioconductor
#'   laod orgdb
#'   GOSemSim::load_OrgDb()
#'
#' @export
get_orgdb <- function(organism) {
  if(is.character(organism)) {
    # convert to scientific organism name
    org_name <- get_organism_name(organism[1])
    # supported organism list
    sup_org  <- Organism.dplyr::supportedOrganisms()
    orgdb    <- sup_org %>%
      dplyr::filter(organism == org_name) %>%
      dplyr::pull(OrgDb) %>%
      unique
    if (is.character(orgdb)) {
      # function, inspired by: GOSemSim:load_OrgDb
      require(orgdb, character.only = TRUE)
      eval(parse(text = orgdb))
    } else {
      sup_org_list <- paste(
        dplyr::pull(sup_org, organism) %>% unique,
        collapse = ", "
      )
      msg <- paste(
        "Using function `Organism.dplyr::supportedOrganisms()`",
        "to list supported organisms:",
        sup_org_list,
        sep = " "
      )
      warning(msg)
      NULL
    }
  } else {
    warning("`organism` require character, failed")
    NULL
  }
}






#' @describeIn get_organism_name Extract the organism name
#'
#' @param x character the names used in other project
#' @param group character The name return, ["organism", "OrgDb", "TxDb"]
#'
#' @export
get_organism_name <- function(x, group = "organism") {
  # common organism
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
        dplyr::pull(sup_org, organism) %>% unique,
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
  ! is.null(get_organism_name(x))
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
is_valid_keys <- function(orgdb, keys, keytype, fks=NULL,
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
#'
#' @export
is_valid_keytype <- function(x, orgdb = NULL, organism = NULL) {
  if(is(x, "character")) {
    # confirm: OrgDb
    if(is(organism, "character")) {
      orgdb <- get_orgdb(organism)
    }
    if(is(orgdb, "OrgDb")) {
      x %in% AnnotationDbi::keytypes(orgdb)
    }
  }
}










#' @describeIn convert_id Convert gene ids between keytypes, using AnnotationDbi
#'
#' @param x gene names
#' @param organism name of the genome, scientific names, eg: Drosophila melanogaster
#'
#' @export
convert_id <- function(x, from_keytype = NULL, to_keytype = "SYMBOL",
                       organism = NULL, na_rm = FALSE) {
  orgdb <- get_orgdb(organism)
  if(! inherits(orgdb, "OrgDb")) {
    msg1 <- "Not a valid orgnism name, (eg: Homo sapiens)"
    stop("'organism' ", msg1)
  }
  # guess keytype
  if(is.null(from_keytype)) {
    from_keytype <- guess_keytype(x, organism)
  }
  kt <- AnnotationDbi::keytypes(orgdb)
  msg2 <- paste(kt, collapse = ", ")
  if(is.null(from_keytype) | ! from_keytype %in% kt) {
    stop("'from_keytype' should be one of: ", msg2)
  }
  if(! all(to_keytype %in% kt)) {
    stop("'to_keytype' should be one of: ", msg2)
  }
  # out <- AnnotationDbi::mapIds(orgdb,
  #                       keys = x,
  #                       column = to_keytype,
  #                       keytype = from_keytype,
  #                       multiVals = "first")
  if(from_keytype %in% to_keytype & length(to_keytype) == 1) {
    msg <- glue::glue("'from_keytype' and 'to_keytype' identical: {from_keytype}")
    message(msg)
    out <- setNames(data.frame(a = x, b = x), nm = c(from_keytype, to_keytype))
  } else {
    out <- AnnotationDbi::select(
      orgdb,
      keys      = x,
      keytype   = from_keytype,
      columns   = c(from_keytype, to_keytype),
      multiVals = "first")
    # failed rows, keys
    i_na <- which(is.na(out[, 2]))
    if(length(i_na)) {
      n_na  <- out[i_na, 1] %>% unique %>% length
      n_pct <- round(n_na / length(x) * 100, 2)
      if(n_na) {
        msg <- glue::glue("{n_na} of {length(x)} ({n_pct}%) of gene IDs failed to convert ids")
        warning(msg)
      }
      if(na_rm) {
        out <- out[-i_na, ]
      }
    }
  }
  out
}



#' @describeIn convert URL to link
#'
#' markdown: [name](url)
#' html: <a href=url target="_blank">name</a>
#'
#' @export
.url_to_link <- function(url, name, style = "markdown") {
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
                       na_rm        = TRUE) %>%
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
                     na_rm        = TRUE) %>%
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
  .url_to_link(url, name, style = style)
}




#' @describeIn is_named_num numeric data, with NAMES, for gsea input
#'
#' @param x string Numeric, with names
#'
#' @export
is_named_num <- function(x) {
  is(x, "numeric") & is(names(x), "character")
}


















