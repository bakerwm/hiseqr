#' Functions for KEGG analysis
#'
#'
#' Usage: KEGGREST functions:
#' keggInfo(), keggList(), keggFind(), keggGet(), keggConv, and keggLink()
#'
#' 1. list all supprted organism
#' keggList("organism")
#'
#' 2. list genes, for one organism
#' keggList("dme")
#'
#' 3. convert ids (idtype)
#' keggConv("ncbi-proteinid", c("hsa:10458", "ece:Z5100"))
#' head(keggConv("eco", "ncbi-geneid")) # all entries
#'
#'
#' @name kegg_utils


#' @description Convert between KEGG id and OrgDb id
#'
#' for most organism:
#' {kegg_code}:{OrgDb_entrezid}
#'
#' for fruitfly:
#' {kegg_code}:{Dmel_}{OrgDb_FLYBASECG}
#' keytype -> FLYBASECG -> dme:Dmel_{CG10000}
#'
#' for prokaryotes
#' do not contains: entrezid, use uniprot instead
#' {up}:{uniprot_id}
#'
#'
#' KEGG supported gene type:
#' ncbi-geneid, ncbi-proteinid, uniprot, flybase (for dme)
#'
#'
#' @param gene_list character geneid
#' @param organism character Name of the organism, support abbr for fruitfly,
#' human, mouse, eg: dm3, dm6, hg38, mouse, mm10, human,
#' @param keytype character or NULL The keytype of gene_list, if NULL,
#' will guess the keytype from keytypes(OrgDb)
#'
#' @export
to_kegg_gene_id <- function(gene_list, organism, keytype = NULL,
                            simplify = TRUE, rm_na = TRUE,
                            return_table = FALSE) {
  #--Check: arg
  if(! is(gene_list, "character")) {
    warning("`gene_list` require character, failed")
    return(NULL)
  }
  #--Check: arg
  organism <- get_organism_name(organism)
  if(! is(organism, "character")) {
    msg <- paste("`organism` not supported:", organism,
                 "; see `Organism.dplyr::supportedOrganisms()` for full list",
                 sep = " ")
    warning(msg)
    return(NULL)
  }
  #--Check: arg
  # convert keytype to "entrezid" (or "flybasecg" for fruit fly)
  if(!is(keytype, "character")) {
    keytype <- guess_keytype(x, organism)
  }
  if(!is_valid_keytype(keytype)) {
    g_str <- paste(gene_list[1:3], collapse = ",")
    message(glue::glue(
      "failed to guess keytype for gene_list: {g_str}"
    ))
    return(NULL)
  }
  # kegg_keytype <- "ENTREZID"
  kegg_keytype <- get_kegg_keytype(organism)
  kegg_code <- get_kegg_code(organism)
  if(keytype == kegg_keytype) {
    message("`gene_list` is KEGG keyType, no need to convert")
    return(gene_list)
  }
  df <- convert_id(gene_list, keytype, kegg_keytype, organism, rm_na = rm_na)
  #-- return data.frame
  # fix Dmel_
  fix1 <- function(x) {
    sapply(x, function(i) {
      if(kegg_keytype == "FLYBASECG") {
        ifelse(is.na(i), NA, paste0("Dmel_", i))
      } else {
        i
      }
    })
  }
  fix2 <- function(x) {
    sapply(x, function(i) {
      ifelse(is.na(i), NA, paste0(kegg_code, ":", i))
    })
  }
  df2 <- df %>%
    dplyr::mutate(across(all_of(kegg_keytype), fix1, .names = "gene_list")) %>%
    dplyr::mutate(across(all_of("gene_list"), fix2, .names = "kegg_gene"))
  #--Output
  if(return_table) {
    if(simplify) {
      out <- dplyr::select(df2, keytype, gene_list)
      if(rm_na) {
        out <- dplyr::filter(out, ! is.na(gene_list))
      }
    } else {
      out <- dplyr::select(df2, keytype, kegg_gene)
      if(rm_na) {
        out <- dplyr::filter(out, ! is.na(kegg_gene))
      }
    }
  } else {
    if(simplify) {
      out <- df2$gene_list
      if(rm_na) {
        out <- purrr::discard(out, is.na)
      }
    } else {
      out <- df2$kegg_gene
      if(rm_na) {
        out <- purrr::discard(out, is.na)
      }
    }
  }
  out
}


#' for clusterProfiler:
#' kegg,
#'
#' for KEGGREST
#' ncbi-geneid
#'
#'
#' for fruit fly: use flybasecg (OrgDb) => ncbi-geneid (kegg in clusterProfiler)
#' for human, mouse: use entrezid (OrgDb) => ncbi-geneid (kegg in clusterProfiler)
#'
#'
#' @param organism character name of the organism
#'
#' @export
get_kegg_keyType <- function(organism) {
  organism <- get_organism_name(organism)
  # ifelse(organism == "Drosophila melanogaster", "FLYBASECG", "ENTREZID")
  ifelse(organism == "Drosophila melanogaster", "ncbi-geneid", "kegg")
  # "ENTREZID"
}


get_kegg_keytype <- function(organism) {
  organism <- get_organism_name(organism)
  # ifelse(organism == "Drosophila melanogaster", "FLYBASECG", "ENTREZID")
  # ifelse(organism == "Drosophila melanogaster", "ncbi-geneid", "kegg")
  "ENTREZID"
}


#' @describeIn kegg_code
#'
#' from: http://www.genome.jp/kegg/catalog/org_list.html
#' using tools: clusterProfiler::search_kegg_organism()
#'
#' using: KEGGREST::keggList("organism")
#'
#' @param x string organism name
#'
#' @import clusterProfiler
#'
#' @export
get_kegg_code <- function(x) {
  organism <- get_organism_name(x) # sci name
  if(!inherits(organism, "character")) {
    warning(glue::glue(
      "x is {x}, is not valid organism name"
    ))
    return(NULL)
  }
  # # Search by KEGGREST::keggList("organism")
  # # up-to 2020-12-20, support 6915 organisms
  # org2 <- KEGGREST::keggList("organism") %>% as.data.frame
  # xs   <- org2 %>%
  #   dplyr::filter(grepl(x, species, ignore.case = TRUE))
  # xa   <- structure(xs$organism, names = xs$species)
  # xai  <- str_similar(x, names(xa))
  # return(xa[xai])
  # search kegg_code by scientific name
  k  <- clusterProfiler::search_kegg_organism("*", by = "scientific_name")
  kc <- setNames(c(k$kegg_code, k$kegg_code),
                 nm = c(k$scientific_name, k$common_name))
  ki <- grep(organism, names(kc), ignore.case = TRUE)
  kk <- kc[ki]
  # output
  if(length(kk) == 1) {
    kk[1]
  } else if(length(kk) > 1) {
    km <- str_similar(x, names(kk), ignore_case = TRUE)
    kk_line <- paste(kk, collapse = ", ")
    msg <- glue::glue("multiple kegg_codes matched: {kk_line}; choose: [{kk[km]}]")
    message(x, ": ", msg)
    kk[km]
  } else {
    msg <- glue::glue("kegg_code not found for: {x}")
    stop(msg)
  }
}





