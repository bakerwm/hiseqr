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









#' @description
#'
#'
#' Convert between KEGG id and OrgDb id
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
                            simplify = TRUE, na_rm = TRUE,
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
  if(! is(keytype, "character")) {
    keytype <- guess_keytype(x, organism)
  }
  # kegg_keytype <- "ENTREZID"
  kegg_keytype <- get_kegg_keytype(organism)
  kegg_code <- get_kegg_code(organism)
  if(keytype == kegg_keytype) {
    message("`gene_list` is KEGG keyType, no need to convert")
    return(gene_list)
  }

  df <- convert_id(gene_list, keytype, kegg_keytype, organism, na_rm = na_rm)
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
      if(na_rm) {
        out <- dplyr::filter(out, ! is.na(gene_list))
      }
    } else {
      out <- dplyr::select(df2, keytype, kegg_gene)
      if(na_rm) {
        out <- dplyr::filter(out, ! is.na(kegg_gene))
      }
    }
  } else {
    if(simplify) {
      out <- df2$gene_list
      if(na_rm) {
        out <- purrr::discard(out, is.na)
      }
    } else {
      out <- df2$kegg_gene
      if(na_rm) {
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
  # Support for common names, genome build, ... for dm, hs, mm
  dm_name <- c("dm", "dm3", "dm6", "fruitfly", "drosophila_melagaster",
               "drosophila melagaster")
  hs_name <- c("hs", "hg19", "hg38", "GRCh37", "GRCh38", "human",
               "homo_sapiens", "homo sapiens")
  mm_name <- c("mm", "mm9", "mm10", "GRCm38", "mouse", "mus_musculus",
               "mus musculus")
  # data.frame: organism, name
  org <- data.frame(
    organism = c(rep("Drosophila melanogaster", times = length(dm_name)),
                 rep("Homo sapiens", times = length(hs_name)),
                 rep("Mus musculus", times = length(mm_name))),
    name     = c(dm_name, hs_name, mm_name))
  # Convert to scientific name
  if(is.character(x)) {
    x <- tolower(x[1])
    if(x %in% org$name) {
      x <- org[org$name == x, "organism"]
    }
  } else {
    warning("`x` require character, failed")
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
  ki <- grep(x, names(kc), ignore.case = TRUE)
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













#' @describeIn prep_kegg_input Prepare data for kegg analysis
#'
#' 1. Guess the keytype of the gene
#' 2. Convert gene to entrezid
#' 3. add fold_change to gene
#' 4. extract orgdb
#'
#' !important: skip orgdb in arg_vars
#' use saveDb() and loadDb() to save/read OrgDb from file
#' Because the object is a reference to a sqlite data base.
#'
#' genes, orgsnism, orgdb, keytype, gsea_gene
#'
#' genes with fold_change, or other ranking
#' @param gene_list gene names
#' @param organism string, eg: Homo sapiens
#' @param ... extra argument
#'
#' @param fold_change numeric, fold change
#' @param keytype character
#'
#' @export
prep_kegg_input <- function(gene_list, organism, ...) {
  #--Default values: BEGIN
  arg_vars <- rlang::list2(
    fold_change = NULL,
    keytype     = NULL
  )
  arg_vars <- purrr::list_modify(arg_vars, !!!list(...))
  #--Default values: END
  if(is(gene_list, "character")) {
    organism  <- get_organism_name(organism)
    if(is(organism, "character")) {
      arg_vars$orgdb  <- get_orgdb(organism)
      if(! is(arg_vars$keytype, "character")) {
        arg_vars$keytype <- guess_keytype(gene_list, organism)
      }
      kegg_code <- get_kegg_code(organism)
      kegg_gene <- to_kegg_gene_id(gene_list, organism, arg_vars$keytype)
      kegg_gsea_gene <- prep_kegg_gsea_input(gene_list, arg_vars$fold_change,
                                             organism,
                                             keytype = arg_vars$keytype)
      # update: kegg_gene
      kegg_keytype <- get_kegg_keytype(organism) # for OrgDb
      kegg_keyType <- get_kegg_keyType(organism) # for clusterProfiler
      #--Output: list
      list(
        gene_list = kegg_gene,
        organism  = organism,
        keytype   = kegg_keytype, #arg_vars$keytype,
        kegg_keyType = kegg_keyType,
        kegg_code = kegg_code,
        kegg_gene = kegg_gene,
        kegg_gsea_gene = kegg_gsea_gene
      )
    } else {
      warning(paste0("`organism` not valid name: ", organism))
      NULL
    }
  } else {
    warning("`x` : expect character, failed")
    NULL
  }
}





#' @describeIn gsea_input Prepare data for GSEA analysis
#' require sorted values (fold_change, ...), with names (gene)
#' require: ENTREZID (FLYBASECG for fruitfly)
#'
#' @param gene_list gene name
#' @param fold_change numeric
#'
#' @export
prep_kegg_gsea_input <- function(gene_list, fold_change,
                                 organism, keytype = NULL,
                                 orgdb = NULL) {
  if(! is(gene_list, "character")) {
    warning("`gene_list` expect characters, failed")
    return(NULL)
  }
  if(! is_named_num(fold_change)) {
    warning("`fold_change` require numeric, with gene_name named")
    return(NULL)
  }
  #--Check: organism, orgdb
  if(is(orgdb, "OrgDb")) {
    organism <- get_orgdb_metadata(orgdb, "ORGANISM")
  }
  organism <- get_organism_name(organism)
  if(! is(organism, "character")) {
    warning("`organism`, `orgdb` failed, either one required",
            "for guessing the keytypes")
    return(NULL)
  }
  kegg_keytype <- get_kegg_keytype(organism)
  if(keytype == kegg_keytype) {
    fc <- fold_change[gene_list]
    fc <- purrr::discard(fc, is.na)
  } else {
    # id to fold_change
    fix1 <- function(x) {
      fold_change[x]
    }
    df <- to_kegg_gene_id(names(fold_change), organism, keytype, return_table = TRUE) %>%
      dplyr::mutate(across(all_of(keytype), fix1, .names = "value")) %>%
      dplyr::select(-1) %>%
      unique
    fold_change2 <- structure(df[, 2], names = df[, 1])
    gene_list2 <- to_kegg_gene_id(gene_list, organism, keytype)
    fc <- fold_change2[gene_list2]
  }
  pct <- round(length(fc) / length(gene_list) * 100, 2)
  if(pct) {
    msg <- glue::glue("{length(fc)} of {length(gene_list)} ",
                      "({pct}%) genes return with fold_change.")
    message(msg)
    sort(fc, decreasing = TRUE)
  } else {
    warning("`gene_list` not mapped in `fold_change`")
    NULL
  }
}





