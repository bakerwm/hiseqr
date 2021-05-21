#' Functions for GO analysis, parsing input args
#'
#' force "entrezid" in analysis pipeline
#'
#' @name go_input


#' @describeIn prep_go_input Prepare data for GO analysis
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
prep_go_input <- function(gene_list, organism, ...) {
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
      kegg_gene <- to_kegg_gene_id(gene_list, organism, arg_vars$keytype, simplify = TRUE)
      gsea_gene <- prep_gsea_input(gene_list, arg_vars$fold_change, organism)
      # kegg type for clusterProfiler
      # kegg, (ncbi-geneid)
      kegg_keytype <- ifelse(organism == "Drosophila melanogaster",
                             "ncbi-geneid", "kegg")
      #--Output: list
      list(
        gene_list    = gene_list,
        organism     = organism,
        orgdb        = arg_vars$orgdb,
        keytype      = arg_vars$keytype,
        gsea_gene    = gsea_gene,
        kegg_code    = kegg_code,
        kegg_gene    = kegg_gene,
        kegg_keytype = kegg_keytype
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
#'
#'
#' @param gene_list gene name
#' @param fold_change numeric
#'
#' @export
prep_gsea_input <- function(gene_list, fold_change,
                            organism = NULL, orgdb = NULL) {
  if(! is(gene_list, "character")) {
    warning("`gene_list` expect characters, failed")
    return(NULL)
  }
  if(! is_named_num(fold_change)) {
    warning("`fold_change` require numeric, with gene_name named")
    return(NULL)
  }
  gene_list_fc <- fold_change[gene_list]
  gene_list_fc <- purrr::discard(gene_list_fc, is.na)
  # hit percentage
  pct <- round(length(gene_list_fc) / length(gene_list) * 100, 2)
  if(pct) {
    msg <- glue::glue("{length(gene_list_fc)} of {length(gene_list)}",
                      "({pct}%) genes return with fold_change.")
    message(msg)
    sort(gene_list_fc, decreasing = TRUE)
  } else {
    # try to gueess the keytype of x, and names(fold_change)
    warning("`gene_list` not found in `fold_change` names ",
            "Auto guessing the keytypes:")
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
    #--Guess keytypes
    ge_keytype <- guess_keytype(gene_list, organism)
    fc_keytype <- guess_keytype(names(fold_change), organism)
    if(is(ge_keytype, "character") & is(fc_keytype, "character")) {
      msg <- paste("gene_list: ", ge_keytype,
                   ", fold_change names: ", fc_keytype,
                   sep = "")
      message(msg)
      #--Check: keytypes
      if(ge_keytype == fc_keytype) {
        warning("keytypes of `gene_list` and `fold_change` names identical",
                "No `gene_list` found in `fold_change`")
        return(NULL)
      }
      #--Convert gene_list to fold_change names
      gene_list2 <- convert_id(gene_list, ge_keytype, fc_keytype, organism)
      fc2        <- fold_change[gene_list2]
      fc2        <- purrr::discard(fc2, is.na)
      pct2       <- round(length(fc2) / length(gene_list) * 100, 2)
      if(pct) {
        msg <- paste(length(fc2), " of ", length(gene_list),
                     " (", pct2, ") genes return with fold_change.")
        message(msg)
        sort(fc2, decreasing = TRUE)
      } else {
        warning("`gene_list` non of the genes get fold_change values")
      }
    } else {
      warning("guessing keytypes failed")
    }
  }
}




# prep_gsea_input <- function(x, fold_change, organism = NULL, orgdb = NULL) {
#   if(! is.character(x)) {
#     warning("'x' require characters")
#     return(NULL)
#   }
#   if(! is_named_num(fold_change)) {
#     warning("'fold_change': Expect numeric vector, failed")
#     return(NULL)
#   }
#   if(inherits(orgdb, "OrgDb")) {
#     organism <- metadata(orgdb) %>%
#       filter(name == "ORGANISM") %>%
#       pull(value)
#   } else if(is_character(organism)) {
#     orgdb <- get_orgdb(organism)
#   } else {
#     warning("organism, orgsb, either one required")
#     return(NULL)
#   }
#
#   #--Convert: Sort items, decreasing -------------------------------------------
#   x <- unique(x)
#   x_keytype   <- guess_keytype(x, organism)
#   fc_gene     <- names(fold_change)
#   fc_keytype  <- guess_keytype(fc_gene, organism)
#   if(! x_keytype == fc_keytype) {
#     x2 <- convert_id(fc_gene,
#                      from_keytype = x_keytype,
#                      to_keytype   = fc_keytype,
#                      organism     = organism,
#                      na_rm        = TRUE) %>%
#       pull(2) %>%
#       unique
#     x_fc <- fold_change[x2]
#   } else {
#     x_fc <- fold_change[x]
#   }
#
#   #--Extract: values from fold_change -------------------------------------------
#   # x_fc <- unique(x_fc)
#   if(length(x_fc)) {
#     pct <- round(length(x_fc) / length(x) * 100, 2)
#     msg <- glue::glue("{length(x_fc)} of {length(x)} ({pct}%) genes get fold_change data.")
#     message(msg)
#   } else {
#     warning("genes not found in fold_change, ...")
#     return(NULL)
#   }
#   sort(x_fc, decreasing = TRUE)
# }











# prep_rnaseq_enrich <- function(x) {
#   if(! is_hiseq_multiple_dir(x)) {
#     stop("`x` expect rnaseq_rx type, failed")
#   }
#   px       <- read_hiseq(x)
#   organism <- px$args$genome
#   # px$args$enrich_dir # to-do
#   enrich_dir <- file.path(x, "enrich")  # to-do, deprecated
#   fix_xls    <- get_fix_xls(x)
#   if(is(fix_xls, "character")) {
#     df <- read_text(fix_xls)
#     # required columns
#     req_cols <- c("Gene", "sig")
#     if(! all(req_cols %in% names(df))) {
#       stop("`x` required columns missing: ", paste(req_cols, collapse = ", "))
#     }
#     # sig genes: up, down, not, sig
#     # files: gene_list.xls, require `Gene`
#     sapply(c("sig", "up", "down", "not"), function(s) {
#       subdir   <- file.path(enrich_dir, s)
#       sig_file <- file.path(subdir, "gene_list.xls")
#       if(! dir.exists(subdir)) {
#         dir.create(subdir, recursive = TRUE, mode = "0755")
#       }
#       if(s == "sig") {
#         s = c("up", "down")
#       }
#       df %>%
#         filter(sig %in% s) %>%
#         write.table(file      = sig_file,
#                     sep       = "\t",
#                     row.names = FALSE,
#                     col.names = TRUE)
#       sig_file
#     })
#   }
# }
#





























