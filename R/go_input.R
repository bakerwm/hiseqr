# Process input genes for GO analysis
#
# force "entrezid" in analysis pipeline
# force "readable=T" for output
#

#' genes
#' genes with foldChange, or other ranking
#'
#' @param x gene names
#' @param organism string, eg: Homo sapiens
#' @param foldChange numeric, fold change
#' @param to_entrezid boolean, whether convert x to entrezid
#'
#' @export
go_input <- function(x, organism, foldChange = NULL, to_entrezid = TRUE) {
  # OrgDb
  orgdb <- get_orgdb(organism)

  # genes
  if(is.character(x)) {
    x_keytype <- get_gene_keytype(x, organism)

    # force to ENTREZID
    if(to_entrezid & ! x_keytype == "ENTREZID") {
      message("Convert genes to ENTREZID")
      x <- AnnotationDbi::mapIds(orgdb,
                                 keys      = x,
                                 column    = "ENTREZID",
                                 keytype   = x_keytype,
                                 multiVals = "first")
      x_keytype = "ENTREZID"
    }
  } else {
    warning(paste0("character input expected, ", class(x), " got"))
    return(NULL)
  }

  # folcChange
  if(is.null(foldChange)) {
    message("foldChange not detected")
  } else if(is_gene_fc(foldChange)){
    # fc, gene name
    fc_gene <- names(foldChange)
    fc_gene_keytype <- get_gene_keytype(fc_gene, organism)

    # check keytype
    if(to_entrezid & ! fc_gene_keytype == "ENTREZID") {
      names(foldChange) <- AnnotationDbi::mapIds(orgdb,
                                                 keys      = names(foldChange),
                                                 column    = "ENTREZID",
                                                 keytype   = fc_gene_keytype,
                                                 multiVals = "first")
      fc_gene_keytype = "ENTREZID"
    }
  } else {
    # error
    warning(paste0("Numeric foldChanges expected, ", class(foldChange), " got"))
    foldChange <- NULL #
  }

  kegg_code <- get_kegg_code(organism) # ! multiple hits?!

  list(orgdb      = orgdb,
       gene       = x,
       keytype    = x_keytype,
       foldChange = foldChange,
       gsea_gene  = gsea_input(x, foldChange),
       kegg_code  = kegg_code) # update foldChange
}




#' prepare input for gsea
#'
#' @param x gene name
#' @param foldChange numeric
#'
#' @export
gsea_input <- function(x, foldChange) {
  if(is_gene_fc(foldChange)) {
    hits <- x[x %in% names(foldChange)]
    if(length(hits) == 0) {
      warning("genes not found in foldChange names")
      return(NULL)
    } else if(length(hits) < length(x)) {
      miss_pct <- 100 - round(length(hits) / length(x) * 100, 1)
      warning(paste0(miss_pct, "% ", "genes missing foldChange value"))
    } else {
      message("Preparing input for GSEA")
    }

    sort(foldChange[hits], decreasing = TRUE)
  } else {
    warning(paste0("Numeric foldChanges expected, ", class(foldChange), " got"))
    return(NULL)
  }
}




