library(hiseqr)
library(dplyr)
library(ggplot2)

# # library(clusterProfiler)

# data(geneList, package="DOSE")
# gene <- names(geneList)[abs(geneList) > 2]
#
# # kk <- enrichKEGG(gene         = gene,
# #                  organism     = 'hsa',
# #                  pvalueCutoff = 0.05)
#
# kk2 <- gseKEGG(geneList     = geneList,
#                organism     = 'hsa',
#                # nPerm        = 1000,
#                minGSSize    = 120,
#                pvalueCutoff = 0.05,
#                verbose      = FALSE)
# head(kk2)


# #
# setwd("~/work/devel_pipeline/hiseq/rnaseq/results/deseq_compare/")
#
# ## GO
# outdir  <- "go_analysis"
# organism = "Drosophila melanogaster"
# orgdb <- get_orgdb("Drosophila melanogaster")
# gene_list <- head(keys(orgdb), 20)
# foldChange = NULL
# readable = TRUE
# pvalueCutoff = 0.5
# qvalueCutoff = 0.5
# foldChange = setNames(sample(-10:10, length(gene_list)), gene_list)
# tmp <- run_go(gene_list, "Drosophila melanogaster", outdir,
#               level = 2,
#               foldChange = foldChange,
#               readable = TRUE,
#               pvalueCutoff = 0.5,
#               qvalueCutoff = 0.5)
#
#
# ## KEGG
# outdir <- "go_analysis"
# organism = "Drosophila melanogaster"
# orgdb <- get_orgdb(organism)
# gene_list <- head(keys(orgdb), 20)
# pvalueCutoff = 0.5
# qvalueCutoff = 0.5
# tmp <- hiseqr::run_kegg(gene_list, organism, outdir,
#                      foldChange = NULL,
#                      pvalueCutoff = 0.5,
#                      qvalueCutoff = 0.5)



# f <- "~/work/devel_pipeline/hiseq/rnaseq/results/hiseq_read1/pe_control.vs.pe_treatment"
# outdir  <- "go_analysis/demo"
#
# p <- get_sig_gene(f)
# gene_list <- p$up$Gene
# foldChange <- setNames(p$up$log2FoldChange, gene_list)
# organism = "Drosophila melanogaster"
# to_entrezid <- TRUE
# pvalueCutoff = 0.5
# qvalueCutoff = 0.5
# tmp <- hiseqr::run_go(gene_list, organism, outdir,
#                       level = 2,
#                       foldChange   = foldChange,
#                       readable     = TRUE,
#                       pvalueCutoff = 0.5,
#                       qvalueCutoff = 0.5,
#                       to_entrezid  = TRUE)

# tmp <- hiseqr::run_kegg(gene_list, organism, outdir,
#                         foldChange = foldChange,
#                         pvalueCutoff = 0.5,
#                         qvalueCutoff = 0.5)


# deseq_dir <- "~/work/devel_pipeline/hiseq/rnaseq/results/hiseq_read1/pe_control.vs.pe_treatment"
# feature    <- "gene"
# ctl_vs_exp <- TRUE
# rnaseq_pipe(deseq_dir, feature, ctl_vs_exp)



##-----------------------------------------------##
# ## test GO - deseq_dir
# setwd("~/work/devel_pipeline/hiseq/rnaseq/results/deseq_compare/")
# deseq_dir <- "~/work/devel_pipeline/hiseq/rnaseq/results/hiseq_read1/pe_control.vs.pe_treatment"
# feature    <- "gene"
# ctl_vs_exp <- TRUE
# go_pipe(deseq_dir, feature = feature, ctl_vs_exp = ctl_vs_exp)

# ## test GO - genes
# outdir  <- "go_analysis/demo"
# deseq_dir <- "~/work/devel_pipeline/hiseq/rnaseq/results/hiseq_read1/pe_control.vs.pe_treatment
# p <- get_sig_gene(deseq_dir)
# gene_list <- p$up$Gene
# foldChange <- setNames(p$up$log2FoldChange, gene_list)
# organism = "Drosophila melanogaster"
# go_pipe(gene_list, organism, outdir, foldChange = foldChange)



# ##-----------------------------------------------##
# ## test RNAseq compare
# dirA <- "~/work/devel_pipeline/hiseq/rnaseq/results/hiseq_read1/pe_control.vs.pe_treatment/"
# dirB <- "~/work/devel_pipeline/hiseq/rnaseq/results/hiseq_read2/pe_control.vs.pe_treatment/"
# feature <- "gene"
# outdir <- "RNAseq_cmp"
# # pd   <- compare_deseq(dirA, dirB, feature)
# compare_deseq_report(dirA, dirB, feature, "RNAseq_cmp")


# ##----------------------------------------------------------------------------##
# ## debug
# setwd("/data/yulab/wangming/work/devel_pipeline/hiseq/rnaseq/debug")
#
# deseq_dir <- "a.vs.b"
# feature   <- "gene"
# ctl_vs_exp <- T
# rnaseq_pipe(deseq_dir, feature, ctl_vs_exp)


deseq_dir <- "/data/yulab/wangming/work/devel_pipeline/hiseq/rnaseq/results/RNAseq/pe_control.vs.pe_treatment"
feature <- "te"
rnaseq_pipe(deseq_dir, feature)




