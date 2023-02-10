
## ---------------------------
## The script process bulk RNA-seq data from Brouiller et al..
##
## The script loads raw and TPM counts of analyzed samples. It then filters matrices to keep relevant genes
## and then performs differential expression analysis using DESeq2.
##
## Author: Pierre-Emmanuel Bonté
##
## ---------------------------

# Loading packages
library(tidyverse)
library(data.table)
library(DESeq2)

# DESeq2 function
process_deseq2_results = function(dds, factor, grp1, grp2, rownames_name =  "gene_name", lfc = FALSE){
  # Without foldchange shrinkage
  if(lfc == FALSE){
    df = DESeq2::results(dds, contrast = c(factor, grp1, grp2) )
    df = as.data.frame(df) %>%
      tibble::rownames_to_column(var = rownames_name) %>%
      dplyr::mutate(cluster = ifelse(log2FoldChange >= 0, yes = grp1, no = grp2),
                    Significant = ifelse(padj <= 0.05 & (!is.na(padj)), yes = "Significant", no = "Not signigficant"),
                    comparison = paste0(grp1, "_vs_", grp2))
    
  # With apeglm foldchange shrinkage
  }else{
    dds <- nbinomWaldTest(dds)
    dds@colData[,factor] = relevel(dds@colData[,factor], ref = grp2)
    dds <- nbinomWaldTest(dds)
    coeff_id = which(resultsNames(dds) %in% paste0(factor, "_", grp1, "_vs_", grp2))
    df = lfcShrink(dds , type="apeglm", coef = coeff_id)
    df = as.data.frame(df) %>%
      rownames_to_column(var = rownames_name) %>%
      dplyr::mutate(cluster = ifelse(log2FoldChange >= 0, yes = grp1, no = grp2),
                    Significant = ifelse(padj <= 0.05 & (!is.na(padj)),  yes = "Significant",  no = "Not signigficant"),
                    comparison = paste0(grp1, "_vs_", grp2))
  }
}


# Loading metadata
metadata <- data.table::fread("github_files/metadata_bulk_samples.tsv")
rownames(metadata) <- metadata$Sample
hg38.genes.annotations = readRDS("github_files//hg38.genes.annotations.rds") %>%
  dplyr::mutate(ENSEMBL = gsub("\\..+", "", Gene_ID)) %>%
  dplyr::select(ENSEMBL, Gene_name)

# Loading data
DC.datasets.raw.counts = readRDS("github_files/DC.datasets.raw.counts.rds")
DC.datasets.TPM.counts = readRDS("github_files/DC.datasets.tpm.counts.rds")

# Reformating expression
DC.datasets.TPM.counts.df = DC.datasets.TPM.counts %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_name") %>%
  tidyr::gather(`Sample`, `TPM`, -gene_name) %>%
  left_join(y = metadata, by = "Sample") 


# Filtering to keep features with at least 1 TPM in all samples from 1 condition
list_genes_1_TPM = DC.datasets.TPM.counts.df %>%
  dplyr::group_by(factor, gene_name) %>%
  dplyr::summarise(above_threshold = all(TPM >= 1) )  %>%
  dplyr::filter(above_threshold == TRUE) %>%
  dplyr::pull(gene_name) %>%
  unique()

DC.datasets.raw.counts.filtered = DC.datasets.raw.counts[list_genes_1_TPM,]
DC.datasets.TPM.counts.filtered = DC.datasets.TPM.counts[list_genes_1_TPM,]

# Differential expression analysis 
# Creating DESeq2 object + performing DEA
dds.factor = DESeqDataSetFromMatrix(countData = DC.datasets.raw.counts.filtered ,
                                    colData = metadata,
                                    design = ~ factor + Donor)
dds.factor = estimateSizeFactors(dds.factor)
dds.factor = DESeq(dds.factor)

# Retrieving DESeq2 normalized counts
DC.datasets.norm.counts.filtered = as.data.frame(counts(dds.factor, normalized = TRUE))


# Retrieving and formating DEA results
# MOCK vs HIV
res_all.pdc_AD8.vs.mock = process_deseq2_results(dds = dds.factor.all, 
                                                 factor = "factor", 
                                                 grp1 = "pDC_AD8", 
                                                 grp2 = "pDC_mock", 
                                                 rownames_name = "gene_name", lfc = TRUE)

res_all.cdc1_AD8.vs.mock = process_deseq2_results(dds = dds.factor.all, 
                                                  factor = "factor", 
                                                  grp1 = "cDC1_AD8", 
                                                  grp2 = "cDC1_mock", 
                                                  rownames_name = "gene_name", lfc = TRUE)

res_all.predc_AD8.vs.mock = process_deseq2_results(dds = dds.factor.all, 
                                                   factor = "factor", 
                                                   grp1 = "AxlDC_AD8", 
                                                   grp2 = "AxlDC_mock", 
                                                   rownames_name = "gene_name", lfc = TRUE)

res_all.predc_AD8vpx.vs.mock = process_deseq2_results(dds = dds.factor.all, 
                                                      factor = "factor", 
                                                      grp1 = "AxlDC_AD8vpx", 
                                                      grp2 = "AxlDC_mock", 
                                                      rownames_name = "gene_name", lfc = TRUE)

res_all.predc_AZTNVP.vs.mock = process_deseq2_results(dds = dds.factor.all, 
                                                      factor = "factor", 
                                                      grp1 = "AxlDC_AZTNVP", 
                                                      grp2 = "AxlDC_mock", 
                                                      rownames_name = "gene_name", lfc = TRUE)


res_all.cdc2_AD8.vs.mock = process_deseq2_results(dds = dds.factor.all,
                                                  factor = "factor",
                                                  grp1 = "cDC2_AD8",
                                                  grp2 = "cDC2_mock",
                                                  rownames_name = "gene_name", lfc = TRUE)

res_all.cdc2_AD8vpx.vs.mock = process_deseq2_results(dds = dds.factor.all, 
                                                     factor = "factor", 
                                                     grp1 = "cDC2_AD8vpx", 
                                                     grp2 = "cDC2_mock", 
                                                     rownames_name = "gene_name", lfc = TRUE)


res_all.cdc2_AZTNVP.vs.mock = process_deseq2_results(dds = dds.factor.all, 
                                                     factor = "factor", 
                                                     grp1 = "cDC2_AZTNVP", 
                                                     grp2 = "cDC2_mock", 
                                                     rownames_name = "gene_name", lfc = TRUE)

# AD8 : DCs pop comparisons
res_all.AD8_predc.vs.pdc = process_deseq2_results(dds = dds.factor.all, 
                                                  factor = "factor", 
                                                  grp1 = "AxlDC_AD8", 
                                                  grp2 = "pDC_AD8", 
                                                  rownames_name = "gene_name", lfc = TRUE)

res_all.AD8_predc.vs.cdc2 = process_deseq2_results(dds = dds.factor.all, 
                                                   factor = "factor", 
                                                   grp1 = "AxlDC_AD8", 
                                                   grp2 = "cDC2_AD8", 
                                                   rownames_name = "gene_name", lfc = TRUE)

res_all.AD8_predc.vs.cdc1 = process_deseq2_results(dds = dds.factor.all, 
                                                   factor = "factor", 
                                                   grp1 = "AxlDC_AD8", 
                                                   grp2 = "cDC1_AD8", 
                                                   rownames_name = "gene_name", lfc = TRUE)


res_all.AD8_pdc.vs.cdc1 = process_deseq2_results(dds = dds.factor.all, 
                                                 factor = "factor", 
                                                 grp1 = "pDC_AD8", 
                                                 grp2 = "cDC1_AD8", 
                                                 rownames_name = "gene_name", lfc = TRUE)


res_all.AD8_pdc.vs.cdc2 = process_deseq2_results(dds = dds.factor.all, 
                                                 factor = "factor", 
                                                 grp1 = "pDC_AD8", 
                                                 grp2 = "cDC2_AD8", 
                                                 rownames_name = "gene_name", lfc = TRUE)


res_all.AD8_cdc1.vs.cdc2 = process_deseq2_results(dds = dds.factor.all, 
                                                  factor = "factor", 
                                                  grp1 = "cDC1_AD8", 
                                                  grp2 = "cDC2_AD8", 
                                                  rownames_name = "gene_name", lfc = TRUE)


# MOCK : DCs pop comparisons
res_all.mock_predc.vs.pdc = process_deseq2_results(dds = dds.factor.all, 
                                                   factor = "factor", 
                                                   grp1 = "AxlDC_mock", 
                                                   grp2 = "pDC_mock", 
                                                   rownames_name = "gene_name", lfc = TRUE)

res_all.mock_predc.vs.cdc2 = process_deseq2_results(dds = dds.factor.all, 
                                                    factor = "factor", 
                                                    grp1 = "AxlDC_mock", 
                                                    grp2 = "cDC2_mock", 
                                                    rownames_name = "gene_name", lfc = TRUE)

res_all.mock_predc.vs.cdc1 = process_deseq2_results(dds = dds.factor.all, 
                                                    factor = "factor", 
                                                    grp1 = "AxlDC_mock", 
                                                    grp2 = "cDC1_mock", 
                                                    rownames_name = "gene_name", lfc = TRUE)


res_all.mock_pdc.vs.cdc1 = process_deseq2_results(dds = dds.factor.all, 
                                                  factor = "factor", 
                                                  grp1 = "pDC_mock", 
                                                  grp2 = "cDC1_mock", 
                                                  rownames_name = "gene_name", lfc = TRUE)


res_all.mock_pdc.vs.cdc2 = process_deseq2_results(dds = dds.factor.all, 
                                                  factor = "factor", 
                                                  grp1 = "pDC_mock", 
                                                  grp2 = "cDC2_mock", 
                                                  rownames_name = "gene_name", lfc = TRUE)


res_all.mock_cdc1.vs.cdc2 = process_deseq2_results(dds = dds.factor.all, 
                                                   factor = "factor", 
                                                   grp1 = "cDC1_mock", 
                                                   grp2 = "cDC2_mock", 
                                                   rownames_name = "gene_name", lfc = TRUE)


# Saving results into lists
list.res = list(res_all.pdc_AD8.vs.mock,
                res_all.cdc1_AD8.vs.mock,
                res_all.predc_AD8.vs.mock,
                res_all.predc_AD8vpx.vs.mock,
                res_all.predc_AZTNVP.vs.mock,
                res_all.cdc2_AD8.vs.mock,
                res_all.cdc2_AD8vpx.vs.mock,
                res_all.cdc2_AZTNVP.vs.mock)

names(list.res) <- c("pdc_AD8.vs.mock",
                     "cdc1_AD8.vs.mock",
                     "AxlDC_AD8.vs.mock",
                     "AxlDC_AD8vpx.vs.mock",
                     "AxlDC_AZTNVP.vs.mock",
                     "cdc2_AD8.vs.mock",
                     "cdc2_AD8vpx.vs.mock",
                     "cdc2_AZTNVP.vs.mock")


list.res2 = list(res_all.AD8_predc.vs.pdc,
                 res_all.AD8_predc.vs.cdc2,
                 res_all.AD8_predc.vs.cdc1,
                 res_all.AD8_pdc.vs.cdc1,
                 res_all.AD8_pdc.vs.cdc2,
                 res_all.AD8_cdc1.vs.cdc2)

names(list.res2) <- c("AD8_AxlDC.vs.pdc",
                      "AD8_AxlDC.vs.cdc2",
                      "AD8_AxlDC.vs.cdc1",
                      "AD8_pdc.vs.cdc1",
                      "AD8_pdc.vs.cdc2",
                      "AD8_cdc1.vs.cdc2")

list.res3 = list(res_all.mock_predc.vs.pdc,
                 res_all.mock_predc.vs.cdc2,
                 res_all.mock_predc.vs.cdc1,
                 res_all.mock_pdc.vs.cdc1,
                 res_all.mock_pdc.vs.cdc2,
                 res_all.mock_cdc1.vs.cdc2)

names(list.res3) <- c("mock_AxlDC.vs.pdc",
                      "mock_AxlDC.vs.cdc2",
                      "mock_AxlDC.vs.cdc1",
                      "mock_pdc.vs.cdc1",
                      "mock_pdc.vs.cdc2",
                      "mock_cdc1.vs.cdc2")


# Filtering and keeping only significant differentially expressed genes
list.res.sig = lapply(list.res, function(x){
  x <- x %>% dplyr::filter(padj <= 0.05 & abs(log2FoldChange) >= 1 )
})

list.res2.sig = lapply(list.res2, function(x){
  x <- x %>% dplyr::filter(padj <= 0.05 & abs(log2FoldChange) >= 1 )
})

list.res3.sig = lapply(list.res3, function(x){
  x <- x %>% dplyr::filter(padj <= 0.05 & abs(log2FoldChange) >= 1 )
})


# Venn diagramm
require(VennDiagram)
venn.diagram(list(AxlDC = list.res.sig$AxlDC_AD8.vs.mock %>%
                    dplyr::filter(log2FoldChange > 1) %>%
                    dplyr::pull(gene_name), 
                  pDC = list.res.sig$pdc_AD8.vs.mock %>%
                    dplyr::filter(log2FoldChange > 1) %>%
                    dplyr::pull(gene_name),
                  cDC2 = list.res.sig$cdc2_AD8.vs.mock %>%
                    dplyr::filter(log2FoldChange > 1) %>%
                    dplyr::pull(gene_name),
                  cDC1 = list.res.sig$cdc1_AD8.vs.mock %>%
                    dplyr::filter(log2FoldChange > 1) %>%
                    dplyr::pull(gene_name)),
             fill = c("#930F00", "#2695FE", "#FA9201", "#9337FF"),
             resolution = 160, 
             disable.logging = TRUE, 
             lwd =0, 
             height = 10,
             width = 10,
             units = "in",
             alpha = c(0.8, 0.8, 0.8, 0.8), 
             cex = 7, 
             fontfamily ="Arial")


