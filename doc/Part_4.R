## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(vic3PCD)
suppressMessages(library(pheatmap))
suppressMessages(library(kableExtra))
suppressMessages(library(DESeq2))
suppressMessages(library(RColorBrewer))

## ----annotation---------------------------------------------------------------
# Load Novel transcripts annotation
annot_novel <- readRDS(system.file("extdata", "GenesTableFull_cp_NOVEL_annotation.rda", package = "vic3PCD"))

# Load default annotation and join it with novel data
annot <- readRDS(system.file("extdata", "GenesTableFull_cp_annotation.rda", package = "vic3PCD")) %>%
  dplyr::bind_rows(annot_novel) %>%
  tibble::remove_rownames()
# rownames(annot) <- annot$gene_id

## ----data, warning=FALSE, message=FALSE, cache=TRUE---------------------------
# Load SE dataset
se <- readRDS(system.file("extdata", "se_vic3_2020.RData", package = "vic3PCD"))

# DESeq dataset
dds <- DESeqDataSet(se, design = ~ condition)

# To make sure we have right category used as reference in the analysis
dds$condition <- relevel(dds$condition, ref = "Control")

# DESeq analysis
dds <- DESeq(dds, test = "Wald", sfType = "poscounts", useT = FALSE, minReplicatesForReplace = 7)

# Filter genes with more than 10 aligned reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

## ----hc, warning=FALSE, fig.width=6, fig.height=6-----------------------------
# Threshold values
val = 1.9
pval = 0.001
clastTree = 6

groups <- clusterGenes(dds = dds,
                        annotationTbl = annot,
                        summarise_clusters = FALSE,
                        value = val,
                        pvalue = pval,
                        cutTree = clastTree,
                        distRows = "euclidean",
                        clusterMethod = "average")



## ----hc table-----------------------------------------------------------------
tbl <- dplyr::arrange(groups,desc(log2FoldChange))

tbl <- tbl[c(1:10), c("gene_id", "proteinId", "X6_3_EP155mix", "X6_4_EP155mix", "X6_5_EP155mix", "Mean", "baseMean", "log2FoldChange", "pvalue")]

kable(tbl, escape = F, linesep = "", booktabs = T, longtable = F, caption = "Hierarchical Clustering (partial)") %>%
  kable_styling(latex_options = c("scale_down", "repeat_header"), bootstrap_options = c("condensed"), font_size = 12, full_width = F) %>%
    column_spec(1, bold = T)

## ----hc sum, warning=FALSE, fig.width=6, fig.height=6-------------------------
groups <- clusterGenes(dds = dds,
                        annotationTbl = annot,
                        summarise_clusters = TRUE,
                        value = val,
                        pvalue = pval,
                        cutTree = clastTree,
                        distRows = "euclidean",
                        clusterMethod = "average")

