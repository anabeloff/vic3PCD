## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include=FALSE-----------------------------------------------------
library(vic3PCD)
suppressMessages(library(DESeq2))
suppressMessages(library(dplyr))
suppressMessages(library(kableExtra))
suppressMessages(library(tibble))

## ----annotation---------------------------------------------------------------
# Load Novel transcripts annotation
annot_novel <- readRDS(system.file("extdata", "GenesTableFull_cp_NOVEL_annotation.rda", package = "vic3PCD"))

# Load default annotation and join it with novel data
annot <- readRDS(system.file("extdata", "GenesTableFull_cp_annotation.rda", package = "vic3PCD")) %>%
  dplyr::bind_rows(annot_novel) %>%
  tibble::remove_rownames()
rownames(annot) <- annot$gene_id

# Random rows
smp = sample(1:length(annot$gene_id), 5)

kable(annot[smp,c("proteinId", "gotrm", "goName", "UP_proteinId", "Protein.names")], escape = F, linesep = "", booktabs = T, longtable = F, caption = "Annotation data set (partial)") %>%
  kable_styling(latex_options = c("scale_down", "repeat_header"), bootstrap_options = c("condensed"), font_size = 12, full_width = F)

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

## ----results dt, fig.width=6, fig.height=4------------------------------------
# DESeq2 results table
res <- DESeq2::results(dds, tidy = F, name = "condition_Barrage_vs_Control")

# MA plot
DESeq2::plotMA(res, ylim=c(-6,12), main = "MA plot")
abline(h=c(-2,2), col="dodgerblue", lwd=2) 


## ----de genes-----------------------------------------------------------------
### Selection criteria
# P-value
pval = 0.001
# LFC
lfc = 1.9

### DE genes table joined with annotation
res_dt <- data.frame(gene_id = row.names(res), res) %>%
  dplyr::mutate(log2FoldChange = round(log2FoldChange, 1)) %>%
  dplyr::filter(pvalue < pval, abs(log2FoldChange) >= lfc) %>%
  dplyr::left_join(annot, by = "gene_id") %>%
  dplyr::arrange(desc(log2FoldChange))


tbl <- res_dt[c(1:10), c("gene_id", "baseMean", "log2FoldChange", "pvalue", "proteinId", "UP_proteinId", "Protein.names")]

kable(tbl, escape = F, linesep = "", booktabs = T, longtable = F, caption = "Annotation data set (partial)") %>%
  kable_styling(latex_options = c("scale_down", "repeat_header"), bootstrap_options = c("condensed"), font_size = 12, full_width = F)

