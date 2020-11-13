## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include=FALSE-----------------------------------------------------
library(vic3PCD)
suppressMessages(library(DESeq2))
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggplot2))
suppressMessages(library(kableExtra))

## ----se, warning=FALSE, message=FALSE-----------------------------------------
se <- readRDS(system.file("extdata", "se_vic3_2020.RData", package = "vic3PCD"))


## ----coldt, warning=FALSE, message=FALSE, echo=FALSE--------------------------
tbl <- as.data.frame(colData(se))

kable(tbl, escape = F, linesep = "", booktabs = T, longtable = F, caption = "Summarized Experiment Metadata") %>%
  kable_styling(latex_options = c("scale_down", "repeat_header"), bootstrap_options = c("condensed"), font_size = 12, full_width = F)

## ----dds, warning=FALSE, message=FALSE----------------------------------------
dds <- DESeqDataSet(se, design = ~ condition)

# Collapse replicas
dds <- collapseReplicates(dds, dds$exset)

# To make sure we have right category used as reference in the analysis
dds$condition <- relevel(dds$condition, ref = "Control")

# DESeq analysis
dds <- DESeq(dds, test = "Wald", sfType = "poscounts", useT = FALSE, minReplicatesForReplace = 7)

# Filter genes with more than 10 aligned reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Normised counts
rld <- rlog(dds, blind = F)

## ----pca----------------------------------------------------------------------
# Manual color
clr = c("dodgerblue3", "coral3")
names(clr) <- unique(colData(rld)[,"condition"])

vic3_plotPCA(counts = assay(rld), group_vector = colData(rld)[,"condition"]) +
  scale_colour_manual(values=clr, name = "Condition:") +
  geom_text(label = colData(rld)[,"exset"], size = 3, colour = "black", fontface = "italic", hjust = 0, nudge_x = 3) +
  xlim(-20,60) +
  ylim(-25,25) 

## ----boxplot------------------------------------------------------------------
vic3_boxplot(counts = assay(rld))

