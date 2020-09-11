#### PCA plot ####

#' PCA plot
#'
#' \code{vic3_plotPCA} returns PCA plot for normalized counts output by DESeq2.
#'
#' The function works similarly to DESeq2 \code{plotPCA}. It accepts count tables which could be accessed with \code{assay} function.
#'
#' @param counts Counts table output \code{assay(x)} or matrix.
#' @param group_vector Column name in DESeq2 \code{colData(x)}.
#' @param ntop Numeric. Number of top genes with most variance.
#' @param returnData Logical (default FALSE). Return PC1 and PC2 data.frame.
#'
#' @examples
#'
#' ## Load SE dataset
#' se <- readRDS(system.file("extdata", "se_vic3_2020.RData", package = "vic3PCD"))
#'
#' ## DESeq dataset
#' dds <- DESeqDataSet(se, design = ~ condition)
#'
#' ## DESeq analysis
#' dds <- DESeq(dds)
#' ## Filter genes with more than 10 aligned reads
#' keep <- rowSums(counts(dds)) >= 10
#' dds <- dds[keep,]
#' ## Normised counts
#' rld <- rlog(dds, blind = F)
#'
#' ## Manual color
#' clr = c("dodgerblue3", "coral3")
#' names(clr) <- unique(colData(rld)[,"condition"])
#'
#' vic3_plotPCA(counts = assay(rld), group_vector = colData(rld)[,"condition"]) +
#' scale_colour_manual(values=clr, name = "Condition:") +
#' geom_text(label = colData(rld)[,"exset"], size = 3, colour = "black", fontface = "italic", hjust = 0, nudge_x = 3) +
#' xlim(-20,60) +
#' ylim(-25,25)

vic3_plotPCA <- function(counts, group_vector, ntop = 500, returnData = FALSE) {

  # calculate row variance
  rv <- rowVars(counts)
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  # subset counts matrix
  counts <- counts[select,]


  # calculate pca
  pca <- prcomp(t(counts))

  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

  # Plot data frame
  pca_dt <- base::data.frame(pca$x)
  pca_dt$group <- group_vector


  # plot
  if (returnData) {
    attr(pca_dt, "percentVar") <- percentVar[1:2]
    return(pca_dt)
  }

  pcapl <- ggplot(data = pca_dt, aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=4) +
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
    coord_fixed() +
    theme_bw()

  return(pcapl)
}



#### BOXPLOT VARIANCE #####
#' Sample variance box plot
#'
#' \code{vic3_boxplot} returns box plot diagram of varince per column in counts matrix.
#'
#' The function accepts matrix or \code{assay(x)} count table.
#'
#' @param counts Counts table output \code{assay(x)} or matrix.
#' @param ntop, numeric Number of top genes with most variance.
#' @examples
#'
#' se <- readRDS(system.file("extdata", "se_vic3_2020.RData", package = "vic3PCD"))
#'
#' vic3_boxplot(assay(se))

vic3_boxplot <- function(counts, ntop = 500) {

  # calculate row variance
  rv <- rowVars(counts)
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  # subset counts matrix
  counts <- counts[select,]

  boxplot(counts, main="", xlab="", ylab="Raw read counts per gene (log10)",axes=FALSE)
  axis(2)
  axis(1,at=c(1:length(colnames(counts))),labels=colnames(counts),las=2,cex.axis=0.8)
}
