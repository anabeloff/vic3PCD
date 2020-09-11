#' Expression data subset
#'
#' \code{exprData}Returns a list made of count matrix and DESeq2 results data.frame.
#'
#' @param value Numeric. LFC threshold to subset dataset based on DESeq2 results values. Taken as absolute value.
#' @param pvalue Numeric. P-value threshhold.
#' @param ds DESeq2 dataset, output of \code{DESeq(x)}.
#' @param shrinkLFC Logical. Apply \code{lfcShrink} function to generate DESeq2 results.
#' @param calcDE Logical. Calcilate differential expression by substracting 'control' from 'expreriment' columns. To work this option requires \code{colData(dds)} to contain column called 'clustering'. This column must tag samples as 'control' or 'experiment'.

exprData <- function(value = NULL,
                     pvalue = 0.1,
                     ds = NA,
                     shrinkLFC = FALSE,
                     calcDE = TRUE) {
  # Select data Table

  #Extracting assay

  if(is.list(ds)) {
    for (u in 1:length(ds)) {
      nt <- normTransform(ds[[u]])
      tb <- SummarizedExperiment::assay(nt)

      if(u == 1) {
        htmdt <- tb
      } else {
        htmdt <- rbind(htmdt, tb)
      }
    }
  } else {
    nt <- normTransform(ds)
    htmdt <- SummarizedExperiment::assay(nt)
  }


  # Normilising controls

  coldt <- SummarizedExperiment::colData(ds)

  if(isTRUE(calcDE)) {
        controls <- base::row.names(coldt[coldt$clustering == "control",])
        htmdt[,c(controls[1])] <- base::rowMeans(htmdt[,controls])

        #Showing Differential expression in Log2fold change scale for individeal columns.
        experiments <- base::row.names(coldt[coldt$clustering == "experiment",])
        htmdt[,experiments] <- htmdt[,experiments] - htmdt[,c(controls[1])]

        htmdt <- base::cbind(htmdt, base::rowMeans(htmdt[,experiments]))
        colnames(htmdt)[length(colnames(htmdt))] <- "Mean"
  }

  #removing inf numbers
  resInf <- results(ds, tidy = T)
  resInf[is.na(resInf$log2FoldChange),] <- 0
  resInf[resInf$log2FoldChange == -Inf,3] <- 0
  resInf[resInf$log2FoldChange == Inf,3] <- 0

  minV <- min(resInf$log2FoldChange)
  message(base::paste("Minimum value is ", minV, sep = ""), appendLF = T)
  maxV <- max(resInf$log2FoldChange)
  message(base::paste("Maximum value is ", maxV, sep = ""), appendLF = T)



  # Shrinkage option
  if(isTRUE(shrinkLFC)) {
    res <- lfcShrink(ds, coef = resultsNames(ds)[2], type = "apeglm")
  } else {
    res <- DESeq2::results(ds)
    res[is.na(res$log2FoldChange), "log2FoldChange"] <- 0
    res[res$log2FoldChange == -Inf,2] <- minV-1
    res[res$log2FoldChange == Inf,2] <- maxV+1
  }



  # Round LFC
  res$log2FoldChange <- round(res$log2FoldChange, 1)

  # Select for significant genes based on p-value
  res <- res[res$pvalue < pvalue,]

  # Select for significant genes based on LFC.
  res <- res[abs(res$log2FoldChange) >= value,]

  # Subset count matrix
  sel <- rownames(res[order(res$log2FoldChange, decreasing = T),])
  htmdt <- htmdt[sel,]


  htmdt <- base::as.matrix(htmdt)

  finalDT <- list(htmdt, res)

  return(finalDT)
}
