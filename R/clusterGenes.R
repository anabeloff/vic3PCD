#################################
### Clustering of DESeq2 data ###
#################################
#' Hierarchical clustering heatmap
#'
#' \code{clusterGenes} Creates a heatmap of differentially expressed genes.
#'
#' The function is useful, but requires many very specific settings. The function estimates differential expression based on \code{results} function from DESeq2.
#' Then it subsets counts matrix (\code{assay(x)}) using threshold values for LFC and then for p-value if provided.
#' Then using data provided in \code{colData(dds)[,"clustering"]} it substracts 'control' samples from 'experiment' generating prelimenary DE values. Then it adds 'means' column containing average values for individual genes.
#' Using sample names specified in \code{colData(dds)[,"clustering"]} as 'experiment' the function then performs hierarchical clustering of differentially expressed genes.
#'
#' @param dds DESeq2 dataset, output of \code{DESeq(x)}.
#' @param value Numeric. LFC threshold to subset dataset based on DESeq2 results values. Taken as absolute value.
#' @param pvalue Numeric. P-value threshhold.
#' @param shrinkLFC Logical. Apply \code{lfcShrink} function to generate DESeq2 results.
#' @param cutTree Integer. Number of clusters to cut tree into.
#' @param distRows Default "euclidean". the distance measure to be used. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring can be given.
#' @param clusterMethod Default "average". Hierarchical clustering the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' @param summarise_clusters Show summarised version of the heatmap. Instead of full heatmap the function will output mean values for each cluster in samples.
#' @param annotationTbl Data frame with 'gene_id' column corresponding to gene IDs in \code{results(dds)} table. The annotation data frame option allows to add additional info to the results output table.
#' @param clusterColumns Character or integer. Specify columns in counts matrix which going to used for clustering. The default setting chooses samples that correspont to \code{colData(dds)[,"clustering"] == 'experiment'}  in 'clustering' column.
#' @param clusterNames Character, optional. Specify names for clusters. Should be the same lenght as 'cutTree'.
#' @param titleExperiment Character. Heatmap main title. If not specified generated automatically.
#' @param test Logical. Output tree only. Can used to determine number of clusters.
#'
#' @examples
#' # Load SE dataset
#' se <- readRDS(system.file("extdata", "se_vic3_2020.RData", package = "vic3PCD"))
#'
#' ## DESeq dataset
#' dds <- DESeqDataSet(se, design = ~ condition)
#' ## To make sure we have right category used as reference in the analysis
#' dds$condition <- relevel(dds$condition, ref = "Control")
#' ## DESeq analysis
#' dds <- DESeq(dds, test = "Wald", sfType = "poscounts", useT = FALSE, minReplicatesForReplace = 7)
#' ## Filter genes with more than 10 aligned reads
#' keep <- rowSums(counts(dds)) >= 10
#' dds <- dds[keep,]
#'
#' ## Annotation Table
#' annot <- readRDS(system.file("extdata", "GenesTableFull_cp_annotation.rda", package = "vic3PCD"))
#'
#' ## Threshold values
#' val = 1.9
#' pval = 0.001
#' clastTree = 6
#'
#' ## Clustering
#' groups <- clusterGenes(dds = dds,
#'                       annotationTbl = annot,
#'                       summarise_clusters = FALSE,
#'                       value = val,
#'                       pvalue = pval,
#'                       cutTree = clastTree,
#'                       distRows = "euclidean",
#'                       clusterMethod = "average")
#######################################

clusterGenes <- function(value = 0,
                         pvalue = 0.1,
                          clusterColumns = NA,
                          summarise_clusters = FALSE,
                          cutTree = NA,
                          distRows = "euclidean",
                          clusterMethod = "average",
                          clusterNames = NULL,
                          annotationTbl = NA,
                          dds = NA,
                          titleExperiment = NA,
                          shrinkLFC = FALSE,
                          test = FALSE) {


  # Default Cluster Names
  if(is.null(clusterNames)) {
    clusterNames <- paste(rep("Cluster", cutTree), c(1:cutTree), sep = "")
  }


  # Default Experiment name
  if(is.na(titleExperiment)) {
    titleExperiment <- base::deparse(base::substitute(dds))
  }
  mainTitle <- paste("Clustering of ", titleExperiment, " Groups: ", cutTree, ".", sep = "")


message("Selecting genes...", appendLF = T)

  # Select data Table

  expressionData <- exprData(value = value, pvalue = pvalue, ds = dds, shrinkLFC = shrinkLFC)
  htmdt <- expressionData[[1]]
  res <- expressionData[[2]]

  #### htmdt selection loop is done
message("done", appendLF = T)

  ## Clustering

  if(is.na(clusterColumns)) {
    coldt <- SummarizedExperiment::colData(dds)
    clusterColumns <- base::row.names(coldt[coldt$clustering == "experiment",])
  }

  dst <- dist(htmdt[,clusterColumns], method = distRows)
  clusters <- hclust(dst, method = clusterMethod)
  clusterCut <- cutree(clusters, cutTree)

  if(test == TRUE) {
    # This part designed to show clustered tree only.
    # Tree is build using the same parameters as in heatmap.
    # It is useful when you need to visualise the tree alone and deside about number of clusters.

    plot(clusters)
    rect.hclust(clusters, k = cutTree, border="blue")
    stop("Choose clusters", call. = FALSE)

  }


message("Creating pheatmap...", appendLF = T)


          ## Heatmap function

          # Annotation for row Clusters
          annot_row <- data.frame(Clusters = factor(clusterCut))
          for(u in c(1:cutTree)) {
            annot_row$Clusters <- sub(paste("^", u, "$", sep = ""), clusterNames[u], annot_row$Clusters)
          }

          # Manual colors for heatmap annotation
          annot <- base::data.frame(SummarizedExperiment::colData(dds))

          #str = RColorBrewer::brewer.pal(n = length(unique(annot$strains)), name ="Greys")
          str = colorRampPalette(brewer.pal(n = 9, name ="Greys"))(length(unique(annot$strains)))
          names(str) <- levels(annot$strains)
          clCut = colorRampPalette(brewer.pal(n = 12, name ="Paired"))(cutTree)
          names(clCut) <- clusterNames
          cond = c("azure3", "azure4")
          names(cond) <- levels(annot$clustering)

          annot_clr <- list(clustering = cond, strains = str, Clusters = clCut)


          # Clusters Summary
          if(summarise_clusters == T) {
            #annot <- annot[annot$clustering == "experiment",]
            annot <- annot[clusterColumns,]

                    main_mt <- data.frame(htmdt[,c(clusterColumns, "Mean")])
                    main_mt$groups <- clusterCut[match(base::row.names(htmdt), names(clusterCut))]
                    for(o in c(1:cutTree)) {
                      main_mt$groups <- gsub(paste("^", o, "$", sep = ""), clusterNames[o], main_mt$groups)
                    }

                    main_mt <- main_mt %>%
                      dplyr::group_by(groups) %>%
                      dplyr::summarise_all(list(mean))

                    rownames(main_mt) <- main_mt$groups


                    annot_row <- data.frame(Clusters = factor(main_mt$groups))
                    base::row.names(annot_row) <- annot_row$Clusters
                    main_mt <- as.matrix(main_mt[,-1])
                    dimnames(main_mt) <- list(annot_row$Clusters, c(base::row.names(annot), "Mean"))

                    main_mt <- main_mt[order(rowMeans(main_mt), decreasing = T),]
                    Variance = round(apply(main_mt, 1, var), digits = 3)
                    main_mt <- cbind(main_mt,Variance)

                    clusters <- F
                    display_numbers <- TRUE

          } else {
                    main_mt <- htmdt[,c(clusterColumns, "Mean")]
                    display_numbers <- FALSE
          }

          # Heat map colors
          minV <- round(min(main_mt))-1
          maxV <- round(max(main_mt))+1
          clr = c(colorRampPalette(rev(brewer.pal(n = 3, name ="PuBu")))(abs(minV)),colorRampPalette(brewer.pal(n = 9, name ="OrRd"))(maxV))

          ## PHEATMAP
          pheatmap(main_mt, show_rownames = F, show_colnames = F, cluster_cols= FALSE, annotation_col = annot[,c("clustering", "strains")], main = mainTitle, border_color = NA,
                   color = clr,
                   breaks = c(minV:-1, 0, 1:maxV),
                   annotation_row = annot_row,
                   annotation_colors = annot_clr,
                   display_numbers = display_numbers, fontsize_number = 10, number_color = "dodgerblue4",
                   cluster_rows = clusters, cutree_rows = cutTree, treeheight_row = 80,
                   cellwidth = 40)

message("done", appendLF = T)



  ressig <- data.frame(htmdt[,c(clusterColumns, "Mean")])

  # Cluster names
  ressig$groups <- clusterCut[match(base::row.names(ressig), names(clusterCut))]

  for(o in c(1:cutTree)) {
    ressig$groups <- gsub(paste("^", o, "$", sep = ""), clusterNames[o], ressig$groups)
  }


  res <- base::as.data.frame(res)
  res$gene_id <- base::row.names(res)
  ressig$gene_id <- base::row.names(ressig)
  ressig <- dplyr::left_join(ressig, res, by = "gene_id")

  # Annotation table join
  if(!is.na(annotationTbl)) {
    colnames(annotationTbl)[1] <- "gene_id"
    ressig <- dplyr::left_join(ressig, annotationTbl, by = "gene_id")
  }

  return(ressig)
}
