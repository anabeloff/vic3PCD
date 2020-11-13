## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----gff db, eval=FALSE-------------------------------------------------------
#  library("Rsamtools")
#  library("GenomicFeatures")
#  library("GenomicAlignments")
#  
#  GFFFILE <- "cp_genome.gff"
#  SPECIES_NAME <- "Cryphonectria parasitica"
#  TXDB_FILE <- "crypa_annotation.sqlite"
#  
#  Format <- "gff3"
#  
#  transdb<-makeTxDbFromGFF(GFFFILE, format = Format, organism = SPECIES_NAME)
#  
#  saveDb(transdb, TXDB_FILE)

## ----se, eval=FALSE-----------------------------------------------------------
#  library("Rsamtools")
#  library("GenomicFeatures")
#  library("GenomicAlignments")
#  
#  SE_NAME <- "se.RData"
#  TXDB_FILE <- "crypa_annotation.sqlite"
#  
#  # Create BAM files list
#  filenames <-list.files(".", recursive=TRUE, pattern=".bam$", full=TRUE)
#  bamfiles <- BamFileList(filenames, yieldSize=200000)
#  
#  #Read annotation DB
#  transdb <- loadDb(TXDB_FILE)
#  
#  genes <- GenomicFeatures::genes(transdb)
#  
#  se <- summarizeOverlaps(features=genes, reads=bamfiles,
#                          mode="Union",
#                          singleEnd=FALSE,
#                          ignore.strand=TRUE,
#                          fragments=TRUE )
#  
#  
#  saveRDS(se, file = SE_NAME)

