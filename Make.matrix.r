library(ArchR)
library(RcppAlgos)
library(parallel)
library(Seurat)
# library(Signac)
library(stringr)
library(future)
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)


# P2G definition cutoffs
corrCutoff <- 0.45       # Default in plotPeak2GeneHeatmap is 0.45
varCutoffATAC <- 0.25   # Default in plotPeak2GeneHeatmap is 0.25
varCutoffRNA <- 0.25    # Default in plotPeak2GeneHeatmap is 0.25

# Coaccessibility cutoffs
coAccCorrCutoff <- 0.5  # Default in getCoAccessibility is 0.5


full.proj=readRDS("...rds")

# Get all peaks
allPeaksGR <- getPeakSet(full.proj)
allPeaksGR$peakName <- (allPeaksGR %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
names(allPeaksGR) <- allPeaksGR$peakName

wd <- "..."
full_p2gGR=readRDS(paste0(wd, "/multilevel_p2gGR.rds")) # NOT merged or correlation filtered
plot_loop_list=readRDS(paste0(wd, "/multilevel_plot_loops.rds"))
full_coaccessibility=readRDS(paste0(wd, "/multilevel_coaccessibility.rds"))

# Get full peak matrix
pSE <- getMatrixFromProject(full.proj, useMatrix="PeakMatrix", binarize=FALSE)
pmat <- assays(pSE)$PeakMatrix
rownames(pmat) <- rowRanges(pSE) %>% {paste0(seqnames(.), "_", start(.), "_", end(.))}
pmat <- pmat[,getCellNames(full.proj)] # Need to get cells in order of ArchR project...

# Get GeneScore matrix
gsmSE <- getMatrixFromProject(full.proj, useMatrix="GeneScoreMatrix", binarize=FALSE)
gsmat <- assays(gsmSE)$GeneScoreMatrix
rownames(gsmat) <- rowData(gsmSE)$name
gsmat <- gsmat[,getCellNames(full.proj)] # Need to get cells in order of ArchR project...

# Get GeneIntegration matrix
gimSE <- getMatrixFromProject(full.proj, useMatrix="GeneIntegrationMatrix", binarize=FALSE)
gimat <- assays(gimSE)$GeneIntegrationMatrix
rownames(gimat) <- rowData(gimSE)$name
gimat <- gimat[,getCellNames(full.proj)] # Need to get cells in order of ArchR project...

saveRDS(pmat,file="...pmat.rds")
saveRDS(gsmat,file="...gsmat.rds")
saveRDS(gimat,file="...gimat.rds")