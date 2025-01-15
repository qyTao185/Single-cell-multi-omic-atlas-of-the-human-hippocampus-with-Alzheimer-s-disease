library(scater)
library(dplyr)
library(patchwork)
library(SingleCellExperiment)
library(scater)
library(ComplexHeatmap)
# library(ConsensusClusterPlus)
library(msigdbr)
library(fgsea)
library(dplyr)
library(tibble)
library(Signac)
library(ggplot2)
library(stringr)
library(EnsDb.Hsapiens.v86)
library(Seurat)
library(ggplot2)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(ArchR)
library(SingleR)
library(viridis)

set.seed(11)
addArchRGenome("hg38") # hg38, mm9, mm10
addArchRThreads(threads = 64)

# Perform LSI reduction and clustering with ATAC data only
#######################################################################
# Add LSI dimreduc
proj.filter2 <- addIterativeLSI(
  ArchRProj = proj.filter2,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 4,
  LSIMethod = 2,
  scaleDims = T,
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2),
    sampleCells = 10000,
    n.start = 10
  ),
  UMAPParams = list(n_neighbors =30,
                    min_dist = 0.3,
                    metric = "cosine",
                    verbose =FALSE),
  varFeatures = 15000,
  dimsToUse = 1:50,
  binarize = T,
  corCutOff = 0.75,
  force = T,
  seed = 11
)

proj.filter2 <- addClusters(
  input = proj.filter2,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "ATAC_clusters",
  resolution = 0.7,
  dimsToUse = 1:50,force = T
)

# Add UMAP based on LSI dims
proj.filter2 <- addUMAP(proj.filter2,nNeighbors = 30,minDist = 0.3,dimsToUse = 1:50,metric = "cosine",force = T,reducedDims = "IterativeLSI")

###################################################################################################
p1 <- plotEmbedding(ArchRProj = proj.filter2, colorBy = "cellColData", name = "sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj.filter2, colorBy = "cellColData", name = "class", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = proj.filter2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

# Add Gene activity matrix using ArchR model
proj.filter2 <- addGeneScoreMatrix(proj.filter2,matrixName = "ArchRGeneScore",force = T)
getAvailableMatrices(proj.filter2)
proj.filter2=addImputeWeights(proj.filter2)

proj.filter2 <- addHarmony(ArchRProj = proj.filter2,reducedDims = "IterativeLSI",name = "Harmony",groupBy = "Sample",force=T)
proj.filter2 <- addClusters(input = proj.filter2,reducedDims = "Harmony",method = "Seurat",name = "Clusters",resolution = 0.8,force=T,seed=1,nOutlier=160)
proj.filter2 <- addUMAP(
    ArchRProj = proj.filter2, 
    reducedDims = "Harmony", 
    name = "UMAPHarmony", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",force=T
)

p1 <- plotEmbedding(ArchRProj = proj.filter2, colorBy = "cellColData", name = "sample", embedding = "UMAPHarmony")
p2 <- plotEmbedding(ArchRProj = proj.filter2, colorBy = "cellColData", name = "class", embedding = "UMAPHarmony")
p3 <- plotEmbedding(ArchRProj = proj.filter2, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p4 <- plotEmbedding(ArchRProj = proj.filter2, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
p5 <- plotEmbedding(ArchRProj = proj.filter2, colorBy = "cellColData", name = "DoubletEnrichment", embedding = "UMAPHarmony")
p6 <- plotEmbedding(ArchRProj = proj.filter2, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAPHarmony")
p7 <- plotEmbedding(ArchRProj = proj.filter2, colorBy = "cellColData", name = "nFrags", embedding = "UMAPHarmony")
plotPDF(p1,p2,p3,p4,p5,p6,p7, name = "Plot-UMAP-Sample-Clusters-Doublet.pdf", ArchRProj = proj.filter2, addDOC = FALSE, width = 6, height = 6)

proj.filter2 <- addTSNE(
    ArchRProj = proj.filter2, 
    reducedDims = "Harmony", 
    name = "TSNEHarmony", 
    perplexity = 30
)
p3 <- plotEmbedding(ArchRProj = proj.filter2, colorBy = "cellColData", name = "Sample", embedding = "TSNEHarmony")
p4 <- plotEmbedding(ArchRProj = proj.filter2, colorBy = "cellColData", name = "Clusters", embedding = "TSNEHarmony")
# ggAlignPlots(p3, p4, type = "h")

p1 <- plotEmbedding(ArchRProj = proj.filter2, colorBy = "cellColData", name = "sample", embedding = "TSNEHarmony")
p2 <- plotEmbedding(ArchRProj = proj.filter2, colorBy = "cellColData", name = "class", embedding = "TSNEHarmony")
p3 <- plotEmbedding(ArchRProj = proj.filter2, colorBy = "cellColData", name = "Sample", embedding = "TSNEHarmony")
p4 <- plotEmbedding(ArchRProj = proj.filter2, colorBy = "cellColData", name = "Clusters", embedding = "TSNEHarmony")
p5 <- plotEmbedding(ArchRProj = proj.filter2, colorBy = "cellColData", name = "DoubletEnrichment", embedding = "TSNEHarmony")
p6 <- plotEmbedding(ArchRProj = proj.filter2, colorBy = "cellColData", name = "TSSEnrichment", embedding = "TSNEHarmony")
p7 <- plotEmbedding(ArchRProj = proj.filter2, colorBy = "cellColData", name = "nFrags", embedding = "TSNEHarmony")
plotPDF(p1,p2,p3,p4,p5,p6,p7, name = "Plot-TSNE-Sample-Clusters-Doublet.pdf", ArchRProj = proj.filter2, addDOC = FALSE, width = 6, height = 6)
# plotPDF(p3,p4, name = "Plot-TSNE-Sample-Clusters-Doublet.pdf", ArchRProj = proj.filter2, addDOC = FALSE, width = 6, height = 6)

saveArchRProject(ArchRProj = proj.filter2, outputDirectory = "Save-proj.filter2", load = FALSE)