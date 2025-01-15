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

scRNA <- readRDS("....rds")
proj.filter4 <- readRDS("....rds")
proj.filter4=addImputeWeights(proj.filter4)

proj.filter5 <- addGeneIntegrationMatrix(
    ArchRProj = proj.filter4, 
    useMatrix = "ArchRGeneScore",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "Harmony",
    seRNA = scRNA,
    addToArrow = FALSE,
    groupRNA = "celltype_merge",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

groupList <- SimpleList()

for (i in names(table(info.data$preClust))){
  
  rna.sub <- scRNA.use[,scRNA.use$celltype_merge == i]
  RNA.cells <- colnames(rna.sub)

  idxSample <- BiocGenerics::which(proj.filter5$predictedGroup_Un == i)
  cellsSample <- proj.filter5$cellNames[idxSample]
  proj.fil <- proj.filter5[cellsSample, ]
  ATAC.cells <- proj.fil$cellNames

  groupList[[i]] <- SimpleList(
    ATAC = ATAC.cells,
    RNA = RNA.cells
  )
}

proj.filter7 <- addGeneIntegrationMatrix(
    ArchRProj = proj.filter5, 
    useMatrix = "ArchRGeneScore",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "Harmony",
    seRNA = scRNA.use,
    addToArrow = TRUE, 
    groupList = groupList,
    groupRNA = "celltype.use",
    nameCell = "predictedCell_Co",
    nameGroup = "predictedGroup_Co",
    nameScore = "predictedScore_Co"
)

p5 <- plotEmbedding(ArchRProj = proj.filter7, colorBy = "cellColData", name = "predictedGroup_Co", embedding = "UMAPHarmony")
p6 <- plotEmbedding(ArchRProj = proj.filter7, colorBy = "cellColData", name = "predictedScore_Co", embedding = "UMAPHarmony")

plotPDF(p5,p6, name = "test.3.pdf", ArchRProj = proj.filter5, addDOC = FALSE, width = 10, height = 10)


proj.filter7 <- addGroupCoverages(ArchRProj = proj.filter7, groupBy = "predictedGroup_Co")
pathToMacs2 <- findMacs2()
proj.filter7 <- addReproduciblePeakSet(
    ArchRProj = proj.filter7, 
    groupBy = "predictedGroup_Co", 
    pathToMacs2 = pathToMacs2
)
proj.filter7 <- addPeakMatrix(proj.filter7,force = T)
proj.filter7 <- addBgdPeaks(proj.filter7)
saveArchRProject(ArchRProj = proj.filter7, outputDirectory = "Save-proj.filter7", load = FALSE)