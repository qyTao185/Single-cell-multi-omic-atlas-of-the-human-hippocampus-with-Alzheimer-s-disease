library(ArchR)
library(parallel)
library(Seurat)
library(Signac)
library(stringr)
library(future)
library(ChIPpeakAnno)
library(stats)
library(BSgenome.Hsapiens.UCSC.hg38)

options(future.globals.maxSize = 50 * 1024^3)
plan("multisession", workers = 22)

set.seed(11)
addArchRGenome("hg38") # hg38, mm9, mm10
addArchRThreads(threads = 64)

proj <- readRDS("...rds")
setwd(".../")

# Add p2g links (no restrictions on FDR, Correlation, Variance cutoff) with raw pvalue
##############################################################################################################
proj.filter1 <- addPeak2GeneLinks(
  ArchRProj = proj ,
  reducedDims = "Harmony",
  useMatrix = "GeneIntegrationMatrix",
  dimsToUse = 1:50,
  scaleDims = NULL,
  corCutOff = 0.75,
  k = 100,
  knnIteration = 500,
  overlapCutoff = 0.8,
  maxDist = 250000,
  scaleTo = 10^4,
  log2Norm = TRUE,
  predictionCutoff = 0.5,
  seed = 1,
  threads = max(floor(getArchRThreads()/2), 1),
  verbose = TRUE,
  logFile = createLogFile("addPeak2GeneLinks")
)

subClusterGroups <- list(
 ...
  )

# subgroups
subClusterCells <- lapply(subClusterGroups, function(x){
  getCellNames(proj)[as.character(proj@cellColData[["group1"]]) %in% x]
  })

for (i in 1:length(subClusterGroups)) {
  proj <- readRDS("...rds")
  proj.filter <- proj[proj$group1 %in% subClusterGroups[i]]
  savepath1=paste0(wd,"/",subClusterGroups[i])
  dir.create(savepath1, showWarnings = FALSE)
  saveArchRProject(ArchRProj = proj.filter, outputDirectory = savepath1, load = TRUE)
  rm(proj)
  gc()
  proj <- readRDS(paste0(savepath1,"/","Save-ArchR-Project.rds"))
  savepath2=paste0(savepath1,"/Output")
  dir.create(savepath2, showWarnings = FALSE)
  
  proj@projectMetadata@listData$outputDirectory=savepath2
  getOutputDirectory(proj)
  set.seed(11)
  addArchRGenome("hg38") # hg38, mm9, mm10
  addArchRThreads(threads = 64)
  
  proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    clusterParams = list(resolution=c(2), sampleCells=30000, maxClusters=6, n.start=10),
    sampleCellsPre = 30000,
    varFeatures = 25000,
    dimsToUse = 1:25,
    force = TRUE
  )
  
  # (Batch correct sample effects)
  proj <- addHarmony(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force = TRUE
  )
  
  # Identify Clusters from Iterative LSI
  proj <- addClusters(
    input = proj,
    reducedDims = "Harmony",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.2,
    force = TRUE
  )
  
  set.seed(1)
  proj <- addUMAP(
    ArchRProj = proj, 
    reducedDims = "Harmony", 
    name = "UMAP", 
    nNeighbors = 35, 
    minDist = 0.4, 
    metric = "cosine",
    force = TRUE
  )
  scRNA <- readRDS("...rds")
  scRNA@meta.data$group1=paste0(scRNA$group,"_",scRNA$celltype.eight)
  scRNA=subset(scRNA,cells=rownames(scRNA@meta.data[which(scRNA$group1==subClusterGroups[i]),]))
  ##########################################################################################
  # Integrate RNA and ATAC data
  ##########################################################################################
  # integration with scRNA:
  set.seed(1)
  proj <- addGeneIntegrationMatrix(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "Harmony",
    seRNA = scRNA, # Can be a seurat object
    addToArrow = TRUE, # add gene expression to Arrow Files (Set to false initially)
    force = TRUE,
    groupRNA = "celltype.use", # used to determine the subgroupings specified in groupList (for constrained integration) Additionally this groupRNA is used for the nameGroup output of this 
    nameCell = "RNA_paired_cell", #Name of column where cell from scRNA is matched to each cell
    nameGroup = "sub.pre", #Name of column where group from scRNA is matched to each cell
    nameScore = "sub.pre.Score" #Name of column where prediction score from scRNA
  )
  proj <- addGroupCoverages(
    ArchRProj=proj, 
    groupBy="predictedGroup_Co", 
    minCells = 50, # The minimum number of cells required in a given cell group to permit insertion coverage file generation. (default = 40)
    force=TRUE
  )
  # Get peaks that were called on this subproject's subclusters from full ArchR project
  full_proj <- readRDS("...rds")
  full_peaks <- getPeakSet(full_proj)
  peaks <- getClusterPeaks(full_proj, clusterNames=unique(proj$group1), peakGR=full_peaks)
  rm(full_proj)
  gc()
  
  # Now add these peaks to the subproject and generate peak matrix
  proj <- addPeakSet(proj, peakSet=peaks, force=TRUE)
  proj <- addPeakMatrix(proj, force=TRUE)
  proj <- addMotifAnnotations(proj, motifSet="cisbp", name="Motif", force=TRUE)
  
  # Calculate coaccessibility
  proj <- addCoAccessibility(
    ArchRProj = proj,
    reducedDims = "Harmony"
  )
  
  # Calculate peak-to-gene links
  proj <- addPeak2GeneLinks(
    ArchRProj = proj,
    reducedDims = "Harmony"
  )
  saveArchRProject(proj)
  rm()
}