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

# Create Arrow and ArchR project
##########################################################################
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = sampleNames,
  minTSS = 0, #Dont set this too high because you can always increase later
  minFrags = 0,
  addTileMat = T,
  addGeneScoreMat = F
)


doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP",useMatrix = "TileMatrix",nTrials=5,LSIMethod = 1,scaleDims = F,
  corCutOff = 0.75,UMAPParams = list(n_neighbors =30, min_dist = 0.3, metric = "cosine", verbose =FALSE),
  dimsToUse = 1:50
)


proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory =".../Output",
  copyArrows = T #This is recommened so that you maintain an unaltered copy for later usage.
)

# Filter out outlier low quality cells and doublets
###############################################################################
# GMM for fragments per cell
library(mclust)
dir.create("QC")
setwd(".../QC")
for (i in sampleNames){
  proj.i <- proj[proj$Sample == i]

  # GMM for fragments per cell
  depth.clust <- Mclust(log10(proj.i$nFrags),G = 2)
  proj.i$depth.cluster <- depth.clust$classification
  proj.i$depth.cluster.uncertainty <- depth.clust$uncertainty

  ggPoint(
    x = log10(proj.i$nFrags),
    y = log10(proj.i$TSSEnrichment+1),
    color = as.character(proj.i$depth.cluster),
    xlabel = "log10(unique fragments)",
    ylabel = "log10(TSS Enrichment+1)"
  ) + ggtitle(paste0("GMM classification:\n",i," log10(fragments)"))
    ggsave(paste0(i,"_depth.pdf"),width = 4,height = 4)

  # GMM for TSS per cell
  TSS.clust <- Mclust(log10(proj.i$TSSEnrichment+1),G = 2)
  proj.i$TSS.cluster <- TSS.clust$classification
  proj.i$TSS.cluster.uncertainty <- TSS.clust$uncertainty

  ggPoint(
    x = log10(proj.i$nFrags),
    y = log10(proj.i$TSSEnrichment+1),
    color = as.character(proj.i$TSS.cluster),
    discrete = T,
    xlabel = "log10(unique fragments)",
    ylabel = "log10(TSS Enrichment+1)"
  ) + ggtitle(paste0("GMM classification:\n",i," TSS Enrichment"))
    ggsave(paste0(i,"_TSS.pdf"),width = 4,height = 4)


  df.TSS <- data.frame(proj.i$cellNames,proj.i$TSS.cluster,proj.i$TSS.cluster.uncertainty,proj.i$TSSEnrichment)
  df.TSS <- dplyr::filter(df.TSS,proj.i.TSS.cluster == "2")
  df.TSS <- dplyr::filter(df.TSS,proj.i.TSS.cluster.uncertainty <= 0.05)
  saveRDS(df.TSS,paste0("df_TSS_",i,".rds"))

  df.depth <- data.frame(proj.i$cellNames,proj.i$depth.cluster,proj.i$depth.cluster.uncertainty,proj.i$nFrags)
  df.depth <- dplyr::filter(df.depth,proj.i.depth.cluster == "2")
  df.depth <- dplyr::filter(df.depth,proj.i.depth.cluster.uncertainty <= 0.05)
  saveRDS(df.depth,paste0("df_depth_",i,".rds"))

  ggPoint(
    x = log10(proj.i$nFrags),
    y = log10(proj.i$TSSEnrichment+1),
    colorDensity = T,
    continuousSet = "sambaNight",
    xlabel = "log10(unique fragments)",
    ylabel = "log10(TSS Enrichment+1)"
  ) +geom_hline(yintercept = log10(min(df.TSS$proj.i.TSSEnrichment)+1),linetype = "dashed")+
    geom_vline(xintercept = min(log10(df.depth$proj.i.nFrags)),linetype = "dashed")+
    ggtitle(paste0("QC thresholds:\n",i))
    ggsave(paste0(i,"_QC.pdf"),width = 4,height = 4)

  ggPoint(
    x = log10(proj.i$nFrags),
    y = log10(proj.i$TSSEnrichment+1),
    color = proj.i$DoubletEnrichment,
    discrete = F,
    continuousSet = "sambaNight",
    xlabel = "log10(unique fragments)",
    ylabel = "log10(TSS Enrichment+1)"
  ) +geom_hline(yintercept = min(log10(df.TSS$proj.i.TSSEnrichment+1)),linetype = "dashed")+
    geom_vline(xintercept = min(log10(df.depth$proj.i.nFrags)),linetype = "dashed")+
    ggtitle(paste0("Doublet Enrichment:\n",i))
    ggsave(paste0(i,"_doublets.pdf"),width = 4,height = 4)

}

i= library_id
library(mclust)
proj.i <- proj[proj$Sample == i]

# GMM for fragments per cell
depth.clust <- Mclust(log10(proj.i$nFrags),G = 2)
proj.i$depth.cluster <- depth.clust$classification
proj.i$depth.cluster.uncertainty <- depth.clust$uncertainty

TSS.clust <- Mclust(log10(proj.i$TSSEnrichment+1),G = 2)
proj.i$TSS.cluster <- TSS.clust$classification
proj.i$TSS.cluster.uncertainty <- TSS.clust$uncertainty

df.TSS <- data.frame(proj.i$cellNames,proj.i$TSS.cluster,proj.i$TSS.cluster.uncertainty,proj.i$TSSEnrichment)
df.TSS$logTSS=log10(df.TSS$proj.i.TSSEnrichment+1)
df.TSS <- dplyr::filter(df.TSS,logTSS >= 0.5)
saveRDS(df.TSS,paste0("df_TSS_",i,".rds"))

df.depth <- data.frame(proj.i$cellNames,proj.i$depth.cluster,proj.i$depth.cluster.uncertainty,proj.i$nFrags)
df.depth$lognFrags=log10(df.depth$proj.i.nFrags)
df.depth <- dplyr::filter(df.depth,lognFrags>= 3.7)
saveRDS(df.depth,paste0("df_depth_",i,".rds"))

ggPoint(
    x = log10(proj.i$nFrags),
    y = log10(proj.i$TSSEnrichment+1),
    colorDensity = T,
    continuousSet = "sambaNight",
    xlabel = "log10(unique fragments)",
    ylabel = "log10(TSS Enrichment+1)"
  ) +geom_hline(yintercept = log10(min(df.TSS$proj.i.TSSEnrichment)+1),linetype = "dashed")+
    geom_vline(xintercept = min(log10(df.depth$proj.i.nFrags)),linetype = "dashed")+
    ggtitle(paste0("QC thresholds:\n",i))
    ggsave(paste0(i,"_QC.pdf"),width = 4,height = 4)

# Filter out low quality cells, and remove doublets
##############################################################################
list.depth <- list.files(pattern = "^df_depth")

df.depth <-  data.frame(cellNames=character(),
                        cluster=character(),
                        cluster.uncertainty=character(),
                        nFrags = character())
for (i in list.depth){
  df <- readRDS(i)
  df <- df[,1:4]
  colnames(df) <- c("cellNames","cluster","cluster.uncertainty","nFrags")
  df.depth <- rbind(df.depth,df)
}

list.TSS <- list.files(pattern = "^df_TSS")

df.TSS <-  data.frame(cellNames=character(),
                      cluster=character(),
                      cluster.uncertainty=character(),
                      TSSEnrichment = character())
for (i in list.TSS){
  df <- readRDS(i)
  df <- df[,1:4]
  colnames(df) <- c("cellNames","cluster","cluster.uncertainty","TSSEnrichment")
  df.TSS <- rbind(df.TSS,df)
}


colnames(df.TSS) <- c("cellNames","TSS.cluster","TSS.cluster.uncertainty","TSSEnrichment")
colnames(df.depth) <- c("cellNames","depth.cluster","depth.cluster.uncertainty","nFrags")

cellsPass <- intersect(df.TSS$cellNames,df.depth$cellNames)

cellsFail <-  proj$cellNames[!(proj$cellNames %in% cellsPass)]

# Screen for high quality barcodes (remove non cellular barcodes)
proj.filter1 <- proj[proj$cellNames %in% cellsPass]


proj.filter2 <- filterDoublets(proj.filter1,filterRatio = 1.5,cutEnrich = 1,cutScore = -Inf)



###############################################################################################################
plotFragmentSizes(proj.filter2)+ggtitle("Fragment Size Histogram")
ggsave("Frags_hist.pdf",width = 6,height = 4)
plotTSSEnrichment(proj.filter2)+ggtitle("TSS Enrichment")
ggsave("TSS.pdf",width = 6,height = 4)