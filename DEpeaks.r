library(ArchR)
library(parallel)
library(Seurat)
# library(Signac)
library(stringr)
library(future)
library(ArchR)
set.seed(11)
addArchRGenome("hg38") # hg38, mm9, mm10
addArchRThreads(threads = 64)
library(ArchR)
library(ChIPpeakAnno)
library(stats)
ArchR::addArchRThreads(threads = 64)

proj <- readRDS("...rds")

celltype.names=c('Astro','Endo','EX_CA','EX_DG','IN','Micro','Oligo','OPC')
addArchRThreads(threads = 64)
DEP2G.list=list()
DEpeak.list=list()
for (i in 1:length(celltype.names)){
    proj$test.group="aa"
    proj$test.group[which(proj$predictedGroup_Un==celltype.names[i] & proj$class=="AD")]="AD"
    proj$test.group[which(proj$predictedGroup_Un==celltype.names[i] & proj$class=="Con")]="Con"
  markersPeaks <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix",
  groupBy = "test.group",
  testMethod = "binomial",binarize=TRUE,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "AD",
  bgdGroups = "Con"
)
    DEpeak.list[[celltype.names[i]]] = markersPeaks
}