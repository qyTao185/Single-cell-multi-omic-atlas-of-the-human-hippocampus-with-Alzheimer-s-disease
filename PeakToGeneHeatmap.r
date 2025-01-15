library(ArchR)
library(parallel)
library(Seurat)
# library(Signac)
library(stringr)
library(future)
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)

subClusterGroups <- c(
 ...
)
# P2G definition cutoffs
corrCutoff <- 0.45       # Default in plotPeak2GeneHeatmap is 0.45
varCutoffATAC <- 0.25   # Default in plotPeak2GeneHeatmap is 0.25
varCutoffRNA <- 0.25    # Default in plotPeak2GeneHeatmap is 0.25

# Coaccessibility cutoffs
coAccCorrCutoff <- 0.5  # Default in getCoAccessibility is 0.5

full.proj <- readRDS("...rds")
wd <- "..."
setwd(wd)
# Get all peaks
allPeaksGR <- getPeakSet(full.proj)
allPeaksGR$peakName <- (allPeaksGR %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
names(allPeaksGR) <- allPeaksGR$peakName

##########################################################################################
# Prepare full-project peak to gene linkages, loops, and coaccessibility (full and subproject links)
##########################################################################################
plot_loop_list <- list()
plot_loop_list[["scalp"]] <- getPeak2GeneLinks(full.proj, corCutOff=corrCutoff, resolution = 100)[[1]]
# coaccessibility_list <- list()
# coAccPeaks <- getCoAccessibility(full.proj, corCutOff=corrCutoff, returnLoops=TRUE)[[1]]
# coAccPeaks$linkName <- (coAccPeaks %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
# coAccPeaks$source <- "scalp"
# coaccessibility_list[["scalp"]] <- coAccPeaks
peak2gene_list <- list()
p2gGR <- getP2G_GR(full.proj, corrCutoff=NULL, varCutoffATAC=-Inf, varCutoffRNA=-Inf, filtNA=FALSE)
p2gGR$source <- "scalp"
peak2gene_list[["scalp"]] <- p2gGR

for(subgroup in subClusterGroups){
  message(sprintf("Reading in subcluster %s", subgroup))
    
  sub_dir <- paste0("...", subgroup,"/Output/Save-ArchR-Project.rds")
  sub_proj <- readRDS(sub_dir)

  subP2G <- getP2G_GR(sub_proj, corrCutoff=NULL, varCutoffATAC=-Inf, varCutoffRNA=-Inf, filtNA=FALSE)
  subP2G$source <- subgroup
  peak2gene_list[[subgroup]] <- subP2G

  plot_loop_list[[subgroup]] <- getPeak2GeneLinks(sub_proj, corCutOff=corrCutoff, resolution = 100)[[1]]

  # Get coaccessibility
  # coAccPeaks <- getCoAccessibility(sub_proj, corCutOff=coAccCorrCutoff, returnLoops=TRUE)[[1]]
  # coAccPeaks$linkName <- (coAccPeaks %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
  # coAccPeaks$source <- subgroup
  # coaccessibility_list[[subgroup]] <- coAccPeaks
}

full_p2gGR <- as(peak2gene_list, "GRangesList") %>% unlist()
# full_coaccessibility <- as(coaccessibility_list, "GRangesList") %>% unlist()

idxATAC <- peak2gene_list[["scalp"]]$idxATAC
names(idxATAC) <- peak2gene_list[["scalp"]]$peakName
full_p2gGR$idxATAC <- idxATAC[full_p2gGR$peakName]

saveRDS(full_p2gGR, file=paste0(wd, "/multilevel_p2gGR.rds")) # NOT merged or correlation filtered
# saveRDS(full_coaccessibility, file=paste0(wd, "/multilevel_coaccessibility.rds"))
saveRDS(plot_loop_list, file=paste0(wd, "/multilevel_plot_loops.rds"))

##########################################################################################
# Upset plot of number of peak to gene links identified per group
##########################################################################################

library(UpSetR)

groups <- unique(full_p2gGR$source)
sub_p2gGR <- full_p2gGR[!is.na(full_p2gGR$Correlation)]

upset_list <- lapply(groups, function(g){
  gr <- sub_p2gGR[sub_p2gGR$source == g &  sub_p2gGR$Correlation > corrCutoff &  sub_p2gGR$VarQATAC > varCutoffATAC &  sub_p2gGR$VarQRNA > varCutoffRNA]
  paste0(gr$peakName, "-", gr$symbol)
  })
names(upset_list) <- groups


plotUpset <- function(plist, main.bar.color="red", keep.order=FALSE){
  # Function to plot upset plots
  upset(
    fromList(plist), 
    sets=names(plist), 
    order.by="freq",
    empty.intersections="off",
    point.size=3, line.size=1.5, matrix.dot.alpha=0.5,
    main.bar.color=main.bar.color,
    keep.order=keep.order
    )
}

pdf(paste0(wd, "/upset_plot_linked_peaks.pdf"), width=10, height=7)
plotUpset(upset_list, main.bar.color="royalblue1", keep.order=FALSE)
dev.off()


library(RColorBrewer)
celltype.names=names(table(full.proj$group1))
cols=c('#7FC97F',"#e066ba","#f3f37a",'#E6AB02',"#084594","#40a2d5","#ee002f","#bb304f",'#7FC97F',"#e066ba","#f3f37a",'#E6AB02',"#084594","#40a2d5","#ee002f","#bb304f")
names(cols)=celltype.names

##########################################################################################
# Plot Peak2Gene heatmap
##########################################################################################

################################
nclust <- 25 
p <- plotPeak2GeneHeatmap(
  full.proj, 
  corCutOff = corrCutoff, 
  groupBy="group1", 
  nPlot = 1000000, returnMatrices=FALSE, 
  k=nclust, seed=1, palGroup=cols
  )
pdf(paste0(wd, sprintf("/peakToGeneHeatmap_LabelClust_k%s.pdf", nclust)), width=16, height=10)
p
dev.off()
################################