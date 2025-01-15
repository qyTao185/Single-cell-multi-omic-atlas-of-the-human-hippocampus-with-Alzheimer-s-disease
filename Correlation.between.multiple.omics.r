library(Seurat)
library(dplyr)
library(ggplot2)

library(RColorBrewer)
library(ComplexHeatmap)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ArchR)
library(parallel)
library(Signac)
library(stringr)
library(future)
library(cowplot)
library(tidyverse)
library(viridis)
theme_set(theme_cowplot())
library(ArchRtoSignac)
library(EnsDb.Hsapiens.v86)

set.seed(11)
addArchRGenome("hg38") # hg38, mm9, mm10
addArchRThreads(threads = 64)
col.eight=c('#7FC97F',"#e066ba","#f3f37a",'#E6AB02',"#084594","#40a2d5","#ee002f","#bb304f")
names(col.eight)=c("Astro","Endo","EX_CA","EX_DG","IN","Micro","Oligo","OPC")
cols = c(brewer.pal(9, "Set1"),brewer.pal(8,"Set2")[1:8],brewer.pal(12,"Paired")[1:12],brewer.pal(8,"Dark2")[1:8],brewer.pal(8,"Accent"))


ATAC.proj <- readRDS("...rds")
scRNA <- readRDS("...rds")
gimat <- readRDS("...gimat.rds")


mycol <- c("#f3f37a",'#E6AB02', "#084594", "#40a2d5", '#7FC97F', "#ee002f", "#bb304f", "#e066ba")
names(mycol) <- my_levels
coluse <- list()
coluse[["typer1"]]=mycol


library(circlize)
breaks <- seq(-2, 2, length.out = length(paletteContinuous("solarExtra")))
color <- circlize::colorRamp2(breaks, paletteContinuous("solarExtra"))

celltype.col <- as.character(sapply(colnames(plot.data), function(x) strsplit(x, ":")[[1]][1]))
celltype.col <- gsub("rna_","",celltype.col)
                                    
celltype.row <- as.character(sapply(rownames(plot.data), function(x) strsplit(x, ":")[[1]][1]))
celltype.row <- gsub("atac_","",celltype.row)
                                    
column_ha = HeatmapAnnotation(df = data.frame(typer1=celltype.col,row.names = colnames(plot.data)),
                            which = "column",
                            gap = unit(1,"mm"),
                            simple_anno_size = unit(0.45, "cm"),
                            col = coluse)
ha = rowAnnotation(typer1 = celltype.row,col = coluse)

                                    pdf("...pdf",width=6,height = 5)
ComplexHeatmap::Heatmap(as.matrix(plot.data),
                   col = color,
                   show_row_names = F,
                   top_annotation = column_ha,
                   left_annotation = ha,
                   use_raster=T,
                   show_row_dend=F,
                   show_column_dend=F,
                   show_column_names = F,
                   column_split=factor(celltype.col,levels=my_levels),
                   row_split=factor(celltype.row,levels=my_levels),
                   cluster_row_slices = FALSE, 
                   cluster_column_slices = FALSE,
                   row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), 
                   border_gp = gpar(col = "white", lty = 1))
dev.off()