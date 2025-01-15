library(Seurat)
library(stringr)
library(future)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggh4x)
library(tidyverse)

scRNA <- readRDS("...rds")
celltype.names=names(table(scRNA$celltype.eight))
meta=scRNA@meta.data
data1.list=list()
for (i in 1:length(celltype.names)){
    meta.sub=meta[meta$celltype.eight==celltype.names[i],]
    data=as.data.frame(table(meta.sub$group))
    colnames(data)=c("group","count")
    data$ID=paste0("total.",celltype.names[i])
    data$celltype=celltype.names[i]
    data1.list[[celltype.names[i]]]=data
}
data2.list=list()
subtype.name=names(table(meta$celltype.use))
subtype.name=subtype.name[-6]
subtype.name=subtype.name[-17]
for (i in 1:length(subtype.name)){
    meta.sub=meta[meta$celltype.use==subtype.name[i],]
    data=as.data.frame(table(meta.sub$group))
    colnames(data)=c("group","count")
    data$ID=paste0("sub.",subtype.name[i])
    #data$celltype=subtype.name[i]
    data2.list[[subtype.name[i]]]=data
}

col=c("#ff736d","#008b93")
names(col)=c("AD","Con")
p1=ggplot(data = data, 
       mapping = aes(x = factor(ID), 
                     y = count, fill=factor(group))) +
  geom_bar(stat = 'identity', position = 'fill') + 
  facet_wrap(~celltype,strip.position = "top")+
  labs(fill="color") + theme_bw()+
  theme(panel.grid =element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(size=1, colour = "black"))+
  theme(panel.grid.major=element_line(colour=NA))+ #
  scale_y_continuous(expand = c(0,0))+#
  scale_fill_manual(values = col)+theme(axis.text.x=element_text(angle=25,size=8,face ="bold",vjust = 0.5, hjust = 0.5))+coord_flip()
pdf(paste0(".../","RNA.Cell.ratio.pdf"),width = 10,height = 15)
print(p1)
dev.off()

proj=readRDS("...rds")
meta=proj@cellColData
celltype.names=names(table(proj$predictedGroup_Un))
data1.list=list()
for (i in 1:length(celltype.names)){
    meta.sub=meta[meta$predictedGroup_Un==celltype.names[i],]
    data=as.data.frame(table(meta.sub$class))
    colnames(data)=c("group","count")
    data$ID=paste0("total.",celltype.names[i])
    data$celltype=celltype.names[i]
    data1.list[[celltype.names[i]]]=data
}
data2.list=list()
subtype.name=names(table(meta$predictedGroup_Co))
subtype.name=subtype.name[-6]
subtype.name=subtype.name[-17]
for (i in 1:length(subtype.name)){
    meta.sub=meta[meta$predictedGroup_Co==subtype.name[i],]
    data=as.data.frame(table(meta.sub$class))
    colnames(data)=c("group","count")
    data$ID=paste0("sub.",subtype.name[i])
    data2.list[[subtype.name[i]]]=data
}

col=c("#ff736d","#008b93")
names(col)=c("AD","Con")
p1=ggplot(data = data, 
       mapping = aes(x = factor(ID), 
                     y = count, fill=factor(group))) +
  geom_bar(stat = 'identity', position = 'fill') + 
  facet_wrap(~celltype,strip.position = "top")+
  labs(fill="color") + theme_bw()+
  theme(panel.grid =element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(size=1, colour = "black"))+
  theme(panel.grid.major=element_line(colour=NA))+ #
  scale_y_continuous(expand = c(0,0))+#
  scale_fill_manual(values = col)+theme(axis.text.x=element_text(angle=25,size=8,face ="bold",vjust = 0.5, hjust = 0.5))+coord_flip()
pdf(paste0(".../","ATAC.Cell.ratio.pdf"),width = 10,height = 15)
print(p1)
dev.off()
