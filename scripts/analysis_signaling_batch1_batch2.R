##########################################################################
##########################################################################
# Project: Natasya's CSD project
# Script purpose: analyze the signaling pathway for BL and CSD scRNA-seq data
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Jan 13 11:01:15 2023
##########################################################################
##########################################################################
rm(list = ls())

library(Seurat)
library(decoupleR)
library(tictoc)
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library("viridis")

version.analysis = '_axolotl_20241216'
resDir = paste0("../results/scRNAseq_signaling.analysis", version.analysis)
RdataDir = paste0(resDir, '/Rdata')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../fromTobie/CSD_batch1_batch2/'

functionDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts'
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_scRNAseq.R')
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_Visium.R')

annot = readRDS(paste0('/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                       'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))

convert_to_geneSymbol = function(gene.ids, annot)
{
  # gene.ids = rownames(mat)
  mm = match(gene.ids, annot$geneID)
  ggs = annot$gene.symbol.toUse[mm]
  
  kk = which(is.na(ggs))
  ggs[kk] = gene.ids[kk]
  #return(make.unique(ggs))
  return(ggs)
  
}

levels = c('CSD_0dpa', 
           'BL_3and5dpa', 'BL_5dpa', 'BL_6dpa', 'BL_7dpa',  'BL_8dpa', 'BL_11dpa', 
           'CSD_3dpa', 'CSD_5dpa', 'CSD_6dpa', 'CSD_7dpa', 'CSD_8dpa', 'CSD_11dpa')

cols = rep(NA, length = length(levels))
names(cols) = levels

cols[1] = 'gray60'
#cols[1:3] = viridis(3)
cols[2:7] = colorRampPalette((brewer.pal(n = 6, name ="Blues")))(6)
cols[8:13] = colorRampPalette((brewer.pal(n = 6, name ="OrRd")))(6)

## gene example of signaling pathways 
sps = readRDS(file = paste0("/groups/tanaka/People/current/jiwang/projects/RA_competence/",
                            '/data/annotations/curated_signaling.pathways_gene.list_v3.rds'))
sps = unique(sps$gene)
sps = toupper(sps)

## TFs 
tfs = readRDS(file = paste0('/groups/tanaka/People/current/jiwang/projects/RA_competence/data',
                            '/annotations/curated_human_TFs_Lambert.rds'))
tfs = unique(tfs$`HGNC symbol`)


########################################################
########################################################
# Section 0: prepare the scRNA-seq samples  
# import the processed scRNA data from Natasia and Tobi 
########################################################
########################################################
aa = readRDS(file = paste0(dataDir, 'BL_SeuratObj.RDS'))
#xx = readRDS(file = paste0("../fromTobie/CSD_batch1/", 'BL_SeuratObj.RDS'))

bb = readRDS(file = paste0(dataDir, 'CSD_SeuratObj.RDS'))

aa$batch = aa$exp
bb$batch = bb$exp

p1 = DimPlot(aa, group.by = 'celltype', label = TRUE, repel = TRUE)
p2 = DimPlot(bb, group.by = 'celltype', label = TRUE, repel = TRUE)

p1 + p2

ggsave(filename = paste0(resDir, '/Tobie_umap.harmony_BL.CSD.split.pdf'), width = 18, height = 6)

xx = merge(aa, bb)

DimPlot(xx, group.by = 'type', label = TRUE, repel = TRUE)

aa = xx
rm(bb)
rm(xx)

aa$condition = aa$orig.ident
aa$condition[which(aa$condition == 'BL_3_5dpa')] = 'BL_3and5dpa'
aa$condition[which(aa$condition == 'CSD_5dpa_A')] = 'CSD_5dpa'
aa$condition[which(aa$condition == 'CSD_5dpa_B')] = 'CSD_5dpa'

aa$sample = aa$condition
aa$sample[grep('BL_', aa$sample)] = 'BL'
aa$sample[grep('CSD_', aa$sample)] = 'CSD'

aa$days = aa$condition
aa$days = gsub('BL_', '', aa$days)
aa$days = gsub('CSD_', '', aa$days)


saveRDS(aa, file = paste0(RdataDir, '/BL.CSD_twoBacthes_merged.rds'))

########################################################
########################################################
# Section I: batch correction for two batches and BL and CSD
# 
########################################################
########################################################
aa = readRDS(file = paste0(RdataDir, '/BL.CSD_twoBacthes_merged.rds'))

aa = NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 8000) # find subset-specific HVGs

aa <- ScaleData(aa)

## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, dims = 1:50, n.neighbors = 30, min.dist = 0.3)

p1 = DimPlot(aa, group.by = 'batch', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'sample', label = TRUE, repel = TRUE)
p3 = DimPlot(aa, group.by = 'condition', label = TRUE, repel = TRUE)

(p1 + p2)/p3 

ggsave(filename = paste0(resDir, '/UMAP_merged.BL.CSD_batch_samples_condition_no.batchCorrection.pdf'), 
       width = 18, height = 12)

DimPlot(aa, group.by = 'celltype', split.by = 'sample', label = TRUE, repel = TRUE)

p1 = DimPlot(aa, group.by = 'batch', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'sample', label = TRUE, repel = TRUE)
p3 = DimPlot(aa, group.by = 'celltype', label = TRUE, repel = TRUE)

(p1 + p2)/p3 

ggsave(filename = paste0(resDir, '/UMAP_merged.BL.CSD_batch_samples_celltypes_no.batchCorrection.pdf'), 
       width = 18, height = 12)

saveRDS(aa, file = paste0(RdataDir, '/BL.CSD_merged_renormalized_8000HVGs.rds'))

##########################################
# make batch correction using harmony 
##########################################
aa = readRDS(file = paste0(RdataDir, '/BL.CSD_merged_renormalized_8000HVGs.rds'))

source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_dataIntegration.R')

aa$dataset = paste0(aa$sample, '_', aa$batch)

xx = IntegrateData_runHarmony(aa, group.by = 'dataset', 
                              nfeatures = 5000, 
                              dims.use = c(1:50),
                              nclust = NULL,
                              redo.normalization.hvg.scale.pca = FALSE)

p1 = DimPlot(xx, group.by = 'batch', label = TRUE, repel = TRUE)
p2 = DimPlot(xx, group.by = 'sample', label = TRUE, repel = TRUE)
p3 = DimPlot(xx, group.by = 'celltype', label = TRUE, repel = TRUE)

(p1 + p2)/p3 

ggsave(filename = paste0(resDir,
                         '/UMAP_merged.BL.CSD_batch_samples_celltypes_Harmony_twoBatch.BL.CSD.pdf'), 
       width = 18, height = 12)

p1 = DimPlot(xx, group.by = 'batch', label = TRUE, repel = TRUE)
p2 = DimPlot(xx, group.by = 'sample', label = TRUE, repel = TRUE)
p3 = DimPlot(xx, group.by = 'condition', label = TRUE, repel = TRUE)

(p1 + p2)/p3 

ggsave(filename = paste0(resDir, 
                         '/UMAP_merged.BL.CSD_batch_samples_condition_Harmony_twoBatch.BL.CSD.pdf'), 
       width = 18, height = 12)


saveRDS(xx, file = paste0(RdataDir, '/BL.CSD_merged_renormalized_8000HVGs_harmony.twoBatches.rds'))


########################################################
########################################################
# Section II: subset relevant cell populations:
# connective tissue, epidermis, macrophage, neutrophils
########################################################
########################################################

aa = readRDS(file = paste0("../results/scRNAseq_signaling.analysis_axolotl_20240116/Rdata", 
                           '/BL.CSD_merged_renormalized_8000HVGs_harmony.twoBatches.rds'))

Idents(aa) = aa$celltype

p1 = DimPlot(aa, group.by = 'batch', label = TRUE, repel = TRUE) 
p2 = DimPlot(aa, group.by = 'sample', label = TRUE, repel = TRUE)
p3 = DimPlot(aa, group.by = 'celltype', label = TRUE, repel = TRUE) + NoLegend()

(p1 + p2)/p3 

# redo clustering 
aa <- FindNeighbors(aa, reduction = "harmony",
                    dims = 1:30,  k.param = 20)
aa <- FindClusters(aa, resolution = 1.0)

p1 = DimPlot(aa, label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'celltype', label = TRUE, repel = TRUE) + NoLegend()
p1 + p2

Idents(aa) = aa$celltype
xx = subset(aa, idents = c('Connective Tissue', 'Epidermis', 'Macrophages', 'Neutrophils'))
DimPlot(xx, group.by = 'celltype',  label = TRUE, repel = TRUE)

aa = xx
rm(xx)

aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 8000) # find subset-specific HVGs

#all.genes <- rownames(aa)
aa <- ScaleData(aa)

## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.3)

DimPlot(aa, group.by = 'celltype',  label = TRUE, repel = TRUE)

aa <- FindNeighbors(aa, dims = 1:30,  k.param = 30)
aa <- FindClusters(aa, resolution = 0.5)

DimPlot(aa, label = TRUE, repel = TRUE)

aa$clusters = aa$seurat_clusters


p1 = DimPlot(aa, group.by = 'batch', label = TRUE, repel = TRUE) 
p2 = DimPlot(aa, group.by = 'sample', label = TRUE, repel = TRUE)
p3 = DimPlot(aa, group.by = 'celltype',  label = TRUE, repel = TRUE)
p4 = DimPlot(aa, group.by = 'clusters',  label = TRUE, repel = TRUE)

(p1 + p2)/(p3 + p4)

ggsave(filename = paste0(resDir, 
                         '/UMAP_merged.BL.CSD_subset_celltypes_newClusters.pdf'), 
       width = 18, height = 12)


saveRDS(aa, file = paste0(RdataDir, '/BL.CSD_merged_subset_CT_MAC_Neu_Epd_umap.rds'))

##########################################
# double check Tobi's CT filtering  
##########################################
Test_Tobi_CT_filtering_steps = FALSE
if(Test_Tobi_CT_filtering_steps){
  
  bb = readRDS(file = paste0(dataDir, 'CSD_SeuratObj.RDS'))
  
  DimPlot(bb, group.by = 'seurat_clusters')
  CSD_int = subset(bb, cells = which(bb$celltype == 'Connective Tissue'))
  
  DimPlot(CSD_int, group.by = 'seurat_clusters')
  
  CSD_CT = subset(bb, idents = c(5,12,13,17,19))
  
  #Cluster cells
  CSD_CT<- FindNeighbors(CSD_CT, dims = 1:100, reduction = "harmony")
  CSD_CT <- FindClusters(CSD_CT, resolution = 1)
  
  CSD_CT <- RunUMAP(CSD_CT, reduction = "harmony", dims = 1:100, min.dist = 0.1)
  
  pdf(paste0(resDir, "/CSD_CT_Harmony_Embedding_PC100_res1.pdf"), 
      width=10,height=10)
  DimPlot(object = CSD_CT, reduction = 'umap',raster = T, pt.size = 4, shuffle= T,label = T)
  DimPlot(object = CSD_CT, reduction = 'umap', raster = T,shuffle= T, pt.size = 4,group.by = "orig.ident")
  DimPlot(object = CSD_CT, reduction = 'umap',raster = T, shuffle= T, pt.size = 4,group.by = "exp")
  FeaturePlot(CSD_CT, raster = T,pt.size = 2, features = c("nFeature_RNA","nCount_RNA","percent.mt","mCherry",
                                                           "eGFP"), order = T, 
              cols = c("grey",rev(viridis_pal(option = "viridis")(12))))
  dev.off()
  
  #Filter low count cluster
  VlnPlot(CSD_CT, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'))

  CSD_CT = subset(CSD_CT, idents = c(1,8,10), invert = T)
  
  #Cluster cells
  CSD_CT<- FindNeighbors(CSD_CT, dims = 1:100, reduction = "harmony")
  CSD_CT <- FindClusters(CSD_CT, resolution = 1)
  
  CSD_CT <- RunUMAP(CSD_CT, reduction = "harmony", dims = 1:100, min.dist = 0.1)
  
  
  pdf(paste0(resDir, "/CSD_CT_Harmony_Embedding_PC100_res1_filtering_1.pdf"), 
      width=10,height=10)
  
  DimPlot(object = CSD_CT, reduction = 'umap',raster = T, pt.size = 4, shuffle= T,label = T)
  DimPlot(object = CSD_CT, reduction = 'umap', raster = T,shuffle= T, pt.size = 4, group.by = "orig.ident")
  DimPlot(object = CSD_CT, reduction = 'umap',raster = T, shuffle= T, pt.size = 4, group.by = "exp")
  FeaturePlot(CSD_CT, raster = T,pt.size = 2, features = c("nFeature_RNA","nCount_RNA","percent.mt","mCherry","eGFP"),
              order = T, cols = c("grey",rev(viridis_pal(option = "viridis")(12))))
  
  dev.off()
  
  #Filter contaminants
  #cluster 12 = muscle sarcolemma
  #cluster 14 = macrophage
  #cluster 11 = high percent.mt
  #cluster 15 = smooth muscle
  
  CSD_CT = subset(CSD_CT, idents = c(11,12,14,15), invert = T)
  
  #Cluster cells
  CSD_CT<- FindNeighbors(CSD_CT, dims = 1:100, reduction = "harmony")
  CSD_CT <- FindClusters(CSD_CT, resolution = 1)
  
  CSD_CT <- RunUMAP(CSD_CT, reduction = "harmony", dims = 1:100, min.dist = 0.1)
  
  
  pdf("CSD_CT_Harmony_Embedding_PC100_res1.pdf",width=10,height=10)
  DimPlot(object = CSD_CT, reduction = 'umap',raster = T, pt.size = 4, shuffle= T,label = T)
  DimPlot(object = CSD_CT, reduction = 'umap', raster = T,shuffle= T, pt.size = 4,group.by = "orig.ident")
  DimPlot(object = CSD_CT, reduction = 'umap',raster = T, shuffle= T, pt.size = 4,group.by = "exp")
  FeaturePlot(CSD_CT, raster = T,pt.size = 2, features = c("nFeature_RNA","nCount_RNA","percent.mt","mCherry","eGFP"),
              order = T, cols = c("grey",rev(viridis_pal(option = "viridis")(12))))
  dev.off()
  
  #Filter low count cluster
  #cluster 1 = low count and blood signal
  CSD_CT = subset(CSD_CT, idents = c(1,5), invert = T)
  
  #Cluster cells
  CSD_CT<- FindNeighbors(CSD_CT, dims = 1:100, reduction = "harmony")
  CSD_CT <- FindClusters(CSD_CT, resolution = 1)
  
  CSD_CT <- RunUMAP(CSD_CT, reduction = "harmony", dims = 1:100, min.dist = 0.1)
  
  CSD_CT@meta.data$time = CSD_CT@meta.data$orig.ident
  CSD_CT@meta.data$time = gsub("dpa","",CSD_CT@meta.data$time)
  CSD_CT@meta.data$time = gsub("_A","",CSD_CT@meta.data$time)
  CSD_CT@meta.data$time = gsub("_B","",CSD_CT@meta.data$time)
  CSD_CT@meta.data$time = gsub("CSD_","",CSD_CT@meta.data$time)
  CSD_CT@meta.data$time = as.numeric(CSD_CT@meta.data$time)
  
  
  CSD_CT@meta.data$injury = "CSD"
  
  #my_spectral_palette = colorRampPalette(colors = grafify:::graf_palettes$okabe_ito[1:8])
  
  pdf("CSD_CT_Harmony_Embedding_PC100_res1.pdf",width=10,height=10)
  DimPlot(object = CSD_CT, reduction = 'umap',raster = T, pt.size = 4, shuffle= T,label = T, 
          cols = my_spectral_palette(length(unique(CSD_CT@active.ident))))
  DimPlot(object = CSD_CT, reduction = 'umap', raster = T,shuffle= T, pt.size = 4,group.by = "orig.ident")
  DimPlot(object = CSD_CT, reduction = 'umap',raster = T, shuffle= T, pt.size = 4,group.by = "exp")
  DimPlot(object = CSD_CT, reduction = 'umap',raster = T, shuffle= T, pt.size = 3,group.by = "time",
          cols = c("grey50",rev(viridis_pal(option = "plasma")(length(unique(CSD_CT@meta.data$time))-1))))
  FeaturePlot(CSD_CT, raster = T,pt.size = 2, features = c("nFeature_RNA","nCount_RNA","percent.mt","mCherry","eGFP"),
              order = T, cols = c("grey",rev(viridis_pal(option = "viridis")(12))))
  dev.off()
  
  gene_ids = c("AMEX60DD048840","AMEX60DD023361","AMEX60DD029426","AMEX60DD055540","AMEX60DD029436", 
               "AMEX60DD037674","AMEX60DD023964","AMEX60DD009987","AMEX60DD030520","AMEX60DD048972",
               "AMEX60DD002658","AMEX60DD012132","AMEX60DD018450","AMEX60DD020580","AMEX60DD052517",
               "AMEX60DD053922","AMEX60DD048332","AMEX60DD052070","AMEX60DD031414","AMEX60DD024964",
               "AMEX60DD045921","AMEX60DD035908","AMEX60DD006619","AMEX60DD013910","AMEX60DD025537")
  
  gene_names =c("KLF5","COL7A1","COL2A1","COL6A1","TWIST3","TNMD","ASPN","MEOX1","MGP",
                "CNMD","OTOS","GREM1","PRRX1","MYH11","ACTA2","TAGLN","COL8A1","C1QB","LGALS3BP",
                "LECT2","S100P","EPCAM","VWF","PLVAP","ALAS2")
  
  gg_Fig <- FeaturePlot(CSD_CT, pt.size = 2,raster = T, features = gene_ids,order = T, slot = "data", repel=T)
  gg_Fig <- lapply( 1:length(gene_ids), function(x) { gg_Fig[[x]] + labs(title=gene_names[x])  & NoAxes()})
  #,cols = c("grey",rev(viridis_pal(option = "mako")(12)))
  
  
  pdf("CSD_CT_Harmony_UMAP_feature_CellAtlas_Marker_raw.pdf",width=16,height=10)
  CombinePlots( gg_Fig )
  dev.off()
  
}


##########################################
# check the condition and time changes for each cell type 
##########################################
aa = readRDS(file = paste0(RdataDir, '/BL.CSD_merged_subset_CT_MAC_Neu_Epd_umap.rds'))

aa$celltype[which(aa$celltype == 'Connective Tissue')] = 'CT'
Idents(aa) = aa$celltype

celltype = 'CT'
sub.obj = subset(x = aa, idents = celltype)

sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = 5000)
sub.obj = ScaleData(sub.obj)
sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE)
ElbowPlot(sub.obj, ndims = 30)

nb.pcs = 30 # nb of pcs depends on the considered clusters or ids
n.neighbors = 30; min.dist = 0.3;
sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", 
                   dims = 1:nb.pcs, 
                   n.neighbors = n.neighbors,
                   min.dist = min.dist)

p1 = DimPlot(sub.obj, group.by = 'batch', label = TRUE, repel = TRUE) 
p2 = DimPlot(sub.obj, group.by = 'sample', label = TRUE, repel = TRUE)
p3 = DimPlot(sub.obj, group.by = 'condition',  label = TRUE, repel = TRUE)
p4 = DimPlot(sub.obj, group.by = 'Phase',  label = TRUE, repel = TRUE)

(p1 + p2)/(p3 + p4)

ggsave(filename = paste0(resDir, 
                         '/UMAP_merged.BL.CSD_subset_celltypes_checkEachcelltypes.aross.condition.time_',
                         'noBatchIntegration_', celltype,  '.pdf'), 
       width = 18, height = 12)


source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_dataIntegration.R')

xx = IntegrateData_runHarmony(sub.obj, group.by = 'dataset', 
                              nfeatures = 5000, 
                              dims.use = c(1:30), 
                              nclust = NULL,
                              redo.normalization.hvg.scale.pca = TRUE)

xx <- RunUMAP(xx, reduction = "harmony", dims = 1:50, n.neighbors = 50, min.dist = 0.1)

p1 = DimPlot(xx, group.by = 'batch', label = TRUE, repel = TRUE) 
p2 = DimPlot(xx, group.by = 'sample', label = TRUE, repel = TRUE)
p3 = DimPlot(xx, group.by = 'condition',  label = TRUE, repel = TRUE)
p4 = DimPlot(xx, group.by = 'Phase',  label = TRUE, repel = TRUE)

(p1 + p2)/(p3 + p4)
ggsave(filename = paste0(resDir, 
                         '/UMAP_merged.BL.CSD_subset_celltypes_checkEachcelltypes.aross.condition.time_',
                         'BatchIntegration_Harmony_', celltype,  '.pdf'), 
       width = 18, height = 12)


xx <- FindNeighbors(xx, dims = 1:50, reduction = "harmony")
xx <- FindClusters(xx, resolution = 1)

pdf(paste0(resDir, "/",  celltype,  "_Harmony_Embedding.pdf"), 
    width=10, height=8)
DimPlot(object = xx, reduction = 'umap',raster = T, pt.size = 4, shuffle= T,label = T)

DimPlot(object = xx, reduction = 'umap', raster = T,shuffle= T, pt.size = 4,group.by = "orig.ident")
DimPlot(object = xx, reduction = 'umap',raster = T, shuffle= T, pt.size = 4,group.by = "exp")
DimPlot(object = xx, reduction = 'umap',raster = T, shuffle= T, pt.size = 3,group.by = "time",
        cols = c("grey50",rev(viridis_pal(option = "plasma")(length(unique(xx@meta.data$time))-1))))
FeaturePlot(xx, raster = T,pt.size = 2, features = c("nFeature_RNA","nCount_RNA","percent.mt","mCherry","eGFP"),
            order = T)

VlnPlot(xx, features = 'nFeature_RNA')
VlnPlot(xx, features = 'nCount_RNA')
VlnPlot(xx, features = 'percent.mt')
VlnPlot(xx, features = 'mCherry')

dev.off()

xx = subset(xx, idents = c(2, 6, 7, 11, 12, 13, 20, 21, 14, 24), invert = T)

xx <- RunUMAP(xx, reduction = "harmony", dims = 1:50,  min.dist = 0.1)

p1 = DimPlot(xx, group.by = 'batch', label = TRUE, repel = TRUE) 
p2 = DimPlot(xx, group.by = 'sample', label = TRUE, repel = TRUE)
p3 = DimPlot(xx, group.by = 'condition',  label = TRUE, repel = TRUE)
p4 = DimPlot(xx, group.by = 'Phase',  label = TRUE, repel = TRUE)

(p1 + p2)/(p3 + p4)

xx <- FindNeighbors(xx, dims = 1:50, reduction = "harmony")
xx <- FindClusters(xx, resolution = 1)


pdf(paste0(resDir, "/",  celltype,  "_Harmony_Embedding_filtered.pdf"), 
    width=10, height=8)
DimPlot(object = xx, reduction = 'umap',raster = T, pt.size = 4, shuffle= T,label = T)

DimPlot(object = xx, reduction = 'umap', raster = T,shuffle= T, pt.size = 4,group.by = "orig.ident")
DimPlot(object = xx, reduction = 'umap',raster = T, shuffle= T, pt.size = 4,group.by = "exp")
DimPlot(object = xx, reduction = 'umap',raster = T, shuffle= T, pt.size = 3,group.by = "time",
        cols = c("grey50",rev(viridis_pal(option = "plasma")(length(unique(xx@meta.data$time))-1))))
FeaturePlot(xx, raster = T,pt.size = 2, features = c("nFeature_RNA","nCount_RNA","percent.mt","mCherry","eGFP"),
            order = T)

VlnPlot(xx, features = 'nFeature_RNA')
VlnPlot(xx, features = 'nCount_RNA')
VlnPlot(xx, features = 'percent.mt')
VlnPlot(xx, features = 'mCherry')

dev.off()

yy = readRDS(file = paste0(dataDir, 'BLct_CSDct_Harmony_SeuratObj.RDS'))
DimPlot(yy, group.by = 'type')

ggs = convert_to_geneSymbol(rownames(xx), annot = annot)
which(ggs == 'PRRX1')
FeaturePlot(xx, features = rownames(xx)[ which(ggs == 'HMGB3')])



##########################################
# epidermis  
##########################################
aa = readRDS(paste0(dataDir, 'BLep_CSDep_Harmony_SeuratObj.RDS'))


p1 = DimPlot(aa, group.by = 'injury')
p2 = DimPlot(object = aa, reduction = 'umap',raster = T, shuffle= T, pt.size = 3,group.by = "time",
        cols = c("grey50",rev(viridis_pal(option = "plasma")(length(unique(xx@meta.data$time))-1))))

aa <- FindNeighbors(aa, dims = 1:50, reduction = "harmony")
aa <- FindClusters(aa, resolution = 0.5)

p3 = DimPlot(aa, group.by = 'seurat_clusters', label = TRUE, repel = TRUE)

(p1 + p2)/p3

ggs = c('AMEX60DD053027', 'AMEX60DD050718', 'AMEX60DD044963',
          'AMEX60DD014606', 'AMEX60DD018809', 'AMEX60DD018809',
          'AMEX60DD055164', 'AMEX60DD055164', 'AMEX60DD022398',
          'AMEX60DD009937')
FeaturePlot(aa, features = ggs)

genes = convert_to_geneSymbol(ggs, annot = annot)




Idents(aa) = aa$seurat_clusters

markers = FindAllMarkers(aa, only.pos = TRUE,  
                         min.pct = 0.1, 
                         logfc.threshold = 0.25)

markers = readRDS(file = paste0(RdataDir, '/Epidermis_markers.rds'))
markers$geneSymbols = convert_to_geneSymbol(markers$gene, annot = annot)
saveRDS(markers, file = paste0(RdataDir, '/Epidermis_markers.rds'))

markers = readRDS(file = paste0(RdataDir, '/Epidermis_markers.rds'))
markers = markers[!is.na(match(markers$geneSymbols, c(tfs, sps))), ]


markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

DoHeatmap(aa, features = top10$gene) + NoLegend()

ggsave(filename = paste0(resDir,  'Epidermis_heatmap_markerGenes_top10.pdf'), 
       width = 12, height = 30)

markers = markers[which(markers$cluster == 4), ]

write.csv2(markers, file = paste0(resDir, '/markerGenes_Blastema_specific_Epidermis.csv'), row.names = TRUE)


features = intersect(markers$gene[which(markers$p_val < 10^-100)], sps)
cat(length(features), ' genes to display \n')

DoHeatmap(aa, features = features) + NoLegend()
ggsave(filename = paste0(saveDir,  'heatmap_markerGenes_signalingGenes_by_', cell_ids, '.pdf'), 
       width = 12, height = 40)

tfs_sel = intersect(tfs, markers$gene[which(markers$p_val < 10^-50)])

DoHeatmap(aa, features = tfs_sel) + NoLegend()
ggsave(filename = paste0(saveDir,  'heatmap_markerGenes_TFs_by_', cell_ids, '.pdf'), 
       width = 12, height = 30)


########################################################
########################################################
# Section II : ligand-receptor anlaysis
# e.g. LIANA and NicheNet with defined subpopulation
# test LIANA and NicheNet using clusters 
########################################################
########################################################
## LR interaction using only the cells from batch1
dataDir = '../results/scRNAseq_signaling.analysis_axolotl_20230308/Rdata/'
aa = readRDS(file = paste0(dataDir, '/BL.CSD_merged_subset_CT_MAC_Neu_Epd_day3_5_8_subtypes_umap.rds'))

##########################################
# test LIANA for all pairs
##########################################
#dataDir = '../results/scRNAseq_signaling.analysis_axolotl_20240116/Rdata/'
#aa = readRDS(file = paste0(dataDir, '/BL.CSD_merged_subset_CT_MAC_Neu_Epd_umap.rds'))

outDir = paste0(resDir, '/LR_analysis_LIANA_mergingCTsubtypes')
additionalLabel = '_fixedCelltypes'

aa$condition = droplevels(aa$condition)

p1 = DimPlot(aa, group.by = 'subtypes', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'condition', label = TRUE, repel = TRUE)
p1 + p2

table(aa$subtypes, aa$condition)


## manually merge subclusters
#aa$subtypes[which(aa$subtypes == 'epidermis_BL_early_2')] = 'epidermis_BL_early'
#aa$subtypes[which(aa$subtypes == 'epidermis_BL_early_1')] = 'epidermis_BL_early'
#aa$subtypes[which(aa$subtypes == 'CT_BL_early_2')] = 'CT_BL_early_1'

DimPlot(aa, group.by = 'subtypes', label = TRUE, repel = TRUE)

ggsave(filename = paste0(resDir,  '/UMAP_subtypes.pdf'), 
       width = 12, height = 8)

p1 = DimPlot(aa, group.by = 'condition', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'celltype', label = TRUE, repel = TRUE)

p1 + p2

ggsave(filename = paste0(resDir,  '/UMAP_condition_celltype.pdf'),
       width = 16, height = 6)

Subset_epidermis = FALSE
if(Subset_epidermis){
  celltype2subset = 'epidermis'
  
  cells = colnames(aa)[which(aa$celltype == celltype2subset)]
  sub.obj = subset(aa, cells = cells) 
  
  sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = 3000)
  sub.obj = ScaleData(sub.obj)
  sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE)
  ElbowPlot(sub.obj, ndims = 30)
  
  nb.pcs = 30 # nb of pcs depends on the considered clusters or ids
  n.neighbors = 30; min.dist = 0.3;
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", 
                     dims = 1:nb.pcs, 
                     n.neighbors = n.neighbors,
                     min.dist = min.dist)
  
  sub.obj$condition = droplevels(sub.obj$condition)
  sub.obj$subtypes[which(sub.obj$condition == 'CSD_3dpa')] = 'epidermis_BL.CSD_early'
  
  
  p1 = DimPlot(sub.obj, group.by = 'condition', label = TRUE, repel = TRUE)
  p2 = DimPlot(sub.obj, group.by = 'subtypes', label = TRUE, repel = TRUE)
  
  p1 + p2
  
  
  ggsave(filename = paste0(outDir, celltype2subset, '_condtion_subtype.pdf'), 
         width = 16, height = 8)
  
  sub.obj <- FindNeighbors(sub.obj, dims = 1:20,  k.param = 30)
  sub.obj <- FindClusters(sub.obj, resolution = 0.3)
  
  p1 = DimPlot(sub.obj, label = TRUE, repel = TRUE)
  p2 = DimPlot(sub.obj, group.by = 'groups', label = TRUE, repel = TRUE)
  
  p1 | p2
  
  saveRDS(sub.obj, file = paste0(RdataDir, '/Early_subset_epidermis.rds'))
  
  
  ggsave(filename = paste0(outDir, celltype2subset, '_groups_subclustering.pdf'), 
         width = 14, height = 6)
  
  ## assign new labels
  sub.obj$subtypes = NA
  sub.obj$subtypes[which(sub.obj$seurat_clusters == '2'| sub.obj$seurat_clusters == '3'|
                           sub.obj$seurat_clusters == '4')] = 'epidermis_BL.CSD_early'
  
  sub.obj$subtypes[which(sub.obj$seurat_clusters == '0')] = 'epidermis_BL_early_1'
  sub.obj$subtypes[which(sub.obj$seurat_clusters == '1')] = 'epidermis_BL_early_2'
  
  DimPlot(sub.obj, group.by = 'subtypes', label = TRUE, repel = TRUE)
  
  ggsave(filename = paste0(outDir, celltype2subset, '_groups_subclustering_manual_labels.pdf'), 
         width = 8, height = 6) 
  
  aa$subtypes[match(colnames(sub.obj), colnames(aa))] = sub.obj$subtypes
  DimPlot(aa, group.by = 'subtypes', label = TRUE, repel = TRUE)
  
  
  
}


## manually merge again the CT
merge_CT_subclusters = FALSE
if(merge_CT_subclusters){
  aa$subtypes[grep('CT_BL_early', aa$subtypes)] = 'CT_BL_early'
  
  aa$subtypes[grep('CT_CSD_early', aa$subtypes)] = 'CT_CSD_early'
  
  DimPlot(aa, group.by = 'subtypes', label = TRUE, repel = TRUE)
  
  ggsave(filename = paste0(resDir,  '/UMAP_subtypes_afterMergingCT.pdf'), 
         width = 12, height = 8)
  
}


aa$subtypes = factor(aa$subtypes, levels = sort(unique(aa$subtypes)))

aa = ScaleData(aa, features = rownames(aa))

Idents(aa) = aa$subtypes
source("functions_ligandReceptor_analysis.R")
aa$celltypes = aa$subtypes

## run LIANA
#source(paste0(functionDir, "/functions_cccInference.R"))
#aa$celltypes = aa$subtypes

liana_test = run_LIANA_defined_celltype(subref = aa, 
                           celltypes = unique(aa$celltypes),
                           additionalLabel = additionalLabel)

#liana_test <- liana_test %>%
#  liana_aggregate(resource = 'Consensus')

liana_test = readRDS(file = paste0(outDir, '/res_lianaTest_Consensus', 
                                   additionalLabel, '.rds'))

liana_test <- liana_test %>%
  liana_aggregate(resource = 'Consensus')

celltypes = unique(aa$celltypes)

receivers = celltypes

df_test = liana_test %>% filter(target %in% receivers & source %in% celltypes) %>% as.data.frame() 

write.table(df_test, file = paste0(outDir, '/res_lianaTest_Consensus', additionalLabel, '.txt'), 
            sep = '\t', quote = FALSE)

saveRDS(liana_test, file = paste0(outDir, '/res_lianaTest_Consensus', 
                                  additionalLabel, '_saved.rds'))


## prepare the LIANA output for circoplot
res = df_test
res = res[, c(1:5, which(colnames(res) == 'natmi.edge_specificity'), 
                which(colnames(res) == 'sca.LRscore'))]
colnames(res)[1:4] = c('sender', 'receiver', 'ligand', 'receptor')


#require(cellcall)
library(SeuratData)
library(Connectome)
library(cowplot)

res = res[order(-res$sca.LRscore), ]

#which(res$ligand == 'GAS6' & res$receptor == 'AXL')
#res[which(res$ligand == 'GAS6' & res$receptor == 'AXL'), ]

colnames(res)[1:2] = c('source', 'target')
res$weight_norm = res$sca.LRscore
res$pair = paste0(res$ligand, ' - ', res$receptor)
res$vector = paste0(res$source, ' - ', res$target)
res$edge = paste0(res$source, ' - ', res$ligand, ' - ', res$receptor, ' - ', res$target)
res$source.ligand = paste0(res$source, ' - ', res$ligand)
res$receptor.target = paste0(res$receptor, ' - ', res$target)

write.table(res, 
            file = paste0(outDir, '/LR_interactions_allPairs_LIANA.txt'), 
            quote = FALSE, row.names = TRUE, col.names = TRUE, sep = '\t')

saveRDS(res, file = paste0(outDir, '/res_lianaTest_for_circosplot.rds'))

##########################################
### plot circosplot
##########################################
#res = readRDS(file = paste0(outDir, '/res_lianaTest_for_circosplot.rds'))
source(paste0(functionDir, '/functions_cccInference.R'))

celltypes = unique(aa$celltypes)
additionalLabel = '_fixedCelltypes'
receivers = celltypes

sender_cells = celltypes[grep('BL', celltypes)]
receiver_cells = sender_cells[grep('CT', sender_cells)]

print(as.character(sender_cells))
print(as.character(receiver_cells))

res = readRDS(file = paste0(outDir, '/res_lianaTest_for_circosplot.rds'))

head(grep('WNT', res$ligand))
res = res[!is.na(match(res$source, sender_cells)) & !is.na(match(res$target, receiver_cells)), ]

res = res[order(res$aggregate_rank), ]

head(grep('WNT', res$ligand))

#cell_color = randomcoloR::distinctColorPalette(length(cells.of.interest))

cells.of.interest = unique(c(res$source, res$target))
print(cells.of.interest)

cell_color = c("#2AD0B7", "#98B304", "#16005e", "#005e45", "#3200F5")
names(cell_color) <- c("mac_BL_early", "neu_BL_early", "epidermis_BL_early",
                       "CT_BL_early", "epidermis_BL.CSD_early")
cell_color = cell_color[match(cells.of.interest, names(cell_color))]

# Mac_BL_early #2AD0B7
# Neu_BL_early #98B304
# Epidermis_BL_early #16005e
# CT_BL_early_1 #005e45
# CT_BL_early_3 #2f7c67
# Epidermis_BL.CSD_early #3200F5

pdfname = paste0(outDir, '/LR_interactions_LIANA_tops_BL_all_v2.pdf')
pdf(pdfname, width=12, height = 8)
for(ntop in c(100, 200, 300))
{
  # ntop = 300
  cat('top LR -- ', ntop, '\n')
  test = res[c(1:ntop), ]
  
  jj = which(test$ligand == 'RGMB'|test$receptor == 'RGMB'|
               test$ligand == "FGFR3")
  if(length(jj) >0){
    test = test[-jj, ]
  }
  
  #test = test[-which(test$ligand == 'SPON1'), ] 
  
  my_CircosPlot(test, 
                weight.attribute = 'weight_norm',
                cols.use = cell_color,
                sources.include = cells.of.interest,
                targets.include = cells.of.interest,
                lab.cex = 0.5,
                title = paste('LR scores top :', ntop))
  
}

dev.off()


sender_cells = celltypes[grep('CSD', celltypes)]
receiver_cells = sender_cells[grep('CT', sender_cells)]

print(as.character(sender_cells))
print(as.character(receiver_cells))

res = readRDS(file = paste0(outDir, '/res_lianaTest_for_circosplot.rds'))

res = res[!is.na(match(res$source, sender_cells)) & !is.na(match(res$target, receiver_cells)), ]

res = res[order(res$aggregate_rank), ]
head(grep('WNT', res$ligand))

cells.of.interest = unique(c(res$source, res$target))
print(cells.of.interest)

cell_color = c("#2AD0B7", "#98B304", "#d693b8", "#3200F5")
names(cell_color) <- c("mac_CSD_early", "neu_CSD_early",  "CT_CSD_early", "epidermis_BL.CSD_early")

# Mac_CSD_early #2AD0B7
# Neu_CSD_early #98B304
# CT_CSD_early_1 #d693b8
# CT_CSD_early_2 #925677
# CT_ CSD_early_3 #61394f
# Epidermis_BL.CSD_early #3200F5

pdfname = paste0(outDir, '/LR_interactions_LIANA_tops_CSD_all_v2.pdf')
pdf(pdfname, width=12, height = 8)

for(ntop in c(100, 200, 300))
{
  # ntop = 200
  cat('top LR -- ', ntop, '\n')
  test = res[c(1:ntop), ]
  
  # test = test[-which(test$ligand == 'SPON1'), ] ## because one gene can't be both ligand and receptor
  # https://github.com/msraredon/Connectome/issues/8
  jj = which(test$ligand == 'RGMB'|test$receptor == 'RGMB'|
               test$ligand == 'APP' | test$ligand == 'APP')
  
  if(length(jj)>0){
    test = test[-jj, ]
  }
  
  #cells.of.interest = unique(c(test$source, test$target))
  #cell_color = randomcoloR::distinctColorPalette(length(cells.of.interest))
  #names(cell_color) <- cells.of.interest
  
  my_CircosPlot(test, 
                weight.attribute = 'weight_norm',
                cols.use = cell_color,
                sources.include = cells.of.interest,
                targets.include = cells.of.interest,
                lab.cex = 0.5,
                title = paste('LR scores top :', ntop))
  
}

dev.off()


##########################################
# circosPlot for BL and CSD-specific LR
##########################################
res = readRDS(file = paste0(outDir, '/res_lianaTest_for_circosplot.rds'))

celltypes = unique(aa$celltypes)
receivers = celltypes

sender_cells = celltypes[grep('BL', celltypes)]
receiver_cells = sender_cells[grep('CT', sender_cells)]

print(as.character(sender_cells))
print(as.character(receiver_cells))

head(grep('WNT', res$ligand))

res_BL = res[!is.na(match(res$source, sender_cells)) & !is.na(match(res$target, receiver_cells)), ]
res_BL = res_BL[order(res_BL$aggregate_rank), ]

head(grep('WNT', res_BL$ligand))

ntop = 200
cat('top LR -- ', ntop, '\n')
test = res_BL[c(1:ntop), ]

jj = which(test$ligand == 'RGMB'|test$receptor == 'RGMB'|
             test$ligand == "FGFR3")
if(length(jj) >0){
  test = test[-jj, ]
}

saveRDS(test, file = paste0(outDir, '/res_lianaTest_for_circosplot_BL_ntop200.rds'))


res = readRDS(file = paste0(outDir, '/res_lianaTest_for_circosplot.rds'))
sender_cells = celltypes[grep('CSD', celltypes)]
receiver_cells = sender_cells[grep('CT', sender_cells)]

print(as.character(sender_cells))
print(as.character(receiver_cells))

res = res[!is.na(match(res$source, sender_cells)) & !is.na(match(res$target, receiver_cells)), ]

res = res[order(res$aggregate_rank), ]
head(grep('WNT', res$ligand))

ntop = 200
cat('top LR -- ', ntop, '\n')
test = res[c(1:ntop), ]

jj = which(test$ligand == 'RGMB'|test$receptor == 'RGMB'|
             test$ligand == 'APP' | test$ligand == 'APP')

if(length(jj)>0){
  test = test[-jj, ]
}
saveRDS(test, file = paste0(outDir, '/res_lianaTest_for_circosplot_CSD_ntop200.rds'))


Idents(aa) = as.factor(aa$celltypes)
source(paste0(functionDir, '/functions_cccInference.R'))

test_bl = readRDS(file = paste0(outDir, '/res_lianaTest_for_circosplot_BL_ntop200.rds'))
test_csd = readRDS(file = paste0(outDir, '/res_lianaTest_for_circosplot_CSD_ntop200.rds'))

ggs = unique(c(test_bl$ligand, test_bl$receptor, test_csd$ligand, test_csd$receptor))
xx = c()
for(n in 1:length(ggs))
{
  xx = c(xx,  unlist(strsplit(as.character(ggs[n]), '_')))
}
ggs = unique(xx)
rm(xx)

ggs = ggs[which(!is.na(match(ggs, rownames(aa))))]

cells_bl = unique(c(test_bl$source, test_bl$target))
cells_csd = unique(c(test_csd$source, test_csd$target))

#rm(test)
ct_bl_csd <- FindMarkers(aa, ident.1 ='CT_BL_early', ident.2 =  "CT_CSD_early",
                         features = ggs, logfc.threshold = 0, min.pct = 0, only.pos = FALSE)
mac_bl_csd = FindMarkers(aa, ident.1 ='mac_BL_early', ident.2 =  "mac_CSD_early",
                         features = ggs, logfc.threshold = 0, min.pct = 0, only.pos = FALSE)
neu_bl_csd = FindMarkers(aa, ident.1 ='neu_BL_early', ident.2 =  "neu_CSD_early",
                         features = ggs, logfc.threshold = 0, min.pct = 0, only.pos = FALSE)
ep_bl_csd = FindMarkers(aa, ident.1 ='epidermis_BL_early', ident.2 =  "epidermis_BL.CSD_early",
                        features = ggs, logfc.threshold = 0, min.pct = 0, only.pos = FALSE)

### compute the enrichment score (ES) for Blastema
test_bl$logFC_ligand = NA
test_bl$logFC_receptor = NA
test_bl$logPct_ligand = NA
test_bl$logPct_receptor = NA
test_bl$es = NA

for(n in 1:nrow(test_bl))
{
  # n = 5
  cat(n, '\n')
  # receptors
  receptors = unlist(strsplit(as.character(test_bl$receptor[n]), '_'))
  #ii2 = which(rownames(ct_bl_csd) == test_bl$receptor[n])
  ii2 = which(!is.na(match(rownames(ct_bl_csd),  receptors)))
  test_bl$logFC_receptor[n] = mean(ct_bl_csd$avg_log2FC[ii2])
  test_bl$logPct_receptor[n] = mean(log2(ct_bl_csd$pct.1[ii2]/ct_bl_csd$pct.2[ii2]))
  
  # ligands
  ligands = unlist(strsplit(as.character(test_bl$ligand[n]), '_'))
  if(test_bl$vector[n] == "CT_BL_early - CT_BL_early"){
    #ii1 = which(rownames(ct_bl_csd) == test_bl$ligand[n])
    ii1 = which(!is.na(match(rownames(ct_bl_csd),  ligands)))
    test_bl$logFC_ligand[n] = mean(ct_bl_csd$avg_log2FC[ii1])
    test_bl$logPct_ligand[n] = mean(log2(ct_bl_csd$pct.1[ii1]/ct_bl_csd$pct.2[ii1]))
  }
  
  if(test_bl$vector[n] == "mac_BL_early - CT_BL_early"){
    ii1 = which(!is.na(match(rownames(mac_bl_csd),  ligands)))
    #ii1 = which(rownames(mac_bl_csd) == test_bl$ligand[n])
    test_bl$logFC_ligand[n] = mean(mac_bl_csd$avg_log2FC[ii1])
    test_bl$logPct_ligand[n] = mean(log2(mac_bl_csd$pct.1[ii1]/mac_bl_csd$pct.2[ii1]))
  }
  
  if(test_bl$vector[n] == "neu_BL_early - CT_BL_early"){
    ii1 = which(!is.na(match(rownames(neu_bl_csd),  ligands)))
    #ii1 = which(rownames(neu_bl_csd) == test_bl$ligand[n])
    test_bl$logFC_ligand[n] = mean(neu_bl_csd$avg_log2FC[ii1])
    test_bl$logPct_ligand[n] = mean(log2(neu_bl_csd$pct.1[ii1]/neu_bl_csd$pct.2[ii1]))
  }
  
  if(test_bl$vector[n] == "epidermis_BL.CSD_early - CT_BL_early"){
    #ii1 = which(rownames(neu_bl_csd) == test_bl$ligand[n])
    test_bl$logFC_ligand[n] = 0
    test_bl$logPct_ligand[n] = log2(1)
  }
  
  if(test_bl$vector[n] == "epidermis_BL_early - CT_BL_early"){
    ii1 = which(!is.na(match(rownames(ep_bl_csd),  ligands)))
    #ii1 = which(rownames(neu_bl_csd) == test_bl$ligand[n])
    test_bl$logFC_ligand[n] = mean(ep_bl_csd$avg_log2FC[ii1])
    test_bl$logPct_ligand[n] = mean(log2(ep_bl_csd$pct.1[ii1]/ep_bl_csd$pct.2[ii1]))
  }
}

test_bl$fc_ligand = 2^test_bl$logFC_ligand
test_bl$fc_receptor = 2^test_bl$logFC_receptor
test_bl$pct_ligand = 2^test_bl$logPct_ligand
test_bl$pct_receptor = 2^test_bl$logPct_receptor

test_bl$es = test_bl$logFC_ligand + test_bl$logFC_receptor + test_bl$logPct_ligand + test_bl$logPct_receptor
test_bl = test_bl[order(-test_bl$es), ]

cells.of.interest = unique(c(test_bl$source, test_bl$target))
print(cells.of.interest)

cell_color = c("#2AD0B7", "#98B304", "#16005e", "#005e45", "#3200F5")
names(cell_color) <- c("mac_BL_early", "neu_BL_early", "epidermis_BL_early",
                       "CT_BL_early", "epidermis_BL.CSD_early")
cell_color = cell_color[match(cells.of.interest, names(cell_color))]

pdfname = paste0(outDir, '/LR_interactions_LIANA_tops_BL_specific_v2.pdf')
pdf(pdfname, width=12, height = 8)
for(es_cut in seq(1, 4, by = 1))
{
  # es_cut = 
  test = test_bl[which(test_bl$es > es_cut), ]
  cat('Es cutoff  -- ', es_cut, '--', nrow(test), ' LR \n')
  
  jj = which(test$ligand == 'RGMB'|test$receptor == 'RGMB'|
               test$ligand == "FGFR3")
  if(length(jj) >0){
    test = test[-jj, ]
  }
  
  #test = test[-which(test$ligand == 'SPON1'), ] 
  
  my_CircosPlot(test, 
                weight.attribute = 'weight_norm',
                cols.use = cell_color,
                sources.include = cells.of.interest,
                targets.include = cells.of.interest,
                lab.cex = 0.5,
                title = paste('ES cutoff :', es_cut))
  
}
dev.off()


### compute the enrichment score (ES) for CSD
test_csd$logFC_ligand = NA
test_csd$logFC_receptor = NA
test_csd$logPct_ligand = NA
test_csd$logPct_receptor = NA
test_csd$es = NA

for(n in 1:nrow(test_csd))
{
  # n = 1
  cat(n, '\n')
  # receptors
  receptors = unlist(strsplit(as.character(test_csd$receptor[n]), '_'))
  #ii2 = which(rownames(ct_bl_csd) == test_csd$receptor[n])
  ii2 = which(!is.na(match(rownames(ct_bl_csd),  receptors)))
  test_csd$logFC_receptor[n] = -mean(ct_bl_csd$avg_log2FC[ii2])
  test_csd$logPct_receptor[n] = -mean(log2(ct_bl_csd$pct.1[ii2]/ct_bl_csd$pct.2[ii2]))
  
  # ligands
  ligands = unlist(strsplit(as.character(test_csd$ligand[n]), '_'))
  if(test_csd$vector[n] == "CT_CSD_early - CT_CSD_early"){
    #ii1 = which(rownames(ct_bl_csd) == test_csd$ligand[n])
    ii1 = which(!is.na(match(rownames(ct_bl_csd),  ligands)))
    test_csd$logFC_ligand[n] = -mean(ct_bl_csd$avg_log2FC[ii1])
    test_csd$logPct_ligand[n] = -mean(log2(ct_bl_csd$pct.1[ii1]/ct_bl_csd$pct.2[ii1]))
  }
  
  if(test_csd$vector[n] == "mac_CSD_early - CT_CSD_early"){
    ii1 = which(!is.na(match(rownames(mac_bl_csd),  ligands)))
    #ii1 = which(rownames(mac_bl_csd) == test_csd$ligand[n])
    test_csd$logFC_ligand[n] = -mean(mac_bl_csd$avg_log2FC[ii1])
    test_csd$logPct_ligand[n] = -mean(log2(mac_bl_csd$pct.1[ii1]/mac_bl_csd$pct.2[ii1]))
  }
  
  if(test_csd$vector[n] == "neu_CSD_early - CT_CSD_early"){
    ii1 = which(!is.na(match(rownames(neu_bl_csd),  ligands)))
    #ii1 = which(rownames(neu_bl_csd) == test_csd$ligand[n])
    test_csd$logFC_ligand[n] = -mean(neu_bl_csd$avg_log2FC[ii1])
    test_csd$logPct_ligand[n] = -mean(log2(neu_bl_csd$pct.1[ii1]/neu_bl_csd$pct.2[ii1]))
  }
  
  if(test_csd$vector[n] == "epidermis_BL.CSD_early - CT_CSD_early"){
    #ii1 = which(rownames(neu_bl_csd) == test_csd$ligand[n])
    test_csd$logFC_ligand[n] = 0
    test_csd$logPct_ligand[n] = log2(1)
  }
  
  # if(test_csd$vector[n] == "epidermis_BL_early - CT_BL_early"){
  #   ii1 = which(!is.na(match(rownames(ep_bl_csd),  ligands)))
  #   #ii1 = which(rownames(neu_bl_csd) == test_csd$ligand[n])
  #   test_csd$logFC_ligand[n] = mean(ep_bl_csd$avg_log2FC[ii1])
  #   test_csd$logPct_ligand[n] = mean(log2(ep_bl_csd$pct.1[ii1]/ep_bl_csd$pct.2[ii1]))
  # }
}

test_csd$fc_ligand = 2^test_csd$logFC_ligand
test_csd$fc_receptor = 2^test_csd$logFC_receptor
test_csd$pct_ligand = 2^test_csd$logPct_ligand
test_csd$pct_receptor = 2^test_csd$logPct_receptor

test_csd$es = test_csd$logFC_ligand + test_csd$logFC_receptor + 
  test_csd$logPct_ligand + test_csd$logPct_receptor
test_csd = test_csd[order(-test_csd$es), ]

cells.of.interest = unique(c(test_csd$source, test_csd$target))
print(cells.of.interest)
cell_color = c("#2AD0B7", "#98B304", "#d693b8", "#3200F5")
names(cell_color) <- c("mac_CSD_early", "neu_CSD_early",  "CT_CSD_early", "epidermis_BL.CSD_early")


pdfname = paste0(outDir, '/LR_interactions_LIANA_tops200_CSD_specific.pdf')
pdf(pdfname, width=12, height = 8)
for(es_cut in seq(1, 4, by = 1))
{
  # es_cut = 
  test = test_csd[which(test_csd$es > es_cut), ]
  cat('Es cutoff  -- ', es_cut, '--', nrow(test), ' LR \n')
  
  #jj = which(test$ligand == 'RGMB'|test$receptor == 'RGMB'|
  #             test$ligand == "FGFR3")
  jj = which(test$ligand == 'RGMB'|test$receptor == 'RGMB'|
               test$ligand == 'APP' | test$ligand == 'APP')
  
  if(length(jj) >0){
    test = test[-jj, ]
  }
  
  #test = test[-which(test$ligand == 'SPON1'), ] 
  my_CircosPlot(test, 
                weight.attribute = 'weight_norm',
                cols.use = cell_color,
                sources.include = cells.of.interest,
                targets.include = cells.of.interest,
                lab.cex = 0.5,
                title = paste('ES cutoff :', es_cut))
  
}

dev.off()


##########################################
##########################################
# Test NicheNet 
##########################################
##########################################
aa = readRDS(file = paste0(RdataDir, '/BL.CSD_merged_subset_CT_MAC_Neu_Epd_day3_5_8_subtypes_umap.rds'))

## manually merge subclusters
aa$subtypes[which(aa$subtypes == 'epidermis_BL_early_2')] = 'epidermis_BL_early'
aa$subtypes[which(aa$subtypes == 'epidermis_BL_early_1')] = 'epidermis_BL_early'
aa$subtypes[which(aa$subtypes == 'CT_BL_early_2')] = 'CT_BL_early_1'

## shorten the subtype names
aa$subtypes[which(aa$subtypes == 'epidermis_BL_early')] = 'epidermis_BL'
aa$subtypes[which(aa$subtypes == 'epidermis_BL.CSD_early')] = 'epidermis_BL.CSD'
aa$subtypes[which(aa$subtypes == 'mac_BL_early')] = 'mac_BL'
aa$subtypes[which(aa$subtypes == 'mac_CSD_early')] = 'mac_CSD'
aa$subtypes[which(aa$subtypes == 'neu_BL_early')] = 'neu_BL'
aa$subtypes[which(aa$subtypes == 'neu_CSD_early')] = 'neu_CSD'

aa$subtypes[which(aa$subtypes == 'CT_BL_early_1')] = 'CT_BL'
aa$subtypes[which(aa$subtypes == 'CT_CSD_early_1')] = 'CT_CSD'

DimPlot(aa, group.by = 'subtypes')

aa$subtypes = factor(aa$subtypes, levels = sort(unique(aa$subtypes)))

aa = ScaleData(aa, features = rownames(aa))

Idents(aa) = aa$subtypes


## test NicheNet
outDir = paste0(resDir, '/LR_analysis_NicheNet_v1')
source("functions_ligandReceptor_analysis.R")


########################################################
########################################################
# Section III: native way of checking gene examples of pathways of interest
# (RA, WNT, BMP, TGF-beta, FGF, Noch, Yap, SHH, PI3K)
########################################################
########################################################

##########################################
# plot gene examples of selected signaling pathways
##########################################
curate.geneList.signalingPathways = FALSE
if(curate.geneList.signalingPathways){
  sps = readRDS(file = paste0('../data/annotations/curated_signaling.pathways_gene.list_v2.rds'))
  #sps = unique(sps$gene)
  xx = read.table('../data/annotations/GO_term_summary_RAR_signalingPathway.txt', header = TRUE, sep = '\t', 
                  row.names = NULL)
  xx = unique(xx[,2])
  xx = data.frame(gene = xx, pathway = rep('RA', length(xx)))
  sps = rbind(sps, xx)
  
  
  xx = read.table('../data/annotations/GO_term_summary_TGFb.txt', header = TRUE, sep = '\t', 
                  row.names = NULL)
  xx = unique(xx[,2])
  xx = data.frame(gene = xx, pathway = rep('TGF_beta', length(xx)))
  sps = rbind(sps, xx)
  
  xx = read.table('../data/annotations/GO_term_summary_NOTCH.txt', header = TRUE, sep = '\t', 
                  row.names = NULL)
  xx = unique(xx[,2])
  xx = data.frame(gene = xx, pathway = rep('NOTCH', length(xx)))
  sps = rbind(sps, xx)
  
  xx = read.table('../data/annotations/GO_term_summary_PI3K.txt', header = TRUE, sep = '\t', 
                  row.names = NULL)
  xx = unique(xx[,2])
  xx = data.frame(gene = xx, pathway = rep('PI3K', length(xx)))
  sps = rbind(sps, xx)
  
  xx = read.table('../data/annotations/GO_term_summary_Hippo.txt', header = TRUE, sep = '\t', 
                  row.names = NULL)
  xx = unique(xx[,2])
  xx = data.frame(gene = xx, pathway = rep('Hippo', length(xx)))
  sps = rbind(sps, xx)
  
  saveRDS(sps, file = paste0('../data/annotations/curated_signaling.pathways_gene.list_v3.rds'))
  
  
}


saveDir = paste0(resDir, '/geneExamples_pathways/')
system(paste0('mkdir -p ', saveDir))

DefaultAssay(aa) = 'RNA'

cell_ids  = 'groups'

## gene example of signaling pathways 
sps = readRDS(file = paste0("/groups/tanaka/People/current/jiwang/projects/RA_competence/",
                            '/data/annotations/curated_signaling.pathways_gene.list_v3.rds'))
sps = unique(sps$gene)
sps = toupper(sps)

## TFs 
tfs = readRDS(file = paste0('/groups/tanaka/People/current/jiwang/projects/RA_competence/data',
                            '/annotations/curated_human_TFs_Lambert.rds'))
tfs = unique(tfs$`HGNC symbol`)
#tfs = as.character(unlist(sapply(tfs, firstup)))

## ligands and receptors from NicheNet
dataPath_nichenet = '/users/jingkui.wang/projects/heart_regeneration/data/NicheNet/'

ligand_target_matrix = readRDS(paste0(dataPath_nichenet,  "ligand_target_matrix.rds"))
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns

# ligand-receptor network, Putative ligand-receptor links from NicheNet
lr_network = readRDS(paste0(dataPath_nichenet, "lr_network.rds"))
lr_network = lr_network %>% mutate(bonafide = ! database %in% c("ppi_prediction","ppi_prediction_go"))
lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% 
  distinct(ligand, receptor, bonafide)

head(lr_network)

ligands = lr_network %>% pull(ligand) %>% unique()
receptors = lr_network %>% pull(receptor) %>% unique()


##########################################
# manual subtyping early BL and CSD clusters
##########################################
markers_saved = c()

aa = readRDS(file = paste0(RdataDir, '/BL.CSD_merged_subset_CT_MAC_Neu_Epd_day3_5_8_subtypes_umap.rds'))
aa$subtypes = factor(aa$subtypes, levels = sort(unique(aa$subtypes)))

aa = ScaleData(aa, features = rownames(aa))
Idents(aa) = aa$subtypes

DimPlot(aa, group.by = 'subtypes', label = TRUE, repel = TRUE)

markers = FindAllMarkers(aa, only.pos = TRUE,  
                         #min.pct = 0.25, 
                         logfc.threshold = 0.25)

markers_saved = unique(c(markers_saved, markers$gene))


## call marker genes
cell_ids = 'manual_subtyping'

markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top10

DoHeatmap(aa, features = top10$gene) + NoLegend()
ggsave(filename = paste0(saveDir,  'heatmap_markerGenes_top20_by_', cell_ids, '.pdf'), 
       width = 12, height = 30)

features = intersect(markers$gene[which(markers$p_val < 10^-100)], sps)
cat(length(features), ' genes to display \n')

DoHeatmap(aa, features = features) + NoLegend()
ggsave(filename = paste0(saveDir,  'heatmap_markerGenes_signalingGenes_by_', cell_ids, '.pdf'), 
       width = 12, height = 40)

tfs_sel = intersect(tfs, markers$gene[which(markers$p_val < 10^-50)])

DoHeatmap(aa, features = tfs_sel) + NoLegend()
ggsave(filename = paste0(saveDir,  'heatmap_markerGenes_TFs_by_', cell_ids, '.pdf'), 
       width = 12, height = 30)


features = intersect(markers$gene[which(markers$p_val < 10^-100)], receptors)
cat(length(features), ' genes to display \n')

DoHeatmap(aa, features = features) + NoLegend()
ggsave(filename = paste0(saveDir,  'heatmap_markerGenes_receptors_by_', cell_ids, '.pdf'), 
       width = 12, height = 40)

features = intersect(markers$gene[which(markers$p_val < 10^-100)], ligands)
cat(length(features), ' genes to display \n')

DoHeatmap(aa, features = features) + NoLegend()
ggsave(filename = paste0(saveDir,  'heatmap_markerGenes_ligands_by_', cell_ids, '.pdf'), 
       width = 12, height = 40)


genes = unique(c(sps, tfs, ligands, receptors))
genes = genes[!is.na(match(genes, rownames(aa)))]
genes = intersect(genes, markers$gene[which(markers$p_val < 10^-50)])

aa$condition = droplevels(aa$condition)
#cols[2:4] = colorRampPalette((brewer.pal(n = 3, name ="Blues")))(3)
#cols[5:8] = colorRampPalette((brewer.pal(n = 4, name ="OrRd")))(4)

cols_sel = c("#9ECAE1", "#3182BD", 
             colorRampPalette((brewer.pal(n = 3, name ="OrRd")))(3))
DimPlot(aa, group.by = 'condition', cols = cols_sel)

p1 = DimPlot(aa, group.by = 'condition', cols = cols_sel, label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'celltype', label = TRUE, repel = TRUE) + NoLegend()
p2 + p1

ggsave(filename = paste0(resDir, '/BL.CSD_merged_subset_CT_MAC_Neu_Epd_day3_5_8_subtypes_conditions.pdf'), 
       width = 12, height = 6)

# feature plots
pdfname = paste0(saveDir, "gene_examples_sps_tfs_ligands_receptors_v2.pdf")
pdf(pdfname, width=16, height = 8)
par(cex =0.7, mar = c(3,0.8,2,5)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)

#for(n in 1:20)
for(n in 1:length(genes))
{
  # n = 1
  gg = genes[n];
  # gg = "Foxa2"
  cat(n, '--', gg, '\n')

  p1 = FeaturePlot(aa, features = gg)
  p2 = VlnPlot(aa, features = gg, group.by = 'subtypes', pt.size = 0.01)  + NoLegend()
    #scale_fill_manual(values=c("blue", 'green', "red",  "blue",
    #                           'red', 'green', 'blue', 'green', 'red', 'green', 'green'))
  #p3 = VlnPlot(aa, features = gg, group.by = 'clusters', pt.size = 0.01) + NoLegend()
  plot(p1 | p2)

}

dev.off()
