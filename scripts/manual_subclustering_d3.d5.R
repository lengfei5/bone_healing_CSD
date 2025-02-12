##########################################################################
##########################################################################
# Project: blastema vs CSD 
# Script purpose: subclustering d3.d5 blastema and CSD
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Mar 13 15:28:52 2023
##########################################################################
##########################################################################

system(paste0('mkdir -p ', resDir, '/manual_annotation/'))

outDir = paste0(resDir, '/manual_annotation/')
cols_sel = cols[c(2, 5, 6)]

library(Seurat)
library(decoupleR)
library(tictoc)
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)

dataDir = '../results/scRNAseq_signaling.analysis_axolotl_20230308/Rdata/'

aa = readRDS(file = paste0(dataDir, '/BL.CSD_merged_subset_CT_MAC_Neu_Epd_umap.rds'))

DimPlot(aa, group.by = 'celltype', label = TRUE, repel = TRUE)

aa$celltype[which(aa$celltype == 'connective tissue')] = 'CT'
Idents(aa) = aa$celltype
aa$condition = factor(aa$condition, levels = levels)

aa$groups = paste0(aa$celltype, '_', aa$condition)

Idents(aa) = aa$condition
xx = subset(aa, idents= c('BL_3_5dpa', 'BL_8dpa', 'CSD_3dpa', 'CSD_5dpa', 'CSD_8dpa'))
counts = xx@assays$RNA@counts

rownames(counts) = convert_to_geneSymbol(rownames(counts), annot = annot)
counts = counts[grep('^AMEX', rownames(counts), invert = TRUE), ]

aa = CreateSeuratObject(counts = counts, assay = 'RNA', meta.data = xx@meta.data) 
aa = NormalizeData(aa, scale.factor = 10000)

# rerun the umap
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000)
aa = ScaleData(aa)
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(aa, ndims = 50)


aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 30, min.dist = 0.1)

Idents(aa) = aa$condition

p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltype')

p1 / p2

saveRDS(aa, file = paste0(RdataDir, '/BL.CSD_merged_subset_CT_MAC_Neu_Epd_day3_5_8.rds'))


########################################################
########################################################
# Section :  subset each cell types 
# 
########################################################
########################################################
aa = readRDS(file = paste0(RdataDir, '/BL.CSD_merged_subset_CT_MAC_Neu_Epd_day3_5_8.rds'))
Idents(aa) = aa$celltype
aa$subtypes = NA

##########################################
# subclustering CT
##########################################
celltype2subset = 'CT'
sub.obj = subset(aa, idents = celltype2subset) 

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

p1 = DimPlot(sub.obj, group.by = 'condition', cols = cols, label = TRUE, repel = TRUE)
p2 = DimPlot(sub.obj, group.by = 'Phase', label = TRUE, repel = TRUE)

p1 + p2

ggsave(filename = paste0(outDir, celltype2subset, '_condtion_cellcycle.pdf'), 
       width = 16, height = 8)

sub.obj <- FindNeighbors(sub.obj, dims = 1:20,  k.param = 30)
sub.obj <- FindClusters(sub.obj, resolution = 0.5)

p1 = DimPlot(sub.obj, label = TRUE, repel = TRUE)
p2 = DimPlot(sub.obj, group.by = 'groups', label = TRUE, repel = TRUE)

p1 | p2

ggsave(filename = paste0(outDir, celltype2subset, '_groups_subclustering.pdf'), 
       width = 14, height = 6)

## assign new labels
sub.obj$subtypes = NA
sub.obj$subtypes[which(sub.obj$sample == 'BL' & sub.obj$seurat_clusters == '1')] = 'CT_BL_early_1'
sub.obj$subtypes[which(sub.obj$sample == 'BL' & sub.obj$seurat_clusters == '2')] = 'CT_BL_early_2'
sub.obj$subtypes[which(sub.obj$sample == 'BL' & sub.obj$seurat_clusters == '3')] = 'CT_BL_early_3'

sub.obj$subtypes[which(sub.obj$sample == 'CSD' & sub.obj$seurat_clusters == '1')] = 'CT_CSD_early_1'
sub.obj$subtypes[which(sub.obj$sample == 'CSD' & sub.obj$seurat_clusters == '0')] = 'CT_CSD_early_2'
sub.obj$subtypes[which(sub.obj$sample == 'CSD' & sub.obj$seurat_clusters == '4')] = 'CT_CSD_early_3'

DimPlot(sub.obj, group.by = 'subtypes', label = TRUE, repel = TRUE)

ggsave(filename = paste0(outDir, celltype2subset, '_groups_subclustering_manual_labels.pdf'), 
       width = 8, height = 6) 

aa$subtypes[match(colnames(sub.obj), colnames(aa))] = sub.obj$subtypes
DimPlot(aa, group.by = 'subtypes', label = TRUE, repel = TRUE)


##########################################
# subclustering epidermis
##########################################
celltype2subset = 'epidermis'

cells = colnames(aa)[which(aa$celltype == celltype2subset & (aa$condition == 'BL_3_5dpa'|
                                                               aa$condition == 'CSD_3dpa'))]
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

p1 = DimPlot(sub.obj, group.by = 'condition', cols = cols, label = TRUE, repel = TRUE)
p2 = DimPlot(sub.obj, group.by = 'Phase', label = TRUE, repel = TRUE)

p1 + p2

ggsave(filename = paste0(outDir, celltype2subset, '_condtion_cellcycle.pdf'), 
       width = 16, height = 8)

sub.obj <- FindNeighbors(sub.obj, dims = 1:20,  k.param = 30)
sub.obj <- FindClusters(sub.obj, resolution = 0.3)

p1 = DimPlot(sub.obj, label = TRUE, repel = TRUE)
p2 = DimPlot(sub.obj, group.by = 'groups', label = TRUE, repel = TRUE)

p1 | p2

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


##########################################
# subclustering macrophage
##########################################
celltype2subset = 'macrophages'

cells = colnames(aa)[which(aa$celltype == celltype2subset & (aa$condition == 'BL_3_5dpa'|
                                                               aa$condition == 'CSD_3dpa'|
                                                               aa$condition == 'CSD_5dpa'))]
sub.obj = subset(aa, cells = cells) 

#sub.obj = subset(aa, idents = celltype2subset) 

sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = 3000)
sub.obj = ScaleData(sub.obj)
sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE)
ElbowPlot(sub.obj, ndims = 30)

nb.pcs = 20 # nb of pcs depends on the considered clusters or ids
n.neighbors = 30; min.dist = 0.3;
sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", 
                   dims = 1:nb.pcs, 
                   n.neighbors = n.neighbors,
                   min.dist = min.dist)

p1 = DimPlot(sub.obj, group.by = 'condition', cols = cols, label = TRUE, repel = TRUE)
p2 = DimPlot(sub.obj, group.by = 'Phase', label = TRUE, repel = TRUE)

p1 + p2

ggsave(filename = paste0(outDir, celltype2subset, '_condtion_cellcycle.pdf'), 
       width = 16, height = 8)

sub.obj <- FindNeighbors(sub.obj, dims = 1:30,  k.param = 30)
sub.obj <- FindClusters(sub.obj, resolution = 0.3)

p1 = DimPlot(sub.obj, label = TRUE, repel = TRUE)
p2 = DimPlot(sub.obj, group.by = 'groups', label = TRUE, repel = TRUE)

p1 | p2

ggsave(filename = paste0(outDir, celltype2subset, '_groups_subclustering.pdf'), 
       width = 14, height = 6)

## assign new labels
sub.obj$subtypes = NA
sub.obj$subtypes[which(sub.obj$seurat_clusters == '1'| sub.obj$seurat_clusters == '4')] = 
  'mac_BL_early'

sub.obj$subtypes[which(sub.obj$seurat_clusters == '0')] = 'mac_CSD_early'
#sub.obj$subtypes[which(sub.obj$seurat_clusters == '1')] = 'mac_CSD_early'

DimPlot(sub.obj, group.by = 'subtypes', label = TRUE, repel = TRUE)

ggsave(filename = paste0(outDir, celltype2subset, '_groups_subclustering_manual_labels.pdf'), 
       width = 8, height = 6) 

aa$subtypes[match(colnames(sub.obj), colnames(aa))] = sub.obj$subtypes
DimPlot(aa, group.by = 'subtypes', label = TRUE, repel = TRUE)


##########################################
# subclustering macrophage
##########################################
celltype2subset = 'neutrophils'

cells = colnames(aa)[which(aa$celltype == celltype2subset & (aa$condition == 'BL_3_5dpa'|
                                                               aa$condition == 'CSD_3dpa'|
                                                               aa$condition == 'CSD_5dpa'))]
sub.obj = subset(aa, cells = cells) 

#sub.obj = subset(aa, idents = celltype2subset) 

sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = 3000)
sub.obj = ScaleData(sub.obj)
sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE)
ElbowPlot(sub.obj, ndims = 30)

nb.pcs = 20 # nb of pcs depends on the considered clusters or ids
n.neighbors = 30; min.dist = 0.3;
sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", 
                   dims = 1:nb.pcs, 
                   n.neighbors = n.neighbors,
                   min.dist = min.dist)

p1 = DimPlot(sub.obj, group.by = 'condition', cols = cols, label = TRUE, repel = TRUE)
p2 = DimPlot(sub.obj, group.by = 'Phase', label = TRUE, repel = TRUE)

p1 + p2

ggsave(filename = paste0(outDir, celltype2subset, '_condtion_cellcycle.pdf'), 
       width = 16, height = 8)

sub.obj <- FindNeighbors(sub.obj, dims = 1:30,  k.param = 30)
sub.obj <- FindClusters(sub.obj, resolution = 0.3)

p1 = DimPlot(sub.obj, label = TRUE, repel = TRUE)
p2 = DimPlot(sub.obj, group.by = 'groups', label = TRUE, repel = TRUE)

p1 | p2

ggsave(filename = paste0(outDir, celltype2subset, '_groups_subclustering.pdf'), 
       width = 14, height = 6)

## assign new labels
sub.obj$subtypes = NA
sub.obj$subtypes[which(sub.obj$seurat_clusters == '1'| sub.obj$seurat_clusters == '0')] = 
  'neu_CSD_early'

sub.obj$subtypes[which(sub.obj$seurat_clusters == '3')] = 'neu_BL_early'
#sub.obj$subtypes[which(sub.obj$seurat_clusters == '1')] = 'mac_CSD_early'

DimPlot(sub.obj, group.by = 'subtypes', label = TRUE, repel = TRUE)

ggsave(filename = paste0(outDir, celltype2subset, '_groups_subclustering_manual_labels.pdf'), 
       width = 8, height = 6) 

aa$subtypes[match(colnames(sub.obj), colnames(aa))] = sub.obj$subtypes
DimPlot(aa, group.by = 'subtypes', label = TRUE, repel = TRUE)

saveRDS(aa, file = paste0(RdataDir, '/BL.CSD_merged_subset_CT_MAC_Neu_Epd_day3_5_8_subtypes.rds'))

##########################################
# umap of cells with new labels 
##########################################
aa = readRDS(file = paste0(RdataDir, '/BL.CSD_merged_subset_CT_MAC_Neu_Epd_day3_5_8_subtypes.rds'))
aa = subset(aa, cells = colnames(aa)[!is.na(aa$subtypes)])

aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000) # find subset-specific HVGs

#all.genes <- rownames(aa)
aa <- ScaleData(aa)

## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.3)

DimPlot(aa, group.by = 'subtypes',  label = TRUE, repel = TRUE)

#aa <- FindNeighbors(aa, dims = 1:30,  k.param = 30)
#aa <- FindClusters(aa, resolution = 0.5)

saveRDS(aa, file = paste0(RdataDir, '/BL.CSD_merged_subset_CT_MAC_Neu_Epd_day3_5_8_subtypes_umap.rds'))

