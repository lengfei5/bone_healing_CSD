##########################################################################
##########################################################################
# Project: Natasia's CSD project
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

version.analysis = '_axolotl_20230308'
resDir = paste0("../results/scRNAseq_signaling.analysis", version.analysis)
RdataDir = paste0(resDir, '/Rdata')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../scRNAseq_dataset/'
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
  return(make.unique(ggs))
  
}

cols = rep(NA, length = 8)
levels = c('CSD_0dpa', 'BL_3_5dpa',  'BL_8dpa', 'BL_11dpa', 
           'CSD_3dpa', 'CSD_5dpa', 'CSD_8dpa', 'CSD_11dpa')

names(cols) = levels

cols[1] = 'gray60'
#cols[1:3] = viridis(3)
cols[2:4] = colorRampPalette((brewer.pal(n = 3, name ="Blues")))(3)
cols[5:8] = colorRampPalette((brewer.pal(n = 4, name ="OrRd")))(4)

########################################################
########################################################
# Section 0: prepare the scRNA-seq samples  
# import the processed scRNA data from Natasia and Tobi 
########################################################
########################################################
aa = readRDS(file = paste0(dataDir, 'BL_SeuratObj.RDS'))
bb = readRDS(file = paste0(dataDir, 'CSD_SeuratObj.RDS'))

p1 = DimPlot(aa, group.by = 'celltype', label = TRUE, repel = TRUE)
p2 = DimPlot(bb, group.by = 'celltype', label = TRUE, repel = TRUE)

p1 + p2

ggsave(filename = paste0(resDir, '/Tobie_umap_BL.CSD.pdf'), width = 18, height = 6)

aa = merge(aa, bb)

rm(bb)

aa$condition = aa$orig.ident
aa$sample = aa$condition
aa$sample[grep('BL_', aa$sample)] = 'BL'
aa$sample[grep('CSD_', aa$sample)] = 'CSD'

aa$days = aa$condition
aa$days = gsub('BL_', '', aa$days)
aa$days = gsub('CSD_', '', aa$days)

aa = NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 8000) # find subset-specific HVGs

all.genes <- rownames(aa)
aa <- ScaleData(aa, features = all.genes)

## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, dims = 1:50, n.neighbors = 30, min.dist = 0.3)

p1 = DimPlot(aa, group.by = 'celltype', label = TRUE, repel = TRUE)

p2 = DimPlot(aa, group.by = 'sample', label = TRUE, repel = TRUE)

p1 + p2

ggsave(filename = paste0(resDir, '/UMAP_merged.BL.CSD_celltypes_samples.pdf'), width = 18, height = 6)

DimPlot(aa, group.by = 'celltype', split.by = 'sample', label = TRUE, repel = TRUE)

ggsave(filename = paste0(resDir, '/UMAP_merged.BL.CSD_celltypes_split.by.sample.pdf'), width = 18, height = 6)

saveRDS(aa, file = paste0(RdataDir, '/BL.CSD_merged_renormalized_8000HVGs_umap.rds'))

##########################################
# subset relevant cell populations:
# connective tissue, epidermis, macrophage, neutrophils
##########################################
aa = readRDS(file = paste0(RdataDir, '/BL.CSD_merged_renormalized_8000HVGs_umap.rds'))

Idents(aa) = aa$celltype

xx = subset(aa, idents = c('connective tissue', 'epidermis', 'macrophages', 'neutrophils'))
DimPlot(xx, group.by = 'celltype',  label = TRUE, repel = TRUE)

aa = xx
rm(xx)

aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 8000) # find subset-specific HVGs

all.genes <- rownames(aa)
aa <- ScaleData(aa, features = all.genes)

## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.3)

DimPlot(aa, group.by = 'celltype',  label = TRUE, repel = TRUE)

aa <- FindNeighbors(aa, dims = 1:30,  k.param = 30)
aa <- FindClusters(aa, resolution = 0.5)

DimPlot(aa, label = TRUE, repel = TRUE)

aa$clusters = aa$seurat_clusters

p1 = DimPlot(aa, group.by = 'celltype',  label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'clusters',  label = TRUE, repel = TRUE)

p1 + p2

ggsave(filename = paste0(resDir, 
                         '/UMAP_merged.BL.CSD_subset_celltypes_newClusters.pdf'), 
       width = 18, height = 6)

DimPlot(aa, group.by = 'celltype', split.by = 'sample', label = TRUE, repel = TRUE)

saveRDS(aa, file = paste0(RdataDir, '/BL.CSD_merged_subset_CT_MAC_Neu_Epd_umap.rds'))

# aa = readRDS(file = paste0(RdataDir, '/BL.CSD_merged_subset_CT_MAC_Neu_Epd_umap.rds'))


DimPlot(aa, group.by = 'condition', cols = cols)


##########################################
# check the condition and time changes for each cell type 
##########################################
aa = readRDS(file = paste0(RdataDir, '/BL.CSD_merged_subset_CT_MAC_Neu_Epd_umap.rds'))

aa$celltype[which(aa$celltype == 'connective tissue')] = 'CT'
Idents(aa) = aa$celltype


#aa$subtypes = NA

for(celltype in unique(aa$celltype))
{
  # celltype = 'neutrophils'
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
  
  p1 = DimPlot(sub.obj, group.by = 'condition', cols = cols, label = TRUE, repel = TRUE)
  p2 = DimPlot(sub.obj, group.by = 'Phase', label = TRUE, repel = TRUE)
  
  p1 + p2
  
  ggsave(filename = paste0(resDir, 
                           '/UMAP_merged.BL.CSD_subset_celltypes_checkEachcelltypes.aross.condition.time',
                           celltype,  '.pdf'), 
         width = 18, height = 6)
  
}


########################################################
########################################################
# Section I: test global pathway analysis 
#  - e.g. decoupleR, Pagoda2 (optional)
# decoupleR original code from https://saezlab.github.io/decoupleR/articles/pw_sc.html
# and https://www.sc-best-practices.org/conditions/gsea_pathway.html
########################################################
########################################################
library(Seurat)
library(decoupleR)
library(tictoc)
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)

aa = readRDS(file = paste0(RdataDir, '/BL.CSD_merged_subset_CT_MAC_Neu_Epd_umap.rds'))

aa$celltype[which(aa$celltype == 'connective tissue')] = 'CT'
Idents(aa) = aa$celltype
aa$condition = factor(aa$condition, levels = levels)

aa$groups = paste0(aa$celltype, '_', aa$condition)

p1 = DimPlot(aa, group.by = 'celltype', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'condition', cols = cols, label = TRUE, repel = TRUE)

p1 + p2

ggsave(filename = paste0(resDir, 
                         '/UMAP_merged.BL.CSD_subset_celltypes_aross.condition.time',
                        'selectedCelltypes.pdf'), 
       width = 18, height = 6)


##########################################
# test progeny 
##########################################
# PROGENy is a comprehensive resource containing a curated collection of pathways and their target genes, 
# with weights for each interaction
net <- get_progeny(organism = 'human', top = 500)
net

# Activity inference with Weighted Mean
# Extract the normalized log-transformed counts
mat <- as.matrix(aa@assays$RNA@data)
ss = colSums(mat)

rownames(mat) = convert_to_geneSymbol(rownames(mat), annot = annot)
mat = mat[grep('^AMEX', rownames(mat), invert = TRUE), ]

# Run wmean
tic()
acts <- run_wmean(mat=mat, net=net, .source='source', .target='target',
                  .mor='weight', times = 100, minsize = 5)
saveRDS(acts, file = paste0(RdataDir, '/res_run_wmean_progeny.rds'))

toc()

acts = readRDS(file = paste0(RdataDir, '/res_run_wmean_progeny.rds'))
acts


## Visualization 
# Extract norm_wmean and store it in pathwayswmean in data
aa[['pathwayswmean']] <- acts %>%
  filter(statistic == 'norm_wmean') %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = aa) <- "pathwayswmean"

# Scale the data
aa <- ScaleData(aa)
aa@assays$pathwayswmean@data <- aa@assays$pathwayswmean@scale.data

p1 <- DimPlot(aa, reduction = "umap", group.by = 'condition', cols = cols, label = TRUE, pt.size = 0.5) + 
  ggtitle('condition')
p2 <- (FeaturePlot(aa, features = c("Trail")) & 
         scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) +
  ggtitle('Trail activity')
p1 | p2

FeaturePlot(aa, features = rownames(aa)) & 
  scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')

ggsave(filename = paste0(resDir, '/progeny_14_signalingPathways.pdf'), 
       width = 30, height = 18)

pdfname = paste0(resDir, "/progyny_14_signalingPathways_split.by.sample.pdf")
pdf(pdfname, width=18, height = 6)
par(cex =1.0, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)

for(n in 1:nrow(aa))
{
  cat(n, '--', rownames(aa)[n], '\n')
  p1 = FeaturePlot(aa, features = rownames(aa)[n], split.by = 'sample') & 
    scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')
  plot(p1)
  
}

dev.off()


## Exploration 
# Extract activities from object as a long dataframe
Idents(aa) = aa$groups
DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE) + NoLegend()

df <- t(as.matrix(aa@assays$pathwayswmean@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(aa)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Transform to wide matrix
top_acts_mat <- df %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

vars = names(table(aa$celltype))
vis = c("BL_11dpa", "BL_8dpa", "BL_3_5dpa", "CSD_0dpa", "CSD_3dpa", "CSD_5dpa", "CSD_8dpa", "CSD_11dpa")
group_sorted = paste(rep(vars, each = length(vis)), vis, sep = "_")

top_acts_mat = top_acts_mat[match(group_sorted, rownames(top_acts_mat)), ]
gap_rows = c(8,16,24)

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 2, length.out=floor(palette_length/2)))


# Plot
pheatmap::pheatmap((top_acts_mat),
                   border_color = NA,
                   cluster_cols = TRUE,
                   cluster_rows = FALSE,
                   gaps_row = gap_rows,
                   color=my_color, 
                   breaks = my_breaks,
                   filename = paste0(resDir, '/progeny_14_signalingPathways_summary_celltype_time.pdf'),
                   width = 6, height = 10)


##########################################
# test TF activity with decoupleR (original code https://saezlab.github.io/decoupleR/articles/tf_sc.html)  
##########################################
library(OmnipathR)

net <- get_dorothea(organism="mouse", levels=c('A', 'B', 'C'))
net

# Extract the normalized log-transformed counts
mat <- as.matrix(aa@assays$RNA@data)

rownames(mat) = convert_to_geneSymbol(rownames(mat), annot = annot)
mat = mat[grep('^AMEX', rownames(mat), invert = TRUE), ]
rownames(mat) = firstup(rownames(mat))

# Run wmean
tic()
acts <- run_wmean(mat=mat, net=net, .source='source', .target='target',
                  .mor='mor', times = 100, minsize = 5)

saveRDS(acts, file = paste0(RdataDir, '/res_run_tfswmean.rds'))

toc()

acts = readRDS(file = paste0(RdataDir, '/res_run_tfswmean.rds'))
acts

# Extract norm_wmean and store it in tfswmean in pbmc
aa[['tfswmean']] <- acts %>%
  filter(statistic == 'norm_wmean') %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = aa) <- "tfswmean"

# Scale the data
aa <- ScaleData(aa)
aa@assays$tfswmean@data <- aa@assays$tfswmean@scale.data


# DimPlot(aa, reduction = "umap", group.by =  "Phase",  label = TRUE, pt.size = 0.5, label.size = 6) + 
#   ggtitle('time points')
# ggsave(filename = paste0(outDir,  '/cellcycle_umap.pdf'), 
#        width = 10, height = 6)
# 
# p0 = DimPlot(aa, reduction = "umap", group.by = 'condition',  label = TRUE, pt.size = 0.5, label.size = 6) + 
#   ggtitle('time points')
# 
# p1 <- DimPlot(aa, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 6) + 
#   NoLegend() + ggtitle('clusters')
# p0 | p1
# 
# ggsave(filename = paste0(outDir,  '/timePoints_clusters.pdf'), 
#        width = 12, height = 6)

DefaultAssay(object = aa) <- "tfswmean"
n_tfs <- 100
# Extract activities from object as a long dataframe
Idents(aa) = aa$clusters

df <- t(as.matrix(aa@assays$tfswmean@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(aa)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Get top tfs with more variable means across clusters
tfs <- df %>%
  group_by(source) %>%
  summarise(std = sd(mean)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)

# Subset long data frame to top tfs and transform to wide matrix
top_acts_mat <- df %>%
  filter(source %in% tfs) %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

# Plot
pheatmap::pheatmap(t(top_acts_mat), border_color = NA, color=my_color, breaks = my_breaks,
         filename = paste0(resDir, '/decoupleR_TF.activity_summary_clusters_top', n_tfs, '.pdf'), 
         width = 6, height = 16) 


Idents(aa) = aa$groups

df <- t(as.matrix(aa@assays$tfswmean@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(aa)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Get top tfs with more variable means across clusters
tfs <- df %>%
  group_by(source) %>%
  summarise(std = sd(mean)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)

# Subset long data frame to top tfs and transform to wide matrix
top_acts_mat <- df %>%
  filter(source %in% tfs) %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

vars = names(table(aa$celltype))
vis = c("BL_11dpa", "BL_8dpa", "BL_3_5dpa", "CSD_0dpa", "CSD_3dpa", "CSD_5dpa", "CSD_8dpa", "CSD_11dpa")
group_sorted = paste(rep(vars, each = length(vis)), vis, sep = "_")

top_acts_mat = top_acts_mat[match(group_sorted, rownames(top_acts_mat)), ]
gap_rows = c(8,16,24)


# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

# Plot
pheatmap::pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks,
                   gaps_row = gap_rows,
                   cluster_cols = TRUE,
                   cluster_rows = FALSE,
                   filename = paste0(resDir, '/decoupleR_TF.activity_summary_groups.pdf'), 
                   width = 20, height = 8) 

outDir = paste0(resDir, '/TF_activity_inferred')
if(!dir.exists(outDir)) dir.create(outDir)

DefaultAssay(object = aa) <- "RNA" 
ggs = convert_to_geneSymbol(rownames(aa), annot = annot)
#for (gene in c('Foxa2', 'Pax6', 'Sox2', 'Sox1'))
#for(gene in colnames(top_acts_mat))
DefaultAssay(object = aa) <- "tfswmean"
tfs = rownames(aa)

for(gene in tfs)
{
  cat('gene -- ', gene, '\n')
  # gene = 'Smad1'
  DefaultAssay(object = aa) <- "tfswmean"
  p2 <- (FeaturePlot(aa, features = gene) & 
           scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) + 
    ggtitle(paste0(gene, ' activity'))
  DefaultAssay(object = aa) <- "RNA" 
  
  feature = rownames(aa)[which(ggs == toupper(gene))]
  if(length(feature) == 1){
    p3 <- FeaturePlot(aa, features = feature) + ggtitle(paste0(gene, ' expression'))
    DefaultAssay(object = aa) <- "tfswmean"
    p2 | p3
    
    ggsave(filename = paste0(outDir, '/decoupleR_TF_activity_expression_', gene, '.pdf'), 
           width = 12, height = 6)
  }else{
    p2 
    ggsave(filename = paste0(outDir, '/decoupleR_TF_activity_expression_', gene, '.pdf'), 
           width = 8, height = 6)
  }
  
}

saveRDS(aa, file = paste0(RdataDir, 
                      'BL.CSD_merged_subset_CT_MAC_Neu_Epd_umap.rds_',
                      'savedTFactivities.decoupleR.rds'))

##########################################
# try ReactomeGSA for global pathway analysis 
# https://bioconductor.org/packages/release/bioc/vignettes/ReactomeGSA/inst/doc/analysing-scRNAseq.html
# 
##########################################
use_ReactomeGSA = FALSE
if(use_ReactomeGSA){
  require(ReactomeGSA)
  require(tictoc)
  
  cell_ids  = 'groups'
  Idents(aa) = aa$groups
  
  # if(cell_ids == 'condition')  {
  #   Idents(aa) = aa$condition
  # }else{
  #   Idents(aa) = aa$clusters
  # }
  # Extract the normalized log-transformed counts
  mat <- aa@assays$RNA@data
  rownames(mat) = convert_to_geneSymbol(rownames(mat), annot = annot)
  xx = CreateSeuratObject(counts = mat, assay = 'RNA', meta.data = aa@meta.data)
  rm(mat)
  
  Idents(xx) = xx$groups
  
  tic()
  gsva_result <- analyse_sc_clusters(object = xx, verbose = TRUE)
  
  toc()
  
  saveRDS(gsva_result, file = paste0(RdataDir, 
                            'BL.CSD_merged_subset_CT_MAC_Neu_Epd_umap.rds_',
                            'gsva_result.rds'))
  rm(xx)
  
  # The resulting object is a standard ReactomeAnalysisResult object.
  gsva_result = readRDS(file = paste0(RdataDir, 
                                      'BL.CSD_merged_subset_CT_MAC_Neu_Epd_umap.rds_',
                                      'gsva_result.rds'))
  
  gsva_result
  
  # pathways returns the pathway-level expression values per cell cluster
  pathway_expression <- pathways(gsva_result)
  
  # simplify the column names by removing the default dataset identifier
  colnames(pathway_expression) <- gsub("\\.Seurat", "", colnames(pathway_expression))
  pathway_expression[1:3, ]
  
  # A simple approach to find the most relevant pathways is to assess the maximum difference 
  # in expression for every pathway:
  
  # find the maximum differently expressed pathway
  max_difference <- do.call(rbind, apply(pathway_expression, 1, function(row) {
    values <- as.numeric(row[2:length(row)])
    return(data.frame(name = row[1], min = min(abs(values)), max = max(abs(values))))
  }))
  
  max_difference$diff <- max_difference$max - max_difference$min
  
  # sort based on the difference
  max_difference <- max_difference[order(max_difference$diff, decreasing = T), ]
  head(max_difference)
  
  plot(max_difference$max, max_difference$diff, cex = 0.5)
  abline(v = 0.2)
  abline(h = 0.02)
  sels = which(max_difference$max>=0.2 & max_difference$diff>0.1)
  
  cat(length(sels), ' pathways passing thresholds \n')
  
  #relevant_pathways <- c("R-HSA-983170", "R-HSA-388841", "R-HSA-2132295", "R-HSA-983705", "R-HSA-5690714")
  jj = intersect(grep('FGF|WNT|BMP|TGF|Retinoic Acid|Hippo|EGFR|NOTCH|GLI|SHH|Hh', max_difference$name),
                 grep('Signaling|signaling|Hippo|GLI|SHH|Hh|ERK', max_difference$name))
  
  max_difference_jj = max_difference[unique(c(jj)), ]
  max_difference = max_difference[unique(c(sels, jj)), ]
  
  gsva_mat = pathway_expression[match(max_difference$name, pathway_expression$Name), ]
  rownames(gsva_mat) = gsva_mat$Name
  gsva_mat = gsva_mat[, -c(1)]
  
  vars = names(table(aa$celltype))
  vis = c("BL_11dpa", "BL_8dpa", "BL_3_5dpa", "CSD_0dpa", "CSD_3dpa", "CSD_5dpa", "CSD_8dpa", "CSD_11dpa")
  group_sorted = paste(rep(vars, each = length(vis)), vis, sep = "_")
  
  gsva_mat = gsva_mat[, match(group_sorted, colnames(gsva_mat))]
  gap_cols = c(8,16,24)
  
  
  # Choose color palette
  palette_length = 100
  my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)
  
  pheatmap::pheatmap(gsva_mat, border_color = NA, color=my_color, 
           #breaks = my_breaks,
           scale = 'row',
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           gaps_col = gap_cols,
           fontsize = 7,
           filename = paste0(resDir, '/ReactomeGSA_global_pathwayActivity_', cell_ids, '_many.pdf'), 
           width = 12, height = 35) 
  
  # selected pathways
  gsva_mat = pathway_expression[match(max_difference_jj$name, pathway_expression$Name), ]
  rownames(gsva_mat) = gsva_mat$Name
  gsva_mat = gsva_mat[, -c(1)]
  
  palette_length = 100
  my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)
  
  gsva_mat = gsva_mat[, match(group_sorted, colnames(gsva_mat))]
  gap_cols = c(8,16,24)
  
  pheatmap::pheatmap(gsva_mat, border_color = NA, color=my_color, 
           #breaks = my_breaks,
           scale = 'row',
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           gaps_col = gap_cols,
           fontsize = 8,
           filename = paste0(resDir, '/ReactomeGSA_global_pathwayActivity_', cell_ids, '_selected.pdf'), 
           width = 12, height = 8) 
  
}

##########################################
# try different method for signaling pathways with decoupleR
# original code from https://www.sc-best-practices.org/conditions/gsea_pathway.html
##########################################
Test_DecoupleR_reactome = FALSE
if(Test_DecoupleR_reactome){
  #reactome = gmt_to_decoupler("c2.cp.reactome.v7.5.1.symbols.gmt")
  msigdb = get_resource("MSigDB")
  
  resource = 'kegg'
  outDir = paste0(resDir, '/decoupler_signaling.activity_singleCell_AUCell/')
  if(!dir.exists(outDir)) dir.create(outDir)
  
  # Get reactome pathways
  #reactome = msigdb %>% filter(collection == 'reactome_pathways')
  #reactome = msigdb %>% filter(collection == 'hallmark')
  reactome = msigdb %>% filter(collection == 'kegg_pathways')
  
  # Filter duplicates
  reactome = reactome %>% dplyr::select(geneset, genesymbol) %>% 
    distinct()
  
  # Filtering genesets to match behaviour of fgsea
  geneset_size = table(as.data.frame(reactome)[, 1])
  geneset_sel = names(geneset_size)[which((geneset_size > 15) & (geneset_size < 500))]
  #res_gsea <- run_fgsea(mat, network, .source='source', .target='target', nproc=1, minsize = 0)
  
  # run aucel for each cell
  mat <- as.matrix(aa@assays$RNA@data)
  
  rownames(mat) = convert_to_geneSymbol(rownames(mat), annot = annot)
  mat = mat[grep('^AMEX', rownames(mat), invert = TRUE), ]
  #rownames(mat) = firstup(rownames(mat))
  
  library(tictoc)
  
  tic() # run aucell for each cell
  aucel = run_aucell(mat = mat, network = reactome, 
                     .source="geneset",
                     .target="genesymbol",
                     minsize = 5
  )
  
  saveRDS(aucel, file = paste0(RdataDir, '/res_run_aucel.rds'))
  
  toc()
  
  aucel = readRDS(file = paste0(RdataDir, '/res_run_aucel.rds'))
  
  
  # Extract norm_wmean and store it in tfswmean in pbmc
  aa[['aucell']] <- aucel %>%
    filter(statistic == 'aucell') %>%
    pivot_wider(id_cols = 'source', names_from = 'condition',
                values_from = 'score') %>%
    column_to_rownames('source') %>%
    Seurat::CreateAssayObject(.)
  
  # Change assay
  DefaultAssay(object = aa) <- "aucell"
  
  # Scale the data
  aa <- ScaleData(aa)
  aa@assays$aucell@data <- aa@assays$aucell@scale.data
  
  Idents(aa) = aa$groups
  DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE)
  
  # FeaturePlot(aa, features = 'HALLMARK-ALLOGRAFT-REJECTION')
  
  df <- t(as.matrix(aa@assays$aucell@data)) %>%
    as.data.frame() %>%
    mutate(cluster = Idents(aa)) %>%
    pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
    group_by(cluster, source) %>%
    summarise(mean = mean(score))
  
  # Get top tfs with more variable means across clusters
  ntop = 120
  top_SP <- df %>%
    group_by(source) %>%
    summarise(std = sd(mean)) %>%
    arrange(-abs(std)) %>%
    #summarise(mean_mean = mean(mean)) %>%
    #arrange(-abs(mean_mean)) %>%
    head(ntop) %>%
    pull(source)
  
  # Subset long data frame to top tfs and transform to wide matrix
  top_aucell_mat <- df %>%
    filter(source %in% top_SP) %>%
    pivot_wider(id_cols = 'cluster', names_from = 'source',
                values_from = 'mean') %>%
    column_to_rownames('cluster') %>%
    as.matrix()
  
  vars = names(table(aa$celltype))
  vis = c("BL_11dpa", "BL_8dpa", "BL_3_5dpa", "CSD_0dpa", "CSD_3dpa", "CSD_5dpa", "CSD_8dpa", "CSD_11dpa")
  group_sorted = paste(rep(vars, each = length(vis)), vis, sep = "_")
  
  top_aucell_mat = top_aucell_mat[match(group_sorted, rownames(top_aucell_mat)), ]
  gap_cols = c(8,16,24)
  
  # Choose color palette
  palette_length = 100
  my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)
  
  my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
                 seq(0.05, 3, length.out=floor(palette_length/2)))
  
  # Plot
  pheatmap(t(top_aucell_mat), border_color = NA, 
           color=my_color, 
           breaks = my_breaks,
           cluster_cols = FALSE,
           gaps_col = gap_cols,
           filename = paste0(resDir, '/decoupleR_Aucell_', resource, '_summary_clusters.pdf'), 
           width = 15, height = 24) 
  
  
  for(sp in rownames(aa))
  {
    cat('signaling  -- ', sp, '\n')
    
    p4 = FeaturePlot(aa, features = sp) + ggtitle(sp)
    
    ggsave(filename = paste0(outDir, '/decoupleR_TF_activity_expression_', resource, '_', sp, '.pdf'), 
           width = 8, height = 6)
    
  }
  
}


########################################################
########################################################
# Section II : ligand-receptor anlaysis
# e.g. LIANA and NicheNet with defined subpopulation
# test LIANA and NicheNet using clusters 
########################################################
########################################################
aa = readRDS(file = paste0(RdataDir, '/BL.CSD_merged_subset_CT_MAC_Neu_Epd_day3_5_8_subtypes_umap.rds'))

## manually merge subclusters
aa$subtypes[which(aa$subtypes == 'epidermis_BL_early_2')] = 'epidermis_BL_early'
aa$subtypes[which(aa$subtypes == 'epidermis_BL_early_1')] = 'epidermis_BL_early'
aa$subtypes[which(aa$subtypes == 'CT_BL_early_2')] = 'CT_BL_early_1'

DimPlot(aa, group.by = 'subtypes')

aa$subtypes = factor(aa$subtypes, levels = sort(unique(aa$subtypes)))

aa = ScaleData(aa, features = rownames(aa))

Idents(aa) = aa$subtypes


## test LIANA
source("functions_ligandReceptor_analysis.R")
aa$celltypes = aa$subtypes


## test NicheNet
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


# 
# tfs = unique(c(tfs, "Zfp703"))
# tfs_sel = intersect(tfs, markers_saved)
# 
# pdfname = paste0(saveDir, "gene_examples_TFs.pdf")
# pdf(pdfname, width=24, height = 6)
# par(cex =0.7, mar = c(3,0.8,2,5)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
# 
# #for(n in 1:20)
# for(n in 1:length(tfs_sel))
# {
#   # n = 1
#   gg = sps_sels[n];
#   # gg = "Foxa2"
#   cat(n, '--', gg, '\n')
#   
#   p1 = FeaturePlot(aa, features = gg)
#   p2 = VlnPlot(aa, features = gg, group.by = 'condition', pt.size = 0.01)  + NoLegend()  
#   #scale_fill_manual(values=c("blue", 'green', "red",  "blue", 
#   #                           'red', 'green', 'blue', 'green', 'red', 'green', 'green'))
#   p3 = VlnPlot(aa, features = gg, group.by = 'clusters', pt.size = 0.01) + NoLegend()
#   plot(p1 | p2| p3)
#   
# }
# 
# dev.off()
# 
