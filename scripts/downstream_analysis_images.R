##########################################################################
##########################################################################
# Project: CSD cartilage closure 
# Script purpose:
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Oct 23 10:43:04 2025
##########################################################################
##########################################################################
outDir = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/image_analysis/axolotl_limb_CSD/",
              "results/downstream_ilastik/")
res = read.csv(file = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/image_analysis/axolotl_limb_CSD/",
                             "results/downstream_ilastik/",
                             "ilastik_segmentation_pct_closure.csv"), 
               header = TRUE, row.names = c(1))


manual = readxl::read_xlsx(path = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/image_analysis/axolotl_limb_CSD/",
                                "AA_LNP-aug2025/",
                                "LNPinj4w-gapMesurement180925.xlsx"),
                           sheet = 1, col_names = TRUE
                           )
manual = data.frame(manual[, c(1:9)])

manual$image = paste0(manual$sample, '-LNP', manual$LNP, '-limb', manual$limb)

mm = match(manual$image, res$image)
res$pct_manual = NA

res$pct_manual[mm] = manual$X.closed/100


## images that need double check
"675930-LNP3-limb2"
"675948-LNP2-limb1"
"675969-LNP3-limb1"
"675982-LNP1-limb1"
"676002-LNP2-limb1"

res[which(res$image == "675930-LNP3-limb2"), ]

res[which(res$image == "675948-LNP2-limb1"), ]

res[which(res$image == "675969-LNP3-limb1"), ]

res[which(res$image == "675982-LNP1-limb1"), ]

res[which(res$image == "676002-LNP2-limb1"), ]


##########################################
# convert the CSD gap into mm 
##########################################
res$csd_gap = res$bone_gap + res$left_cuttingPlan + res$right_cuttingPlan
#res$csd_gap = res$csd_gap * 1.5287730727470145*0.001 * 2.2
res$csd_gap_mm = res$csd_gap * 4.317044*0.001

write.csv2(res, file = paste0(outDir, 'ilastik_segmentation_pct_closure_boneGap_mm.csv'), row.names = TRUE, 
           quote = FALSE)




##########################################
# for 2nd run, merge the manual cropping and automatic running 
##########################################
outDir = paste0("/groups/tanaka/People/current/jiwang/projects/image_analysis/",
                #"axolotl_limb_CSD/export_aug2026/4w_LNPtreated_07082026/results")
                #"axolotl_limb_CSD/export_aug2026/4w_LNP_AA_19052026/results")
                "axolotl_limb_CSD/export_aug2026/4w_LNPtreated_gapMeasure_100926/results")

res = read.csv2(file = paste0(outDir, "/ilastik_segmentation_pct_closure_addedManualCropping_manualCorrect.csv"), 
               header = TRUE, row.names = c(1))

files = list.files(path = outDir, 
                   pattern = '*.csv', full.names = TRUE)

ff = basename(files)
files = files[which(basename(files) != 'ilastik_segmentation_pct_closure_addedManualCropping_manualCorrect.csv')]

for(n in 1:length(files))
{
  # n = 14
  x = read.csv(file = files[n], header = TRUE, row.names = c(1))
  x$image = gsub('_Probabilities.tif','', x$image)
  if(nrow(x) > 1) x = x[nrow(x), ]
  
  kk = which(res$image == x$image)
  if(length(kk) == 1 & nrow(x) == 1){
    cat('-- Correct result of image : ', x$image, '--\n')
    res[kk, ] = x
  }else{
    cat('Error in image -- ', n, '--', x$image, '\n')
  }
}

write.csv2(res, file = paste0(outDir, 
                              '/ilastik_segmentation_pct_closure_addedManualCropping_manualCorrect_v2.csv'), 
           row.names = TRUE, quote = FALSE)



########################################################
########################################################
# Section II: process and analyze the HCR signal
# 
########################################################
########################################################
rm(list = ls())

library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library("viridis")
library(sf)

project_to_line <- function(points, A, B) {
  points <- as.matrix(points)
  A <- as.numeric(A)
  B <- as.numeric(B)
  
  direction <- B - A
  denominator <- sum(direction^2)
  
  if (denominator == 0) {
    stop("A and B must be different points.")
  }
  
  # Subtract A from every point
  relative <- sweep(points, 2, A, FUN = "-")
  
  # Position of each projection along line A -> B
  t <- as.vector(relative %*% direction / denominator)
  
  # Projected x-y coordinates
  projected <- sweep(
    outer(t, direction),
    2,
    A,
    FUN = "+"
  )
  
  colnames(projected) <- c("x_projected", "y_projected")
  
  data.frame(
    points,
    t = t,
    projected
  )
}


inputDir = paste0("/groups/tanaka/People/current/jiwang/projects/image_analysis/axolotl_limb_CSD/", 
                  "HCR_Kazald1/images_processed/cellmasks_cellpose_imageJmacro/features_collected")

outDir = paste0("/groups/tanaka/People/current/jiwang/projects/image_analysis/axolotl_limb_CSD/HCR_Kazald1/",
                "results/imageJcellpose_mask_v2/")

if(!dir.exists(outDir)){
  dir.create(outDir)
}

files = list.files(path = inputDir, 
                   pattern = '*.csv', full.names = TRUE)

if(length(grep('removedLabelIDs', files)) > 0){
  files = files[grep('removedLabelIDs', files, invert = TRUE)] 
}

ff = basename(files)
ff = gsub('_parameterCollection_HCRintensity.csv', '', ff)
ff = sapply(ff, function(x){unlist(strsplit(x, '-'))[1]})


crops = data.frame(readxl::read_xlsx(path = paste0("/groups/tanaka/People/current/jiwang/projects/image_analysis/axolotl_limb_CSD/HCR_Kazald1/", 
                                                   'coordinates_HCR_230926.xlsx'), sheet = 1,
                          col_names = TRUE))

names = crops$file
names = sapply(names, function(x){x = unlist(strsplit(x, '_')); paste0(x[c(1:(length(x)-2))], 
                                                                       collapse = '_')})

names[1] = "11d_LNP3_20x_limb1"
names[5] ="11d_LNP4_limb3"
names[which(names == "11d_LNP4_limb3")] = "11d_LNP4_20x_limb3"
names[which(names == "11dpi_LNP1_20x")] = "11dpi_LNP1_20x_an1"

mm = match(ff, names)

ff[which(is.na(match(ff, names)))]
names[which(is.na(match(names, ff)))]

crops$oldFile = crops$file
crops$file[mm] = ff

area_cutoff = 3.8
prrx_cutoff = 3.0
dapi_cutoff = 3.0

keep = c()

for(n in 1:length(files))
#for(n in 1:2)
{
  # n = 14
  
  file_name = gsub('_parameterCollection_HCRintensity.csv', '', basename(files[n]))
  file_name = unlist(strsplit(file_name, '-'))[1]
  
  cat(n, ' -- ', file_name, '\n')
  
  x = read.csv(file = files[n], header = TRUE, row.names = c(1))
  
  cat(' nb of cells detected before filtering :', nrow(x), '\n')
  
  x$area_mask = log10(x$area_mask)
  x$intensity_hcr = log10(x$intensity_mean_hcr) 
  x$intensity_prrx = log10(x$intensity_mean_cellmarker)
  x$intensity_gfp = log10(x$intensity_mean_lnpGFP)
  x_lims = max(x$centroid.1_mask)
  y_lims = max(x$centroid.0_mask)
  
  kk = which(crops$file == file_name)
  cat(kk, 'th line of file in crop coordinates \n')
  
  Discard_blood_cells = FALSE
  if(n == 14){
    blood = read.csv(file = paste0(inputDir, '/11dpi_LNP2_an3_20x-AiryStitch-14_removedLabelIDs.csv'),
                     header = FALSE)
    
    mm = match(x$label_mask, blood$V1)
    
    x = x[which(is.na(mm)), ]
    
  }
  
  pdfname = paste0(outDir, file_name, 'cell_selection_HCR_quanty_projection.pdf')
  pdf(pdfname, width=12, height = 8)
  
  #plot(x$eccentricity_mask, x$solidity_mask)
  plot(x$area_mask, x$eccentricity_mask, cex = 0.5)
  #area_cutoff = crops$area_mask_cutoff[kk]
  abline(v = area_cutoff, col = 'red')
  abline(v = 4, col = 'red')
  
  # filter cells with size
  #sels = which(x$area_mask > area_cutoff)
  
  sels = which(x$area_mask > area_cutoff)
  #sels = which(x$area_mask > 4 & x$eccentricity_mask < 0.8)
  
  cat(' nb of cells after size filtering :', length(sels), '\n')
  x = x[sels, ]
  
  
  p1 = ggplot() +
    geom_point(data = data.frame(x), mapping = aes(x = centroid.1_mask, y = centroid.0_mask, 
                                       fill = intensity_hcr), 
               shape = 21, size = 2) +
    #scale_fill_viridis_c() + 
    scale_fill_viridis_c(option = "magma") + 
    theme_bw() + 
    ggtitle('cells after size filtering')
  
  plot(p1)
  
  # p2 = ggplot() +
  #   geom_point(data = data.frame(x), mapping = aes(x = centroid.1_mask, y = centroid.0_mask, 
  #                                                  fill = intensity_gfp), 
  #              shape = 21, size = 2) +
  #   #scale_fill_viridis_c() + 
  #   scale_fill_viridis_c(option = "magma") + 
  #   theme_bw() + 
  #   ggtitle('cells after size filtering')
  # 
  # p1 + p2
  # 
  
  plot(x$intensity_prrx, log10(x$intensity_mean_dapi), cex = 0.2)
  #prrx_cutoff = crops$Prrx_intensity_cutoff[kk]
  #dapi_cutoff = crops$Dapi_intensity_cutoff[kk]
  abline(v = prrx_cutoff, col = 'red')
  abline(h = dapi_cutoff, col = 'red')
  
  # plot(x$intensity_gfp, x$intensity_prrx,  cex = 0.2)
  # #prrx_cutoff = crops$Prrx_intensity_cutoff[kk]
  # #dapi_cutoff = crops$Dapi_intensity_cutoff[kk]
  # abline(h = prrx_cutoff, col = 'red')
  # abline(h = dapi_cutoff, col = 'red')
  
  sels = which(x$intensity_prrx > prrx_cutoff & log10(x$intensity_mean_dapi) > dapi_cutoff)
  # sels = which(x$intensity_prrx > prrx_cutoff & 
  #                log10(x$intensity_mean_dapi) > 3.2 &
  #                x$intensity_gfp > 3.2)
  cat(' nb of cells detected after Prxx and dapi filtering  :', length(sels), '\n')
  
  x = x[sels, ]
  
  
  coord_xy = ceiling(as.numeric(crops[kk, c(2:9)]))
  coord_x = coord_xy[c(1, 3, 5, 7)]
  coord_y = coord_xy[c(2, 4, 6, 8)]
  
  cat(coord_x, '\n')
  cat(coord_y, '\n')
  
  # coord_x[2] = 4000
  # coord_y[2] = 7500
  # coord_y[4] = 10000
  
  # select points inside a rectangle defined by four corners from GPT
 
  #resize_factor = 1
  # Four rectangle corners; order does not matter
  corners <- data.frame(
    x = as.numeric(coord_x[c(1,2, 4, 3)]),
    y = as.numeric(coord_y[c(1,2, 4, 3)])
  )
  
  # corners <- data.frame(
  #   x = as.numeric(c(coord_x[c(1,2)], 30000, coord_x[c(4, 3)])),
  #   y = as.numeric(c(coord_y[c(1,2)], 8000, coord_y[c(4, 3)]))
  # )
  # 
  
  p2 = ggplot() +
    geom_point(data = data.frame(x), mapping = aes(x = centroid.1_mask, y = centroid.0_mask, 
                                                   fill = intensity_hcr), 
               shape = 21, size = 2) +
    #scale_fill_viridis_c() + 
    xlim(0, x_lims) + 
    ylim(0, y_lims) +
    scale_fill_viridis_c(option = "magma") + 
    theme_bw() +
    ggtitle('cells after Prxx and dapi filtering') +
    geom_polygon(data = corners, aes(x = x, y = y),
                 fill = "orange", alpha = 0.3, color = "magenta")
  
  plot(p2)
  
  # Order corners around their centroid
  center_x <- mean(corners$x)
  center_y <- mean(corners$y)
  
  angles <- atan2(
    corners$y - center_y,
    corners$x - center_x
  )
  
  ordered_corners <- corners[order(angles), ]
  
  # Close the polygon by repeating the first corner
  polygon_coordinates <- rbind(
    as.matrix(ordered_corners[, c("x", "y")]),
    as.matrix(ordered_corners[1, c("x", "y")])
  )
  
  # Create an sf polygon
  rectangle <- st_sfc(st_polygon(list(polygon_coordinates)))
  
  # Convert points to sf
  points <- data.frame(
    id = x$label_mask,
    x = x$centroid.1_mask,
    y = x$centroid.0_mask
  )
  points_sf <- st_as_sf(
    points,
    coords = c("x", "y"),
    remove = FALSE
  )
  
  inside <- lengths(st_intersects(points_sf, rectangle)) > 0
  
  # Select points
  selected_points <- points[inside, ]
  
  index_fibs = match(selected_points$id, x$label_mask)
  
  x = x[index_fibs, ]
  
  cat(' nb of cells detected within the cropping rectangle  :', nrow(x), '\n')
  
  
  p3 = ggplot() +
    geom_point(data = data.frame(x), mapping = aes(x = centroid.1_mask, y = centroid.0_mask, 
                                                   fill = intensity_hcr), 
               shape = 21, size = 2) +
    #scale_fill_viridis_c() + 
    xlim(0, x_lims) + 
    ylim(0, y_lims) +
    scale_fill_viridis_c(option = "magma") + 
    theme_bw() +
    ggtitle('cells within the cropping region') +
    geom_polygon(data = corners, aes(x = x, y = y),
                 fill = "orange", alpha = 0.3, color = "magenta")
  
  plot(p3)
  
  ## project all points to the line defined by the two centers
  points <- data.frame(
    x = x$centroid.1_mask,
    y = x$centroid.0_mask
  )
  
  A <- c(mean(as.numeric(coord_x[1:2])), mean(as.numeric(coord_y[1:2])))
  B <- c(mean(as.numeric(coord_x[3:4])), mean(as.numeric(coord_y[3:4])))
  
  result <- project_to_line(points, A, B)
  #result
  
  x$project = result$t
  #plot(x$project, x$intensity_hcr)
  
  p4 = ggplot(data.frame(x), aes(project, intensity_hcr)) +
    geom_point() +
    geom_smooth(span = 0.6, formula = 'y ~ x', method = 'loess') +
    theme_classic() +
    ggtitle('all cells within regions')
  
  plot(p4)
  
  xx = x[which(x$project > 0.1 & x$project < 0.9), ]
  
  ggplot(data.frame(xx), aes(project, intensity_hcr)) +
    geom_point() +
    geom_smooth(span = 0.6, formula = 'y ~ x', method = 'loess') +
    theme_classic() + 
    xlim(0, 1) + 
    ggtitle('cells without borders')
  
  
  dev.off()
  
  x$sample = file_name
  
  keep = rbind(keep, x[, match(c('sample', 'project', 'intensity_hcr'), colnames(x))])
  
}

write.csv2(keep, file = paste0(outDir, 'image_hcrIntensity_gapProjection_all.csv'), 
          row.names = FALSE)


# xx = keep
# 
# keep = read.csv2(file = paste0(outDir, 'image_hcrIntensity_gapProjection.csv'), header = TRUE)
# 
# keep = keep[which(keep$sample != "11dpi_LNP2_an2_20x"), ]
# 
# keep = rbind(keep, xx)


xx = keep[grep('_LNP2_|_LNP3_|_LNP4_', keep$sample), ]
xx$LNP = sapply(xx$sample, function(x){unlist(strsplit(x, '_'))[2]})

ggplot(data.frame(xx), aes(x = project, y = intensity_hcr, group = sample, colour = LNP)) +
  #geom_point(size = 0.01)  + 
  theme_classic() +  
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  labs( x = 'Relative Position in CSD gap', y = 'HCR intensity' ) +
  xlim(0, 1) +
  geom_smooth(span = 0.6, formula = 'y ~ x', method = 'loess', se = FALSE) 
  
ggsave(filename = paste0(outDir, 'hcrIntensity_gapProjection_LNP234.pdf'), 
       width = 8, height = 4)


xx = keep
xx$LNP = sapply(xx$sample, function(x){unlist(strsplit(x, '_'))[2]})

ggplot(data.frame(xx), aes(x = project, y = intensity_hcr, group = sample, colour = LNP)) +
  #geom_point(size = 0.01)  + 
  theme_classic() +  
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  labs( x = 'Relative Position in CSD gap', y = 'HCR intensity' ) +
  xlim(0, 1) +
  geom_smooth(span = 0.6, formula = 'y ~ x', method = 'loess', se = FALSE) 

ggsave(filename = paste0(outDir, 'hcrIntensity_gapProjection_LNPall.pdf'), 
       width = 8, height = 4)




########################################################
########################################################
# Section III: some statistical test
# 
########################################################
########################################################
rm(list = ls())
setwd("/groups/tanaka/People/current/jiwang/projects/bone_healing_CSD/scripts")

library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)

##########################################
# first table test
##########################################
res = read.csv2(file = paste0("../tables/gap-closure_4w_AA_ALL_new-old_quantified_ilastik.csv"),
               header = TRUE)

df_long <- pivot_longer(
  res,
  cols = everything(),             # Convert all columns except gene
  names_to = "condition",   # Column names become conditions
  values_to = "proportion",  # Values become expression
  values_drop_na = TRUE
)

ggplot(df_long, aes(x = condition, y = proportion)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
  theme_classic() +
  labs(x = NULL, y = "Value")

#xx = df_long[which(df_long$condition == "Wnt3a.Gfp"| df_long$condition == 'Wnt3a.Wnt5a.Gfp'), ]

cc = unique(df_long$condition)

compares = data.frame(c1 = c(rep("Wnt3a.Lef1.Gfp", 4),   "Wnt3a.Gfp"),
                      c2 = c("Gfp", "Wnt5a.Gfp", "Tgfb.Inhba.Gfp", "Wnt3a.Gfp", "Gfp"), 
                      pval = NA)

for(n in 1:nrow(compares))
{
  # n = 1
  test = wilcox.test(x = df_long$proportion[which(df_long$condition ==compares$c1[n])], 
              y = df_long$proportion[which(df_long$condition ==compares$c2[n])], 
              exact = FALSE, 
              alternative = c('greater'))
 
  compares$pval[n] = test$p.value 
}

x = df_long$proportion[which(df_long$condition =="Wnt3a.Lef1.Gfp")]
y = df_long$proportion[which(df_long$condition =="Wnt3a.Gfp")]

ks.test(x = x, y = y, alternative = c('less'))
kSamples::ad.test(x, y)
twosamples::cvm_test(x, y)

# for some reason, the coin package makes you make the 'group' column/variable a factor
df = df_long[which(df_long$condition == "Wnt3a.Lef1.Gfp" | df_long$condition == "Wnt3a.Gfp"), ]
df$condition = factor(df$condition, levels = c("Wnt3a.Lef1.Gfp", "Wnt3a.Gfp"))

library(coin)
independence_test(proportion ~ condition, data = df, alternative = "greater")

write.csv2(compares, file = "../tables/gap-closure_4w_AA_ALL_new-old_quantified_ilastik_wilcoxTest.csv")

##########################################
# second test
##########################################
res = read.csv2(file = paste0("../tables/2wEdU_of_DAPI_LNPtreated.csv"),
                header = TRUE)

df_long <- pivot_longer(
  res,
  cols = everything(),             # Convert all columns except gene
  names_to = "condition",   # Column names become conditions
  values_to = "proportion",  # Values become expression
  values_drop_na = TRUE
)

ggplot(df_long, aes(x = condition, y = proportion)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
  theme_classic() +
  labs(x = NULL, y = "Value")

#xx = df_long[which(df_long$condition == "Wnt3a.Gfp"| df_long$condition == 'Wnt3a.Wnt5a.Gfp'), ]

cc = unique(df_long$condition)

compares = data.frame(c1 = c(rep("blastema", 8),  rep("X.W3aGFP", 8), rep("X.W3aL1GFP", 8)), 
                      c2 = c(cc[which(cc != "blastema")], cc[which(cc != "X.W3aGFP")], 
                             cc[which(cc != "X.W3aL1GFP")]), 
                      pval = NA)

for(n in 1:nrow(compares))
{
  # n = 1
  test = wilcox.test(x = df_long$proportion[which(df_long$condition ==compares$c1[n])], 
                     y = df_long$proportion[which(df_long$condition ==compares$c2[n])], 
                     exact = FALSE, 
                     alternative = c("two.sided"))
  
  compares$pval[n] = test$p.value 
}

ks.test(x = df_long$proportion[which(df_long$condition =="Wnt3a.Lef1.Gfp")], 
        y = df_long$proportion[which(df_long$condition =="Wnt3a.Gfp")]
        #alternative = c('greater')
)

write.csv2(compares, file = "../tables/2wEdU_of_DAPI_LNPtreated_wilcoxTest.csv")


##########################################
# third test 
##########################################
res = readxl::read_xlsx(path = paste0("../tables/11092026EdUofCherry_v2.xlsx"),
                col_names = TRUE, sheet = 1)

res = res[which(res$proportions != ' ' & !is.na(res$proportions)), ]

res = res[which(res$sample != 'intact'), ]

res$sample = gsub(' ', '_', res$sample)
res$sample = gsub('blastema', 'bl', res$sample)
res$sample = gsub('csd', 'CSD_', res$sample)

ggplot(res, aes(x=as.factor(proportions), fill=as.factor(sample) )) +
  geom_bar() +
  #geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
  theme_classic() +
  labs(x = NULL, y = "Value")

#xx = df_long[which(df_long$condition == "Wnt3a.Gfp"| df_long$condition == 'Wnt3a.Wnt5a.Gfp'), ]

res$condition = factor(res$sample)
cc = unique(res$condition)

df_long = res
df_long$proportion = df_long$proportions

compares = data.frame(c1 = c("bl_1w", "bl_2w", "bl_3w"), 
                      c2 = c("CSD_1w", "CSD_2w", "CSD_3w"), 
                      pval = NA)

for(n in 1:nrow(compares))
{
  # n = 1
  test = wilcox.test(x = df_long$proportion[which(df_long$condition ==compares$c1[n])], 
                     y = df_long$proportion[which(df_long$condition ==compares$c2[n])], 
                     exact = FALSE, 
                     alternative = c("two.sided"))
  
  compares$pval[n] = test$p.value 
}

ks.test(x = df_long$proportion[which(df_long$condition =="Wnt3a.Lef1.Gfp")], 
        y = df_long$proportion[which(df_long$condition =="Wnt3a.Gfp")]
        #alternative = c('greater')
)

write.csv2(compares, file = "../tables/11092026EdUofCherry_wilcoxTest.csv")

