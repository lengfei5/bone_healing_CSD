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
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library("viridis")

outDir = "/groups/tanaka/People/current/jiwang/projects/image_analysis/axolotl_limb_CSD/HCR_Kazald1/"

files = list.files(path = paste0(outDir, 'results'), 
                   pattern = '*.csv', full.names = TRUE)

ff = basename(files)

crops = read.csv2(file = paste0(outDir, 'HCR_crop_coordinates.csv'), header = TRUE)

for(n in 1:length(files))
{
  # n = 1
  cat(n, ' -- ', basename(files[n]), '\n')
  
  file_name = gsub('_parameterCollection_HCRintensity.csv', '', basename(files[n]))
  
  x = read.csv(file = files[n], header = TRUE, row.names = c(1))
  
  cat(' nb of cells detected before filtering :', nrow(x), '\n')
  
  pdfname = paste0(outDir, '/results/', file_name, 'cell_selection_HCR_quanty_projection.pdf')
  pdf(pdfname, width=12, height = 8)
  
  
  x$area_mask = log10(x$area_mask)
  
  plot(x$eccentricity_mask, x$solidity_mask)
  plot(log10(x$area_mask), x$eccentricity_mask)
  abline(v = 3, col = 'red')
  
  # filter cells with size
  sels = which(x$area_mask > 10^3)
  cat(' nb of cells after size filtering :', length(sels), '\n')
  x = x[sels, ]
  
  x$intensity_hcr = log10(x$intensity_mean_hcr) 
  x$intensity_prrx = log10(x$intensity_mean_cellmarker)
  
  ggplot() +
    geom_point(data = data.frame(x), mapping = aes(x = centroid.1_mask, y = centroid.0_mask, 
                                       fill = intensity_hcr), 
               shape = 21, size = 2) +
    #scale_fill_viridis_c() + 
    scale_fill_viridis_c(option = "magma") + 
    theme_bw() + 
    ggtitle('cells after size filtering')
  
  
  x_lims = max(x$centroid.1_mask)
  y_lims = max(x$centroid.0_mask)
  
  plot(x$intensity_prrx, log10(x$intensity_mean_dapi), cex = 0.2)
  prrx_cutoff = 3.0
  dapi_cutoff = 3.0
  abline(v = prrx_cutoff, col = 'red')
  abline(h = dapi_cutoff, col = 'red')
  
  sels = which(x$intensity_prrx > prrx_cutoff & log10(x$intensity_mean_dapi) > dapi_cutoff)
  cat(' nb of cells detected after Prxx and dapi filtering  :', length(sels), '\n')
  
  x = x[sels, ]
  
  
  kk = grep(gsub("_parameterCollection_HCRintensity.csv", '', basename(files[n])), crops$file)
  coord_xy = ceiling(as.numeric(crops[kk, c(2:ncol(crops))])/0.061810168997668995)
  coord_x = coord_xy[c(1, 3, 5, 7)]
  coord_y = coord_xy[c(2, 4, 6, 8)]
  cat(coord_x)
  cat(coord_y, '\n')
  
  
  # select points inside a rectangle defined by four corners from GPT
  library(sf)
  
  # Four rectangle corners; order does not matter
  corners <- data.frame(
    x = as.numeric(coord_x[c(1,2,4,3)]*3),
    y = as.numeric(coord_y[c(1,2,4,3)]*3)
  )
  
  ggplot() +
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
                 fill = "orange", alpha = 0.3, color = "black")
  
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
  
  
  ggplot() +
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
  
  
  ## project all points to the line defined by the two centers
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
  
  ggplot(data.frame(x), aes(project, intensity_hcr)) +
    geom_point() +
    geom_smooth(span = 0.5, formula = 'y ~ x') +
    theme_classic() +
    ggtitle('all cells within regions')
  
  
  xx = x[which(x$project > 0.1 & x$project < 0.9), ]
  ggplot(data.frame(xx), aes(project, intensity_hcr)) +
    geom_point() +
    geom_smooth(span = 0.5, formula = 'y ~ x') +
    theme_classic() +
    ggtitle('cells without borders')
  
  dev.off()
  
}



