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
                "axolotl_limb_CSD/export_aug2026/4w_LNPtreated_07082026/results")
                #"axolotl_limb_CSD/export_aug2026/4w_LNP_AA_19052026/results")

res = read.csv(file = paste0(outDir, "/ilastik_segmentation_pct_closure.csv"), 
               header = TRUE, row.names = c(1))

files = list.files(path = outDir, 
                   pattern = '*.csv', full.names = TRUE)
ff = basename(files)
files = files[which(basename(files) != 'ilastik_segmentation_pct_closure.csv')]

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

write.csv2(res, file = paste0(outDir, '/ilastik_segmentation_pct_closure_addedManualCropping.csv'), 
           row.names = TRUE, quote = FALSE)

