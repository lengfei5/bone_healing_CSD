run_LIANA_defined_celltype = function(subref,
                                      celltypes,
                                      receiver_cells = NULL,
                                      ntop = 100, 
                                      additionalLabel = '_fixedCelltypes')
{
  require(liana)
  require(Seurat)
  require(scater)
  require(scran)
  # subref = aa; additionalLabel = '_fixedCelltypes'; celltypes = unique(aa$celltypes)
  
  system(paste0('mkdir -p ', outDir))
  # source('functions_scRNAseq.R') 
  
  sce <- as.SingleCellExperiment(subref)
  colLabels(sce) = as.factor(sce$celltypes)
  rownames(sce) = toupper(rownames(sce))
  
  ave.counts <- calculateAverage(sce, assay.type = "counts")
  
  #hist(log10(ave.counts), breaks=100, main="", col="grey80",
  #     xlab=expression(Log[10]~"average count"))
  
  num.cells <- nexprs(sce, byrow=TRUE)
  #smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells",
  #              xlab=expression(Log[10]~"average count"))
  
  # detected in >= 5 cells, ave.counts >=5 but not too high
  genes.to.keep <- num.cells > 5 & ave.counts >= 10^-4  & ave.counts <10^3  
  summary(genes.to.keep)
  
  sce <- sce[genes.to.keep, ]
  
  ## run the liana wrap function by specifying resource and methods
  # Resource currently included in OmniPathR (and hence `liana`) include:
  show_resources()
  # Resource currently included in OmniPathR (and hence `liana`) include:
  show_methods()
  
  liana_test <- liana_wrap(sce,  
                           method = c("natmi", "connectome", "logfc", "sca", "cytotalk"  
                                      #'cellphonedb'
                           ),
                           resource = c("Consensus", 'CellPhoneDB', "OmniPath", "LRdb",
                                        "CellChatDB",  "CellTalkDB"), 
                           assay.type = "logcounts", 
                           idents_col = 'celltypes')
  
  # Liana returns a list of results, each element of which corresponds to a method
  liana_test %>% glimpse
  
  # We can aggregate these results into a tibble with consensus ranks
  saveRDS(liana_test, file = paste0(outDir, '/res_lianaTest_Consensus', 
                                    additionalLabel, '.rds'))
  
  return(liana_test)
  
  
  Test = FALSE
  if(Test){
    #if(is.na(receiver_cells)){ # loop over all cell type candidates
    #  celltypes = as.character(celltypes)
    #  receiver_cells = as.character(celltypes)
    #}
    
    ntop = 200
    manually_specifying_sender_receiver = FALSE
    if(manually_specifying_sender_receiver){
      
      
      
      for(ntop in c(100, 200, 500)){
        
        # ntop = 100
        liana_test %>%
          liana_dotplot(source_groups = sender_cells,
                        target_groups = receiver_cells,
                        ntop = ntop)
        ggsave(filename = paste0(outDir, '/liana_LR_prediction_BL_early', 
                                 additionalLabel, 
                                 #'_receiverCells.', receiver_cells[m], 
                                 '_ntop.', ntop, '.pdf'), 
               width = 3*length(sender_cells)*length(receiver_cells), 
               height = 0.25*ntop, limitsize = FALSE)
        
      }
      
      pdfname = paste0(outDir, '/liana_BL_celltype_communication_freqHeatmap', 
                       additionalLabel, '.pdf')
      pdf(pdfname, width=20, height = 8)
      
      liana_trunc <- liana_test %>% 
        filter(source %in% sender_cells) %>%
        filter(target %in% receiver_cells) %>%
        # only keep interactions concordant between methods
        filter(aggregate_rank <= 0.01) # this can be FDR-corr if n is too high
      
      heat_freq(liana_trunc)
      
      dev.off()
      
      sender_cells = c("CT_CSD_early_1", "mac_CSD_early",  "epidermis_BL.CSD_early",  "neu_CSD_early")
      receiver_cells = sender_cells
      
      for(ntop in c(100, 200, 500)){
        
        # ntop = 100
        liana_test %>%
          liana_dotplot(source_groups = sender_cells,
                        target_groups = receiver_cells,
                        ntop = ntop)
        ggsave(filename = paste0(outDir, '/liana_LR_prediction_CSD_early', 
                                 additionalLabel, 
                                 #'_receiverCells.', receiver_cells[m], 
                                 '_ntop.', ntop, '.pdf'), 
               width = 3*length(sender_cells)*length(receiver_cells), 
               height = 0.25*ntop, limitsize = FALSE)
        
      }
      
      pdfname = paste0(outDir, '/liana_CSDcelltype_communication_freqHeatmap', 
                       additionalLabel, '.pdf')
      pdf(pdfname, width=20, height = 8)
      
      liana_trunc <- liana_test %>% 
        filter(source %in% sender_cells) %>%
        filter(target %in% receiver_cells) %>%
        # only keep interactions concordant between methods
        filter(aggregate_rank <= 0.01) # this can be FDR-corr if n is too high
      
      heat_freq(liana_trunc)
      
      dev.off()
      
      
    }
    
  }
    
  
  
}

########################################################
########################################################
# Section : # NichetNet analysis using differential niche analysis
# 
########################################################
########################################################
Run_NicheNet_Seurat = function()
{
  library(pryr) # monitor the memory usage
  require(ggplot2)
  library(nichenetr)
  library(Seurat) # please update to Seurat V4
  library(tidyverse)
  library(circlize)
  library(RColorBrewer)
  require(scran)
  require(scater)
  library(nichenetr)
  library(Seurat)
  library(tidyverse)
  library(knitr)
  library(tidyr)
  library(reshape2)
  library(gridExtra)
  library(ComplexHeatmap)
  library(clusterProfiler)
  library(plotly)
  library(ggalluvial)
  library(stringr)
  
  options(stringsAsFactors = F)
  
  dataPath_nichenet = '/users/jingkui.wang/projects/heart_regeneration/data/NicheNet/'
  
  system(paste0('mkdir -p ', outDir))
  
  mem_used()
  
  ##########################################
  # load processed scRNA-seq and visium data
  ##########################################
  refs = aa
  
  refs$celltypes = refs$subtypes
  
  table(refs$celltypes)
  Idents(refs) = as.factor(refs$celltypes)
  
  ntop = c(50, 100, 200)
  
  ########################################################
  ########################################################
  # Section : test default NicheNet with Seurat object 
  # 
  ########################################################
  ########################################################
  
  ##########################################
  # step 0):  load the Nichenet data
  ##########################################
  # NicheNet’s ligand-target prior model
  ligand_target_matrix = readRDS(paste0(dataPath_nichenet,  "ligand_target_matrix.rds"))
  ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
  
  # ligand-receptor network, Putative ligand-receptor links from NicheNet
  lr_network = readRDS(paste0(dataPath_nichenet, "lr_network.rds"))
  lr_network = lr_network %>% mutate(bonafide = ! database %in% c("ppi_prediction","ppi_prediction_go"))
  lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% 
    distinct(ligand, receptor, bonafide)
  
  # If wanted, users can remove ligand-receptor interactions that were predicted based on 
  # protein-protein interactions and 
  # only keep ligand-receptor interactions that are described in curated databases.
  # lr_network = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
  head(lr_network)
  
  # ## weighted integrated networks 
  weighted_networks = readRDS(paste0(dataPath_nichenet,  "weighted_networks.rds"))
  head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
  head(weighted_networks$gr) # interactions and their weights in the gene regulatory network
  
  weighted_networks_lr = weighted_networks$lr_sig %>% 
    dplyr::rename(ligand = from, receptor = to) %>%
    inner_join(lr_network %>% 
                 distinct(ligand,receptor), by = c("ligand","receptor"))
  
  
  ## specificy manually the receiver in two niches, no time points
  timepoint = 'early'
  celltypes = c('CT_BL', 'epidermis_BL', 'mac_BL', 'neu_BL')  
  celltypes_ctl = c('CT_CSD', 'epidermis_BL.CSD', 'mac_CSD', 'neu_CSD')
  
  receivers = celltypes
  receivers_ctl = celltypes_ctl
  
  # cat(n, '-- ', timepoint, '\n')
  cat('selected cell types : \n ', celltypes, '\n')
  cat('controled cell types : \n ', celltypes_ctl, '\n')
  
  cat(' -- receiver cells  :  ', receivers, '\n')
  cat(' -- receiver cell of control niche :  ', receivers_ctl, '\n')
  
  ##########################################
  # from here the arguments are 
  # - BZ cell types
  # - RZ cell types
  # - receiver cell of BZ niche
  # - receiver cell of RZ niche
  ##########################################
  celltypes_sel = unique(c(celltypes, receivers, 
                           celltypes_ctl, receivers_ctl)
  )
  
  Idents(refs) = as.factor(refs$celltypes)
  subref = subset(refs, cells = colnames(refs)[!is.na(match(refs$celltypes, celltypes_sel))])
  subref$celltypes = droplevels(as.factor(subref$celltypes))
  table(subref$celltypes)
  
  rm(celltypes_sel)
  
  cat('celltype to consider -- ', names(table(subref$celltypes)), '\n')
  table(subref$celltypes)
  
  set.seed(0)
  subref = subset(x = subref, downsample = 3000) # downsample the CM and EC for the sake of speed
  table(subref$celltypes)
  Idents(subref) = subref$celltypes
  
  ##########################################
  # Gene filtering, in particular filtering gene without UMI counts or lowly expressed genes
  ##########################################
  Gene.filtering.preprocessing = TRUE
  if(Gene.filtering.preprocessing){
    sce <- as.SingleCellExperiment(subref)
    ggs = rownames(sce)
    
    colLabels(sce) = as.factor(sce$celltypes)
    rownames(sce) = get_geneName(rownames(sce))
    
    ave.counts <- calculateAverage(sce, assay.type = "counts")
    
    #hist(log10(ave.counts), breaks=100, main="", col="grey80",
    #     xlab=expression(Log[10]~"average count"))
    
    num.cells <- nexprs(sce, byrow=TRUE)
    #smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells",
    #              xlab=expression(Log[10]~"average count"))
    
    # detected in >= 10 cells, ave.counts cutoff are both very lose
    genes.to.keep <- num.cells > 5 & ave.counts >= 10^-5  & ave.counts <10^4  
    summary(genes.to.keep)
    sce <- sce[genes.to.keep, ]
    ggs = ggs[genes.to.keep]
    
    rownames(sce) = make.names(rownames(sce), unique = TRUE)
    
    geneDup = data.frame(gene = ggs, gg.uniq = rownames(sce), stringsAsFactors = FALSE)
    geneDup$geneSymbol = get_geneName(geneDup$gene)
    
    kk = which(geneDup$geneSymbol != geneDup$gg.uniq)
    
    subref = as.Seurat(sce, counts = "counts", data = "logcounts")
    rm(sce)
    Idents(subref) = as.factor(subref$celltypes)
    
    # saveRDS(geneDup, paste0(outDir, '/geneSymbol_duplication_inLRanalysis.rds'))
    
  }
  
  ##########################################
  # double check the subref before running nichenet
  ##########################################
  #subref = readRDS(file = paste0(outDir, '/seuratObject_snRNAseq_subset_for_NicheNet.rds'))
  #subref = subset(x = subref, downsample = 1000) 
  table(subref$celltypes)
  subref$celltype = subref$celltypes
  seurat_obj = SetIdent(subref, value = "celltype")
  
  table(seurat_obj$celltypes)
  #seurat_obj$celltype = as.factor(seurat_obj$celltypes)
  Idents(seurat_obj) = seurat_obj$celltype
  rm(subref)
  
  ##########################################
  # Step 1:  Define the niches/microenvironments of interest
  ##########################################
  seurat_obj$celltype = as.factor(seurat_obj$celltype)
  table(seurat_obj$celltype)
  
  cat('cell types : ', celltypes, '\n')
  cat('cell types in control niches : ', celltypes_ctl, '\n')
  
  # receiver_cells = 'CM_BZ'
  cat('--- receiver cells in main niches: ', receivers, ' ---\n')
  cat('--- receiver cells in control niches: ', receivers_ctl, ' ---\n')
  
  lfc_cutoff = 0.15
  expression_pct = 0.10
  top_n_target = 500
  
  for(n in 1:length(receivers))
  {
    # n = 1
    receiver = receivers[n] 
    receiver_ctl = receivers_ctl[n]
    
    cat('-- define niches with receiver ', receiver, ' and control ', receiver_ctl, '--\n')
    
    # define the niches
    niches = list(
      "sample_niche" = list(
        "sender" = celltypes, # including autocrine
        "receiver" = receiver),
      "control_niche" = list(
        "sender" = celltypes_ctl, 
        "receiver" = receiver_ctl)
    )
    
    res = prioritize_ligands_between_niches(seurat_obj, niches,
                                            expression_pct = expression_pct,
                                            lfc_cutoff = lfc_cutoff,
                                            top_n_target = top_n_target)
    
    saveRDS(res, file = paste0(outDir, '/nichenet_prioritization_tables_output',
                               '_receiver.', receiver,
                               '_receiver.control.', receiver_ctl, 
                               '_timepoint.', timepoint, '.rds'))
    
    prioritization_tables = res[[1]]
    output = res[[2]]
    
    ##########################################
    # 8). Visualization of the Differential NicheNet output 
    ##########################################
    top_ligand_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
      dplyr::select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% 
      group_by(ligand) %>% 
      top_n(1, prioritization_score) %>% 
      ungroup() %>% 
      dplyr::select(ligand, receptor, niche) %>% 
      dplyr::rename(top_niche = niche)
    
    top_ligand_receptor_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
      dplyr::select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% 
      group_by(ligand, receptor) %>% 
      top_n(1, prioritization_score) %>% 
      ungroup() %>% 
      dplyr::select(ligand, receptor, niche) %>% 
      dplyr::rename(top_niche = niche)
    
    for(ntop in c(50, 100, 200))
    {
      # ntop = 50
      cat('ntop -- ', ntop, '\n')
      ligand_prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
        dplyr::select(niche, sender, receiver, ligand, prioritization_score) %>% 
        group_by(ligand, niche) %>% top_n(1, prioritization_score) %>% 
        ungroup() %>% distinct() %>% 
        inner_join(top_ligand_niche_df) %>% 
        filter(niche == top_niche) %>% 
        group_by(niche) %>% 
        top_n(n = ntop, prioritization_score) %>% 
        ungroup() # get the top50 ligands per niche
      
      receiver_oi = receiver
      
      filtered_ligands = ligand_prioritized_tbl_oi %>% 
        filter(receiver == receiver_oi) %>% 
        pull(ligand) %>% unique()
      
      prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
        filter(ligand %in% filtered_ligands) %>% 
        dplyr::select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% 
        distinct() %>% 
        inner_join(top_ligand_receptor_niche_df) %>% 
        group_by(ligand) %>% 
        filter(receiver == receiver_oi) %>% 
        top_n(2, prioritization_score) %>% ungroup() 
      
      lfc_plot = make_ligand_receptor_lfc_plot(receiver_oi, 
                                               prioritized_tbl_oi, 
                                               prioritization_tables$prioritization_tbl_ligand_receptor, 
                                               plot_legend = FALSE, 
                                               heights = NULL, widths = NULL)
      lfc_plot
      
      ggsave(paste0(outDir, '/Ligand_receptors_LFC', 
                    '_receiver.', receiver,
                    '_receiver.control..', receiver_ctl, 
                    '_timepoint.', timepoint,
                    '_ntop.', ntop, '.pdf'), 
             width = 12, height = 12*ntop/50, limitsize = FALSE)
      
      
      exprs_activity_target_plot = 
        make_ligand_activity_target_exprs_plot(receiver_oi, 
                                               prioritized_tbl_oi,  
                                               prioritization_tables$prioritization_tbl_ligand_receptor,  
                                               prioritization_tables$prioritization_tbl_ligand_target, 
                                               output$exprs_tbl_ligand,  
                                               output$exprs_tbl_target, 
                                               lfc_cutoff, 
                                               ligand_target_matrix, 
                                               plot_legend = FALSE, 
                                               heights = NULL, widths = NULL)
      
      exprs_activity_target_plot$combined_plot
      ggsave(paste0(outDir, '/Combined_plots_ligand_noFiltering',
                    '_receiver.', receiver,
                    '_receiver.control.', receiver_ctl, 
                    '_timepoint.', timepoint,
                    '_ntop.', ntop, '.pdf'), 
             width = 60, height = 12*ntop/50, limitsize = FALSE)
    }
    
  }
  
}


prioritize_ligands_between_niches = function(seurat_obj, 
                                             niches,
                                             assay_oi = "RNA",
                                             expression_pct = 0.10,
                                             lfc_cutoff = 0.15,
                                             include_spatial_info_sender = FALSE,
                                             include_spatial_info_receiver = FALSE,
                                             top_n_target = 250
)
{
  ##########################################
  # step 0): Nichenet data loaded already 
  ##########################################
  # NicheNet’s ligand-target prior model
  ligand_target_matrix = readRDS(paste0(dataPath_nichenet,  "ligand_target_matrix.rds"))
  ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
  
  # ligand-receptor network, Putative ligand-receptor links from NicheNet
  lr_network = readRDS(paste0(dataPath_nichenet, "lr_network.rds"))
  lr_network = lr_network %>% mutate(bonafide = ! database %in% c("ppi_prediction","ppi_prediction_go"))
  lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% 
    distinct(ligand, receptor, bonafide)
  
  # If wanted, users can remove ligand-receptor interactions that were predicted based on 
  # protein-protein interactions and 
  # only keep ligand-receptor interactions that are described in curated databases.
  # lr_network = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
  head(lr_network)
  
  # ## weighted integrated networks 
  weighted_networks = readRDS(paste0(dataPath_nichenet,  "weighted_networks.rds"))
  head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
  head(weighted_networks$gr) # interactions and their weights in the gene regulatory network
  
  weighted_networks_lr = weighted_networks$lr_sig %>% 
    dplyr::rename(ligand = from, receptor = to) %>%
    inner_join(lr_network %>% 
                 distinct(ligand,receptor), by = c("ligand","receptor"))
  
  
  ##########################################
  # step 2. Calculate differential expression between the niches 
  # issue solved because find markers for all genes, should run only for receptors and ligands
  ##########################################
  # only ligands important for sender cell types
  DE_sender = calculate_niche_de(seurat_obj = seurat_obj %>% subset(features = lr_network$ligand %>% unique()), 
                                 niches = niches, 
                                 type = "sender", 
                                 assay_oi = assay_oi)
  
  # only receptors now, later on: DE analysis to find targets
  DE_receiver = calculate_niche_de(seurat_obj = seurat_obj %>% subset(features = lr_network$receptor %>% unique()), 
                                   niches = niches, 
                                   type = "receiver", 
                                   assay_oi = assay_oi)
  
  DE_sender = DE_sender %>% 
    mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), 
                               ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))
  DE_receiver = DE_receiver %>% 
    mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), 
                               ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))
  
  # expected percentage of cells expressing the genes 
  DE_sender_processed = process_niche_de(DE_table = DE_sender, 
                                         niches = niches, 
                                         expression_pct = expression_pct, 
                                         type = "sender")
  DE_receiver_processed = process_niche_de(DE_table = DE_receiver, 
                                           niches = niches, 
                                           expression_pct = expression_pct, 
                                           type = "receiver")
  
  # Combine sender-receiver DE based on L-R pairs:
  # As mentioned above, DE of ligands from one sender cell type is determined be calculating DE 
  # between that cell type, and all the sender cell types of the other niche. 
  # To summarize the DE of ligands of that cell type we have several options: 
  # we could take the average LFC, but also the minimum LFC compared to the other niche. 
  # We recommend using the minimum LFC, because this is the strongest specificity measure of ligand expression, 
  # because a high min LFC means that a ligand is more strongly expressed in the cell type of niche 1 
  # compared to all cell types of niche 2 (in contrast to a high average LFC, 
  # which does not exclude that one or more cell types in niche 2 also strongly express that ligand).
  specificity_score_LR_pairs = "min_lfc"
  DE_sender_receiver = combine_sender_receiver_de(DE_sender_processed, 
                                                  DE_receiver_processed, 
                                                  lr_network, 
                                                  specificity_score = specificity_score_LR_pairs)
  
  
  ##########################################
  # Step 3) Optional: Calculate differential expression between the different spatial regions 
  ##########################################
  # if not spatial info to include: put this to false 
  
  # user adaptation required on own dataset
  spatial_info = tibble(celltype_region_oi = "CAF_High", 
                        celltype_other_region = "myofibroblast_High", 
                        niche =  "pEMT_High_niche", celltype_type = "sender") 
  specificity_score_spatial = "lfc"
  
  # this is how this should be defined if you don't have spatial info (mock spatial info)
  if(include_spatial_info_sender == FALSE & include_spatial_info_receiver == FALSE){
    spatial_info = tibble(celltype_region_oi = NA, celltype_other_region = NA) %>% 
      mutate(niche =  niches %>% names() %>% head(1), celltype_type = "sender")
  } 
  
  
  if(include_spatial_info_sender == TRUE){
    sender_spatial_DE = calculate_spatial_DE(seurat_obj = seurat_obj %>% 
                                               subset(features = lr_network$ligand %>% unique()), 
                                             spatial_info = spatial_info %>% filter(celltype_type == "sender"))
    sender_spatial_DE_processed = process_spatial_de(DE_table = sender_spatial_DE, 
                                                     type = "sender", 
                                                     lr_network = lr_network, 
                                                     expression_pct = expression_pct, 
                                                     specificity_score = specificity_score_spatial)
    
    # add a neutral spatial score for sender celltypes in which the spatial is not known / not of importance
    sender_spatial_DE_others = get_non_spatial_de(niches = niches, 
                                                  spatial_info = spatial_info, 
                                                  type = "sender", 
                                                  lr_network = lr_network)
    sender_spatial_DE_processed = sender_spatial_DE_processed %>% bind_rows(sender_spatial_DE_others)
    
    sender_spatial_DE_processed = sender_spatial_DE_processed %>% 
      mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))
    
  }else{
    #add a neutral spatial score for all sender celltypes (for none of them, spatial is relevant in this case)
    sender_spatial_DE_processed = get_non_spatial_de(niches = niches, spatial_info = spatial_info, 
                                                     type = "sender", lr_network = lr_network)
    sender_spatial_DE_processed = sender_spatial_DE_processed %>% 
      mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))  
    
  }
  
  if(include_spatial_info_receiver == TRUE){
    receiver_spatial_DE = calculate_spatial_DE(seurat_obj = seurat_obj %>% 
                                                 subset(features = lr_network$receptor %>% unique()), 
                                               spatial_info = spatial_info %>% 
                                                 filter(celltype_type == "receiver"))
    receiver_spatial_DE_processed = process_spatial_de(DE_table = receiver_spatial_DE, 
                                                       type = "receiver", 
                                                       lr_network = lr_network, 
                                                       expression_pct = expression_pct, 
                                                       specificity_score = specificity_score_spatial)
    
    # add a neutral spatial score for receiver celltypes in which the spatial is not known / not of importance
    receiver_spatial_DE_others = get_non_spatial_de(niches = niches, 
                                                    spatial_info = spatial_info, 
                                                    type = "receiver", 
                                                    lr_network = lr_network)
    receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% bind_rows(receiver_spatial_DE_others)
    
    receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% 
      mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
    
  }else{
    # # add a neutral spatial score for all receiver celltypes (for none of them, spatial is relevant in this case)
    receiver_spatial_DE_processed = get_non_spatial_de(niches = niches, 
                                                       spatial_info = spatial_info, 
                                                       type = "receiver", 
                                                       lr_network = lr_network)
    receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% 
      mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
  }
  
  
  ##########################################
  # Step 4). Calculate ligand activities and infer active ligand-target links
  ##########################################
  ## first define the target genes by DE
  # It is always useful to check the number of genes in the geneset before doing the ligand activity analysis. 
  # We recommend having between 20 and 1000 genes in the geneset of interest, 
  # and a background of at least 5000 genes for a proper ligand activity analysis. 
  # If you retrieve too many DE genes, it is recommended to use a higher lfc_cutoff threshold. 
  # We recommend using a cutoff of 0.15 if you have > 2 receiver cells/niches to compare and 
  # use the min_lfc as specificity score. 
  # If you have only 2 receivers/niche, we recommend using a higher threshold (such as using 0.25). 
  # If you have single-cell data like Smart-seq2 with high sequencing depth, 
  #we recommend to also use higher threshold.
  # Smart-seq2 data and only 2 niches to compare, 
  # so we will use a stronger LFC threshold to keep less DE genes, but more trustworthy ones.
  # lfc_cutoff = 0.15 # recommended for 10x as min_lfc cutoff. 
  specificity_score_targets = "min_lfc"
  
  # here DE all genes in the receiver cells, not only for ligand and receptors as before
  cat('find DE targets in receiver cell, may take some time \n')
  DE_receiver_targets = calculate_niche_de_targets(seurat_obj = seurat_obj, 
                                                   niches = niches, 
                                                   lfc_cutoff = lfc_cutoff, 
                                                   expression_pct = expression_pct, 
                                                   assay_oi = assay_oi)
  
  DE_receiver_processed_targets = process_receiver_target_de(DE_receiver_targets = DE_receiver_targets, 
                                                             niches = niches, 
                                                             expression_pct = expression_pct, 
                                                             specificity_score = specificity_score_targets)
  
  
  background = DE_receiver_processed_targets  %>% pull(target) %>% unique()
  geneset_niche1 = DE_receiver_processed_targets %>% 
    filter(receiver == niches[[1]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & 
             target_present == 1) %>% pull(target) %>% unique()
  geneset_niche2 = DE_receiver_processed_targets %>% 
    filter(receiver == niches[[2]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & 
             target_present == 1) %>% pull(target) %>% unique()
  
  # Good idea to check which genes will be left out of the ligand activity analysis 
  # (=when not present in the rownames of the ligand-target matrix).
  # If many genes are left out, this might point to some issue in the gene naming 
  # (eg gene aliases and old gene symbols, bad human-mouse mapping)
  geneset_niche1 %>% setdiff(rownames(ligand_target_matrix))
  geneset_niche2 %>% setdiff(rownames(ligand_target_matrix))
  
  length(geneset_niche1)
  length(geneset_niche2)
  
  niche_geneset_list = list(
    "BZ_niche" = list(
      "receiver" = niches[[1]]$receiver,
      "geneset" = geneset_niche1,
      "background" = background),
    "noBZ_niche" = list(
      "receiver" = niches[[2]]$receiver,
      "geneset" = geneset_niche2 ,
      "background" = background)
  )
  
  # get_ligand_activites_targets calling the function predict_ligand_activities for each niche
  ligand_activities_targets = get_ligand_activities_targets(niche_geneset_list = niche_geneset_list, 
                                                            ligand_target_matrix = ligand_target_matrix, 
                                                            top_n_target = top_n_target)
  
  head(ligand_activities_targets)
  ligand_activities_targets %>% arrange(-activity) %>% filter(receiver %in% niches$BZ_niche$receiver) %>%
    filter(ligand == 'ITGA4')
  ligand_activities_targets %>% arrange(-activity) %>% filter(receiver %in% niches$RZ_niche$receiver) %>%
    filter(ligand == 'ITGA4')
  
  ##########################################
  # step 5. Calculate (scaled) expression of ligands, receptors and targets 
  # across cell types of interest (log expression values and expression fractions)
  # we will calculate average (scaled) expression, and fraction of expression, of 
  # ligands, receptors, and target genes across all cell types of interest
  ##########################################
  features_oi = union(lr_network$ligand, lr_network$receptor) %>% 
    union(ligand_activities_targets$target) %>% setdiff(NA)
  
  require(dplyr)
  # save dotplot 
  dotplot = suppressWarnings(Seurat::DotPlot(seurat_obj %>% subset(idents = niches %>% unlist() %>% unique()), 
                                             features = features_oi, assay = assay_oi))
  
  exprs_tbl = dotplot$data %>% as_tibble()
  exprs_tbl = exprs_tbl %>% 
    dplyr::rename(celltype = id, gene = features.plot, expression = avg.exp, 
                  expression_scaled = avg.exp.scaled, fraction = pct.exp) %>%
    mutate(fraction = fraction/100) %>% as_tibble() %>% 
    dplyr::select(celltype, gene, expression, expression_scaled, fraction) %>% 
    distinct() %>% arrange(gene) %>% mutate(gene = as.character(gene))
  
  exprs_tbl_ligand = exprs_tbl %>% filter(gene %in% lr_network$ligand) %>% 
    dplyr::rename(sender = celltype, ligand = gene, ligand_expression = expression, 
                  ligand_expression_scaled = expression_scaled, ligand_fraction = fraction)
  
  exprs_tbl_receptor = exprs_tbl %>% 
    filter(gene %in% lr_network$receptor) %>% 
    dplyr::rename(receiver = celltype, receptor = gene, receptor_expression = expression, 
                  receptor_expression_scaled = expression_scaled, receptor_fraction = fraction)
  exprs_tbl_target = exprs_tbl %>% 
    filter(gene %in% ligand_activities_targets$target) %>% 
    dplyr::rename(receiver = celltype, target = gene, target_expression = expression, 
                  target_expression_scaled = expression_scaled, target_fraction = fraction)
  
  exprs_tbl_ligand = exprs_tbl_ligand %>% 
    mutate(scaled_ligand_expression_scaled = scale_quantile_adapted(ligand_expression_scaled)) %>% 
    mutate(ligand_fraction_adapted = ligand_fraction) %>% 
    mutate_cond(ligand_fraction >= expression_pct, ligand_fraction_adapted = expression_pct)  %>% 
    mutate(scaled_ligand_fraction_adapted = scale_quantile_adapted(ligand_fraction_adapted))
  
  exprs_tbl_receptor = exprs_tbl_receptor %>% 
    mutate(scaled_receptor_expression_scaled = scale_quantile_adapted(receptor_expression_scaled)) %>% 
    mutate(receptor_fraction_adapted = receptor_fraction) %>% 
    mutate_cond(receptor_fraction >= expression_pct, receptor_fraction_adapted = expression_pct) %>% 
    mutate(scaled_receptor_fraction_adapted = scale_quantile_adapted(receptor_fraction_adapted))
  
  ##########################################
  # step 6. Expression fraction and receptor
  # score ligand-receptor interactions based on expression strength of the receptor, 
  # in such a way that we give higher scores to the most strongly expressed receptor of a certain ligand, 
  # in a certain celltype. 
  # This will not effect the rank of individual ligands later on, 
  # but will help in prioritizing the most important receptors per ligand 
  # (next to other factors regarding the receptor - see later).
  ##########################################
  exprs_sender_receiver = lr_network %>% 
    inner_join(exprs_tbl_ligand, by = c("ligand")) %>% 
    inner_join(exprs_tbl_receptor, by = c("receptor")) %>% 
    inner_join(DE_sender_receiver %>% distinct(niche, sender, receiver))
  
  ligand_scaled_receptor_expression_fraction_df = exprs_sender_receiver %>% 
    group_by(ligand, receiver) %>% 
    mutate(rank_receptor_expression = dense_rank(receptor_expression), 
           rank_receptor_fraction  = dense_rank(receptor_fraction)) %>% 
    mutate(ligand_scaled_receptor_expression_fraction = 0.5* ((rank_receptor_fraction / max(rank_receptor_fraction)) + 
                                                                ((rank_receptor_expression / max(rank_receptor_expression))) ))  %>% 
    distinct(ligand, receptor, receiver, ligand_scaled_receptor_expression_fraction, bonafide) %>% 
    distinct() %>% ungroup() 
  
  ##########################################
  # step 7. Prioritization of ligand-receptor and ligand-target links 
  ##########################################
  prioritizing_weights = c(
    "scaled_ligand_score" = 5, # # niche-specific expression of the ligand: Recommended 5 (between 1-5)
    # scaled ligand expression in one sender compared to the other cell types in the dataset
    "scaled_ligand_expression_scaled" = 1,
    "ligand_fraction" = 1, # Ligands expressed in a smaller fraction of cells of cell type than cutoff (default: 0.10)
    "scaled_ligand_score_spatial" = 0, 
    # receptor DE score: niche-specific expression, Recommended 0.5 (>=0.5 and lower than "scaled_ligand_score")
    "scaled_receptor_score" = 1, 
    "scaled_receptor_expression_scaled" = 0.5, # Recommended weight: 0.5
    # Receptors that are expressed in a smaller fraction of cells of a cell type than exprs_cutoff(default: 0.10) 
    # will get a lower ranking, proportional to their fraction, Recommended weight: 1. 
    "receptor_fraction" = 1, 
    # this factor let us give higher weights to the most highly expressed receptor of a ligand in the receiver.
    # Recommended value: 1 (minimum: 0.5)
    "ligand_scaled_receptor_expression_fraction" = 1,
    "scaled_receptor_score_spatial" = 0,
    # Absolute ligand activity: to further prioritize ligand-receptor pairs based on their predicted effect of 
    # the ligand-receptor interaction on the gene expression in the receiver cell type - 
    # prioritizing_weights argument: "scaled_activity". Recommended weight: 0, 
    # unless absolute enrichment of target genes is of specific interest.
    "scaled_activity" = 0, 
    # Normalized ligand activity: to further prioritize ligand-receptor pairs based on their predicted effect of
    # the ligand-receptor interaction on the gene expression in the receiver cell type, Recommended weight: >=1.
    "scaled_activity_normalized" = 4,
    "bona_fide" = 1)
  
  output = list(DE_sender_receiver = DE_sender_receiver, 
                ligand_scaled_receptor_expression_fraction_df = ligand_scaled_receptor_expression_fraction_df, 
                sender_spatial_DE_processed = sender_spatial_DE_processed, 
                receiver_spatial_DE_processed = receiver_spatial_DE_processed,
                ligand_activities_targets = ligand_activities_targets, 
                DE_receiver_processed_targets = DE_receiver_processed_targets, 
                exprs_tbl_ligand = exprs_tbl_ligand,  
                exprs_tbl_receptor = exprs_tbl_receptor, 
                exprs_tbl_target = exprs_tbl_target
  )
  
  prioritization_tables = get_prioritization_tables(output, prioritizing_weights)
  
  prioritization_tables$prioritization_tbl_ligand_receptor %>% 
    filter(receiver == niches[[1]]$receiver) %>% 
    head(10)
  
  prioritization_tables$prioritization_tbl_ligand_target %>% 
    filter(receiver == niches[[2]]$receiver) %>% 
    head(10)
  
  Test_by_plotting = FALSE
  if(Test_by_plotting){
    
    top_ligand_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
      dplyr::select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% 
      group_by(ligand) %>% 
      top_n(1, prioritization_score) %>% 
      ungroup() %>% 
      dplyr::select(ligand, receptor, niche) %>% 
      dplyr::rename(top_niche = niche)
    
    top_ligand_receptor_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
      dplyr::select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% 
      group_by(ligand, receptor) %>% 
      top_n(1, prioritization_score) %>% 
      ungroup() %>% 
      dplyr::select(ligand, receptor, niche) %>% 
      dplyr::rename(top_niche = niche)
    
    ntop = 50
    cat('ntop -- ', ntop, '\n')
    
    xx =  prioritization_tables$prioritization_tbl_ligand_receptor %>% 
      filter(sender == 'RBC') %>% 
      dplyr::group_by(ligand_receptor, niche) %>%
      filter(ligand_receptor %in% c('FGF17--FGFR1', 'ITGA4--VCAM1', 'PTDSS1--ERBB2', 'WNT7B--LRP5')) %>% 
      data.frame()
    
    ligand_prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
      dplyr::select(niche, sender, receiver, ligand, prioritization_score) %>% 
      group_by(ligand, niche) %>% top_n(1, prioritization_score) %>% 
      ungroup() %>% distinct() %>% 
      inner_join(top_ligand_niche_df) %>% 
      filter(niche == top_niche) %>% 
      group_by(niche) %>% 
      top_n(n = ntop, prioritization_score) %>% 
      ungroup() # get the top50 ligands per niche
    
    #  select top ligand-receptor pairs for cell population of interest
    # (here, we will take the top 2 scoring receptors per prioritized ligand)
    receiver_oi = receiver
    filtered_ligands = ligand_prioritized_tbl_oi %>% 
      filter(receiver == receiver_oi) %>% 
      pull(ligand) %>% unique()
    
    prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
      filter(ligand %in% filtered_ligands) %>% 
      dplyr::select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% 
      distinct() %>% 
      inner_join(top_ligand_receptor_niche_df) %>% 
      group_by(ligand_receptor) %>%
      #group_by(ligand) %>% 
      #filter(receiver == receiver_oi) %>% 
      top_n(2, prioritization_score) %>% ungroup() 
    
    
    lfc_plot = make_ligand_receptor_lfc_plot(receiver_oi, 
                                             prioritized_tbl_oi, 
                                             prioritization_tables$prioritization_tbl_ligand_receptor, 
                                             plot_legend = FALSE, 
                                             heights = NULL, widths = NULL)
    lfc_plot
    
    
  }
  
  return(list(prioritization_tables, output))
  
  
}


make_ligand_receptor_lfc_plot_customized = function(receiver_oi, 
                                                    prioritized_tbl_oi, 
                                                    prioritization_tbl_ligand_receptor, 
                                                    plot_legend = TRUE, 
                                                    heights = NULL, 
                                                    widths = NULL)
{
  # prioritization_tbl_ligand_receptor = prioritization_tables$prioritization_tbl_ligand_receptor
  
  filtered_ligand_receptors = prioritized_tbl_oi %>% pull(ligand_receptor) %>% unique()
  
  ordered_ligand_receptors = prioritization_tbl_ligand_receptor %>% 
    filter(ligand_receptor %in% filtered_ligand_receptors) %>% 
    select(niche, sender, ligand, receptor, ligand_receptor, prioritization_score) %>% 
    distinct() %>% 
    group_by(ligand_receptor) %>% 
    summarise(prioritization_score = max(prioritization_score)) %>% 
    inner_join(prioritization_tbl_ligand_receptor %>% 
                 select(niche, sender, ligand, receptor, ligand_receptor, prioritization_score) %>% 
                 distinct()) %>% 
    arrange(sender, prioritization_score)
  ordered_ligand_receptors_max_ligand_score = prioritization_tbl_ligand_receptor %>% 
    filter(ligand_receptor %in% filtered_ligand_receptors) %>% 
    select(niche, sender, ligand, prioritization_score) %>% 
    distinct() %>% 
    group_by(ligand) %>% 
    summarise(prioritization_score_ligand = max(prioritization_score)) %>% 
    inner_join(prioritization_tbl_ligand_receptor %>% 
                 select(niche, sender, ligand, prioritization_score) %>% 
                 distinct()) %>% 
    arrange(sender, prioritization_score_ligand) %>% distinct()
  
  ordered_ligand_receptors = ordered_ligand_receptors %>% 
    inner_join(ordered_ligand_receptors_max_ligand_score) %>% 
    arrange(sender, prioritization_score_ligand, prioritization_score)
  ordered_ligand_receptors = ordered_ligand_receptors %>% 
    mutate(ligand_receptor_ordered = factor(ligand_receptor, ordered = T, 
                                            levels = unique(ordered_ligand_receptors$ligand_receptor))) %>% 
    distinct(ligand_receptor, ligand, receptor, ligand_receptor_ordered, niche) %>% 
    dplyr::rename(niche_prior = niche)
  
  plot_data = prioritization_tbl_ligand_receptor %>% inner_join(ordered_ligand_receptors)
  
  p_lig_lfc = plot_data %>% 
    mutate(scores_test = abs(scaled_activity_normalized) * abs(ligand_score)) %>%
    #filter(receiver %in% niches[[1]]$receiver) %>%
    #filter(niche %in% names(niches)[1]) %>% 
    #ggplot(aes(sender, ligand_receptor_ordered, fill = ligand_score)) +
    ggplot(aes(sender, ligand_receptor_ordered, fill = scores_test)) +
    geom_tile(color = "black") +
    facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.75, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + 
    labs(fill = "Ligand:\nmin LFC vs\nother niches") + 
    ylab("Prioritized Ligand-Receptor pairs") + xlab("Ligand LFC\n in Sender")
  
  max_lfc = max(abs(plot_data$ligand_score) %>% max(), abs(plot_data$receptor_score) %>% max())
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "OrRd"), 
                                           #%>% 
                                           #   rev(),
                                           # values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  
                                           #limits = c(-1*max_lfc, max_lfc),
                                           limits = c(0, max_lfc),
  )
  
  p_lig_lfc = p_lig_lfc + custom_scale_fill
  p_lig_lfc
  
  p_lig_lfc_ctl = plot_data %>% 
    filter(receiver %in% niches[[2]]$receiver) %>%
    filter(niche %in% names(niches)[2]) %>% 
    ggplot(aes(sender, ligand_receptor_ordered, fill = ligand_score)) +
    geom_tile(color = "black") +
    facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.75, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + 
    labs(fill = "Ligand:\nmin LFC vs\nother niches")  + 
    ylab("Prioritized Ligand-Receptor pairs") + xlab("Ligand LFC\n in Sender") +
    custom_scale_fill
  
  design = "A#B"
  p_L_niches = patchwork::wrap_plots(A = p_lig_lfc, 
                                     B = p_lig_lfc_ctl + ylab(""), 
                                     nrow = 1, guides = "collect", 
                                     design = design, 
                                     widths = c(plot_data$sender %>% unique() %>% length(), 1, 
                                                plot_data$receiver %>% unique() %>% length() +0.5))
  p_L_niches
  
  p_rec_lfc = plot_data %>%
    ggplot(aes(receiver, ligand_receptor_ordered, fill = receptor_score)) +
    geom_tile(color = "black") +
    facet_grid(~receiver, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 0),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.75, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Receptor:\nmin LFC vs\nother niches")  + xlab("Receptor LFC\n in Receiver")
  max_lfc = max(abs(plot_data$ligand_score) %>% max(), abs(plot_data$receptor_score) %>% max())
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))
  
  p_rec_lfc = p_rec_lfc + custom_scale_fill
  p_rec_lfc
  
  design = "A#B"
  p_LR_pair = patchwork::wrap_plots(A = p_lig_lfc, B = p_rec_lfc + ylab(""), nrow = 1, guides = "collect", design = design, widths = c(plot_data$sender %>% unique() %>% length(), 1 ,plot_data$receiver %>% unique() %>% length() +0.5))
  p_LR_pair
  
  
  
}







