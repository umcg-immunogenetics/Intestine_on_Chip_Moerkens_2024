#### Data quality control and filtering

# locations of the object
seurat_objects_loc <- #path to folder where the intestine-on-chip data object is saved
# read the data object that is demultiplexed, filtered for doublets assigned by SoupOrCell, assigned to donor information based on genotypes and contains relevant metadata  
seurat_object_assigned <- paste(seurat_objects_loc, 'intestine_on_chip_annotated.rds', sep = '')
ooac <- readRDS(seurat_object_assigned)

#libraries
library(Seurat)
library(ggplot2)

#plot threshold for percentage of mitochondrial genes per cell
VlnPlot(ooac, features = "percent.mt") & geom_hline(yintercept = 20)

#discard dead cells based on percentage mitochondrial genes, here 20% or higher is discarded
ooac <- subset(ooac, subset = percent.mt < 20)

#discard cluster 7, the cluster with dying cells and MT-gene enrichment
ooac <- ooac[, !(ooac@meta.data[['seurat_clusters']] %in% c(7))]

#plot threshold for number of reads (transcripts) per cell
VlnPlot(ooac, features = "nCount_RNA") & geom_hline(yintercept = 100000)

#discard doublets cells based on number of transcripts, here 100.000 transcripts per cell or higher is discarded
ooac <- subset(ooac, subset = nCount_RNA < 100000)

#plot threshold for number of genes per cell
VlnPlot(ooac, features = "nFeature_RNA") & geom_hline(yintercept = 2000)

#discard debris and dead cells based on number of genes per cell, here 2000 genes per cell or lower is discarded
ooac <- subset(ooac, subset = nFeature_RNA > 2000)

# do normal normalization
ooac <- NormalizeData(ooac)
#Run FindVariableFeatures
ooac <- FindVariableFeatures(ooac, selection.method = "vst", nfeatures = 2000)
# do SCT
ooac <- SCTransform(ooac)
# perform PCA
ooac <- RunPCA(ooac)
# do UMAP 2d reduction
ooac <- RunUMAP(ooac, assay = 'SCT', reduction = 'pca', dims = 1:30)
# do neighbours
ooac <- FindNeighbors(ooac, assay = 'SCT', reduction = 'pca', dims = 1:30)
# do clustering
ooac <- FindClusters(ooac, resolution = 1.2)


#### Remove doublets based on scDblFinder

#Libraries
library(Seurat)
library(tidyverse)
library(scDblFinder)
library(scater)

# Functions
findDoublets <- function(lane, obj) {
  start.time <- Sys.time()
  
  lane_obj <- subset(obj, subset = lane == lane)
  lane_obj <- as.SingleCellExperiment(lane_obj)
  doublets <- recoverDoublets(lane_obj, doublets = lane_obj$known_doublet, samples = table(lane_obj$sample_id))
  
  doublets <- as.data.frame(doublets)
  doublets$barcode <- lane_obj$barcode_lane
  doublets$lane <- lane_obj$lane
  
  return(doublets)
  
  end.time() <- Sys.time()
  time.taken <- end.time - start.time
  
  print(paste0('Lane ', lane, 'done. Time taken: ', time.taken))
}

# locations of data objects
seurat_objects_loc <- #path to folder where the intestine-on-chip data object is saved
# location of the data object containing the different-genotype doublets assiged by SoupOrCell (before trimming) 
seurat_object_demultiplexed <- paste(seurat_objects_loc, 'intestine_on_chip_demultiplexed.rds', sep = '')

#read your Seurat object before trimming (with SoupOrCell assigned doublets)
obj <- readRDS(seurat_object_demultiplexed)

# Getting the metadata from your object
metadata <- obj@meta.data

# Copying your column soup_assignment 
metadata$sample_id <- metadata$soup_assignment

# Adding the column with TRUE or FALSE for doublets and doublets samples for NA 
metadata <- metadata %>%
  mutate(known_doublet = if_else(soup_status == 'doublet', TRUE, FALSE)) %>%
  mutate(sample_id = replace(sample_id, soup_status != 'singlet', NA))

# Save the metadata back to the object
obj@meta.data <- metadata

# This line will select the lanes you have in your dataset 
lanes <- unique(obj$lane)

# Now we apply the function that we defined above “findDoublets” to each lane in “lanes” in our obj
all_doublets <- map(lanes, ~findDoublets(lane = .x, obj = obj))

# The all_doublets objects is a list with a dataframe per lane, the line below turns it into a single dataframe
all_doublets <- bind_rows(all_doublets)

# Retrieve the barcodes of the doublets (both intra-sample/same-genotype (predicted == TRUE) and inter-sample/different-genotype (known == TRUE)
all_doublets_true <- all_doublets %>%
  filter(known == TRUE | predicted == TRUE) %>% pull(barcode)

#Rename the ‘barcode’ column to ‘barcode_lane’ in the all_doublets objects
#Rename the other columns to the preferred name to add as metadata to the Seurat object
names(all_doublets)[names(all_doublets) == "barcode"] <- "barcode_lane"
names(all_doublets)[names(all_doublets) == "known"] <- "doublet_association_known"
names(all_doublets)[names(all_doublets) == "predicted"] <- "doublet_association_predicted"
names(all_doublets)[names(all_doublets) == "proportion"] <- "doublet_association_proportion"

#get all the unique barcodes of the predicted doublets, removing the duplicates
all_doublets_unique <- unique(all_doublets)

#Function to add specific columns from all_doublets_unique to the Seurat object
add_doublet_association_assignments <- function(seurat_object, doublet_output){
  # add the lane+barcode as rownames
  rownames(doublet_output) <- doublet_output[['barcode_lane']]
  # remove the column we don't want
  doublet_output[['lane']] <- NULL
  doublet_output[['barcode_lane']] <- NULL
  # add the doublet prepend
  # colnames(doublet_output) <- paste('scr', colnames(doublet_output), sep = '_')
  # now add each column to the object
  for(column in colnames(doublet_output)){
    seurat_object <- AddMetaData(seurat_object, doublet_output[column])
  }
  return(seurat_object)
}

# add the doublet association assignments to the Seurat object that was generated after previous trimming steps 
ooac <- add_doublet_association_assignments(ooac, all_doublets_unique)

# verify that the doublets from ‘all_doublets’ object are the same as doublets in the new metadata column 
DimPlot(ooac, reduction = 'umap', group.by = 'doublet_association_predicted', cols = c('TRUE' = 'red', 'FALSE' = 'grey'))
DimPlot(ooac, reduction = 'umap', cells.highlight = all_doublets_true, label = T, repel = T)

# remove the doublets called by doublet association script
ooac <- ooac[, ooac@meta.data[['doublet_association_predicted']] == 'FALSE']

#after discarding cells from object, re-do the clustering
# do normal normalization
ooac <- NormalizeData(ooac)
#Run FindVariableFeatures
ooac <- FindVariableFeatures(ooac, selection.method = "vst", nfeatures = 2000)
# do SCT, default assay of command is ‘RNA’, new assay is ‘SCT’. Sets the default assay to SCT.
ooac <- SCTransform(ooac)
# perform PCA
ooac <- RunPCA(ooac)
# do UMAP 2d reduction
ooac <- RunUMAP(ooac, assay = 'SCT', reduction = 'pca', dims = 1:30)
# do neighbours
ooac <- FindNeighbors(ooac, assay = 'SCT', reduction = 'pca', dims = 1:30)
# do clustering
ooac <- FindClusters(ooac, resolution = 1.2)


#### Correct for the variable cell_line using the reciprocal PCA method and reduce the resolution for the clusters and annotate the clusters and compartments

#libraries
library(Seurat)
library(SeuratData)
library(patchwork)
library(ggplot2)
library(viridisLite)
library(viridis)
library(RColorBrewer)
library(cowplot)
library(glmGamPoi)

#set seed to ensure reproducible data
set.seed(6)

# split the dataset into a list of seurat objects 
cell_line.list <- SplitObject(ooac, split.by = "cell_line")

# normalize and identify variable features for each dataset independently
cell_line.list <- lapply(X = cell_line.list, FUN = SCTransform)

# select features that are repeatedly variable across datasets for integration run PCA on each dataset using these features
features <- SelectIntegrationFeatures(object.list = cell_line.list, nfeatures = 3000)
cell_line.list <- PrepSCTIntegration(object.list = cell_line.list, anchor.features = features)
cell_line.list <- lapply(X = cell_line.list, FUN = RunPCA, features = features)

#identify anchors 
cell_line.anchors <- FindIntegrationAnchors(object.list = cell_line.list, normalization.method = "SCT", anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 7)

# create an 'integrated' data assay
cell_line.combined.sct <- IntegrateData(anchorset = cell_line.anchors, normalization.method = "SCT", dims = 1:30)

# Re-do normalization, FindVariableFeatures and ScaleData in RNA assay
cell_line.combined.sct <- NormalizeData(cell_line.combined.sct, assay = "RNA")
cell_line.combined.sct <- FindVariableFeatures(cell_line.combined.sct, assay = "RNA", selection.method = "vst", nfeatures = 2000)
cell_line.combined.sct <- ScaleData(cell_line.combined.sct, assay = 'RNA')

#repeat SCTransform, use v2 as it is suited for DE gene analysis.
cell_line.combined.sct <- SCTransform(cell_line.combined.sct, assay = "RNA", vst.flavor = "v2")

#Perform the clustering using the ‘integrated’ assay
DefaultAssay(cell_line.combined.sct) <- 'integrated'

#re-do clustering on the integrated assay for annotation and visulation of clusters
# perform PCA
cell_line.combined.sct <- RunPCA(cell_line.combined.sct)
# do UMAP 2d reduction
cell_line.combined.sct <- RunUMAP(cell_line.combined.sct, assay = 'integrated', reduction = 'pca', dims = 1:30)
# do neighbours
cell_line.combined.sct <- FindNeighbors(cell_line.combined.sct, assay = 'integrated', reduction = 'pca', dims = 1:30)
# do clustering
cell_line.combined.sct <- FindClusters(cell_line.combined.sct, resolution = 1.2)

# change name of Seurat object to ooac
ooac <- cell_line.combined.sct


#### Perform cell type and compartment annotation

# Alter resolution
ooac <- FindClusters(ooac, resolution = 0.34)
DimPlot(ooac, reduction = "umap", label = TRUE)

# add cell type annotation  
ooac@meta.data$annotation_res0.34 <- 'none'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.34 == 0, 'annotation_res0.34'] <- 'enterocyte_mature_1'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.34 == 1, 'annotation_res0.34'] <- 'paneth'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.34 == 2, 'annotation_res0.34'] <- 'enterocyte_mature_2'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.34 == 3, 'annotation_res0.34'] <- 'transit_amplifying_stem_cell'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.34 == 4, 'annotation_res0.34'] <- 'enterocyte_progenitor'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.34 == 5, 'annotation_res0.34'] <- 'mesenchymal_dividing'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.34 == 6, 'annotation_res0.34'] <- 'epithelial_IFNy'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.34 == 7, 'annotation_res0.34'] <- 'mesenchymal_fibroblast'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.34 == 8, 'annotation_res0.34'] <- 'epithelial_NTS'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.34 == 9, 'annotation_res0.34'] <- 'mesenchymal_neuronal'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.34 == 10, 'annotation_res0.34'] <- 'VIM_EPCAM_1'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.34 == 11, 'annotation_res0.34'] <- 'goblet'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.34 == 12, 'annotation_res0.34'] <- 'enteroendocrine'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.34 == 13, 'annotation_res0.34'] <- 'epithelial_IFNb'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.34 == 14, 'annotation_res0.34'] <- 'mesenchymal_WNT4'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.34 == 15, 'annotation_res0.34'] <- 'VIM_EPCAM_2'

# improve and clarify cell type names
ooac@meta.data$annotation_res0.34_new2 <- 'none'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'enterocyte_mature_1', 'annotation_res0.34_new2'] <- 'Enterocyte type 1'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'paneth', 'annotation_res0.34_new2'] <- 'Paneth-like cell'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'enterocyte_mature_2', 'annotation_res0.34_new2'] <- 'Enterocyte type 2'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'transit_amplifying_stem_cell', 'annotation_res0.34_new2'] <- 'TA/stem cell'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'enterocyte_progenitor', 'annotation_res0.34_new2'] <- 'Enterocyte progenitor'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'mesenchymal_dividing', 'annotation_res0.34_new2'] <- 'Dividing mesenchymal/neural cell'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'epithelial_IFNy', 'annotation_res0.34_new2'] <- 'IFNγ-responding epithelial cell'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'mesenchymal_fibroblast', 'annotation_res0.34_new2'] <- 'Myofibroblast'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'epithelial_NTS', 'annotation_res0.34_new2'] <- 'Mesenchymal-like epithelial precursor'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'mesenchymal_neuronal', 'annotation_res0.34_new2'] <- 'Neuron'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'VIM_EPCAM_1', 'annotation_res0.34_new2'] <- 'Mesenchymal-like epithelial cell type 1'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'goblet', 'annotation_res0.34_new2'] <- 'Goblet cell'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'enteroendocrine', 'annotation_res0.34_new2'] <- 'Enteroendocrine cell'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'epithelial_IFNb', 'annotation_res0.34_new2'] <- 'IFNβ-responding epithelial cell'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'mesenchymal_WNT4', 'annotation_res0.34_new2'] <- 'WNT4-positive neural cell'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'VIM_EPCAM_2', 'annotation_res0.34_new2'] <- 'Mesenchymal-like epithelial cell type 2'

# simplify cell type names and clustering
ooac@meta.data$annotation_res0.34_simple2 <- 'none'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'enterocyte_mature_1', 'annotation_res0.34_simple2'] <- 'Enterocyte'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'paneth', 'annotation_res0.34_simple2'] <- 'Paneth-like cell'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'enterocyte_mature_2', 'annotation_res0.34_simple2'] <- 'Enterocyte'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'transit_amplifying_stem_cell', 'annotation_res0.34_simple2'] <- 'TA/stem cell'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'enterocyte_progenitor', 'annotation_res0.34_simple2'] <- 'Enterocyte'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'mesenchymal_dividing', 'annotation_res0.34_simple2'] <- 'Dividing mesenchymal/neural cell'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'epithelial_IFNy', 'annotation_res0.34_simple2'] <- 'IFNγ-responding epithelial cell'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'mesenchymal_fibroblast', 'annotation_res0.34_simple2'] <- 'Myofibroblast'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'epithelial_NTS', 'annotation_res0.34_simple2'] <- 'Mesenchymal-like epithelial cell'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'mesenchymal_neuronal', 'annotation_res0.34_simple2'] <- 'Neuron'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'VIM_EPCAM_1', 'annotation_res0.34_simple2'] <- 'Mesenchymal-like epithelial cell'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'goblet', 'annotation_res0.34_simple2'] <- 'Goblet cell'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'enteroendocrine', 'annotation_res0.34_simple2'] <- 'Enteroendocrine cell'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'epithelial_IFNb', 'annotation_res0.34_simple2'] <- 'IFNβ-responding epithelial cell'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'mesenchymal_WNT4', 'annotation_res0.34_simple2'] <- 'WNT4-positive neural cell'
ooac@meta.data[ooac@meta.data$annotation_res0.34 == 'VIM_EPCAM_2', 'annotation_res0.34_simple2'] <- 'Mesenchymal-like epithelial cell'

# improve condition names
ooac@meta.data$condition_new <- 'none'
ooac@meta.data[ooac@meta.data$condition == 'EM/EM', 'condition_new'] <- 'EM'
ooac@meta.data[ooac@meta.data$condition == 'EM/DM', 'condition_new'] <- 'EM-DM'
ooac@meta.data[ooac@meta.data$condition == 'DM/DM', 'condition_new'] <- 'DM'
ooac@meta.data[ooac@meta.data$condition == 'EM/DM+IFNb', 'condition_new'] <- 'EM-DM + IFN-β'
ooac@meta.data[ooac@meta.data$condition == 'EM/DM+IFNy', 'condition_new'] <- 'EM-DM + IFN-γ'

#Function to assign cellular compartments
get_average_expression_per_group <- function(seurat_object, metadata_column, genes, use_sct=F){
  # initialize the table
  expression_table <- NULL
  # check each group
  for(group in unique(seurat_object@meta.data[[metadata_column]])){
    # we won't do NA, obviously
    if(!is.na(group)){
      # subset to that group
      seurat_object_group <- seurat_object[, !is.na(seurat_object@meta.data[[metadata_column]]) & seurat_object@meta.data[[metadata_column]] == group]
      # put per gene in a list
      exp_per_list <- list()
      # check each gene
      for(gene in genes){
        # get the mean expression
        mean_expression <- NULL
        # depending on the assay
        if(use_sct){
          mean_expression <- mean(as.vector(unlist(seurat_object_group@assays$SCT@counts[gene, ])))
        }
        else{
          mean_expression <- mean(as.vector(unlist(seurat_object_group@assays$RNA@data[gene, ])))
        }
        # put in the list
        exp_per_list[[gene]] <- mean_expression
      }
      # turn into a row
      exp_row <- data.frame(exp_per_list)
      # add the group
      exp_row[['group']] <- group
      # set the correct order
      exp_row <- exp_row[, c('group', genes)]
      # add to the big matrix
      if(is.null(expression_table)){
        expression_table <- exp_row
      }
      else{
        expression_table <- rbind(expression_table, exp_row)
      }
    }
  }
  return(expression_table)
}


# get the mean expression of marker genes per cluster
avg_exp_cluster <- get_average_expression_per_group(ooac, 'seurat_clusters', c('EPCAM', 'VIM'))
# use same assay as the avg exp one
DefaultAssay(ooac) <- 'RNA'
# set the threshold for the markers
ooac@meta.data$epcam_positive <- NA
ooac@meta.data[ooac@meta.data[['seurat_clusters']] %in% avg_exp_cluster[avg_exp_cluster[['EPCAM']] >= 1.0, 'group'], 'epcam_positive'] <- T
ooac@meta.data[ooac@meta.data[['seurat_clusters']] %in% avg_exp_cluster[avg_exp_cluster[['EPCAM']] < 1.0, 'group'], 'epcam_positive'] <- F
ooac@meta.data$vim_positive <- NA
ooac@meta.data[ooac@meta.data[['seurat_clusters']] %in% avg_exp_cluster[avg_exp_cluster[['VIM']] >= 2.7, 'group'], 'vim_positive'] <- T
ooac@meta.data[ooac@meta.data[['seurat_clusters']] %in% avg_exp_cluster[avg_exp_cluster[['VIM']] < 2.7, 'group'], 'vim_positive'] <- F

# add the compartment based on marker expression
ooac@meta.data$compartment <- 'none'
ooac@meta.data[ooac@meta.data$epcam_positive == T, 'compartment'] <- 'epithelial'
ooac@meta.data[ooac@meta.data$vim_positive == T, 'compartment'] <- 'mesenchymal'

# improve compartment names
ooac@meta.data$compartment_new2 <- 'none'
ooac@meta.data[ooac@meta.data$compartment == 'epithelial', 'compartment_new'] <- 'Epithelial cell'
ooac@meta.data[ooac@meta.data$compartment == 'mesenchymal', 'compartment_new'] <- 'Mesenchymal/Neural cell'

#set Default assay back to ‘SCT’ for further analysis
DefaultAssay(ooac) <- 'SCT'

# Assign the cellular detection rate (CDR)
cdr <- scale(colSums(GetAssayData(ooac, assay = "SCT", slot = "data") > 0))
ooac$CDR <- cdr


#### Subset the data object in object containing the media conditions (EM, EM-DM, DM) and and object containing the interferon stimulation conditions (EM-DM, EM-DM+IFNβ, EM-DM+IFNγ)

#generate ‘media’ object
ooac_media <- ooac[, ooac@meta.data[['condition']] %in% c('EM/EM', 'EM/DM', 'DM/DM')]

#Remove the 1 cell remaining in epithelial_IFNb cluster (from EM-DM data)
ooac_media <- ooac_media[, !(ooac_media@meta.data[['annotation_res0.34']] %in% c('epithelial_IFNb'))]

#generate ‘IFN’ object
ooac_IFN <- ooac[, ooac@meta.data[['condition']] %in% c('EM/DM', 'EM/DM+IFNb', 'EM/DM+IFNy')]

#These objects were saved ('intestine_on_chip_media.rds' and 'intestine_on_chip_IFN.rds') and used to generate the Figures as outlined in the other scripts

