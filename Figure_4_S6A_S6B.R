###########################################################################################################################
### Figure 4 ##############################################################################################################
###########################################################################################################################

#libraries
library(Seurat)
library(ggplot2)
library(SeuratData)
library(dplyr)
library(speckle)
library(limma)
library(scater)
library(patchwork)
library(edgeR)
library(statmod)
library(tidyverse)
library(viridisLite)
library(viridis)
library(RColorBrewer)

#Download the Gut Cell Atlas data here: https://www.gutcellatlas.org
#Use: 'Full single cell RNA-seq dataset of 428K intestinal cells from fetal, pediatric, adult donors, and up to 11 intestinal regions.'
#Described in Elmentaite, R., Kumasaka, N., Roberts, K. et al. Cells of the human intestinal tract mapped across space and time. Nature 597, 250–255 (2021). https://doi.org/10.1038/s41586-021-03852-1
objects_loc <- #path to folder where the Gut Cell Atlas data is saved
object_loc_all_sct_clus <- paste(objects_loc, 'elmentaite_2021_all_sct_clus.rds', sep = '')

# read complete 'Elmentaite object' (how we have referred to the Gut Cell Atlas dataset)
elmentaite_2021 <- readRDS(object_loc_all_sct_clus)
DefaultAssay(elmentaite_2021) <- 'SCT'

#Subset the data to remove IBD samples 
elmentaite_2021$Diagnosis <- as.character(elmentaite_2021$Diagnosis)
elmentaite_2021 <- elmentaite_2021[, !(elmentaite_2021@meta.data[['Diagnosis']] %in% c('Pediatric Crohn Disease'))]

# locations of intestine-on-chip objects
seurat_objects_loc <- #path to folder where the intestine-on-chip data of the media conditions is saved
seurat_object_media <- paste(seurat_objects_loc, 'intestine_on_chip_media.rds', sep = '')
# read the intestine-on-chip object
ooac <- readRDS(seurat_object_media)
#Set the default assay to SCT
DefaultAssay(ooac) <- 'SCT'

# find transfer anchors
anchors <- FindTransferAnchors(
  reference = elmentaite_2021,
  query = ooac,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:50
)
# for some reason the UMAP did not stick
#ooac <- RunUMAP(ooac, dims=1:30, reduction='pca', return.model = T)
elmentaite_2021 <- RunUMAP(elmentaite_2021, dims=1:30, reduction='pca', return.model = T)
# do the reference mapping, generating a new object 'ooac_elmentaite_2021'
ooac_elmentaite_2021 <- MapQuery(
  anchorset = anchors,
  query = ooac,
  reference = elmentaite_2021,
  refdata = list(
    region.elmentaite.2021 = "Region",
    celltype.elmentaite.2021 = "cell_type",
    category.elmentaite.2021 = "category",
    age.group.elmentaite.2021 = 'Age_group',
    region.code.elmentaite.2021 = 'Region.code'
  ),
  reference.reduction = "pca",
  reduction.model = "umap"
)

#transfer the predictions from the ooac_elmentaite_2021 object to ooac object
ooac <- AddMetaData(ooac, ooac_elmentaite_2021[, colnames(ooac)]@meta.data['predicted.region.elmentaite.2021'], 'predicted.region.elmentaite.2021.full')
ooac <- AddMetaData(ooac, ooac_elmentaite_2021[, colnames(ooac)]@meta.data['predicted.region.elmentaite.2021.score'], 'predicted.region.elmentaite.2021.full.score')
ooac <- AddMetaData(ooac, ooac_elmentaite_2021[, colnames(ooac)]@meta.data['predicted.celltype.elmentaite.2021'], 'predicted.celltype.elmentaite.2021.full')
ooac <- AddMetaData(ooac, ooac_elmentaite_2021[, colnames(ooac)]@meta.data['predicted.celltype.elmentaite.2021.score'], 'predicted.celltype.elmentaite.2021.full.score')
ooac <- AddMetaData(ooac, ooac_elmentaite_2021[, colnames(ooac)]@meta.data['predicted.category.elmentaite.2021'], 'predicted.category.elmentaite.2021.full')
ooac <- AddMetaData(ooac, ooac_elmentaite_2021[, colnames(ooac)]@meta.data['predicted.category.elmentaite.2021.score'], 'predicted.category.elmentaite.2021.full.score')
ooac <- AddMetaData(ooac, ooac_elmentaite_2021[, colnames(ooac)]@meta.data['predicted.age.group.elmentaite.2021'], 'predicted.age.group.elmentaite.2021.full')
ooac <- AddMetaData(ooac, ooac_elmentaite_2021[, colnames(ooac)]@meta.data['predicted.age.group.elmentaite.2021.score'], 'predicted.age.group.elmentaite.2021.full.score')
ooac <- AddMetaData(ooac, ooac_elmentaite_2021[, colnames(ooac)]@meta.data['predicted.region.code.elmentaite.2021'], 'predicted.region.code.elmentaite.2021.full')
ooac <- AddMetaData(ooac, ooac_elmentaite_2021[, colnames(ooac)]@meta.data['predicted.region.code.elmentaite.2021.score'], 'predicted.region.code.elmentaite.2021.full.score')

#transfer reductions
ooac@reductions$ref.umap.elmentaite.2021.all <- ooac_elmentaite_2021[, colnames(ooac)]@reductions$ref.umap

#Plot the predicted variables
order <- c('EM',
           'EM-DM',
           'DM')

ooac$condition_new <- factor(ooac$condition_new, levels = order)

#Figure 4B: Region
ooac@meta.data[ooac@meta.data$predicted.region.elmentaite.2021.full == 'LargeInt', 'predicted.region.elmentaite.2021.full'] <- 'Large intestine'
ooac@meta.data[ooac@meta.data$predicted.region.elmentaite.2021.full == 'SmallInt', 'predicted.region.elmentaite.2021.full'] <- 'Small intestine'

order <- c('Small intestine',
           'Large intestine')

ooac$predicted.region.elmentaite.2021.full <- as.character(ooac$predicted.region.elmentaite.2021.full)
ooac$predicted.region.elmentaite.2021.full <- factor(ooac$predicted.region.elmentaite.2021.full, levels = order)

DimPlot(ooac, reduction = 'umap', group.by = 'predicted.region.elmentaite.2021.full', cols=alpha(c("#1B9E77", "#D95F02"))) +
  labs(title="Predicted intestinal region")+
  theme(plot.title = element_text(hjust = 0.5, size=15), axis.title=element_text(size=12), axis.text=element_text(size=10))

#Figure 4A: Age group
ooac@meta.data[ooac@meta.data$predicted.age.group.elmentaite.2021.full == 'First trim', 'predicted.age.group.elmentaite.2021.full'] <- 'First trimester'
ooac@meta.data[ooac@meta.data$predicted.age.group.elmentaite.2021.full == 'Second trim', 'predicted.age.group.elmentaite.2021.full'] <- 'Second trimester'

order <- c('First trimester',
           'Second trimester',
           'Pediatric',
           'Adult')

ooac$predicted.age.group.elmentaite.2021.full <- factor(ooac$predicted.age.group.elmentaite.2021.full, levels = order)

DimPlot(ooac, reduction = 'umap', group.by = 'predicted.age.group.elmentaite.2021.full', cols=alpha(c("#8AB5D3","orange", "#8CC34B",  "#FB8072"))) +
  labs(title="Predicted age")+
  theme(plot.title = element_text(hjust = 0.5, size=15), axis.title=element_text(size=12), axis.text=element_text(size=10))

#Figure S6B: Category
DimPlot(ooac, reduction = 'umap', group.by = 'predicted.category.elmentaite.2021.full', cols=alpha(c("#FC8D62", "#8DA0CB", "#8CC34B"))) +
  labs(title="Predicted cell type")+
  theme(plot.title = element_text(hjust = 0.5, size=15), axis.title=element_text(size=12), axis.text=element_text(size=10))

#Figure 4C: Cell type
#Simplify predicted cell type classification based on the description in the publication of the Gut Cell Atlas (Elmentaite, R., Kumasaka, N., Roberts, K. et al. Cells of the human intestinal tract mapped across space and time. Nature 597, 250–255 (2021). https://doi.org/10.1038/s41586-021-03852-1)
ooac$predicted.celltype.elmentaite.2021.full.2 = ooac$predicted.celltype.elmentaite.2021.full
ooac$predicted.celltype.elmentaite.2021.full.2[ooac$predicted.celltype.elmentaite.2021.full.2 == 'D cells (SST+)'] <- 'Enteroendocrine cell'
ooac$predicted.celltype.elmentaite.2021.full.2[ooac$predicted.celltype.elmentaite.2021.full.2 == 'EC cells (TAC1+)'] <- 'Enteroendocrine cell'
ooac$predicted.celltype.elmentaite.2021.full.2[ooac$predicted.celltype.elmentaite.2021.full.2 == 'EECs'] <- 'Enteroendocrine cell'
ooac$predicted.celltype.elmentaite.2021.full.2[ooac$predicted.celltype.elmentaite.2021.full.2 == 'I cells (CCK+)'] <- 'Enteroendocrine cell'
ooac$predicted.celltype.elmentaite.2021.full.2[ooac$predicted.celltype.elmentaite.2021.full.2 == 'M/X cells (MLN/GHRL+)'] <- 'Enteroendocrine cell'
ooac$predicted.celltype.elmentaite.2021.full.2[ooac$predicted.celltype.elmentaite.2021.full.2 == 'Progenitor (NEUROG3+)'] <- 'Enteroendocrine cell'
#specify proximal and distal progenitors as 'fetal'
ooac$predicted.celltype.elmentaite.2021.full.2[ooac$predicted.celltype.elmentaite.2021.full.2 == 'Proximal progenitor'] <- 'Fetal proximal progenitor'
ooac$predicted.celltype.elmentaite.2021.full.2[ooac$predicted.celltype.elmentaite.2021.full.2 == 'Distal progenitor'] <- 'Fetal distal progenitor'
#BEST4+ epithelial cells 
ooac$predicted.celltype.elmentaite.2021.full.2[ooac$predicted.celltype.elmentaite.2021.full.2 == 'BEST4+ epithelial'] <- 'BEST4+ epithelial cell'
#Paneth 
ooac$predicted.celltype.elmentaite.2021.full.2[ooac$predicted.celltype.elmentaite.2021.full.2 == 'Paneth'] <- 'Paneth cell'
#group Branch A1 (iMN), Branch A2 (IPAN/IN), Branch B1 (eMN) as 'neuronal branches', note only fetal samples have neuronal compartment in elmentaite dataset
ooac$predicted.celltype.elmentaite.2021.full.2[ooac$predicted.celltype.elmentaite.2021.full.2 == 'Branch A1 (iMN)'] <- 'Enteric neuron'
ooac$predicted.celltype.elmentaite.2021.full.2[ooac$predicted.celltype.elmentaite.2021.full.2 == 'Branch A2 (IPAN/IN)'] <- 'Enteric neuron'
ooac$predicted.celltype.elmentaite.2021.full.2[ooac$predicted.celltype.elmentaite.2021.full.2 == 'Branch B1 (eMN)'] <- 'Enteric neuron'
#group Stromal 1, 2, 3 as 'stromal subtypes'
ooac$predicted.celltype.elmentaite.2021.full.2[ooac$predicted.celltype.elmentaite.2021.full.2 == 'Stromal 1 (ADAMDEC1+)'] <- 'Stromal subtype'
ooac$predicted.celltype.elmentaite.2021.full.2[ooac$predicted.celltype.elmentaite.2021.full.2 == 'Stromal 2 (CH25H+)'] <- 'Stromal subtype'
ooac$predicted.celltype.elmentaite.2021.full.2[ooac$predicted.celltype.elmentaite.2021.full.2 == 'Stromal 2 (NPY+)'] <- 'Stromal subtype'
ooac$predicted.celltype.elmentaite.2021.full.2[ooac$predicted.celltype.elmentaite.2021.full.2 == 'Stromal 3 (C7+)'] <- 'Stromal subtype'
ooac$predicted.celltype.elmentaite.2021.full.2[ooac$predicted.celltype.elmentaite.2021.full.2 == 'Stromal 3 (KCNN3+)'] <- 'Stromal subtype'
#write ICC and SMC in full, singular form of stem cell
ooac$predicted.celltype.elmentaite.2021.full.2[ooac$predicted.celltype.elmentaite.2021.full.2 == 'ICC'] <- 'Interstitial cell of Cajal'
ooac$predicted.celltype.elmentaite.2021.full.2[ooac$predicted.celltype.elmentaite.2021.full.2 == 'SMC (PLPP2+)'] <- 'Smooth muscle cell'
ooac$predicted.celltype.elmentaite.2021.full.2[ooac$predicted.celltype.elmentaite.2021.full.2 == 'Stem cells'] <- 'Stem cell'
#Follicular dendritic cells
ooac$predicted.celltype.elmentaite.2021.full.2[ooac$predicted.celltype.elmentaite.2021.full.2 == 'FDC'] <- 'Lymph node fibroblast'
#Mesothelium
ooac$predicted.celltype.elmentaite.2021.full.2[ooac$predicted.celltype.elmentaite.2021.full.2 == 'Mesothelium (RGS5+)'] <- 'Mesothelium'
#Mesoderm
ooac$predicted.celltype.elmentaite.2021.full.2[ooac$predicted.celltype.elmentaite.2021.full.2 == 'Mesoderm 1 (HAND1+)'] <- 'Mesoderm'
#Capital letters
ooac$predicted.celltype.elmentaite.2021.full.2[ooac$predicted.celltype.elmentaite.2021.full.2 == 'cycling ENCC/glia'] <- 'Cycling ENCC/glia'
ooac$predicted.celltype.elmentaite.2021.full.2[ooac$predicted.celltype.elmentaite.2021.full.2 == 'cycling neuroblast'] <- 'Cycling neuroblast'
ooac$predicted.celltype.elmentaite.2021.full.2[ooac$predicted.celltype.elmentaite.2021.full.2 == 'cycling stromal'] <- 'Cycling stromal cell'
ooac$predicted.celltype.elmentaite.2021.full.2[ooac$predicted.celltype.elmentaite.2021.full.2 == 'ENCC/glia Progenitor'] <- 'ENCC/glia progenitor'

#Plot cell type annotation
order <- c('Stem cell',
           'TA',
           'Enterocyte',
           'Colonocyte',
           'Paneth cell',
           'Goblet cell',
           'Enteroendocrine cell',
           'Microfold cell',
           'Fetal proximal progenitor',
           'Fetal distal progenitor',
           'Cycling ENCC/glia',
           'Cycling stromal cell',
           'Cycling neuroblast',
           'Mesoderm',
           'Stromal subtype',
           'Smooth muscle cell',
           'Enteric neuron',
           'ENCC/glia progenitor',
           'Neuroblast',
           'Interstitial cell of Cajal',
           'BEST4+ epithelial cell',
           'Lymph node fibroblast',
           'Mesothelium')

ooac$predicted.celltype.elmentaite.2021.full.2 <- as.character(ooac$predicted.celltype.elmentaite.2021.full.2)
ooac$predicted.celltype.elmentaite.2021.full.2 <- factor(ooac$predicted.celltype.elmentaite.2021.full.2, levels = order)

#Only include cell types with more than 10 cells predicted - exclude BEST4+ epithelial cell, Lymph node fibroblast and Mesothelium
DimPlot(ooac, reduction = 'umap', group.by = 'predicted.celltype.elmentaite.2021.full.2',
        cols = alpha(c('Cycling stromal cell' = '#1B9E77', 'Cycling neuroblast' = '#FF7F00', 'Cycling ENCC/glia' = '#A6CEE3', 'Mesoderm' = '#E7298A', 'Stromal subtype' = '#8CBF30', 'Smooth muscle cell' = '#E6AB02', 'Enteric neuron' = '#FF9999', 'ENCC/glia progenitor' = '#AA0000', 'Interstitial cell of Cajal' = '#9467BD', 'Neuroblast' = '#1F78B4'))) +
  labs(title="Predicted mesenchymal and neural subtypes")+
  theme(plot.title = element_text(hjust = 0.5, size=13), axis.title=element_text(size=12), axis.text=element_text(size=10))

DimPlot(ooac, reduction = 'umap', group.by = 'predicted.celltype.elmentaite.2021.full.2',
        cols = alpha(c('Stem cell' = '#A6CEE3', 'TA' = '#1F78B4', 'Enterocyte' = '#B2DF8A', 'Colonocyte' = '#33A02C', 'Paneth cell' = '#E31A1C', 'Microfold cell' = '#FF7F00',  'Goblet cell' = '#FF8080', 'Enteroendocrine cell' = '#FDBF6F', 'Fetal distal progenitor' = '#6A3D9A', 'Fetal proximal progenitor' = '#CAB2D6'))) +
  labs(title="Predicted epithelial subtypes")+
  theme(plot.title = element_text(hjust = 0.5, size=13), axis.title=element_text(size=12), axis.text=element_text(size=10))

#Figure S6A: Human fetal endoderm atlas projection
#Download the required data to run the scoreHIO function here: https://github.com/Camp-Lab/scoreHIO
library(uwot)
library(RANN)
library(quadprog)
library(dplyr)
library(beanplot)
library(scoreHIO)

ooac <- SetIdent(ooac, value = "annotation_res0.34_new2")

ooac_fetal_endoderm_atlas <- score_fidelity(
  que_obj = ooac,
  organ_ref_dir = #path to the folder where the Human Fetal Endoderm Atlas data is saved,
  group_by = "orig.ident"
)

#Transfer prediction from ooac_fetal_endoderm_atlas to ooac object
ooac <- AddMetaData(ooac, ooac_fetal_endoderm_atlas[, colnames(ooac)]@meta.data['Mapped_fetal_organ'], 'Mapped_fetal_organ')

DimPlot(ooac, reduction = 'umap', group.by = 'Mapped_fetal_organ', cols=alpha(c('Esophagus'="#6CAEDF", 'Intestine'="#976BB6", 'Lung'="#85C567", 'Pancreas'="#E66161", 'Stomach'="#FFA500"))) +
  labs(title="Predicted fetal endodermal organ")+
  theme(plot.title = element_text(hjust = 0.5, size=14), axis.title=element_text(size=12), axis.text=element_text(size=10))

#Generate barplots of the predicted variables
n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
my_cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

ooac@meta.data <- ooac@meta.data %>%
  mutate(condition_cell_line = paste(condition, cell_line, sep = '_'))

ooac_epithelial <- ooac[, ooac@meta.data[['compartment_new2']] %in% c('Epithelial cell')]
ooac_mesenchymal <- ooac[, ooac@meta.data[['compartment_new2']] %in% c('Mesenchymal/Neural cell')]

order_fetal_organ1 <- c('Stomach',
                        'Pancreas',
                        'Lung',
                        'Intestine')

order_fetal_organ2 <- c('Stomach',
                        'Pancreas',
                        'Lung',
                        'Intestine',
                        'Esophagus')

order_category <- c('Neuronal',
                    'Mesenchymal',
                    'Epithelial')

order_region <- c('Large intestine',
                  'Small intestine')

order_age <- c('Adult',
               'Pediatric',
               'Second trimester',
               'First trimester')


par(mar=c(5.1, 4.1, 4.1, 15.1), xpd=TRUE)

#Figure S6A: Fetal organ, separate epithelial and mesenchymal/neural cells
props <- propeller(clusters=ooac_epithelial$Mapped_fetal_organ, sample=ooac_epithelial$condition_cell_line, group=ooac_epithelial$condition, transform = 'logit')
props_subset <- subset(props, select=c('PropMean.EM.EM', 'PropMean.EM.DM', 'PropMean.DM.DM'))
props_subset <- props_subset[match(order_fetal_organ1, rownames(props_subset)), ]
props_subset_M <- data.matrix(props_subset)
colnames(props_subset_M) <- c('EM','EM-DM','DM')
barplot(props_subset_M, legend=FALSE, ylab="Proportions", col =alpha(c("#FFA500", "#E66161", "#85C567", "#976BB6")), cex.lab=1.2, cex.names=1.2, args.legend = list(x = "topright", inset = c(-0.3, 0.0), cex = 0.9), main="Predicted fetal organ - epithelial cells")

props <- propeller(clusters=ooac_mesenchymal$Mapped_fetal_organ, sample=ooac_mesenchymal$condition_cell_line, group=ooac_mesenchymal$condition, transform = 'logit')
props_subset <- subset(props, select=c('PropMean.EM.EM', 'PropMean.EM.DM', 'PropMean.DM.DM'))
props_subset <- props_subset[match(order_fetal_organ2, rownames(props_subset)), ]
props_subset_M <- data.matrix(props_subset)
colnames(props_subset_M) <- c('EM','EM-DM','DM')
barplot(props_subset_M, legend=FALSE, ylab="Proportions", col =alpha(c("#FFA500","#E66161","#85C567","#976BB6","#6CAEDF")), cex.lab=1.2, cex.names=1.2, args.legend = list(x = "topright", inset = c(-0.3, 0.0), cex = 0.9), main="Predicted fetal organ - mesenchymal and neural cells")

#Figure S6B: category, separate epithelial and mesenchymal/neural cells
props <- propeller(clusters=ooac_epithelial$predicted.category.elmentaite.2021.full, sample=ooac_epithelial$condition_cell_line, group=ooac_epithelial$condition, transform = 'logit')
props_subset <- subset(props, select=c('PropMean.EM.EM', 'PropMean.EM.DM', 'PropMean.DM.DM'))
props_subset <- props_subset[match(order_category, rownames(props_subset)), ]
props_subset_M <- data.matrix(props_subset)
colnames(props_subset_M) <- c('EM','EM-DM','DM')
barplot(props_subset_M, legend=FALSE, ylab="Proportions", col = alpha(c("#8CC34B", "#8DA0CB", "#FC8D62")), cex.lab=1.2, cex.names=1.2, args.legend = list(x = "topright", inset = c(-0.3, 0.0), cex = 0.9), main="Predicted cell type - epithelial cells")

props <- propeller(clusters=ooac_mesenchymal$predicted.category.elmentaite.2021.full, sample=ooac_mesenchymal$condition_cell_line, group=ooac_mesenchymal$condition, transform = 'logit')
props_subset <- subset(props, select=c('PropMean.EM.EM', 'PropMean.EM.DM', 'PropMean.DM.DM'))
props_subset <- props_subset[match(order_category, rownames(props_subset)), ]
props_subset_M <- data.matrix(props_subset)
colnames(props_subset_M) <- c('EM','EM-DM','DM')
barplot(props_subset_M, legend=FALSE, ylab="Proportions", col = alpha(c("#8CC34B", "#8DA0CB", "#FC8D62")), cex.lab=1.2, cex.names=1.2, args.legend = list(x = "topright", inset = c(-0.3, 0.0), cex = 0.9), main="Predicted cell type - mesenchymal and neural cells")

par(mar=c(5.1, 4.1, 4.1, 5.1), xpd=TRUE)
#Figure 4B: Region
props <- propeller(clusters=ooac$predicted.region.elmentaite.2021.full, sample=ooac$condition_cell_line, group=ooac$condition, transform = 'logit')
props_subset <- subset(props, select=c('PropMean.EM.EM', 'PropMean.EM.DM', 'PropMean.DM.DM'))
props_subset <- props_subset[match(order_region, rownames(props_subset)), ]
props_subset_M <- data.matrix(props_subset)
colnames(props_subset_M) <- c('EM','EM-DM','DM')
barplot(props_subset_M, legend=FALSE, ylab="Proportions", col = alpha(c("#D95F02", "#1B9E77")), cex.lab=1.2, cex.names=1.2, args.legend = list(x = "bottom", inset = c(0.0, -0.3), cex = 0.9), main="Predicted intestinal region") 

#Figure 4A: Age group
props <- propeller(clusters=ooac$predicted.age.group.elmentaite.2021.full, sample=ooac$condition_cell_line, group=ooac$condition, transform = 'logit')
props_subset <- subset(props, select=c('PropMean.EM.EM', 'PropMean.EM.DM', 'PropMean.DM.DM'))
props_subset <- props_subset[match(order_age, rownames(props_subset)), ]
props_subset_M <- data.matrix(props_subset)
colnames(props_subset_M) <- c('EM','EM-DM','DM')
barplot(props_subset_M, legend=FALSE, ylab="Proportions", col = alpha(c("#FB8072", "#8CC34B", "orange", "#8AB5D3")), cex.lab=1.2, cex.names=1.2, args.legend = list(x = "bottom", inset = c(0.0, -0.3), cex = 0.9), main="Predicted age") 

#Figure 4D and 4I: Correlation plot: comparing intestine-on-chip and Elmentaite dataset
#Subset Gut Cell Atlas object for the epithelial, mesenchymal, neural category
elmentaite_2021 <- elmentaite_2021[, elmentaite_2021@meta.data[['category']] %in% c('Epithelial', 'Mesenchymal', 'Neuronal')]
#Set category as character
elmentaite_2021$category <- as.character(elmentaite_2021$category)

#Normalize data using the RNA assay
elmentaite_2021 <- NormalizeData(elmentaite_2021, assay = "RNA")

#set these factor vectors to character vectors 
elmentaite_2021$cell_type <- as.character(elmentaite_2021$cell_type)
elmentaite_2021$category <- as.character(elmentaite_2021$category)
elmentaite_2021$Age_group <- as.character(elmentaite_2021$Age_group)
elmentaite_2021$Region <- as.character(elmentaite_2021$Region)

#Rename
elmentaite_2021$Age_group[elmentaite_2021$Age_group == 'First trim'] <- 'First trimester'
elmentaite_2021$Age_group[elmentaite_2021$Age_group == 'Second trim'] <- 'Second trimester'
elmentaite_2021$Region[elmentaite_2021$Region == 'LargeInt'] <- 'Large intestine'
elmentaite_2021$Region[elmentaite_2021$Region == 'SmallInt'] <- 'Small intestine'

#Simplify cell type classification based on the description in the publication of the Gut Cell Atlas (Elmentaite, R., Kumasaka, N., Roberts, K. et al. Cells of the human intestinal tract mapped across space and time. Nature 597, 250–255 (2021). https://doi.org/10.1038/s41586-021-03852-1)
elmentaite_2021$cell_type_2 = elmentaite_2021$cell_type
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'D cells (SST+)'] <- 'Enteroendocrine cell'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'EC cells (NPW+)'] <- 'Enteroendocrine cell'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'EC cells (TAC1+)'] <- 'Enteroendocrine cell'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'EECs'] <- 'Enteroendocrine cell'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'I cells (CCK+)'] <- 'Enteroendocrine cell'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'K cells (GIP+)'] <- 'Enteroendocrine cell'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'L cells (PYY+)'] <- 'Enteroendocrine cell'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'M/X cells (MLN/GHRL+)'] <- 'Enteroendocrine cell'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'N cells (NTS+)'] <- 'Enteroendocrine cell'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'β cells (INS+)'] <- 'Enteroendocrine cell'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Progenitor (NEUROG3+)'] <- 'Enteroendocrine cell'
#group proximal and distal progenitors as 'fetal progenitors'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Proximal progenitor'] <- 'Fetal proximal progenitor'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Distal progenitor'] <- 'Fetal distal progenitor'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'CLDN10+ cells'] <- 'Fetal CLDN10+ cell'
#group BEST2+ goblet cells and goblet cells
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'BEST2+ Goblet cell'] <- 'Goblet cell'
#group BEST4+ epithelial cells and enterocytes
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'BEST4+ epithelial'] <- 'BEST4+ epithelial cell'
#singular form of stem cell
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Stem cells'] <- 'Stem cell'
#Tuft cell
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Tuft'] <- 'Tuft cell'
#Paneth cell
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Paneth'] <- 'Paneth cell'
#group Branch A1 (iMN), Branch A2 (IPAN/IN), Branch B1 (eMN) as 'neuronal branches'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Branch A1 (iMN)'] <- 'Enteric neuron'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Branch A2 (IPAN/IN)'] <- 'Enteric neuron'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Branch A3 (IPAN/IN)'] <- 'Enteric neuron'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Branch A4 (IN)'] <- 'Enteric neuron'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Branch B1 (eMN)'] <- 'Enteric neuron'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Branch B2 (eMN)'] <- 'Enteric neuron'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Branch B3 (IPAN)'] <- 'Enteric neuron'
#group glial cells
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Glia 1 (DHH+)'] <- 'Glial cell'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Glia 2 (ELN+)'] <- 'Glial cell'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Glia 3 (BCAN+)'] <- 'Glial cell'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Differentiating glia'] <- 'Glial cell'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Adult Glia'] <- 'Glial cell'
#group Stromal 1, 2, 3, 4 as 'stromal subtypes'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Stromal 1 (ADAMDEC1+)'] <- 'Stromal subtype'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Stromal 1 (CCL11+)'] <- 'Stromal subtype'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Stromal 2 (CH25H+)'] <- 'Stromal subtype'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Stromal 2 (NPY+)'] <- 'Stromal subtype'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Stromal 3 (C7+)'] <- 'Stromal subtype'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Stromal 3 (KCNN3+)'] <- 'Stromal subtype'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Transitional Stromal 3 (C3+)'] <- 'Stromal subtype'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Stromal 4 (MMP1+)'] <- 'Stromal subtype'
#write ICC in full
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'ICC'] <- 'Interstitial cell of Cajal'
#group smooth muscle cell populations 
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'SMC (PLPP2+)'] <- 'Smooth muscle cell'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'SMC (PART1/CAPN3+)'] <- 'Smooth muscle cell'
#group mesoderm populations 
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Mesoderm 1 (HAND1+)'] <- 'Mesoderm'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Mesoderm 2 (ZEB2+)'] <- 'Mesoderm'
#group mesothelium populations 
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Mesothelium'] <- 'Mesothelium'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Mesothelium (PRG4+)'] <- 'Mesothelium'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Mesothelium (RGS5+)'] <- 'Mesothelium'
#group (mesenteric) lymph node immune-organizing fibroblasts
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'mLN Stroma (FMO2+)'] <- 'Lymph node fibroblast'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'T reticular'] <- 'Lymph node fibroblast'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'mLTo'] <- 'Lymph node fibroblast'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'FDC'] <- 'Lymph node fibroblast'
#group myofibroblasts
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'myofibroblast'] <- 'Myofibroblast'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'myofibroblast (RSPO2+)'] <- 'Myofibroblast'
#group pericytes
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'angiogenic pericyte'] <- 'Pericyte'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Contractile pericyte (PLN+)'] <- 'Pericyte'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Immature pericyte'] <- 'Pericyte'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'Pericyte'] <- 'Pericyte'
#Capital letters
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'cycling ENCC/glia'] <- 'Cycling ENCC/glia'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'cycling neuroblast'] <- 'Cycling neuroblast'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'cycling stromal'] <- 'Cycling stromal cell'
elmentaite_2021$cell_type_2[elmentaite_2021$cell_type_2 == 'ENCC/glia Progenitor'] <- 'ENCC/glia progenitor'

#Age_group_2 annotation, fetal samples (of first and second trimester) that are not enriched for Epcam are pooled as 'fetal'
elmentaite_2021$Age_group_2 = elmentaite_2021$Age_group
elmentaite_2021@meta.data[elmentaite_2021@meta.data$Sample.name == 'BRC2026', 'Age_group_2'] <- 'Fetal'
elmentaite_2021@meta.data[elmentaite_2021@meta.data$Sample.name == 'BRC2029', 'Age_group_2'] <- 'Fetal'
elmentaite_2021@meta.data[elmentaite_2021@meta.data$Sample.name == 'BRC2258', 'Age_group_2'] <- 'Fetal'
elmentaite_2021@meta.data[elmentaite_2021@meta.data$Sample.name == 'BRC2259', 'Age_group_2'] <- 'Fetal'
elmentaite_2021@meta.data[elmentaite_2021@meta.data$Age_group == 'Second trimester', 'Age_group_2'] <- 'Fetal'

#subset Elmentaite dataset
#Dataset for Figure 4D 
elmentaite_2021_SI <- elmentaite_2021[, elmentaite_2021@meta.data[['Region']] %in% c('Small intestine')]
#Dataset for Figure 4I (other fetal and pediatric samples were enriched for Epcam, leading to skewed average gene expression, so are not used in this analysis)
elmentaite_2021_SI_age <- elmentaite_2021_SI[, elmentaite_2021_SI@meta.data[['Age_group_2']] %in% c('Adult', 'Fetal')]

#perform subsetting intestine-on-chip object
ooac_epithelial <- ooac[, ooac@meta.data[['compartment_new']] %in% c('Epithelial cell')]
ooac_mesenchymal <- ooac[, ooac@meta.data[['compartment_new']] %in% c('Mesenchymal cell')]

#Generate the DEG list
#Figure 4D:
ooac <- SetIdent(ooac, value = "annotation_res0.34_simple2")
cluster.markers_clusters_simple2 <- FindAllMarkers(ooac, min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST", latent.vars = c("cell_line", "CDR"))
DE_genes <- cluster.markers_clusters_simple2
#libraries
library(org.Hs.eg.db)
library(msigdbr)
#Function: Create column in DE gene object that indicates the direction of log2FC
add.down_or_up.column <- function(df){
  df$Direction <- "Upregulated"
  df$Direction[df$avg_log2FC < 0] <- "Downregulated"
  df$Direction <- as.factor(df$Direction)
  df$Direction <- factor(df$Direction, levels = c("Upregulated", "Downregulated"))
  return(df)
}
#Function: Create column in DE gene object that indicates the absolute value of log2FC, to be able to set threshold for both up- and downregulated genes
add.abs.log2FC <- function(df){
  df$abs_log2FC <- abs(df$avg_log2FC)
  return(df)
}
#Duplicate the ‘gene’ column and name ‘SYMBOL’
DE_genes$SYMBOL <- DE_genes$gene
#Add the up or down annotation for each log2FC
DE_genes <- add.down_or_up.column(DE_genes)
#Add column with the absolute value of log2FC to include up- and down-regulated genes:
DE_genes <- add.abs.log2FC(DE_genes)
#Filter the significant genes. Use p < 0.05 or p < 0.01. 
filtered_DE_genes <- DE_genes[DE_genes$p_val_adj < 0.01,]
filtered_log2FC_DE_genes <- filtered_DE_genes[filtered_DE_genes$avg_log2FC > 0.5,]

#Figure 4I:
ooac <- SetIdent(ooac, value = "compartment_new")
compartments <- list(complete = NULL)
#Function1 for EM vs EM-DM
DEA.subset.sc1 <- function(subset_test){
  subset_result <- FindMarkers(ooac, ident.1 = "EM", ident.2 = "EM-DM", group.by = "condition_new", subset.ident = subset_test, test.use = "MAST", min.pct = 0.1, logfc.threshold = 0.25, latent.vars = c("cell_line", "CDR"))
  return(subset_result)
}
cluster.markers_complete_EMvsEM_DM <- DEA.subset.sc1(compartments$complete)
#Function2 for EM vs DM
DEA.subset.sc2 <- function(subset_test){
  subset_result <- FindMarkers(ooac, ident.1 = "EM", ident.2 = "DM", group.by = "condition_new", subset.ident = subset_test, test.use = "MAST", min.pct = 0.1, logfc.threshold = 0.25, latent.vars = c("cell_line", "CDR"))
  return(subset_result)
}
cluster.markers_complete_EMvsDM <- DEA.subset.sc2(compartments$complete)
#Function3 for EM-DM vs DM
DEA.subset.sc3 <- function(subset_test){
  subset_result <- FindMarkers(ooac, ident.1 = "EM-DM", ident.2 = "DM", group.by = "condition_new", subset.ident = subset_test, test.use = "MAST", min.pct = 0.1, logfc.threshold = 0.25, latent.vars = c("cell_line", "CDR"))
  return(subset_result)
}
cluster.markers_complete_EM_DMvsDM <- DEA.subset.sc3(compartments$complete)
cluster.markers_complete <- list(cluster.markers_complete_EMvsEM_DM, cluster.markers_complete_EMvsDM, cluster.markers_complete_EM_DMvsDM)
names(cluster.markers_complete) <- c("EM vs EM-DM", "EM vs DM", "EM-DM vs DM")
DE_genes <- cluster.markers_complete
#Functions, same as for Figure 4D
EM_EMvsEM_DM <- as.data.frame(DE_genes[['EM vs EM-DM']])
EM_EMvsDM_DM <- as.data.frame(DE_genes[['EM vs DM']])
EM_DMvsDM_DM <- as.data.frame(DE_genes[['EM-DM vs DM']])
#Generate a column named 'gene' from rownames
EM_EMvsEM_DM$gene <- rownames(EM_EMvsEM_DM)
EM_EMvsDM_DM$gene <- rownames(EM_EMvsDM_DM)
EM_DMvsDM_DM$gene <- rownames(EM_DMvsDM_DM)
#Add column 'comparison' to indicate the DE test that is performed to the data
EM_EMvsEM_DM <- EM_EMvsEM_DM %>% add_column(comparison = "EM vs EM-DM")
EM_EMvsDM_DM <- EM_EMvsDM_DM %>% add_column(comparison = "EM vs DM")
EM_DMvsDM_DM <- EM_DMvsDM_DM %>% add_column(comparison = "EM-DM vs DM")
#merge the 3 DE lists into one DE_genes object
DE_genes <- rbind(EM_EMvsEM_DM, EM_EMvsDM_DM, EM_DMvsDM_DM)
#Duplicate the ‘gene’ column and name ‘SYMBOL’
DE_genes$SYMBOL <- DE_genes$gene
#Add the up or down annotation for each log2FC
DE_genes <- add.down_or_up.column(DE_genes)
#Add column with the absolute value of log2FC to include up- and down-regulated genes:
DE_genes <- add.abs.log2FC(DE_genes)
#Filter the significant genes. Use p < 0.05 or p < 0.01. 
filtered_DE_genes <- DE_genes[DE_genes$p_val_adj < 0.01,]
#Filter based on absoluted log2FC to include both up and downregulated genes, as for the different condition comparisons EM vs EM-DM this means including both genes upregulated in EM as well as EM-DM
filtered_log2FC_DE_genes <- filtered_DE_genes[filtered_DE_genes$abs_log2FC > 0.5,]

#Generate average gene expression matrix for intestine-on-chip data
#Change idents to cell type clusters: for Figure 4D use "annotation_res0.34_simple2", for Figure 4I use "condition_new"
ooac <- SetIdent(ooac, value = "annotation_res0.34_simple2")
ooac_av_ex <- AverageExpression(ooac, assays = 'RNA')
ooac_av_ex <- as.data.frame(ooac_av_ex)
#subset genes in gene list of interest
#Use the 'filtered_log2FC_DE_genes' object for the relevant Figure (4D or 4I)
ooac_av_ex <- as.matrix(ooac_av_ex[rownames(ooac_av_ex) %in% filtered_log2FC_DE_genes$gene, ])
ooac_av_ex <- t(scale(t(ooac_av_ex)))
#Format sample names
colnames(ooac_av_ex) <- gsub("\\.","-",colnames(ooac_av_ex))
colnames(ooac_av_ex) <- gsub("RNA-","",colnames(ooac_av_ex))

#Generate average gene expression matrix for Elmentaite dataset
#Change idents to cell type clusters: for Figure 4D use "cell_type_2", for Figure 4I use "Age_group_2"
#For Figure 4I use object 'elmentaite_2021_SI_age'
elmentaite_2021_SI <- SetIdent(elmentaite_2021_SI, value = "cell_type_2")
ooac_av_ex_el <- AverageExpression(elmentaite_2021_SI, assays = 'RNA')
ooac_av_ex_el <- as.data.frame(ooac_av_ex_el)
#subset genes in gene list of interest
#Use the 'filtered_log2FC_DE_genes' object for the relevant Figure (4D or 4I)
ooac_av_ex_el <- as.matrix(ooac_av_ex_el[rownames(ooac_av_ex_el) %in% filtered_log2FC_DE_genes$gene, ])
ooac_av_ex_el <- t(scale(t(ooac_av_ex_el)))
#Format sample names
colnames(ooac_av_ex_el) <- gsub("\\.","/",colnames(ooac_av_ex_el))
colnames(ooac_av_ex_el) <- gsub("RNA/"," ",colnames(ooac_av_ex_el))

#combine the average expression matrix of intestine-on-chip data and elmentaite data based on rownames (i.e. selected gene list)
#very important to set the 'sort=FALSE' to preserve the order of the genes, as this will influence the clustering in Heatmap
ooac_av_ex_merged <- merge(ooac_av_ex, ooac_av_ex_el, by=0, all=TRUE, sort=FALSE)
rownames(ooac_av_ex_merged) <- ooac_av_ex_merged$Row.names
ooac_av_ex_merged <- ooac_av_ex_merged[,-1]

#rename samples in expression matrix
colnames(ooac_av_ex_merged) <- gsub("\\/","-",colnames(ooac_av_ex_merged))
colnames(ooac_av_ex_merged) <- gsub("\\-"," ",colnames(ooac_av_ex_merged))

# Reordering cell types to order used for visualisation
#Figure 4I:
ooac_av_ex_merged <- ooac_av_ex_merged[, c('EM',
                                           'EM DM',
                                           'DM',
                                           ' Fetal',
                                           ' Adult')]

#Figure 4D:
#The reference data cell types have a 'space' before the name
ooac_av_ex_merged <- ooac_av_ex_merged[, c('TA stem cell',
                                           'Enterocyte',
                                           'Paneth like cell',
                                           'Goblet cell',
                                           'Enteroendocrine cell',
                                           'Mesenchymal like epithelial cell',
                                           'Dividing mesenchymal neural cell',
                                           'Myofibroblast',
                                           'Neuron',
                                           'WNT4 positive neural cell',
                                           ' Stem cell',
                                           ' TA',
                                           ' Enterocyte',
                                           ' Colonocyte',
                                           ' BEST4  epithelial cell',
                                           ' Paneth cell',
                                           ' Goblet cell',
                                           ' Enteroendocrine cell',
                                           ' Microfold cell',
                                           ' Tuft cell',
                                           ' Fetal proximal progenitor',
                                           ' Fetal distal progenitor',
                                           ' Fetal CLDN10  cell',
                                           ' Cycling stromal cell',
                                           ' Stromal subtype',
                                           ' Myofibroblast',
                                           ' Smooth muscle cell',
                                           ' Pericyte',
                                           ' Interstitial cell of Cajal',
                                           ' Mesoderm',
                                           ' Mesothelium',
                                           ' Lymph node fibroblast',
                                           ' Cycling neuroblast',
                                           ' Neuroblast',
                                           ' Enteric neuron',
                                           ' Cycling ENCC glia',
                                           ' ENCC glia progenitor',
                                           ' Glial cell')]

#Calculate correlation values and plot
library(Hmisc)
library(corrplot)
res2 <- rcorr(as.matrix(ooac_av_ex_merged))
#set p-values from diagonal from 'NA' to '0' to correspond to corrplot() output
diag(res2$P) <- 0
#Figure 4D: subset the correlations and p-values so that the intestine-on-chip cell types are on y-axis and Gut Cell Atlas cell types on the x-axis
corr_values <- res2$r
corr_values <- corr_values[1:10, 11:38]
p_values <- res2$P
p_values <- p_values[1:10, 11:38]
#Figure 4I: subset the correlations and p-values so that the intestine-on-chip data is on y-axis and Gut Cell Atlas data is the x-axis
corr_values <- res2$r
corr_values <- corr_values[1:3, 1:5]
p_values <- res2$P
p_values <- p_values[1:3, 1:5]
#Correlation plot (insignificant correlations are left blank)
corrplot(corr_values, type="full", order="original", tl.cex=1.2,
         col = rev(COL2("RdBu", n = 200)),
         p.mat = p_values, sig.level = 0.01, insig = "blank", 
         tl.col = "black") 


#Generate the composition barplots for Figure 4E-H
#use the 'elmentaite_2021' (subsetted for epithelial, mesenchymal, neural category) and 'ooac' objects as generated above (with the simplified cell type annotations).
#set these factor vectors to character vectors
elmentaite_2021$cell_type <- as.character(elmentaite_2021$cell_type)
elmentaite_2021$category <- as.character(elmentaite_2021$category)
elmentaite_2021$Age_group <- as.character(elmentaite_2021$Age_group)
elmentaite_2021$Region <- as.character(elmentaite_2021$Region)
elmentaite_2021$Sample.name <- as.character(elmentaite_2021$Sample.name)
elmentaite_2021$Region.code <- as.character(elmentaite_2021$Region.code)

#for Figure 4H use elmentaite_2021$cell_type_3 and ooac$annotation_res0.34_simple3 to synchronize intestine-on-chip and elmentaite annotation:
elmentaite_2021$cell_type_3 = elmentaite_2021$cell_type_2
elmentaite_2021$cell_type_3[elmentaite_2021$cell_type_3 == 'Paneth cell'] <- 'Paneth(-like) cell'
elmentaite_2021$cell_type_3[elmentaite_2021$cell_type_3 == 'TA'] <- 'TA/stem cell'
elmentaite_2021$cell_type_3[elmentaite_2021$cell_type_3 == 'Stem cell'] <- 'TA/stem cell'
ooac$annotation_res0.34_simple3 = ooac$annotation_res0.34_simple2
ooac$annotation_res0.34_simple3[ooac$annotation_res0.34_simple3 == 'Paneth-like cell'] <- 'Paneth(-like) cell'

#combine the samples that were seperated based on EPCAM-selection
elmentaite_2021$Sample.name_2 = elmentaite_2021$Sample.name
elmentaite_2021$Sample.name_2[elmentaite_2021$Sample.name_2 == 'T036NEG'] <- 'T036'
elmentaite_2021$Sample.name_2[elmentaite_2021$Sample.name_2 == 'T036POS'] <- 'T036'
elmentaite_2021$Sample.name_2[elmentaite_2021$Sample.name_2 == 'T110NEG'] <- 'T110'
elmentaite_2021$Sample.name_2[elmentaite_2021$Sample.name_2 == 'T110POS'] <- 'T110'
#combine the samples of the two ileal locations, just as was done with duodenal locations
elmentaite_2021$Region.code_2 = elmentaite_2021$Region.code
elmentaite_2021$Region.code_2[elmentaite_2021$Region.code_2 == 'ILE1'] <- 'ILE'
elmentaite_2021$Region.code_2[elmentaite_2021$Region.code_2 == 'ILE2'] <- 'ILE'
#remove certain individuals that have skewed cell type ratios
elmentaite_2021 <- elmentaite_2021[, !(elmentaite_2021@meta.data[['Sample.name']] %in% c('A32 (411C)'))]

#subsetting
elmentaite_2021_epi <- elmentaite_2021[, elmentaite_2021@meta.data[['category']] %in% c('Epithelial')]
elmentaite_2021_mes_neur <- elmentaite_2021[, elmentaite_2021@meta.data[['category']] %in% c('Mesenchymal', 'Neuronal')]
ooac_epithelial <- ooac[, ooac@meta.data[['compartment_new']] %in% c('Epithelial cell')]
ooac_mesenchymal <- ooac[, ooac@meta.data[['compartment_new']] %in% c('Mesenchymal cell')]
ooac_EEC <- ooac[, ooac@meta.data[['predicted.celltype.elmentaite.2021.full.2']] %in% c('Enteroendocrine cell')]

#generate color palette
n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
my_cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#Figure 4F: compartment/category barplot celltype
#'Age_group_2' annotation does not include the pediatric and fetal first trimester (except BRC2026, BRC2029, BRC2258, BRC2259) samples that have been enriched for Epcam-positive cells
#synchronize category naming
ooac$compartment_new3 = ooac$compartment_new2
ooac$compartment_new3[ooac$compartment_new3 == 'Epithelial cell'] <- 'Epithelial'
ooac$compartment_new3[ooac$compartment_new3 == 'Mesenchymal/Neural cell'] <- 'Mesenchymal/Neural'
elmentaite_2021$category_2 = elmentaite_2021$category
elmentaite_2021$category_2[elmentaite_2021$category_2 == 'Neuronal'] <- 'Neural'
order_category <- c('Neural',
                    'Mesenchymal',
                    'Mesenchymal/Neural',
                    'Epithelial')
#elmentaite_2021 barplot
elmentaite_2021@meta.data <- elmentaite_2021@meta.data %>%
  mutate(Region_age_group = paste(Region, Age_group_2, sep = '_'))
elmentaite_2021@meta.data <- elmentaite_2021@meta.data %>%
  mutate(Region_sample_name = paste(Region.code_2, Sample.name_2, sep = '_'))
props <- propeller(clusters=elmentaite_2021$category_2, sample=elmentaite_2021$Region_sample_name, group=elmentaite_2021$Region_age_group, transform = 'logit')
props_subset <- subset(props, select=c('PropMean.SmallInt_Adult', 'PropMean.SmallInt_Fetal'))
props_subset_M <- data.matrix(props_subset)
barplot(props_subset_M, legend=TRUE, ylab="Proportions", col = alpha(my_cols[55:70]), args.legend = list(x = "topright", inset = c(-0.02, 0), cex = 0.8))
#ooac barplot
ooac@meta.data <- ooac@meta.data %>% 
  mutate(condition_cell_line = paste(condition, cell_line, sep = '_'))
props_2 <- propeller(clusters=ooac$compartment_new3, sample=ooac$condition_cell_line, group=ooac$condition, transform = 'logit')
props_subset_2 <- subset(props_2, select=c('PropMean.EM.EM', 'PropMean.EM.DM', 'PropMean.DM.DM'))
props_subset_2_M <- data.matrix(props_subset_2)
barplot(props_subset_2_M, legend=TRUE, ylab="Proportions", col = alpha(my_cols[55:70]), args.legend = list(x = "topright", inset = c(-0.02, 0), cex = 0.8))
#merge together
props_3 <- merge(props, props_2, by = 'row.names', all = TRUE)
props_3[is.na(props_3)] <- 0
rownames(props_3) <- props_3$Row.names
#plot the combined proportions in barplot
props_subset <- subset(props_3, select=c('PropMean.EM.EM', 'PropMean.EM.DM', 'PropMean.DM.DM', 'PropMean.SmallInt_Fetal', 'PropMean.SmallInt_Adult'))
props_subset <- props_subset[match(order_category, rownames(props_subset)), ]
props_subset_M <- data.matrix(props_subset)
colnames(props_subset_M) <- c('EM','EM-DM','DM', 'Fetal', 'Adult')
barplot(props_subset_M, legend=TRUE, cex.names=1.0, ylab="Proportions", cex.lab=1.2, col = alpha(my_cols[59:56]), args.legend = list(x = "bottomleft", inset = c(0.0, -0.7), cex = 1.0), main="Epithelial to mesenchymal/neural ratio")

#Figure 4E: predicted EEC subtypes
ooac_EEC@meta.data <- ooac_EEC@meta.data %>% 
  mutate(condition_cell_line = paste(condition, cell_line, sep = '_'))
props_2 <- propeller(clusters=ooac_EEC$predicted.celltype.elmentaite.2021.full, sample=ooac_EEC$condition_cell_line, group=ooac_EEC$condition, transform = 'logit')
props_subset_2 <- subset(props_2, select=c('PropMean.EM.EM', 'PropMean.EM.DM', 'PropMean.DM.DM'))
props_subset_2_M <- data.matrix(props_subset_2)
colnames(props_subset_2_M) <- c('EM','EM-DM','DM')
barplot(props_subset_2_M, legend=TRUE, cex.lab=1.2, cex.names=1.2, ylab="Proportions", col = alpha(my_cols[64:75]), args.legend = list(x = "bottomleft", inset = c(0.0, -1.05), cex = 1.3), main="Predicted enteroendocrine subtypes")


#Figure 4H: epithelial composition
order_epithelial <- c('BEST4+ epithelial cell',
                      'Tuft cell',
                      'Microfold cell',
                      'Mesenchymal-like epithelial cell',
                      'Enteroendocrine cell',
                      'Goblet cell',
                      'Paneth(-like) cell',
                      'Enterocyte',
                      'TA/stem cell')
#elmentaite_2021 barplot
elmentaite_2021_epi@meta.data <- elmentaite_2021_epi@meta.data %>%
  mutate(Region_age_group = paste(Region, Age_group, sep = '_'))
elmentaite_2021_epi@meta.data <- elmentaite_2021_epi@meta.data %>%
  mutate(Region_sample_name = paste(Region.code_2, Sample.name_2, sep = '_'))
props <- propeller(clusters=elmentaite_2021_epi$cell_type_3, sample=elmentaite_2021_epi$Region_sample_name, group=elmentaite_2021_epi$Region_age_group, transform = 'logit')
#Since I only use adult small intestine as comparison, remove all other annotations:
props <- props[!(props$PropMean.SmallInt_Adult<=0.0001),]
props_subset <- subset(props, select=c('PropMean.SmallInt_Adult'))
props_subset_M <- data.matrix(props_subset)
barplot(props_subset_M, legend=TRUE, ylab="Proportions", col = alpha(my_cols[55:70]), args.legend = list(x = "topright", inset = c(-0.02, 0), cex = 0.8))
#ooac barplot
ooac_epithelial@meta.data <- ooac_epithelial@meta.data %>% 
  mutate(condition_cell_line = paste(condition, cell_line, sep = '_'))
props_2 <- propeller(clusters=ooac_epithelial$annotation_res0.34_simple3, sample=ooac_epithelial$condition_cell_line, group=ooac_epithelial$condition, transform = 'logit')
props_subset_2 <- subset(props_2, select=c('PropMean.EM.EM', 'PropMean.EM.DM', 'PropMean.DM.DM'))
props_subset_2_M <- data.matrix(props_subset_2)
barplot(props_subset_2_M, legend=TRUE, ylab="Proportions", col = alpha(my_cols[55:70]), args.legend = list(x = "topright", inset = c(-0.02, 0), cex = 0.8))
#merge together
props_3 <- merge(props, props_2, by = 'row.names', all = TRUE)
props_3[is.na(props_3)] <- 0
rownames(props_3) <- props_3$Row.names
#plot the combined proportions in barplot
props_subset <- subset(props_3, select=c('PropMean.EM.EM', 'PropMean.EM.DM', 'PropMean.DM.DM', 'PropMean.SmallInt_Adult'))
props_subset <- props_subset[match(order_epithelial, rownames(props_subset)), ]
props_subset_M <- data.matrix(props_subset)
colnames(props_subset_M) <- c('EM','EM-DM','DM', 'Adult small intestine')
barplot(props_subset_M, legend=TRUE, cex.lab=1.2, cex.names=1.0, ylab="Proportions", col = alpha(c("dark green","blue","#D2691E","#725191","orange","#E31A1C","#FB9A99","#A8D58D","#A6CEE3")), args.legend = list(x = "bottomleft", inset = c(0.0, -1.0), cex = 0.9), main="Epithelial composition")

#Figure 4G: Cell cycle phase in epithelial cells
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ooac <- CellCycleScoring(ooac, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
elmentaite_2021_epi <- CellCycleScoring(elmentaite_2021_epi, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#Add column with cell cycle phase change to 'dividing', 'non-dividing' annotation
ooac$dividing = ooac$Phase
ooac$dividing[ooac$dividing == 'G2M'] <- 'Dividing'
ooac$dividing[ooac$dividing == 'S'] <- 'Dividing'
ooac$dividing[ooac$dividing == 'G1'] <- 'Non-dividing'
elmentaite_2021_epi$dividing = elmentaite_2021_epi$Phase
elmentaite_2021_epi$dividing[elmentaite_2021_epi$dividing == 'G2M'] <- 'Dividing'
elmentaite_2021_epi$dividing[elmentaite_2021_epi$dividing == 'S'] <- 'Dividing'
elmentaite_2021_epi$dividing[elmentaite_2021_epi$dividing == 'G1'] <- 'Non-dividing'
ooac_epithelial <- ooac[, ooac@meta.data[['compartment_new']] %in% c('Epithelial cell')]
order_epithelial <- c('Non-dividing',
                      'Dividing')
#elmentaite_2021 barplot
elmentaite_2021_epi@meta.data <- elmentaite_2021_epi@meta.data %>%
  mutate(Region_age_group = paste(Region, Age_group, sep = '_'))
elmentaite_2021_epi@meta.data <- elmentaite_2021_epi@meta.data %>%
  mutate(Region_sample_name = paste(Region.code_2, Sample.name_2, sep = '_'))
props <- propeller(clusters=elmentaite_2021_epi$dividing, sample=elmentaite_2021_epi$Region_sample_name, group=elmentaite_2021_epi$Region_age_group, transform = 'logit')
props_subset <- subset(props, select=c('PropMean.SmallInt_Adult', 'PropMean.SmallInt_Pediatric', 'PropMean.SmallInt_Second.trim', 'PropMean.SmallInt_First.trim'))
props_subset_M <- data.matrix(props_subset)
barplot(props_subset_M, legend=TRUE, ylab="Proportions", col = alpha(my_cols[55:70]), args.legend = list(x = "topright", inset = c(-0.02, 0), cex = 0.8))
#ooac barplot
ooac_epithelial@meta.data <- ooac_epithelial@meta.data %>% 
  mutate(condition_cell_line = paste(condition, cell_line, sep = '_'))
props_2 <- propeller(clusters=ooac_epithelial$dividing, sample=ooac_epithelial$condition_cell_line, group=ooac_epithelial$condition, transform = 'logit')
props_subset_2 <- subset(props_2, select=c('PropMean.EM.EM', 'PropMean.EM.DM', 'PropMean.DM.DM'))
props_subset_2_M <- data.matrix(props_subset_2)
barplot(props_subset_2_M, legend=TRUE, ylab="Proportions", col = alpha(my_cols[55:70]), args.legend = list(x = "topright", inset = c(-0.02, 0), cex = 0.8))
#merge together
props_3 <- merge(props, props_2, by = 'row.names', all = TRUE)
props_3[is.na(props_3)] <- 0
rownames(props_3) <- props_3$Row.names
props_subset <- subset(props_3, select=c('PropMean.EM.EM', 'PropMean.EM.DM', 'PropMean.DM.DM', 'PropMean.SmallInt_First.trim', 'PropMean.SmallInt_Second.trim', 'PropMean.SmallInt_Pediatric', 'PropMean.SmallInt_Adult'))
props_subset <- props_subset[match(order_epithelial, rownames(props_subset)), ]
props_subset_M <- data.matrix(props_subset)
colnames(props_subset_M) <- c('EM','EM-DM','DM', 'First trimester', 'Second trimester', 'Pediatric', 'Adult')
barplot(props_subset_M, legend=TRUE, las=2, cex.names=1.0, cex.lab=1.2,  ylab="Proportions", col = alpha(my_cols[55:80]), args.legend = list(x = "topright", inset = c(-0.2, 0), cex = 1.0), main="Cell division in epithelial cells")


