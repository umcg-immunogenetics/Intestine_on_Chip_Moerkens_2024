###########################################################################################################################
### Figure 3 ##############################################################################################################
###########################################################################################################################

#libraries
library(Seurat)
library(SeuratData)
library(dplyr)
library(ggplot2)
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

# locations of the outputs
seurat_objects_loc <- #path to where the seurat objects are saved
# Seurat object of the media conditions EM, EM-DM, DM
seurat_object_media <- paste(seurat_objects_loc, 'intestine_on_chip_media.rds', sep = '')
# read object
ooac <- readRDS(seurat_object_media)
#Set the default assay to SCT
DefaultAssay(ooac) <- 'SCT'

#Set order of the media conditions 
order <- c('EM',
           'EM-DM',
           'DM')

ooac$condition_new <- factor(ooac$condition_new, levels = order)

#Define a list of colors to be used in visualisations
n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
my_cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#Figure 3A: Visualize the epithelial and mesenchymal/neural compartment
DimPlot(ooac, reduction = 'umap', group.by = 'compartment_new2', cols=my_cols[56:57]) +
  labs(title="Media conditions combined")+
  theme(plot.title = element_text(hjust = 0.5, size=12), axis.title=element_text(size=10), axis.text=element_text(size=10))

DimPlot(ooac, reduction = 'umap', group.by = 'compartment_new2', split.by='condition_new',  cols=my_cols[56:57]) +
  labs(title="Celltype composition")+
  theme(plot.title = element_text(hjust = 0.5, size=15), axis.title=element_text(size=12), axis.text=element_text(size=10), legend.position = "bottom")

#Figure S4 A, B: Visualize gene expression levels of canonical genes
FeaturePlot(ooac, ncol=3, features = c('EPCAM', 'VIM'), split.by='condition_new')
FeaturePlot(ooac, ncol=3, features = c('LGR5','SMOC2','MKI67','PCNA', 'FABP2','RBP2', 'MT1E', 'LYZ','PRSS2','MUC2', 'CHGA', 'COL1A2','DCN', 'TAGLN', 'ELAVL3', 'ELAVL4', 'NTN1', 'WNT4'))


#Figure 3A: Plot cell type annotation
order <- c('TA/stem cell',
           'Enterocyte progenitor',
           'Enterocyte type 1',
           'Enterocyte type 2',
           'Paneth-like cell',
           'Goblet cell',
           'Enteroendocrine cell',
           'Mesenchymal-like epithelial precursor',
           'Mesenchymal-like epithelial cell type 1',
           'Mesenchymal-like epithelial cell type 2',
           'Dividing mesenchymal/neural cell',
           'Myofibroblast',
           'Neuron',
           'WNT4-positive neural cell')

ooac$annotation_res0.34_new2 <- as.character(ooac$annotation_res0.34_new2)
ooac$annotation_res0.34_new2 <- factor(ooac$annotation_res0.34_new2, levels = order)

#14 colors
dark2 <- brewer.pal(6, "Dark2")
paired <- brewer.pal(10, "Paired")
colors <- c(paired, dark2)

DimPlot(ooac, reduction = 'umap', group.by = 'annotation_res0.34_new2') +
  labs(title="Media conditions combined")+
  theme(plot.title = element_text(hjust = 0.5, size=15), axis.title=element_text(size=12), axis.text=element_text(size=10))

DimPlot(ooac, reduction = 'umap', group.by = 'annotation_res0.34_new2', split.by ='condition_new', cols=colors) +
  labs(title="Epithelial and mesenchymal subtype composition")+
  theme(plot.title = element_text(hjust = 0.5, size=12), axis.title=element_text(size=10), axis.text=element_text(size=10), legend.text = element_text(size = 10))


#Figure 3B: Dotplot cell type markers
order <- c('WNT4-positive neural cell',
           'Neuron',
           'Myofibroblast',
           'Dividing mesenchymal/neural cell',
           'Mesenchymal-like epithelial cell type 2',
           'Mesenchymal-like epithelial cell type 1',
           'Mesenchymal-like epithelial precursor',
           'Enteroendocrine cell',
           'Goblet cell',
           'Paneth-like cell',
           'Enterocyte type 2',
           'Enterocyte type 1',
           'Enterocyte progenitor',
           'TA/stem cell')

ooac$annotation_res0.34_new2 <- as.character(ooac$annotation_res0.34_new2)
ooac$annotation_res0.34_new2 <- factor(ooac$annotation_res0.34_new2, levels = order)

General_markers <- c('EPCAM', 'CDX2', 'CDH1', 'VIM')
Stem_markers <- c('SMOC2', 'LGR5')
Proliferation_markers <- c('MKI67', 'TOP2A', 'PCNA', 'CENPF')
Enterocyte_markers  <- c('FGB', 'FGG', 'FABP2', 'RBP2', 'APOA4', 'CYP3A5', 'ALDOB', 'ANPEP', 'SI', 'SLC2A2', 'SLC39A4', 'MT1E', 'MT1G', 'MT1H', 'MT2A')
Paneth_markers <- c('LYZ', 'PRSS2', 'PLA2G2A')
Goblet_markers <- c('MUC2', 'TFF3', 'ZG16', 'CLCA1')
EEC_markers <- c('CHGA', 'NEUROD1', 'MLN', 'GHRL', 'GCG')
Fibroblast_markers <- c('COL1A1', 'COL1A2', 'DCN', 'TAGLN', 'ACTA2', 'ACTG2', 'MYL9')
Neuronal_markers <- c('ELAVL3', 'ELAVL4', 'GAP43', 'CRABP1', 'ONECUT2', 'STMN2', 'TUBB2B', 'NCAM1', 'NTN1', 'SLIT2', 'NTRK2', 'WNT4')

features <- list("General" = General_markers,"Stem cell"= Stem_markers, "Proliferation" = Proliferation_markers, "Enterocyte" = Enterocyte_markers, "Paneth" = Paneth_markers,"Goblet" =Goblet_markers, "EEC" = EEC_markers, "Myofibroblast" = Fibroblast_markers, "Neural" = Neuronal_markers)

DotPlot(object = ooac, features=features, group.by = "annotation_res0.34_new2") + theme(axis.text.x = element_text(angle = 40, hjust = 1), axis.text=element_text(size=12))



#Figure 3C: Barplot cell type composition
n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
my_cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

ooac@meta.data <- ooac@meta.data %>%
  mutate(condition_cell_line = paste(condition, cell_line, sep = '_'))

ooac_epithelial <- ooac[, ooac@meta.data[['compartment']] %in% c('epithelial')]
ooac_mesenchymal <- ooac[, ooac@meta.data[['compartment']] %in% c('mesenchymal')]

order <- c('Mesenchymal/Neural cell',
           'Epithelial cell')
ooac$compartment_new2 <- as.character(ooac$compartment_new2)
ooac$compartment_new2 <- factor(ooac$compartment_new2, levels = order)

order_epithelial <- c('Mesenchymal-like epithelial cell',
                      'Enteroendocrine cell',
                      'Goblet cell',
                      'Paneth-like cell',
                      'Enterocyte',
                      'TA/stem cell')

order_mesenchymal <- c('WNT4-positive neural cell',
                       'Neuron',
                       'Myofibroblast',
                       'Dividing mesenchymal/neural cell')


par(mar=c(8.1, 4.1, 4.1, 2.1), xpd=TRUE)

#Cellular compartment
props <- propeller(clusters=ooac$compartment_new2, sample=ooac$condition_cell_line, group=ooac$condition, transform = 'logit')
props_subset <- subset(props, select=c('PropMean.EM.EM', 'PropMean.EM.DM', 'PropMean.DM.DM'))
props_subset_M <- data.matrix(props_subset)
colnames(props_subset_M) <- c('EM','EM-DM','DM')
barplot(props_subset_M, legend=TRUE, ylab="Proportions", col = alpha(my_cols[57:56]), cex.lab=1.2, args.legend = list(x = "bottom", inset = c(0.0, -0.4), cex = 1.0), main="Cell type composition")

par(mar=c(10.1, 4.1, 4.1, 2.1), xpd=TRUE)

#Epithelial compartment
props <- propeller(clusters=ooac_mesenchymal$annotation_res0.34_simple2, sample=ooac_mesenchymal$condition_cell_line, group=ooac_mesenchymal$condition, transform = 'logit')
props_subset <- subset(props, select=c('PropMean.EM.EM', 'PropMean.EM.DM', 'PropMean.DM.DM'))
props_subset <- props_subset[match(order_mesenchymal, rownames(props_subset)), ]
props_subset_M <- data.matrix(props_subset)
colnames(props_subset_M) <- c('EM','EM-DM','DM')
barplot(props_subset_M, legend=TRUE, ylab="Proportions", col = alpha(my_cols[58:55]), cex.lab=1.2, args.legend = list(x = "bottom", inset = c(0.0, -0.5), cex = 0.9), main="Mesenchymal and neural composition")

#Mesenchymal/neural compartment
props <- propeller(clusters=ooac_epithelial$annotation_res0.34_simple2, sample=ooac_epithelial$condition_cell_line, group=ooac_epithelial$condition, transform = 'logit')
props_subset <- subset(props, select=c('PropMean.EM.EM', 'PropMean.EM.DM', 'PropMean.DM.DM'))
props_subset <- props_subset[match(order_epithelial, rownames(props_subset)), ]
props_subset_M <- data.matrix(props_subset)
colnames(props_subset_M) <- c('EM','EM-DM','DM')
barplot(props_subset_M, legend=TRUE, ylab="Proportions", col = alpha(c("#725191","orange","#E31A1C","#FB9A99","#A8D58D","#A6CEE3")), cex.lab=1.2, args.legend = list(x = "bottom", inset = c(0.0, -0.7), cex = 0.9), main="Epithelial composition")

#set the graph margins back to normal
par(mar=c(5.1, 4.1, 4.1, 2.1))


#Figure S5B: Celltype variablity among biological replicates
#Select 14 colors
dark2 <- brewer.pal(6, "Dark2")
paired <- brewer.pal(10, "Paired")
colors <- c(paired, dark2)
#composition plot
summary <- ooac@meta.data %>% 
  group_by(condition_new, cell_line, annotation_res0.34_new2) %>% 
  summarise(count = n()) %>%
  mutate(total.cells = sum(count)) %>%
  mutate(pct = count/total.cells)
summary$condition_new <- as.factor(summary$condition_new)
# Reordering cell types
order <- c('TA/stem cell',
           'Enterocyte progenitor',
           'Enterocyte type 1',
           'Enterocyte type 2',
           'Paneth-like cell',
           'Goblet cell',
           'Enteroendocrine cell',
           'Mesenchymal-like epithelial precursor',
           'Mesenchymal-like epithelial cell type 1',
           'Mesenchymal-like epithelial cell type 2',
           'Dividing mesenchymal/neural cell',
           'Myofibroblast',
           'Neuron',
           'WNT4-positive neural cell')
summary$annotation_res0.34_new2 <- factor(summary$annotation_res0.34_new2, levels = order)
#Plot compartment
order <- c('EM',
           'EM-DM',
           'DM')
summary$condition_new <- factor(summary$condition_new, levels = order)
#Rename 
summary$cell_line[summary$cell_line == "Geni002"] <- "Donor 1"
summary$cell_line[summary$cell_line == "Geni007"] <- "Donor 2"
summary$cell_line[summary$cell_line == "Geni012"] <- "Donor 3"
colnames(summary)[colnames(summary) == "annotation_res0.34_new2"] <- "Celltypes"
#Plot
ggplot(summary, aes(fill=Celltypes, y=pct, x=cell_line)) + 
  geom_bar(position="fill", stat="identity") +
  facet_wrap(vars(condition_new)) +
  theme_classic() +
  scale_fill_manual(values=colors) +
  labs(title = 'Celltype composition variability among donors',
       y = 'Proportion',
       x = '')
#repeat for 'compartment' by replacing annotation_res0.34_new with compartment_new2 when generating the 'summary' object
#or for epithelial and mesenchymal object separately
#colors for epithelial subtypes: 'paired', mesenchymal subtypes: 'dark2', compartment:c("#FC8D62", "#8DA0CB")

#Figure S5C: Cell cycle plot
n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
my_cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

ooac@meta.data <- ooac@meta.data %>%
  mutate(condition_cell_line = paste(condition, cell_line, sep = '_'))

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ooac <- CellCycleScoring(ooac, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(ooac, reduction = 'umap', group.by = 'Phase')

ooac_epithelial <- ooac[, ooac@meta.data[['compartment']] %in% c('epithelial')]
ooac_mesenchymal <- ooac[, ooac@meta.data[['compartment']] %in% c('mesenchymal')]

order_cycle <- c('G1',
                 'G2M',
                 'S')

par(mar=c(5.1, 4.1, 4.1, 4.1), xpd=TRUE)

#Epithelial compartment 
props <- propeller(clusters=ooac_epithelial$Phase, sample=ooac_epithelial$condition_cell_line, group=ooac_epithelial$condition, transform = 'logit')
props_subset <- subset(props, select=c('PropMean.EM.EM', 'PropMean.EM.DM', 'PropMean.DM.DM'))
props_subset <- props_subset[match(order_cycle, rownames(props_subset)), ]
props_subset_M <- data.matrix(props_subset)
colnames(props_subset_M) <- c('EM','EM-DM','DM')
barplot(props_subset_M, legend=TRUE, ylab="Proportions", col = alpha(my_cols[55:70]), cex.lab=1.2, args.legend = list(x = "topright", inset = c(-0.25, 0.0), cex = 0.8), main="Epithelial cells: Cell cycle phase")

#Mesenchymal/neural compartment 
props <- propeller(clusters=ooac_mesenchymal$Phase, sample=ooac_mesenchymal$condition_cell_line, group=ooac_mesenchymal$condition, transform = 'logit')
props_subset <- subset(props, select=c('PropMean.EM.EM', 'PropMean.EM.DM', 'PropMean.DM.DM'))
props_subset <- props_subset[match(order_cycle, rownames(props_subset)), ]
props_subset_M <- data.matrix(props_subset)
colnames(props_subset_M) <- c('EM','EM-DM','DM')
barplot(props_subset_M, legend=TRUE, ylab="Proportions", col = alpha(my_cols[55:70]), cex.lab=1.2, args.legend = list(x = "topright", inset = c(-0.25, 0.0), cex = 0.8), main="Mesenchymal cells: Cell cycle phase")

#TA/stem cell cluster
ooac_transit_amplifying_stem_cell <- ooac[, ooac@meta.data[['annotation_res0.34']] %in% c('transit_amplifying_stem_cell')]
props <- propeller(clusters=ooac_transit_amplifying_stem_cell$Phase, sample=ooac_transit_amplifying_stem_cell$condition_cell_line, group=ooac_transit_amplifying_stem_cell$condition, transform = 'logit')
props_subset <- subset(props, select=c('PropMean.EM.EM', 'PropMean.EM.DM', 'PropMean.DM.DM'))
props_subset <- props_subset[match(order_cycle, rownames(props_subset)), ]
props_subset_M <- data.matrix(props_subset)
colnames(props_subset_M) <- c('EM','EM-DM','DM')
barplot(props_subset_M, legend=TRUE, ylab="Proportions", col = alpha(my_cols[55:70]), cex.lab=1.2, args.legend = list(x = "topright", inset = c(-0.25, 0.0), cex = 0.8), main="TA/stem cells: Cell cycle phase")

#Cell cycle phase Dimplot in corresponding colors:
order_phase <- c('S',
                 'G2M',
                 'G1')
ooac$Phase <- as.character(ooac$Phase)
ooac$Phase <- factor(ooac$Phase, levels = order_phase)
DimPlot(ooac, reduction = 'umap', group.by = 'Phase', cols = alpha(c("#8DA0CB", "#FC8D62", "#66C2A5"))) + labs(title="Media conditions combined")
DimPlot(ooac, reduction = 'umap', group.by = 'Phase', split.by='condition_new', cols = alpha(c("#8DA0CB", "#FC8D62", "#66C2A5"))) + labs(title="Cell cycle phase")

#Figure S6C: EEC subtypes- gene expression
#Subset the EEC cells for EM-DM and DM condition
ooac_EMDM <- ooac[, ooac@meta.data[['condition_new']] %in% c('EM-DM')]
ooac_DM <- ooac[, ooac@meta.data[['condition_new']] %in% c('DM')]
ooac_EMDM_EEC <- ooac_EMDM[, ooac_EMDM@meta.data[['annotation_res0.34_new']] %in% c('Enteroendocrine cell')]
ooac_DM_EEC <- ooac_DM[, ooac_DM@meta.data[['annotation_res0.34_new']] %in% c('Enteroendocrine cell')]
#Plot gene expression of intestinal hormones.  
VlnPlot(ooac_EMDM_EEC, ncol=5, features = c("CHGA", 'NEUROD1', 'MLN', 'GHRL', 'SST', 'INS', 'GCG', 'NTS', 'SCT', 'TPH1', slot="counts"))
VlnPlot(ooac_DM_EEC, ncol=5, features = c("CHGA", 'NEUROD1', 'MLN', 'GHRL', 'SST', 'INS', 'GCG', 'NTS', 'SCT', 'TPH1', slot="counts"))


#Figure 3D: DE analysis between selected cell types in different media conditions
#Check that the clusters splitted by condition contain 30 or more cells, otherwise exclude.
table(ooac@meta.data[, c('condition_new', 'annotation_res0.34_new')])

ooac <- SetIdent(ooac, value = "annotation_res0.34_new")

#Generate a list of cell types to be tested
#EM vs EM-DM
cluster_selection_EMvsEM_DM <- list(enterocytes = c("Enterocyte type 1", "Enterocyte type 2", "Enterocyte progenitor"), 
                                    paneth = "Paneth cell",
                                    transit_amplifying_stem_cell = "TA/stem cell",
                                    mesenchymal_fibroblast = "Fibroblast-like cell",
                                    mesenchymal_neuronal = "Neuronal-like cell"
)

#EM vs DM
cluster_selection_EMvsDM <- list(enterocytes = c("Enterocyte type 1", "Enterocyte type 2", "Enterocyte progenitor"), 
                                 paneth = "Paneth cell",
                                 transit_amplifying_stem_cell = "TA/stem cell",
                                 mesenchymal_fibroblast = "Fibroblast-like cell",
                                 mesenchymal_neuronal = "Neuronal-like cell"
)

#EM-DM vs DM
cluster_selection_EM_DMvsDM <- list(enterocytes = c("Enterocyte type 1", "Enterocyte type 2", "Enterocyte progenitor"), 
                                    paneth = "Paneth cell",
                                    transit_amplifying_stem_cell = "TA/stem cell",
                                    mesenchymal_fibroblast = "Fibroblast-like cell",
                                    mesenchymal_neuronal = "Neuronal-like cell"
)

#Function1 for EM vs EM-DM
DEA.models.sc1 <- function(model_test){
  model_result <- FindMarkers(ooac, ident.1 = "EM", ident.2 = "EM-DM", group.by = "condition_new", subset.ident = model_test, test.use = "MAST", min.pct = 0.1, logfc.threshold = 0.25, latent.vars = c("cell_line", "CDR"))
  return(model_result)
}

cluster.markers_EMvsEM_DM <- lapply(cluster_selection_EMvsEM_DM, DEA.models.sc1)

#Function2 for EM vs DM
DEA.models.sc2 <- function(model_test){
  model_result <- FindMarkers(ooac, ident.1 = "EM", ident.2 = "DM", group.by = "condition_new", subset.ident = model_test, test.use = "MAST", min.pct = 0.1, logfc.threshold = 0.25, latent.vars = c("cell_line", "CDR"))
  return(model_result)
}

cluster.markers_EMvsDM <- lapply(cluster_selection_EMvsDM, DEA.models.sc2)

#Function3 for EM-DM vs DM
DEA.models.sc3 <- function(model_test){
  model_result <- FindMarkers(ooac, ident.1 = "EM-DM", ident.2 = "DM", group.by = "condition_new", subset.ident = model_test, test.use = "MAST", min.pct = 0.1, logfc.threshold = 0.25, latent.vars = c("cell_line", "CDR"))
  return(model_result)
}

cluster.markers_EM_DMvsDM <- lapply(cluster_selection_EM_DMvsDM, DEA.models.sc3)

DE_genes_EM_EMvsEM_DM <- cluster.markers_EMvsEM_DM
DE_genes_EM_EMvsDM_DM <- cluster.markers_EMvsDM
DE_genes_EM_DMvsDM_DM <- cluster.markers_EM_DMvsDM

#Subset for specific cluster of interest
ooac_transit_amplifying_stem_cell <- ooac[, ooac@meta.data[['annotation_res0.34']] %in% c('transit_amplifying_stem_cell')]
ooac_paneth <- ooac[, ooac@meta.data[['annotation_res0.34']] %in% c('paneth')]
ooac_myofibroblast <- ooac[, ooac@meta.data[['annotation_res0.34']] %in% c('mesenchymal_fibroblast')]
ooac_neuron <- ooac[, ooac@meta.data[['annotation_res0.34']] %in% c('mesenchymal_neuronal')]
ooac_enterocytes <- ooac[, ooac@meta.data[['annotation_res0.34']] %in% c('enterocyte_mature_1', 'enterocyte_mature_2', 'enterocyte_progenitor')]

#libraries
library(tidyverse)
library(clusterProfiler)
library(ReactomePA)
library(enrichplot)
library(org.Hs.eg.db)
library(msigdbr)

#functions
#Create column in DE gene object that indicates the direction of log2FC
add.down_or_up.column <- function(df){
  df$Direction <- "Upregulated"
  df$Direction[df$avg_log2FC < 0] <- "Downregulated"
  df$Direction <- as.factor(df$Direction)
  df$Direction <- factor(df$Direction, levels = c("Upregulated", "Downregulated"))
  return(df)
}

#Create column in DE gene object that indicates the absolute value of log2FC, to be able to set threshold for both up- and downregulated genes
add.abs.log2FC <- function(df){
  df$abs_log2FC <- abs(df$avg_log2FC)
  return(df)
}

#Function to add column with ENSEMBL and ENTREZID
SYMBOLToEntrezENSEMBLGene <- function(f) {
  df <- f
  gene <- df$SYMBOL
  gene.df <- bitr(gene, fromType = "SYMBOL",
                  toType = c("ENTREZID"),
                  OrgDb = org.Hs.eg.db)
  df <- dplyr::full_join(df, gene.df, by = "SYMBOL")
  return(df)
}

#Select the DE genes of interest: here as example performed using enterocytes clusters, can be repeated for the other cell types
EM_EMvsEM_DM <- as.data.frame(DE_genes_EM_EMvsEM_DM[['enterocytes']])
EM_EMvsDM_DM <- as.data.frame(DE_genes_EM_EMvsDM_DM[['enterocytes']])
EM_DMvsDM_DM <- as.data.frame(DE_genes_EM_DMvsDM_DM[['enterocytes']])

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

#Filter the significant genes.  
filtered_DE_genes <- DE_genes[DE_genes$p_val_adj < 0.01,]

#Filter based on absoluted log2FC to include both up and downregulated genes, as for the different condition comparisons EM_EMvsEM_DM this means including both genes upregulated in EM_EM as well as EM_DM
filtered_log2FC_DE_genes <- filtered_DE_genes[filtered_DE_genes$abs_log2FC > 0.5,]


# Generate metadata matrix
metadata <- ooac_enterocytes@meta.data
short.metadata <- metadata[!duplicated(metadata$condition_new),]
rownames(short.metadata) <- short.metadata$condition_new

# Generate average gene expression matrix for DEGs
ooac_enterocytes <- SetIdent(ooac_enterocytes, value = "condition_new")
ooac_av_ex <- AverageExpression(ooac_enterocytes, assays = 'SCT')
ooac_av_ex <- as.data.frame(ooac_av_ex)
ooac_av_ex <- as.matrix(ooac_av_ex[rownames(ooac_av_ex) %in% filtered_log2FC_DE_genes$gene, ])
ooac_av_ex <- t(scale(t(ooac_av_ex)))

#Format samples names
colnames(ooac_av_ex) <- gsub("\\.","-",colnames(ooac_av_ex))
colnames(ooac_av_ex) <- gsub("SCT-","",colnames(ooac_av_ex))
#set the sample names of expression matrix in the same order as sample names in metadata matrix
ooac_av_ex <- ooac_av_ex[ , rownames(short.metadata)]
#extract the metadata column you want to use to order the heatmap
condition <- short.metadata$condition_new

### Extract different clusters 
#calculate the distance between each gene in the matrix
gene_dist <- dist(ooac_av_ex)

#perform hierarchical clustering using hclust()
gene_hclust <- hclust(gene_dist, method = "complete")

# The default `plot()` function can be used to produce a simple dendrogram
plot(gene_hclust, labels = FALSE)

#cut the dendogram to specify the number of cluster you want
cutree <- cutree(gene_hclust, k = 3)

#convert vector into tibble
gene_cluster_pseudobulk <- cutree(gene_hclust, k = 3) %>% 
  # turn the named vector into a tibble
  enframe() %>% 
  # rename some of the columns
  rename(name = 'gene', value = 'cluster')

head(gene_cluster_pseudobulk)

#Add the EntrezID and ENSEMBL ID
gene_cluster_pseudobulk <- gene_cluster_pseudobulk %>% rename(gene = 'SYMBOL')
gene_cluster_pseudobulk <- SYMBOLToEntrezENSEMBLGene(gene_cluster_pseudobulk)
gene_cluster_pseudobulk_unique = gene_cluster_pseudobulk[!duplicated(gene_cluster_pseudobulk$SYMBOL),]

#Pathway enrichment Gene Ontology: Biological process
DEG.bp <- compareCluster(ENTREZID~cluster, data=gene_cluster_pseudobulk, fun="enrichGO", OrgDb = org.Hs.eg.db, ont = "BP", readable=T) 
DEG.bp@compareClusterResult$minuslog10p.adjust <- -log10(DEG.bp@compareClusterResult$p.adjust)

#Visualize the top 5 pathways per cluster
dotplot(DEG.bp, 
        x="cluster", 
        showCategory=5, 
        color = "minuslog10p.adjust", 
        font.size = 12, label_format = 50) + 
  scale_colour_gradient(low="#D0E0F3", high="#257EE4")+
  labs(colour = "-log10(adj. P-value)") 


#Display average gene expression as Heatmap
library(ComplexHeatmap)

## make the black color map to 0. the yellow map to highest and the purle map to the lowest
col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#FF00FF", "black", "#FFFF00"))

# Reordering cell types to order used for graph
order <- c('EM',
           'EM-DM',
           'DM')

short.metadata$condition_new <- factor(short.metadata$condition_new, levels = order)

#check number of genes per cluster
table(gene_cluster_pseudobulk_unique$cluster)

#select the genes to highlight
#enterocytes
genes_enterocytes <- c("CYP3A5","ANPEP", "SI","RBP2", "RBP4", "FABP1", "FABP6","APOA2", "APOB", "APOA4", "SLC2A2", "SLC39A14", "SLC5A12", "SLC15A1", "SLC9A3R1", "MT1E", "MT1F", "MT1G", "MT1H", "MT1X", "MT2A")
#EMT paneth/TA/enterocyte consensus - all related to 'ECM proteoglycans' and 'Extracellular matrix organization' pathways
genes_EMT <- c("VCAN", "SPARC", "LAMB1", 'COL18A1', 'COL4A2', 'LAMA1', 'VIM')
#M phase TA/stem cells (shared genes between Reactome output 'M phase' and Seurat's g2m.genes)
genes_m_phase <- c("CDK1",   "UBE2C",  "NDC80",   "CENPF",  "SMC4",   "CCNB2", "CDC20",  "KIF23",  "CENPE",  "NEK2")
#Translation TA/stem cells: Translation initiation factors and ribosomal proteins
genes_translation <- c("EIF2S3", "EIF3E", "RPL12", "EIF4B", "RPS2", "RPL26", "RPS4Y1", "EIF3D")
#Metallothionein genes TA/stem cells
genes_mt <- c('MT2A', 'MT1F', 'MT1G')
#Myofibroblast: all genes in pathway 'smooth muscle contraction' and all 'collagen genes'
genes_myofibroblast <- c("CALM2",  "CALD1",  "TPM2",  "MYL6",  "TPM1","MYL12B", "TPM4", "MYH10", "COL3A1", "COL5A2", "COL8A1", "COL5A1")
#Translation, ribosomal proteins and translation initiation factors
genes_translation_myo <- c("RPL11",  "RPS8",   "RPL5",   "RPS7",   "RPS27A", "RPL32", 'EIF3E', 'EIF4B', 'EIF3L', "EIF3F")
#neuron
genes_neuron <- c("NTNG1", "MAP2", "PAK3", "MAP6", "TANC2", "MYCBP2", "IGF1R", "MAPT", "KIDINS220", "NEDD4L", "ROBO1")
#Paneth markers 
genes_paneth <- c("MUC13", "TFF3",  "TFF1", "INAVA", 'RBP4', 'LGR4', 'JAG1')

#Paneth-like cells
genes_to_show <-c(genes_EMT, genes_paneth)
#Enterocytes
genes_to_show <-c(genes_enterocytes, genes_EMT)
#TA/stem cells 
genes_to_show <-c(genes_EMT, genes_m_phase, genes_translation, genes_mt)
#Myofibroblasts
genes_to_show <-c(genes_myofibroblast, genes_translation_myo)
#Neurons
genes_to_show <-c(genes_neuron, genes_translation_myo)

ha = rowAnnotation(foo = anno_mark(at = which(rownames(ooac_av_ex) %in% genes_to_show),
                                   labels = rownames(ooac_av_ex)[rownames(ooac_av_ex)%in%genes_to_show]))

#complex
ht <- Heatmap(ooac_av_ex, name = "Expression",  
              column_split = factor(short.metadata$condition_new),
              row_split = factor(gene_cluster_pseudobulk_unique$cluster),
              cluster_columns = FALSE,
              right_annotation = ha,
              show_row_names = FALSE,
              show_column_dend = FALSE,
              cluster_column_slices = TRUE,
              column_title_gp = gpar(fontsize = 14),
              column_gap = unit(1.0, "mm"),
              cluster_rows = TRUE,
              show_row_dend = TRUE,
              col = col_fun,
              row_names_gp = gpar(fontsize = 6),
              column_title_rot = 0,
              top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))),
              show_column_names = FALSE,
              use_raster = TRUE,
              raster_quality = 5)


#Figure S4C: DE analysis and pathway enrichment between cell types
#Use the same functions as described above for Figure 3D
#Change idents to cell type clusters:
ooac <- SetIdent(ooac, value = "annotation_res0.34_new2")
#DE analysis
cluster.markers_clusters_new2 <- FindAllMarkers(ooac, min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST", latent.vars = c("cell_line", "CDR"))

DE_genes <- cluster.markers_clusters_new2
#Duplicate the ‘gene’ column and name ‘SYMBOL’
DE_genes$SYMBOL <- DE_genes$gene
#Add the up or down annotation for each log2FC
DE_genes <- add.down_or_up.column(DE_genes)
#Add column with the absolute value of log2FC to include up- and down-regulated genes:
DE_genes <- add.abs.log2FC(DE_genes)
#Filter the significant genes.  
filtered_DE_genes <- DE_genes[DE_genes$p_val_adj < 0.01,]
filtered_log2FC_DE_genes <- filtered_DE_genes[filtered_DE_genes$avg_log2FC > 1.0,]

#Add the columns with ENTREZID, there is always a certain percentage that it can’t map
filtered_log2FC_DE_genes <- SYMBOLToEntrezENSEMBLGene(filtered_log2FC_DE_genes)

#Pathway enrichment using Gene Ontology: Biological process
DEG.bp <- compareCluster(ENTREZID~cluster+Direction, data=filtered_log2FC_DE_genes, fun="enrichGO", OrgDb = org.Hs.eg.db, ont = "BP", readable=T) 
DEG.bp@compareClusterResult$minuslog10p.adjust <- -log10(DEG.bp@compareClusterResult$p.adjust)

#Plot the top 5 pathways per cell type in a dotplot
#Order the cell types
order <- c('TA/stem cell',
           'Enterocyte progenitor',
           'Enterocyte type 1',
           'Enterocyte type 2',
           'Paneth-like cell',
           'Goblet cell',
           'Enteroendocrine cell',
           'Mesenchymal-like epithelial precursor',
           'Mesenchymal-like epithelial cell type 1',
           'Mesenchymal-like epithelial cell type 2',
           'Dividing mesenchymal/neural cell',
           'Myofibroblast',
           'Neuron',
           'WNT4-positive neural cell')

DEG.bp@compareClusterResult$cluster <- factor(DEG.bp@compareClusterResult$cluster, levels = order)

dotplot(DEG.bp, 
        x="cluster", 
        showCategory=5, 
        color = "minuslog10p.adjust", 
        font.size = 9, label_format = 50) + 
  scale_colour_gradient(low="#D0E0F3", high="#257EE4")+
  labs(colour = "-log10(adj. P-value)") +
  theme(axis.text.x = element_text(size = 14, angle = 50, hjust=1.0, color = 'black'), axis.text.y = element_text(size = 12)) 

