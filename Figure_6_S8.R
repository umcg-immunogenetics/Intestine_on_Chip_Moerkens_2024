###########################################################################################################################
### Figure 6 ##############################################################################################################
###########################################################################################################################

#libraries
library(Seurat)
library(SeuratData)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(speckle)
library(limma)
library(scater)
library(patchwork)
library(edgeR)
library(statmod)
library(viridisLite)
library(viridis)
library(RColorBrewer)
library(clusterProfiler)
library(ReactomePA)
library(enrichplot)
library(org.Hs.eg.db)
library(msigdbr)
library(ComplexHeatmap)

# locations of the objects
seurat_objects_loc <- #path to the folder where the intestine-on-chip data of the interferon stimulations is saved
seurat_object_IFN <- paste(seurat_objects_loc, 'intestine_on_chip_IFN.rds', sep = '')

# read object
ooac <- readRDS(seurat_object_IFN)
#Set the default assay to SCT
DefaultAssay(ooac) <- 'SCT'

#Subset ooac in different compartments
ooac_epithelial <- ooac[, ooac@meta.data[['compartment_new']] %in% c('Epithelial cell')]
ooac_mesenchymal <- ooac[, ooac@meta.data[['compartment_new']] %in% c('Mesenchymal cell')]

#Figure 6A and S8B: DE analysis between stimulation and control
ooac <- SetIdent(ooac, value = "compartment_new")
#Generate a list of all different elements to be tested. 
compartments <- list(complete = NULL, 
                     epithelial = "Epithelial cell",
                     mesenchymal = "Mesenchymal cell"
)
#Function1 for EM-DM + IFN-γ vs EM-DM
DEA.subset.sc1 <- function(subset_test){
  subset_result <- FindMarkers(ooac, ident.1 = "EM-DM + IFN-γ", ident.2 = "EM-DM", group.by = "condition_new", subset.ident = subset_test, test.use = "MAST", min.pct = 0.1, logfc.threshold = 0.25, latent.vars = c("cell_line", "CDR"))
  return(subset_result)
}
cluster.markers_complete_EM_DM_IFNyvsEM_DM <- DEA.subset.sc1(compartments$complete)
cluster.markers_epithelial_EM_DM_IFNyvsEM_DM <- DEA.subset.sc1(compartments$epithelial)
cluster.markers_mesenchymal_EM_DM_IFNyvsEM_DM <- DEA.subset.sc1(compartments$mesenchymal)
#Function2 for EM-DM + IFN-β vs EM-DM
DEA.subset.sc2 <- function(subset_test){
  subset_result <- FindMarkers(ooac, ident.1 = "EM-DM + IFN-β", ident.2 = "EM-DM", group.by = "condition_new", subset.ident = subset_test, test.use = "MAST", min.pct = 0.1, logfc.threshold = 0.25, latent.vars = c("cell_line", "CDR"))
  return(subset_result)
}
cluster.markers_complete_EM_DM_IFNbvsEM_DM <- DEA.subset.sc2(compartments$complete)
cluster.markers_epithelial_EM_DM_IFNbvsEM_DM <- DEA.subset.sc2(compartments$epithelial)
cluster.markers_mesenchymal_EM_DM_IFNbvsEM_DM <- DEA.subset.sc2(compartments$mesenchymal)
#Function3 for EM-DM + IFN-γ vs EM-DM + IFN-β
DEA.subset.sc3 <- function(subset_test){
  subset_result <- FindMarkers(ooac, ident.1 = "EM-DM + IFN-γ", ident.2 = "EM-DM + IFN-β", group.by = "condition_new", subset.ident = subset_test, test.use = "MAST", min.pct = 0.1, logfc.threshold = 0.25, latent.vars = c("cell_line", "CDR"))
  return(subset_result)
}
cluster.markers_complete_EM_DM_IFNyvsEM_DM_IFNb <- DEA.subset.sc3(compartments$complete)
cluster.markers_epithelial_EM_DM_IFNyvsEM_DM_IFNb <- DEA.subset.sc3(compartments$epithelial)
cluster.markers_mesenchymal_EM_DM_IFNyvsEM_DM_IFNb <- DEA.subset.sc3(compartments$mesenchymal)
#Create lists for each tested compartment combination
cluster.markers_complete <- list(cluster.markers_complete_EM_DM_IFNyvsEM_DM, cluster.markers_complete_EM_DM_IFNbvsEM_DM, cluster.markers_complete_EM_DM_IFNyvsEM_DM_IFNb)
cluster.markers_epithelial <- list(cluster.markers_epithelial_EM_DM_IFNyvsEM_DM, cluster.markers_epithelial_EM_DM_IFNbvsEM_DM, cluster.markers_epithelial_EM_DM_IFNyvsEM_DM_IFNb)
cluster.markers_mesenchymal <- list(cluster.markers_mesenchymal_EM_DM_IFNyvsEM_DM, cluster.markers_mesenchymal_EM_DM_IFNbvsEM_DM, cluster.markers_mesenchymal_EM_DM_IFNyvsEM_DM_IFNb)
#Naming the items in the lists
names(cluster.markers_complete) <- c("IFN-γ vs control", "IFN-β vs control", "IFN-γ vs IFN-β")
names(cluster.markers_epithelial) <- c("IFN-γ vs control", "IFN-β vs control", "IFN-γ vs IFN-β")
names(cluster.markers_mesenchymal) <- c("IFN-γ vs control", "IFN-β vs control", "IFN-γ vs IFN-β")

#test the DE genes per cellular compartment sequentially
DE_genes <- cluster.markers_complete
DE_genes <- cluster.markers_epithelial
DE_genes <- cluster.markers_mesenchymal

#Functions
add.down_or_up.column <- function(df){
  df$Direction <- "Upregulated"
  df$Direction[df$avg_log2FC < 0] <- "Downregulated"
  df$Direction <- as.factor(df$Direction)
  df$Direction <- factor(df$Direction, levels = c("Upregulated", "Downregulated"))
  return(df)
}
add.abs.log2FC <- function(df){
  df$abs_log2FC <- abs(df$avg_log2FC)
  return(df)
}

#Extract DE genes per comparison
IFNyvscontrol <- as.data.frame(DE_genes[['IFN-γ vs control']])
IFNbvscontrol <- as.data.frame(DE_genes[['IFN-β vs control']])
IFNyvsIFNb <- as.data.frame(DE_genes[['IFN-γ vs IFN-β']])
#Generate a column named 'gene' from rownames
IFNyvscontrol$gene <- rownames(IFNyvscontrol)
IFNbvscontrol$gene <- rownames(IFNbvscontrol)
IFNyvsIFNb$gene <- rownames(IFNyvsIFNb)
#Add column 'comparison' to indicate the DE test that is performed to the data
IFNyvscontrol <- IFNyvscontrol %>% add_column(comparison = "IFN-γ vs control")
IFNbvscontrol <- IFNbvscontrol %>% add_column(comparison = "IFN-β vs control")
IFNyvsIFNb <- IFNyvsIFNb %>% add_column(comparison = "IFN-γ vs IFN-β")
#Duplicate the ‘gene’ column and name ‘SYMBOL’
IFNyvscontrol$SYMBOL <- IFNyvscontrol$gene
IFNbvscontrol$SYMBOL <- IFNbvscontrol$gene
IFNyvsIFNb$SYMBOL <- IFNyvsIFNb$gene
#Add the up or down annotation for each log2FC
IFNyvscontrol <- add.down_or_up.column(IFNyvscontrol)
IFNbvscontrol <- add.down_or_up.column(IFNbvscontrol)
IFNyvsIFNb <- add.down_or_up.column(IFNyvsIFNb)
#Add column with the absolute value of log2FC to include up- and down-regulated genes:
IFNyvscontrol <- add.abs.log2FC(IFNyvscontrol)
IFNbvscontrol <- add.abs.log2FC(IFNbvscontrol)
IFNyvsIFNb <- add.abs.log2FC(IFNyvsIFNb)
#Filter the significant genes.
filtered_IFNyvscontrol <- IFNyvscontrol[IFNyvscontrol$p_val_adj < 0.01,]
filtered_IFNbvscontrol <- IFNbvscontrol[IFNbvscontrol$p_val_adj < 0.01,]
filtered_IFNyvsIFNb <- IFNyvsIFNb[IFNyvsIFNb$p_val_adj < 0.01,]
#Filter based on log2FC - average for IFNy/IFNbvscontrol and absolute for IFNyvsIFNb
filtered_log2FC_IFNyvscontrol <- filtered_IFNyvscontrol[filtered_IFNyvscontrol$avg_log2FC > 0.3,]
filtered_log2FC_IFNbvscontrol <- filtered_IFNbvscontrol[filtered_IFNbvscontrol$avg_log2FC > 0.3,]
filtered_log2FC_IFNyvsIFNb <- filtered_IFNyvsIFNb[filtered_IFNyvsIFNb$abs_log2FC > 0.3,]
#Filter top 15 highest log2FC of IFNy and IFNb list
top15_IFNy <- filtered_log2FC_IFNyvscontrol %>% top_n(n = 15, wt = avg_log2FC)
top15_IFNb <- filtered_log2FC_IFNbvscontrol %>% top_n(n = 15, wt = avg_log2FC)
top15 <- rbind(top15_IFNy, top15_IFNb)
gene_list_15 <- top15$gene
gene_list_15_unique <- unique(gene_list_15)

#Figure 6A and S8B: DEGs upon interferon stimulation. 
#Change 'ooac' into 'ooac_epithelial' or 'ooac_mesenchymal' when plotting compartment-specific DE genes
Dotplot <- DotPlot(ooac, features = gene_list_15_unique, group.by = "condition_new") + coord_flip()
Dotplot + theme(axis.text.x = element_text(angle = 60, hjust = 1))

#Pathway enrichment
filtered_log2FC_DE_genes <- rbind(filtered_log2FC_IFNyvscontrol, filtered_log2FC_IFNbvscontrol)
#Function
SYMBOLToEntrezENSEMBLGene <- function(f) {
  df <- f
  gene <- df$SYMBOL
  gene.df <- bitr(gene, fromType = "SYMBOL",
                  toType = c("ENTREZID"),
                  OrgDb = org.Hs.eg.db)
  df <- dplyr::full_join(df, gene.df, by = "SYMBOL")
  return(df)
}
filtered_log2FC_DE_genes <- SYMBOLToEntrezENSEMBLGene(filtered_log2FC_DE_genes)
#REACTOME pathway enrichment
DE_pathways <- compareCluster(ENTREZID~comparison, data=filtered_log2FC_DE_genes, fun="enrichPathway", readable=T)
DE_pathways@compareClusterResult$minuslog10p.adjust <- -log10(DE_pathways@compareClusterResult$p.adjust)
#Figure 6A: Plot the top 5 pathways per stimulation in dotplot
dotplot(DE_pathways, 
        x="comparison", 
        showCategory=5, 
        color = "minuslog10p.adjust", 
        font.size = 9, label_format = 50) +
  scale_colour_gradient(low="#D0E0F3", high="#257EE4")+
  labs(colour = "-log10(adj. P-value)") + 
  theme(axis.text.x = element_text(size = 11, angle = 50, hjust=1.0, color = 'black'))

#Figure 6B: DEGs overlap with CeD and IBD-specific DEGs
#Filter based on log2FC - selection was done using absolute log2FC values as we don't have directionality information for the CeD and IBD datasets
#Use DEGs based on complete set (not epithelial or mesenchymal compartment)
filtered_log2FC_IFNyvscontrol <- filtered_IFNyvscontrol[filtered_IFNyvscontrol$abs_log2FC > 0.3,]
filtered_log2FC_IFNbvscontrol <- filtered_IFNbvscontrol[filtered_IFNbvscontrol$abs_log2FC > 0.3,]
#Load data of CeD: RNAseq whole biopsy, Loberman-Nachum et al. (https://pubmed.ncbi.nlm.nih.gov/31700112/), 'Supplementary Dataset 1' Core 878 genes
#Already selected based on FDR and log2FC (878 genes)
CeDvsCTRL <- #read the DEG data
#Load data of IBD (Ulcerative colitis): scRNAseq epithelial layer, Parikh et al. (https://pubmed.ncbi.nlm.nih.gov/30814735/), 'Supplementary Data' Inflamed versus Normal
#Already selected based on FDR (1147 genes)
IBDvsCTRL <- #read the DEG data
#Euler diagram IFN-y, repeat for IFN-b
list <- list(CeDvsCTRL = CeDvsCTRL$gene,
              IBDvsCTRL = IBDvsCTRL$gene,
              IFNyvsCTRL = filtered_log2FC_IFNyvscontrol$SYMBOL)
m.genes <- list_to_matrix(list)
fit <- euler(m.genes)
colnames(m.genes)
rownames(fit)
p.eu <- plot(fit,
             quantities = TRUE,
             edges = T,
             fills = list(fill = c("lavenderblush2", 
                                   "lightblue2", 
                                   "lightsalmon")),
             legend = list(side = "right"),
             labels = F,
             main = "DEGs")
#Test significance 
#Total number of genes in SCT-data slot: 24890 (universe for Fisher's exact test)
#IFNb vs CeD and IBD
fisher.test(x = matrix(c(43,1828,35,22984), nrow = 2), alternative = "greater")
#IFNy vs CeD and IBD
fisher.test(x = matrix(c(93,1778,60,22959), nrow = 2), alternative = "greater")
p_vals <- #vector of p-values
p.adjust(p_vals, method='bonferroni')


#Select genes that can function as specific identifiers for IFNy and IFNb responding cells using above code (top DEGs-under Figure 6A and S8B) for 'complete', 'epithelial' and 'mesenchymal/neural' DE genes
#IFNb: IFI6, IFIT1, MX1
#IFNy: IRF1, APOL6, TAP1

#Verify that the selected genes are present in the list filtered_log2FC_IFNyvsIFNb, and thus different between and specific for IFNy and IFNb
View(filtered_log2FC_IFNyvsIFNb)
#complete: check
#epithelial: check 
#mesenchymal: check

#Generate gene expression plots of the IFNb candidate genes 
FeaturePlot(ooac, features = c('IFI6', 'IFIT1', 'MX1'), ncol=3)
FeaturePlot(ooac[, ooac@meta.data[['condition_new']] %in% c('EM-DM + IFN-β')], features = c('IFI6', 'IFIT1', 'MX1'), ncol=3)

#Generate gene expression plots of the IFNy candidate genes
FeaturePlot(ooac, features = c('IRF1', 'APOL6', 'TAP1'), ncol=3)
FeaturePlot(ooac[, ooac@meta.data[['condition_new']] %in% c('EM-DM + IFN-γ')], features = c('IRF1', 'APOL6', 'TAP1'), ncol=3)

#Figure S8E: 
#IFN type I receptor expression
FeaturePlot(ooac, features = c('IFNAR1', 'IFNAR2'), ncol=2) 
#IFN type II receptor expression
FeaturePlot(ooac, features = c('IFNGR1', 'IFNGR2'), ncol=2) 

#Annotate interferon-responding cells based on gene expression of selected DEGs within the IFNb and IFNy-specific clusters
ooac@meta.data$IFIT1_positive <- ooac@assays$SCT@data['IFIT1', ] > 1
plot_grid(FeaturePlot(ooac, slot = 'data', features=c('IFIT1')), DimPlot(ooac, group.by = 'IFIT1_positive') + scale_color_manual(values = c('gray', 'red')))
ooac@meta.data$MX1_positive <- ooac@assays$SCT@data['MX1', ] > 1.5
plot_grid(FeaturePlot(ooac, slot = 'data', features=c('MX1')), DimPlot(ooac, group.by = 'MX1_positive') + scale_color_manual(values = c('gray', 'red')))
ooac@meta.data$IFI6_positive <- ooac@assays$SCT@data['IFI6', ] > 1
plot_grid(FeaturePlot(ooac, slot = 'data', features=c('IFI6')), DimPlot(ooac, group.by = 'IFI6_positive') + scale_color_manual(values = c('gray', 'red')))
plot_grid(DimPlot(ooac, group.by = 'MX1_positive') + scale_color_manual(values = c('gray', 'red')), DimPlot(ooac, group.by = 'IFIT1_positive') + scale_color_manual(values = c('gray', 'red')), DimPlot(ooac, group.by = 'IFI6_positive') + scale_color_manual(values = c('gray', 'red')))
ooac@meta.data$IRF1_positive <- ooac@assays$SCT@data['IRF1', ] > 2
plot_grid(FeaturePlot(ooac, slot = 'data', features=c('IRF1')), DimPlot(ooac, group.by = 'IRF1_positive') + scale_color_manual(values = c('gray', 'red')))
ooac@meta.data$APOL6_positive <- ooac@assays$SCT@data['APOL6', ] > 1
plot_grid(FeaturePlot(ooac, slot = 'data', features=c('APOL6')), DimPlot(ooac, group.by = 'APOL6_positive') + scale_color_manual(values = c('gray', 'red')))
ooac@meta.data$TAP1_positive <- ooac@assays$SCT@data['TAP1', ] > 1
plot_grid(FeaturePlot(ooac, slot = 'data', features=c('TAP1')), DimPlot(ooac, group.by = 'TAP1_positive') + scale_color_manual(values = c('gray', 'red')))
plot_grid(DimPlot(ooac, group.by = 'TAP1_positive') + scale_color_manual(values = c('gray', 'red')), DimPlot(ooac, group.by = 'IRF1_positive') + scale_color_manual(values = c('gray', 'red')), DimPlot(ooac, group.by = 'APOL6_positive') + scale_color_manual(values = c('gray', 'red')))

#Figure S8C: Visualize all marker gene positive cells
plot_grid(DimPlot(ooac, group.by = 'MX1_positive') + scale_color_manual(values = c('gray', 'red')) + NoLegend(), DimPlot(ooac, group.by = 'IFIT1_positive') + scale_color_manual(values = c('gray', 'red')) + NoLegend(), DimPlot(ooac, group.by = 'IFI6_positive') + scale_color_manual(values = c('gray', 'red')) + NoLegend(),
          DimPlot(ooac, group.by = 'TAP1_positive') + scale_color_manual(values = c('gray', 'red')) + NoLegend(), DimPlot(ooac, group.by = 'IRF1_positive') + scale_color_manual(values = c('gray', 'red')) + NoLegend(), DimPlot(ooac, group.by = 'APOL6_positive') + scale_color_manual(values = c('gray', 'red')) + NoLegend(), ncol=3)

#generate new annotation column 'stimulation'
ooac@meta.data$stimulation <- 'none'
ooac@meta.data[ooac@meta.data$condition_new == 'EM-DM + IFN-β', 'stimulation'] <- 'IFN-β non-responding cell'
ooac@meta.data[ooac@meta.data$condition_new == 'EM-DM + IFN-γ', 'stimulation'] <- 'IFN-γ non-responding cell'
ooac@meta.data[ooac@meta.data$IFIT1_positive == T &
                 ooac@meta.data$MX1_positive == T &
                 ooac@meta.data$IFI6_positive == T &
                 ooac@meta.data$condition_new == 'EM-DM + IFN-β', 'stimulation'] <- 'IFN-β responding cell'
ooac@meta.data[ooac@meta.data$TAP1_positive == T &
                 ooac@meta.data$IRF1_positive == T &
                 ooac@meta.data$APOL6_positive == T &
                 ooac@meta.data$condition_new == 'EM-DM + IFN-γ', 'stimulation'] <- 'IFN-γ responding cell'
ooac@meta.data[ooac@meta.data$condition_new == 'EM-DM', 'stimulation'] <- 'Control'

DimPlot(ooac, reduction = 'umap', group.by = 'stimulation', split.by='condition_new', cols = c('IFN-β responding cell' = '#7CAE00', 'IFN-γ responding cell' = '#F8766D'))
DimPlot(ooac, reduction = 'umap', group.by = 'stimulation', cols = c('IFN-β responding cell' = 'red'))
DimPlot(ooac, reduction = 'umap', group.by = 'stimulation', cols = c('IFN-γ responding cell' = 'red'))

#Figure 6C: generate Dimplot with IFN-responding cells colored and rest gray
DimPlot(ooac, reduction = 'umap', group.by = 'stimulation', split.by='condition_new',
        cols = c('IFN-γ responding cell' = '#F8766D', 'IFN-β responding cell' = '#7CAE00')) +
  labs(title="Interferon-responding cells")+
  theme(plot.title = element_text(hjust = 0.5, size=15), axis.title=element_text(size=12), axis.text=element_text(size=10))

# Perform integration correcting for cell line and stimulation
set.seed(6)
# Generate column with cell_line and condition together
ooac@meta.data <- ooac@meta.data %>%
  mutate(cell_line_condition = paste(cell_line, condition, sep = '_'))
cell_line.list <- SplitObject(ooac, split.by = "cell_line_condition")
cell_line.list <- lapply(X = cell_line.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = cell_line.list, nfeatures = 3000)
cell_line.list <- PrepSCTIntegration(object.list = cell_line.list, anchor.features = features)
cell_line.list <- lapply(X = cell_line.list, FUN = RunPCA, features = features)
cell_line.anchors <- FindIntegrationAnchors(object.list = cell_line.list, normalization.method = "SCT", anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 7)
cell_line.combined.sct <- IntegrateData(anchorset = cell_line.anchors, normalization.method = "SCT", dims = 1:30)
cell_line.combined.sct <- SCTransform(cell_line.combined.sct, assay = "RNA", vst.flavor = "v2")
DefaultAssay(cell_line.combined.sct) <- 'integrated'
cell_line.combined.sct <- RunPCA(cell_line.combined.sct)
cell_line.combined.sct <- RunUMAP(cell_line.combined.sct, assay = 'integrated', reduction = 'pca', dims = 1:30)
cell_line.combined.sct <- FindNeighbors(cell_line.combined.sct, assay = 'integrated', reduction = 'pca', dims = 1:30)
cell_line.combined.sct <- FindClusters(cell_line.combined.sct, resolution = 0.4)

#Figure S8F:
DimPlot(cell_line.combined.sct, reduction = 'umap', group.by = 'stimulation', cols = c('IFN-β responding cell' = '#7CAE00', 'IFN-γ responding cell' = '#F8766D'))
#save as svg file with width=5.8, height=3.5

# change name of Seurat object to ooac
ooac <- cell_line.combined.sct

#Verify active assay to be ‘integrated’
ooac@active.assay

# add annotation based on new clustering and original cell type annotation 
ooac@meta.data$annotation_res0.4 <- 'none'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.4 == 0, 'annotation_res0.4'] <- 'Enterocyte progenitor'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.4 == 1, 'annotation_res0.4'] <- 'Enterocyte type 1 and 2'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.4 == 2, 'annotation_res0.4'] <- 'Enterocyte type 1 and 2'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.4 == 3, 'annotation_res0.4'] <- 'Paneth-like cell'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.4 == 4, 'annotation_res0.4'] <- 'TA/stem cell'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.4 == 5, 'annotation_res0.4'] <- 'Paneth-like cell'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.4 == 6, 'annotation_res0.4'] <- 'Enterocyte type 1 and 2'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.4 == 7, 'annotation_res0.4'] <- 'Mesenchymal-like epithelial cell'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.4 == 8, 'annotation_res0.4'] <- 'Dividing mesenchymal/neural cell'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.4 == 9, 'annotation_res0.4'] <- 'Myofibroblast'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.4 == 10, 'annotation_res0.4'] <- 'Enteroendocrine cell'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.4 == 11, 'annotation_res0.4'] <- 'Goblet cell'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.4 == 12, 'annotation_res0.4'] <- 'Myofibroblast'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.4 == 13, 'annotation_res0.4'] <- 'Enterocyte progenitor'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.4 == 14, 'annotation_res0.4'] <- 'Neuron'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.4 == 15, 'annotation_res0.4'] <- 'Mesenchymal-like epithelial cell'
ooac@meta.data[ooac@meta.data$integrated_snn_res.0.4 == 16, 'annotation_res0.4'] <- 'WNT4-positive neural cell'

DimPlot(ooac, reduction = "umap", group.by = "annotation_res0.4", label=TRUE)

#Re-run the epithelial and mesenchymal compartment annotation with the newly assigned clusters
#Function
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

avg_exp_cluster <- get_average_expression_per_group(ooac, 'seurat_clusters', c('EPCAM', 'VIM'))
DefaultAssay(ooac) <- 'RNA'
ooac@meta.data$epcam_positive <- NA
ooac@meta.data[ooac@meta.data[['seurat_clusters']] %in% avg_exp_cluster[avg_exp_cluster[['EPCAM']] >= 1.0, 'group'], 'epcam_positive'] <- T
ooac@meta.data[ooac@meta.data[['seurat_clusters']] %in% avg_exp_cluster[avg_exp_cluster[['EPCAM']] < 1.0, 'group'], 'epcam_positive'] <- F
ooac@meta.data$vim_positive <- NA
ooac@meta.data[ooac@meta.data[['seurat_clusters']] %in% avg_exp_cluster[avg_exp_cluster[['VIM']] >= 2.7, 'group'], 'vim_positive'] <- T
ooac@meta.data[ooac@meta.data[['seurat_clusters']] %in% avg_exp_cluster[avg_exp_cluster[['VIM']] < 2.7, 'group'], 'vim_positive'] <- F
# add the compartment
ooac@meta.data$compartment <- 'none'
ooac@meta.data[ooac@meta.data$epcam_positive == T, 'compartment'] <- 'epithelial'
ooac@meta.data[ooac@meta.data$vim_positive == T, 'compartment'] <- 'mesenchymal'
#change compartment names
ooac@meta.data$compartment_new <- 'none'
ooac@meta.data[ooac@meta.data$compartment == 'epithelial', 'compartment_new'] <- 'Epithelial cell'
ooac@meta.data[ooac@meta.data$compartment == 'mesenchymal', 'compartment_new'] <- 'Mesenchymal cell'
#generate comparment_new2
ooac@meta.data$compartment_new2 <- 'none'
ooac@meta.data[ooac@meta.data$compartment_new == 'Epithelial cell', 'compartment_new2'] <- 'Epithelial cell'
ooac@meta.data[ooac@meta.data$compartment_new == 'Mesenchymal cell', 'compartment_new2'] <- 'Mesenchymal/Neural cell'
#Neuron has compartment 'none', add manually to mesenchymal/neural compartment
ooac@meta.data[ooac@meta.data$annotation_res0.4 == 'Neuron', 'compartment_new2'] <- 'Mesenchymal/Neural cell'
ooac@meta.data[ooac@meta.data$annotation_res0.4 == 'Neuron', 'compartment_new'] <- 'Mesenchymal cell'

DefaultAssay(ooac) <- 'SCT'

#Figure 6D-F: Composition of interferon-responding cells
ooac$annotation_res0.4 <- as.character(ooac$annotation_res0.4)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ooac <- CellCycleScoring(ooac, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(ooac, reduction = 'umap', group.by = 'Phase')

n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
my_cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

ooac@meta.data <- ooac@meta.data %>%
  mutate(stimulation_cell_line = paste(stimulation, cell_line, sep = '_'))

ooac_epithelial <- ooac[, ooac@meta.data[['compartment_new']] %in% c('Epithelial cell')]
ooac_mesenchymal <- ooac[, ooac@meta.data[['compartment_new']] %in% c('Mesenchymal cell')]

#Plot cell type annotation
par(mar=c(5.1, 4.1, 4.1, 15.1), xpd=TRUE)

#order
order_epithelial <- c('Mesenchymal-like epithelial cell',
                      'Enteroendocrine cell',
                      'Goblet cell',
                      'Paneth-like cell',
                      'Enterocyte type 1 and 2',
                      'Enterocyte progenitor',
                      'TA/stem cell')

order_mesenchymal <- c('WNT4-positive neural cell',
                       'Neuron',
                       'Myofibroblast',
                       'Dividing mesenchymal/neural cell')

#Figure 6E: Composition plot of interferon-responding cell types
#Manually add IFNb and IFNy labels to plot
props <- propeller(clusters=ooac_epithelial$annotation_res0.4, sample=ooac_epithelial$stimulation_cell_line, group=ooac_epithelial$stimulation, transform = 'logit')
props_subset <- subset(props, select=c('PropMean.Control', 'PropMean.IFN.β.responding.cell', 'PropMean.IFN.β.non.responding.cell', 'PropMean.IFN.γ.responding.cell', 'PropMean.IFN.γ.non.responding.cell'))
props_subset <- props_subset[match(order_epithelial, rownames(props_subset)), ]
props_subset_M <- data.matrix(props_subset)
colnames(props_subset_M) <- c('EM-DM', 'Responding', 'Non-responding', 'Responding', 'Non-responding')
barplot(props_subset_M, legend=TRUE, ylab="Proportions", col = alpha(c("#725191","orange","#E31A1C","#FB9A99","#A8D58D", "#1F78B4", "#A6CEE3")), cex.lab=1.2, args.legend = list(x = "topright", inset = c(-0.4, 0.0), cex = 0.9), main="Epithelial composition")

#Figure 6E 
#Manually add IFNb and IFNy labels to plot
props <- propeller(clusters=ooac_mesenchymal$annotation_res0.4, sample=ooac_mesenchymal$stimulation_cell_line, group=ooac_mesenchymal$stimulation, transform = 'logit')
props_subset <- subset(props, select=c('PropMean.Control', 'PropMean.IFN.β.responding.cell', 'PropMean.IFN.β.non.responding.cell', 'PropMean.IFN.γ.responding.cell', 'PropMean.IFN.γ.non.responding.cell'))
props_subset <- props_subset[match(order_mesenchymal, rownames(props_subset)), ]
props_subset_M <- data.matrix(props_subset)
colnames(props_subset_M) <- c('EM-DM', 'Responding', 'Non-responding', 'Responding', 'Non-responding')
barplot(props_subset_M, legend=TRUE, ylab="Proportions", col = alpha(my_cols[58:55]), cex.lab=1.2, args.legend = list(x = "topright", inset = c(-0.4, 0.0), cex = 0.9), main="Mesenchymal/neural composition")

#Figuer 6F: Composition plot cell cycle phase
#Manually add IFNb and IFNy labels to plot
props <- propeller(clusters=ooac_epithelial$Phase, sample=ooac_epithelial$stimulation_cell_line, group=ooac_epithelial$stimulation, transform = 'logit')
props_subset <- subset(props, select=c('PropMean.Control', 'PropMean.IFN.β.responding.cell', 'PropMean.IFN.β.non.responding.cell', 'PropMean.IFN.γ.responding.cell', 'PropMean.IFN.γ.non.responding.cell'))
props_subset <- props_subset[order(row.names(props_subset)), ]
props_subset_M <- data.matrix(props_subset)
colnames(props_subset_M) <- c('EM-DM', 'Responding', 'Non-responding', 'Responding', 'Non-responding')
barplot(props_subset_M, legend=TRUE, ylab="Proportions", col = alpha(my_cols[55:70]), cex.lab=1.2, args.legend = list(x = "topright", inset = c(-0.15, 0.0), cex = 0.9), main="Cell cycle phase: Epithelial cells")

#Figure 6F 
#Manually add IFNb and IFNy labels to plot
props <- propeller(clusters=ooac_mesenchymal$Phase, sample=ooac_mesenchymal$stimulation_cell_line, group=ooac_mesenchymal$stimulation, transform = 'logit')
props_subset <- subset(props, select=c('PropMean.Control', 'PropMean.IFN.β.responding.cell', 'PropMean.IFN.β.non.responding.cell', 'PropMean.IFN.γ.responding.cell', 'PropMean.IFN.γ.non.responding.cell'))
props_subset <- props_subset[order(row.names(props_subset)), ]
props_subset_M <- data.matrix(props_subset)
colnames(props_subset_M) <- c('EM-DM', 'Responding', 'Non-responding', 'Responding', 'Non-responding')
barplot(props_subset_M, legend=TRUE, ylab="Proportions", col = alpha(my_cols[55:70]), cex.lab=1.2, args.legend = list(x = "topright", inset = c(-0.15, 0.0), cex = 0.9), main="Cell cycle phase: Mesenchymal/neural cells")

#Figure S8D:
#Dimplot compartment
DimPlot(ooac, reduction = 'umap', group.by = 'compartment_new2', cols=my_cols[56:57]) +
  labs(title="Stimulation conditions combined")+
  theme(plot.title = element_text(hjust = 0.5, size=15), axis.title=element_text(size=12), axis.text=element_text(size=12))
#Dimplot annotation_res0.34_new2
dark2 <- brewer.pal(4, "Dark2")
paired <- brewer.pal(10, "Paired")
colors <- c(paired, '#FFD633', "#CC9600", dark2)
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
           'IFNβ-responding epithelial cell',
           'IFNγ-responding epithelial cell',
           'Dividing mesenchymal/neural cell',
           'Myofibroblast',
           'Neuron',
           'WNT4-positive neural cell')

ooac$annotation_res0.34_new2 <- as.character(ooac$annotation_res0.34_new2)
ooac$annotation_res0.34_new2 <- factor(ooac$annotation_res0.34_new2, levels = order)
DimPlot(ooac, reduction = 'umap', group.by = 'annotation_res0.34_new2', cols=colors) +
  labs(title="Cellular subtype composition")+
  theme(plot.title = element_text(hjust = 0.5, size=12), axis.title=element_text(size=10), axis.text=element_text(size=10), legend.text = element_text(size = 10))

#Figure 6D: Calculate the percentage of responding and non-responding cells per stimulation (epithelial and mesenchymal/neural cells seperately)
ooac_IFNb_IFNy <- ooac[, ooac@meta.data[['condition_new']] %in% c('EM-DM + IFN-β', 'EM-DM + IFN-γ')]
summary <- ooac_IFNb_IFNy@meta.data %>% 
  group_by(compartment_new, condition_new, cell_line, stimulation) %>% 
  summarise(count = n()) %>%
  mutate(total.cells = sum(count)) %>%
  mutate(pct = count/total.cells)

p <- ggplot(summary, aes(x = compartment_new, y = pct, fill = stimulation)) +
  geom_boxplot(outlier.shape = NA, color = '#564862', linewidth = 0.5) +
  geom_point(size = 1, color = '#564862', position = position_jitterdodge(jitter.width = 0)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 11, angle = 90, hjust = 1, vjust = 0.5, color = 'black'),
        axis.text.y = element_text(size = 11, color = 'black'),
        axis.title.y = element_text(size = 11, color = 'black'),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual('Condition', values = c("#1B9E77", "#D95F02", "#7570B3", "#66C2A5", "#FC8D62", "#8DA0CB")) +
  labs(title = 'Proportion of responding cells',
       y = 'Proportion',
       x = '')
p

#Display in barplot
ooac_IFNb_IFNy@meta.data$cytokine <- 'none'
ooac_IFNb_IFNy@meta.data[ooac_IFNb_IFNy@meta.data$condition_new == 'EM-DM + IFN-β', 'cytokine'] <- 'IFN-β'
ooac_IFNb_IFNy@meta.data[ooac_IFNb_IFNy@meta.data$condition_new == 'EM-DM + IFN-γ', 'cytokine'] <- 'IFN-γ'

ooac_IFNb_IFNy@meta.data$responding <- 'none'
ooac_IFNb_IFNy@meta.data[ooac_IFNb_IFNy@meta.data$stimulation == 'IFN-β non-responding cell', 'responding'] <- 'Non-responding'
ooac_IFNb_IFNy@meta.data[ooac_IFNb_IFNy@meta.data$stimulation == 'IFN-β responding cell', 'responding'] <- 'Responding'
ooac_IFNb_IFNy@meta.data[ooac_IFNb_IFNy@meta.data$stimulation == 'IFN-γ non-responding cell', 'responding'] <- 'Non-responding'
ooac_IFNb_IFNy@meta.data[ooac_IFNb_IFNy@meta.data$stimulation == 'IFN-γ responding cell', 'responding'] <- 'Responding'

ooac_IFNb_IFNy@meta.data <- ooac_IFNb_IFNy@meta.data %>%
  mutate(compartment_cytokine = paste(compartment_new, cytokine, sep = '_'))

ooac_IFNb_IFNy@meta.data <- ooac_IFNb_IFNy@meta.data %>%
  mutate(compartment_cytokine_cell_line = paste(compartment_cytokine, cell_line, sep = '_'))

#Figure 6D
#Manually add Epithelial and mesenchymal labels to plot
par(mar=c(5.1, 4.1, 4.1, 15.1), xpd=TRUE)
props <- propeller(clusters=ooac_IFNb_IFNy$responding, sample=ooac_IFNb_IFNy$compartment_cytokine_cell_line, group=ooac_IFNb_IFNy$compartment_cytokine, transform = 'logit')
props_subset <- subset(props, select=c('PropMean.Epithelial.cell_IFN.β', 'PropMean.Epithelial.cell_IFN.γ', 'PropMean.Mesenchymal.cell_IFN.β', 'PropMean.Mesenchymal.cell_IFN.γ'))
props_subset_M <- data.matrix(props_subset)
colnames(props_subset_M) <- c('IFN-β', 'IFN-γ', 'IFN-β', 'IFN-γ')
barplot(props_subset_M, legend=TRUE, ylab="Proportions", col = alpha(my_cols[58:70]), cex.lab=1.2, args.legend = list(x = "topright", inset = c(-0.35, 0.0), cex = 0.9), main="Proportion of responding cells")
