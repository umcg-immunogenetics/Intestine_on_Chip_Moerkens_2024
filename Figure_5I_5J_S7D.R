###########################################################################################################################
### Figure 5 ##############################################################################################################
###########################################################################################################################

#Figure 5I, 5J and S7D: Cell-cell communication

#libraries
library(Seurat)
library(SeuratData)
library(dplyr)
library(tidyverse)
library(magrittr)
library(liana)
library(circlize)
library(purrr)

# locations of the objects
seurat_objects_loc <- #path to the folder where the seurat object of the media conditions is saved
seurat_object_media <- paste(seurat_objects_loc, 'intestine_on_chip_media.rds', sep = '')

# read object
ooac <- readRDS(seurat_object_media)

# set ident that you want to use for cell-cell communication testing
ooac <- SetIdent(ooac, value = "annotation_res0.34_simple2")

#generate individual condition objects 
ooac_EM <- ooac[, ooac@meta.data[['condition_new']] %in% c('EM')]
ooac_EM_DM <- ooac[, ooac@meta.data[['condition_new']] %in% c('EM-DM')]
ooac_DM <- ooac[, ooac@meta.data[['condition_new']] %in% c('DM')]

#run CellPhoneDBv2 using the Liana wrapper
cpdb <- liana_wrap(ooac,
                   method = 'cellphonedb',
                   resource = 'CellPhoneDB')

#Filter the liana outputs based on p-value
cpdb <- cpdb %>% filter(pvalue <= 0.01) 

# extract pathway annotation per ligand from CellChatDB 
library(CellChat)
#Extract human database
CellChatDB <- CellChatDB.human 
#Extract the pathway and ligand gene information
CellChatDB_pathway = data.frame(CellChatDB$interaction$ligand, CellChatDB$interaction$pathway_name)
#Remove duplicate rows
CellChatDB_pathway <- CellChatDB_pathway %>% distinct()
#Some ligands are assigned to multiple pathways, which creates duplicate values after merging with liana output
#Verify presence of these genes in the liana output dataframe
CellChatDB_pathway$CellChatDB.interaction.ligand[duplicated(CellChatDB_pathway$CellChatDB.interaction.ligand)]
#Result: CDH1 and GCG are present in cpdb liana output, choose one pathway for each of these genes
CellChatDB_pathway <- CellChatDB_pathway[-(442),] # removing row number
CellChatDB_pathway <- CellChatDB_pathway[-(335),] # removing row number

#If genes are not annotated that clearly belong to one of the pathways, indicate as vector below
setdiff(cpdb$ligand, CellChatDB_pathway$CellChatDB.interaction.ligand)
vec1 <- c('COL5A2', 'COLLAGEN')
vec2 <- c('COL27A1', 'COLLAGEN')
vec3 <- c('COL18A1', 'COLLAGEN')
vec4 <- c('COL17A1', 'COLLAGEN')
vec5 <- c('COL14A1', 'COLLAGEN')
vec6 <- c('FGFR4', 'FGF')
vec7 <- c('FGFR2', 'FGF')
vec8 <- c('FGFR3', 'FGF')
vec9 <- c('CXCL17', 'CXCL')
vec10 <- c('BMP3', 'BMP')
vec11 <- c('EPHB6', 'EPHB')
vec12 <- c('COL21A1', 'COLLAGEN')
vec13 <- c('COL16A1', 'COLLAGEN')
vec14 <- c('COL12A1', 'COLLAGEN')
vec15 <- c('COL26A1', 'COLLAGEN')
vec16 <- c('COL7A1', 'COLLAGEN')
vec17 <- c('COL5A1', 'COLLAGEN')
vec18 <- c('COL3A1', 'COLLAGEN')
vec19 <- c('COL13A1', 'COLLAGEN')
vec20 <- c('COL11A1', 'COLLAGEN')
vec21 <- c('COL8A1', 'COLLAGEN')
vec22 <- c('RSPO3', 'WNT')
vec23 <- c('CADM4', 'CADM')
vec24 <- c('CEACAM6', 'CEACAM')
vec25 <- c('LTB', 'LT')
vec26 <- c('PLXNB1', 'SEMA')
CellChatDB_pathway <- rbind(CellChatDB_pathway, vec1, vec2, vec3, vec4, vec5, vec6, vec7, vec8, vec9, vec10, vec11, vec12, vec13, vec14, vec15, vec16, vec17, vec18, vec19, vec20, vec21, vec22, vec23, vec24, vec25, vec26)

#Add pathway annotation to the liana object
CellChatDB_pathway$ligand <- CellChatDB_pathway$CellChatDB.interaction.ligand
cpdb_pathway <- merge(cpdb, CellChatDB_pathway, by = "ligand", all.x=TRUE)
names(cpdb_pathway)[names(cpdb_pathway) == "CellChatDB.interaction.pathway_name"] <- "pathway_name"

# Reordering cell types
cpdb_pathway$source <- as.factor(cpdb_pathway$source)
cpdb_pathway$target <- as.factor(cpdb_pathway$target)
order <- c('TA/stem cell',
           'Enterocyte',
           'Paneth-like cell',
           'Goblet cell',
           'Enteroendocrine cell',
           'Mesenchymal-like epithelial cell',
           'Dividing mesenchymal/neural cell',
           'Myofibroblast',
           'Neuron',
           'WNT4-positive neural cell')
cpdb_pathway$source <- factor(cpdb_pathway$source, levels = order)
cpdb_pathway$target <- factor(cpdb_pathway$target, levels = order)

#Create objects for the pathways with the most interactions
cpdb_BMP <- cpdb_pathway[cpdb_pathway$pathway_name %in% c('BMP'),]
#Non canonical and canonical WNT pathways were grouped
cpdb_WNT <- cpdb_pathway[cpdb_pathway$pathway_name %in% c('WNT', 'ncWNT'),]
cpdb_COLLAGEN <- cpdb_pathway[cpdb_pathway$pathway_name %in% c('COLLAGEN'),]
cpdb_NOTCH <- cpdb_pathway[cpdb_pathway$pathway_name %in% c('NOTCH'),]
cpdb_LAMININ <- cpdb_pathway[cpdb_pathway$pathway_name %in% c('LAMININ'),]
cpdb_FGF <- cpdb_pathway[cpdb_pathway$pathway_name %in% c('FGF'),]
#All identified EPH pahtways were grouped
cpdb_EPH <- cpdb_pathway[cpdb_pathway$pathway_name %in% c('EPHA', 'EPHB'),]
#All identified SEMA pathways were grouped
cpdb_SEMA <- cpdb_pathway[cpdb_pathway$pathway_name %in% c('SEMA3', 'SEMA4'),]

#Figure 5J and S7D: Heatmap of interaction frequency between cell types
heat_freq(cpdb_BMP, font_size=12)

#Figure 5I: Interaction barplot grouped by pathways
#Run CellPhoneDBv2 per medium condition
liana_function <- function(object_test){
  liana_output <- liana_wrap(object_test, method = 'cellphonedb', resource = 'CellPhoneDB')
  filtered_liana_output <- liana_output %>% filter(pvalue <= 0.01)
  cpdb_pathway <- merge(filtered_liana_output, CellChatDB_pathway, by = "ligand", all.x=TRUE)
  names(cpdb_pathway)[names(cpdb_pathway) == "CellChatDB.interaction.pathway_name"] <- "pathway_name"
  cpdb_pathway$source <- as.factor(cpdb_pathway$source)
  cpdb_pathway$target <- as.factor(cpdb_pathway$target)
  cpdb_pathway$source <- factor(cpdb_pathway$source, levels = order)
  cpdb_pathway$target <- factor(cpdb_pathway$target, levels = order)
  return(cpdb_pathway)
}

cpdb_EM <- liana_function(ooac_EM)
cpdb_EM_DM <- liana_function(ooac_EM_DM)
cpdb_DM <- liana_function(ooac_DM)

#apply same pathway grouping as above
cpdb_EM$pathway_name_2 <- cpdb_EM$pathway_name
cpdb_EM$pathway_name_2[cpdb_EM$pathway_name_2 == 'ncWNT'] <- 'WNT'
cpdb_EM$pathway_name_2[cpdb_EM$pathway_name_2 == 'EPHA'] <- 'EPH'
cpdb_EM$pathway_name_2[cpdb_EM$pathway_name_2 == 'EPHB'] <- 'EPH'
cpdb_EM$pathway_name_2[cpdb_EM$pathway_name_2 == 'SEMA3'] <- 'SEMA'
cpdb_EM$pathway_name_2[cpdb_EM$pathway_name_2 == 'SEMA4'] <- 'SEMA'
cpdb_EM$pathway_name_2[is.na(cpdb_EM$pathway_name_2)] <- "Other"

cpdb_EM_DM$pathway_name_2 <- cpdb_EM_DM$pathway_name
cpdb_EM_DM$pathway_name_2[cpdb_EM_DM$pathway_name_2 == 'ncWNT'] <- 'WNT'
cpdb_EM_DM$pathway_name_2[cpdb_EM_DM$pathway_name_2 == 'EPHA'] <- 'EPH'
cpdb_EM_DM$pathway_name_2[cpdb_EM_DM$pathway_name_2 == 'EPHB'] <- 'EPH'
cpdb_EM_DM$pathway_name_2[cpdb_EM_DM$pathway_name_2 == 'SEMA3'] <- 'SEMA'
cpdb_EM_DM$pathway_name_2[cpdb_EM_DM$pathway_name_2 == 'SEMA4'] <- 'SEMA'
cpdb_EM_DM$pathway_name_2[is.na(cpdb_EM_DM$pathway_name_2)] <- "Other"

cpdb_DM$pathway_name_2 <- cpdb_DM$pathway_name
cpdb_DM$pathway_name_2[cpdb_DM$pathway_name_2 == 'ncWNT'] <- 'WNT'
cpdb_DM$pathway_name_2[cpdb_DM$pathway_name_2 == 'EPHA'] <- 'EPH'
cpdb_DM$pathway_name_2[cpdb_DM$pathway_name_2 == 'EPHB'] <- 'EPH'
cpdb_DM$pathway_name_2[cpdb_DM$pathway_name_2 == 'SEMA3'] <- 'SEMA'
cpdb_DM$pathway_name_2[cpdb_DM$pathway_name_2 == 'SEMA4'] <- 'SEMA'
cpdb_DM$pathway_name_2[is.na(cpdb_DM$pathway_name_2)] <- "Other"

#Calculate percentage of interactions for each pathway
summary_EM <- cpdb_EM %>% 
  group_by(pathway_name_2) %>% 
  summarise(count = n()) %>%
  mutate(total.interactions = sum(count)) %>%
  mutate(pct = (count/total.interactions)*100)
summary_EM_DM <- cpdb_EM_DM %>% 
  group_by(pathway_name_2) %>% 
  summarise(count = n()) %>%
  mutate(total.interactions = sum(count)) %>%
  mutate(pct = (count/total.interactions)*100)
summary_DM <- cpdb_DM %>% 
  group_by(pathway_name_2) %>% 
  summarise(count = n()) %>%
  mutate(total.interactions = sum(count)) %>%
  mutate(pct = (count/total.interactions)*100)

#Get the pathway names that occur fewer than 2%
summary_EM <- summary_EM %>% mutate(
  pathway_name_2 = case_when(
    pct >= 2 ~ "Other",
    TRUE ~ pathway_name_2
  )
)

summary_EM_DM <- summary_EM_DM %>% mutate(
  pathway_name_2 = case_when(
    pct >= 2 ~ "Other",
    TRUE ~ pathway_name_2
  )
)

summary_DM <- summary_DM %>% mutate(
  pathway_name_2 = case_when(
    pct >= 2 ~ "Other",
    TRUE ~ pathway_name_2
  )
)

#Group the pathways that occur fewer than 2% in all conditions as 'Other'
removed_pathways <- c('AGT', 'APP', 'AVP', 'BAFF', 'CADM','CD137', 'CCL', 'CDH','CEACAM', 'CSF', 'CXCL', 'DESMOSOME', 'EDA', 'EGF', 'ENHO', 'GHRELIN', 'GIPR', 'GUCA', 'GALANIN', 'GAS', 'GRN', 'INSULIN', 'KIT', 'LIFR', 'MIF', 'MK', 'NCAM', 'NPR2', 'NECTIN', 'NODAL', 'NPR1', 'NRG', 'OXT', 'PDGF', 'PROS', 'PSAP', 'PTN', 'PACAP', 'SCT', 'SPP1', 'SOMATOSTATIN', 'TENASCIN', 'THBS', 'TRAIL', 'TWEAK', 'VISFATIN', 'VEGI', 'VTN')
cpdb_EM$pathway_name_2[cpdb_EM$pathway_name_2 %in% removed_pathways] <- 'Other'
cpdb_EM_DM$pathway_name_2[cpdb_EM_DM$pathway_name_2 %in% removed_pathways] <- 'Other'
cpdb_DM$pathway_name_2[cpdb_DM$pathway_name_2 %in% removed_pathways] <- 'Other'

# Generate new percentages for pathway abundance
summary_EM <- cpdb_EM %>% 
  group_by(pathway_name_2) %>% 
  summarise(count = n()) %>%
  mutate(total.interactions = sum(count)) %>%
  mutate(pct = (count/total.interactions)*100)
summary_EM_DM <- cpdb_EM_DM %>% 
  group_by(pathway_name_2) %>% 
  summarise(count = n()) %>%
  mutate(total.interactions = sum(count)) %>%
  mutate(pct = (count/total.interactions)*100)
summary_DM <- cpdb_DM %>% 
  group_by(pathway_name_2) %>% 
  summarise(count = n()) %>%
  mutate(total.interactions = sum(count)) %>%
  mutate(pct = (count/total.interactions)*100)

#Merge the summary files into one file
summary_EM <- as.data.frame(summary_EM)
rownames(summary_EM) <- summary_EM$pathway_name_2
summary_EM = subset(summary_EM, select = -c(count,total.interactions, pathway_name_2) )
colnames(summary_EM)[which(names(summary_EM) == "pct")] <- "EM"

summary_EM_DM <- as.data.frame(summary_EM_DM)
rownames(summary_EM_DM) <- summary_EM_DM$pathway_name_2
summary_EM_DM = subset(summary_EM_DM, select = -c(count,total.interactions, pathway_name_2) )
colnames(summary_EM_DM)[which(names(summary_EM_DM) == "pct")] <- "EM-DM"

summary_DM <- as.data.frame(summary_DM)
rownames(summary_DM) <- summary_DM$pathway_name_2
summary_DM = subset(summary_DM, select = -c(count,total.interactions, pathway_name_2) )
colnames(summary_DM)[which(names(summary_DM) == "pct")] <- "DM"

#very important to set the 'sort=FALSE' to preserve the order of the genes, as this will influence the clustering in Heatmap
summary_merged <- merge(summary_EM, summary_EM_DM, by=0, all=TRUE, sort=FALSE)
rownames(summary_merged) <- summary_merged$Row.names
summary_merged <- summary_merged[,-1]
summary_merged <- merge(summary_merged, summary_DM, by=0, all=TRUE, sort=FALSE)
rownames(summary_merged) <- summary_merged$Row.names
summary_merged <- summary_merged[,-1]

#Plot the percentages
summary_merged <- summary_merged[order(summary_merged$EM),]
summary_merged <- as.matrix(summary_merged)
n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
my_cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=TRUE)
barplot(summary_merged, col = alpha(my_cols[55:70]) , border="white", xlab="Condition", ylab="Proportion of interactions (%)", legend=TRUE, args.legend = list(x = "topright", inset = c(-0.3, 0.0), cex = 1.0), main="Interaction pathways")

