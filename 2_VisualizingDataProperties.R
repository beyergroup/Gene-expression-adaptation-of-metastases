## Code from Luise Nagel

library(ggplot2)
library(ggpubr)
library(Seurat)
library(patchwork)

# Loading and sorting filtered and preprocessed data 
mySE <- readRDS("../../Data/AllDonors_Filtered/SE_filteredAndPreprocessed.rds")

mySE_UMAP <- as.data.frame(mySE@colData)
mySE_UMAP$UMAP1 <- as.numeric(Embeddings(mySE, reduction = "umap")[, "UMAP_1"])
mySE_UMAP$UMAP2 <- as.numeric(Embeddings(mySE, reduction = "umap")[, "UMAP_2"])
mySE_UMAP$SeuratCluster <- as.factor(mySE$SeuratCluster)


# Export data
## Extended Data Table 1
extDataTab2 <- mySE_UMAP
extDataTab2 <- extDataTab1[, c(1,2,3,4,5,8,11,12,7,9,6,15,13,14)]
write.csv(extDataTab1, paste0("../../Paper/ExtendedDataTables/ExtendedDataTable2_", myDate, ".csv"))

# Plotting overview of cell types, donors etc
PlotData <- mySE@colData
PlotData$TissueSum <- sapply(PlotData$Tissue, \(x){
    if(x == "acHT"|x == "scHT"){
        "cHT"
    }
    else if(x == "scPT"){
        "cPT"
    }
    else{
        x
    }
})

## Figure 1a: # of cells per donors per origin tissue
plotDonorTissue <- ggplot(mapping = aes(x = PlotData$Donor, color = PlotData$TissueSum, 
                                        fill = PlotData$TissueSum)) +
    geom_bar() +
    scale_color_viridis_d(option = "D", begin = 0, end = 0.9) +
    scale_fill_viridis_d(option = "D", begin = 0, end = 0.9) +
    xlab("Donor") +
    ylab("Number of cells") +
    theme_bw(base_size = 11) 


## Extended Data Figure 1a: # of cells per donors per cell type
ggplot(mapping = aes(x = PlotData$Celltypes, color = PlotData$Donor, 
                                           fill = PlotData$Donor)) +
    geom_bar() +
    scale_color_manual(values = c("#000000", "#701e81", "#bd3976", "#f7725b", "#fecc8f")) +
    scale_fill_manual(values = c("#000000", "#701e81", "#bd3976", "#f7725b", "#fecc8f")) +
    xlab("Cell types") +
    ylab("Number of cells") +
    scale_x_discrete(guide = guide_axis(angle = 70)) +
    theme_bw(base_size = 11) 


## Extended Data Figure 1c: # of cells per donors per seurat cluster 
ggplot(mapping = aes(x = PlotData$SeuratCluster, color = PlotData$Donor, 
                     fill = PlotData$Donor)) +
    geom_bar() +
    scale_color_manual(values = c("#000000", "#701e81", "#bd3976", "#f7725b", "#fecc8f")) +
    scale_fill_manual(values = c("#000000", "#701e81", "#bd3976", "#f7725b", "#fecc8f")) +
    xlab("Cluster numbers") +
    ylab("Number of cells") +
    theme_bw(base_size = 15) 

# Plotting UMAPS
## Coloring cell type
mySE_UMAP$Celltype4PaperVisual <- mySE_UMAP$Celltype 
mySE_UMAP$Celltype4PaperVisual[mySE_UMAP$Celltype4PaperVisual == "ENS"] <- "Enteric nervous system cells"
mySE_UMAP$Celltype4PaperVisual[mySE_UMAP$Celltype4PaperVisual == "LiverEpithelial"] <- "Liver epithelial cells"
mySE_UMAP$Celltype4PaperVisual[mySE_UMAP$Celltype4PaperVisual == "Epithelium"] <- "Liver epithelial cells"
mySE_UMAP$Celltype4PaperVisual[mySE_UMAP$Celltype4PaperVisual == "ColonEpithelial"] <- "Colon epithelial cells"
mySE_UMAP$Celltype4PaperVisual[mySE_UMAP$Celltype4PaperVisual == "Endothelium"] <- "Endothelium cells"
mySE_UMAP$Celltype4PaperVisual[mySE_UMAP$Celltype4PaperVisual == "BCell"] <- "B cells"
mySE_UMAP$Celltype4PaperVisual[mySE_UMAP$Celltype4PaperVisual == "TCell"] <- "T cells"
mySE_UMAP$Celltype4PaperVisual[mySE_UMAP$Celltype4PaperVisual == "Myofibroblasts"] <- "Fibroblasts"
mySE_UMAP$Celltype4PaperVisual[mySE_UMAP$Celltype4PaperVisual == "Stroma/CAF"] <- "Fibroblasts"
mySE_UMAP$Celltype4PaperVisual[mySE_UMAP$Celltype4PaperVisual == "Tuft"] <- "Tuft cells"
mySE_UMAP$Celltype4PaperVisual[mySE_UMAP$Celltype4PaperVisual == "PrimaryTumor"] <- "Primary colon tumor"
mySE_UMAP$Celltype4PaperVisual[mySE_UMAP$Celltype4PaperVisual == "Metastasis"] <- "Liver metastasis"


## Figure 1c: UMAP with cell types colored
ggplot(mySE_UMAP, aes(x = UMAP1, y = UMAP2, color = Celltype4PaperVisual)) +
    geom_point(alpha = 0.1) +
    scale_color_manual(values = c("#a3c7ec", "#197a21", "#41197a", "#29cba3",
                                  "#8621a3", "#7ad151", "#c11509", "#1d558e",
                                  "#ec5813", "#3886d6", "#123559", "#a3a3a3")) +    
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    theme_bw(base_size = 17) +
    theme(aspect.ratio = 1)


##  Figure 1d: UMAP colored by sample tissues
ggplot(mySE_UMAP, aes(x = UMAP1, y = UMAP2, color = TissueSum)) +
    geom_point(alpha = 0.2) +
    scale_color_viridis_d(begin = 0, end = 0.9, option = "D") +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    theme_bw(base_size = 17) +
    theme(aspect.ratio = 1)


## Extended Data Figure 1b: UMAP colored by donor
ggplot(mySE_UMAP, aes(x = UMAP1, y = UMAP2, color = Donor, alpha = Donor)) +
    geom_point() +
    scale_color_manual(values = c("#000000", "#701e81", "#bd3976", "#f7725b", "#fecc8f")) +
    scale_alpha_manual(values = c(0.9, 0.7, 0.5, 0.3, 0.1)) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    theme_bw(base_size = 17) +
    theme(aspect.ratio = 1)


## Extended Data Figure 1d: UMAP colored by seurat cluster
ggplot(mySE_UMAP, aes(x = UMAP1, y = UMAP2, color = SeuratCluster)) +
    geom_point(alpha = 0.3) +
    scale_color_viridis_d(option = "B", begin = 0, end = 1) +
    #scale_color_viridis_c() +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    theme_bw(base_size = 17) 




