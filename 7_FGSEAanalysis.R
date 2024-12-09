# Code from Luise Nagel

## Summarizing data.frames from earlier analysis required (see previous R scripts)

detach("package:Seurat", unload = T) ## unload Seurat so patchwork can work
library(cluster)
library(factoextra)
library(fgsea)
library(ggplot2)
library(ggpubr)
library(MASS)
library(patchwork)
library(qusage)
library(reshape2)

# Loading data
pathwaysC2CP <- read.gmt("../../CanonicalPathways/C2_CP_030123/c2.cp.v2022.1.Hs.symbols.gmt")

## Keeping only pathways from specific datab bases
pathwaysC2CP <- pathwaysC2CP[sapply(gsub("_.*", "", names(pathwaysC2CP)), \(x){
    x %in% c("BIOCARTA", "KEGG", "PID", "REACTOME")
})]


# Performing FGSEA on our scRNA-seq data
## Perform donor-wise FGSEA on MT & PT
FGSEA_MTPT_MTPTonly_D1D1 <- setNames(sumDE_MTPT$MTPT_D1D1, row.names(sumDE_MTPT))
FGSEA_MTPT_MTPTonly_D1D1 <- fgseaMultilevel(pathways = pathwaysC2CP, minSize = 10,
                                            stats = FGSEA_MTPT_MTPTonly_D1D1, eps = 0)

FGSEA_MTPT_MTPTonly_D5D5 <- setNames(sumDE_MTPT$MTPT_D5D5, row.names(sumDE_MTPT))
FGSEA_MTPT_MTPTonly_D5D5 <- fgseaMultilevel(pathways = pathwaysC2CP, minSize = 10,
                                            stats = FGSEA_MTPT_MTPTonly_D5D5, eps = 0)

sumFGSEA_MTPT <-  data.frame(MTPT_D1D1_NES = FGSEA_MTPT_MTPTonly_D1D1$NES,
                             MTPT_D5D5_NES = FGSEA_MTPT_MTPTonly_D5D5$NES,
                             MTPT_D1D1_Padj = FGSEA_MTPT_MTPTonly_D1D1$padj,
                             MTPT_D5D5_Padj = FGSEA_MTPT_MTPTonly_D5D5$padj,
                             row.names = FGSEA_MTPT_MTPTonly_D1D1$pathway)

sumFGSEA_MTPT$MTPT_specific <- apply(sumFGSEA_MTPT[, 1:2] > 0, 1, \(x){
    x <- paste(as.numeric(x),collapse="")
    if(x == "11"){
        "MT"
    }
    else if(x == "00"){
        "PT"
    }
    else{
        FALSE
    }
})

sumFGSEA_MTPT$meanNES <- rowMeans(sumFGSEA_MTPT[, c("MTPT_D1D1_NES", "MTPT_D5D5_NES")])
sumFGSEA_MTPT$meanPadj <- rowMeans(sumFGSEA_MTPT[, c("MTPT_D1D1_Padj", "MTPT_D5D5_Padj")])

## Plotting example pathways
### Figure 2c: Pathways higher in PT
examplePathwaysPT <- sumFGSEA_MTPT[c("REACTOME_RESPIRATORY_ELECTRON_TRANSPORT", 
                                     "REACTOME_MITOCHONDRIAL_BIOGENESIS",
                                     "REACTOME_SIGNALING_BY_WNT",
                                     "REACTOME_INTERLEUKIN_12_SIGNALING",
                                     "REACTOME_TRANSCRIPTIONAL_REGULATION_BY_RUNX2"),]
examplePathwaysPT$Pathways <- c("Respiratory electron transport", 
                                "Mitochondrial biogenesis",
                                "Signaling by WNT",
                                "Interleukin 12 signaling",
                                "Transkriptional regulation by RUNX2")

examplePathwaysPT <- examplePathwaysPT[order(examplePathwaysPT$meanNES, decreasing = T),]
examplePathwaysPT$Pathways <- stringr::str_wrap(examplePathwaysPT$Pathways, 20)
examplePathwaysPT$Pathways <- factor(examplePathwaysPT$Pathways, levels = examplePathwaysPT$Pathways)

ggplot(examplePathwaysPT, aes(y = meanNES, x = Pathways)) +
    geom_bar(stat = "identity", color = "#440154", fill = "#440154") +
    coord_flip()  +
    ylim(-2.5, 0) +
    ylab("mean NES") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")


### Figure 2c: Pathways higher in MT
examplePathwaysMT <- sumFGSEA_MTPT[c("PID_REG_GR_PATHWAY", 
                                     "KEGG_MAPK_SIGNALING_PATHWAY",
                                     "BIOCARTA_P53HYPOXIA_PATHWAY",
                                     "BIOCARTA_PPARA_PATHWAY",
                                     "REACTOME_GLYCOLYSIS"),]
examplePathwaysMT$Pathways <- c("REG GR pathway", 
                                "MAPK signaling pathway",
                                "P53 hypoxia pathway",
                                "PPAR-alpha pathway",
                                "Glycolysis")

examplePathwaysMT <- examplePathwaysMT[order(examplePathwaysMT$meanNES, decreasing = F),]
examplePathwaysMT$Pathways <- stringr::str_wrap(examplePathwaysMT$Pathways, 20)
examplePathwaysMT$Pathways <- factor(examplePathwaysMT$Pathways, levels = examplePathwaysMT$Pathways)

ggplot(examplePathwaysMT, aes(y = meanNES, x = Pathways)) +
    geom_bar(stat = "identity", color = "#7ad151", fill = "#7ad151") +
    coord_flip()  +
    ylim(0, 2.5) +
    ylab("mean NES") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")



## Perform donor-wise FGSEA on MTLE & PTCE
FGSEA_MTPT_D1D1 <- setNames(sumDE_LECE_MTPT$MTPT_D1D1, row.names(sumDE_LECE_MTPT))
FGSEA_MTPT_D1D1 <- fgseaMultilevel(pathways = pathwaysC2CP, minSize = 10,
                                   stats = FGSEA_MTPT_D1D1, eps = 0)

FGSEA_MTPT_D5D5 <- setNames(sumDE_LECE_MTPT$MTPT_D5D5, row.names(sumDE_LECE_MTPT))
FGSEA_MTPT_D5D5 <- fgseaMultilevel(pathways = pathwaysC2CP, minSize = 10,
                                   stats = FGSEA_MTPT_D5D5, eps = 0)

FGSEA_LECE_D1D1 <- setNames(sumDE_LECE_MTPT$LECE_D1D1, row.names(sumDE_LECE_MTPT))
FGSEA_LECE_D1D1 <- fgseaMultilevel(pathways = pathwaysC2CP, minSize = 10,
                                   stats = FGSEA_LECE_D1D1, eps = 0)

FGSEA_LECE_D3D5 <- setNames(sumDE_LECE_MTPT$LECE_D3D5, row.names(sumDE_LECE_MTPT))
FGSEA_LECE_D3D5 <- fgseaMultilevel(pathways = pathwaysC2CP, minSize = 10,
                                   stats = FGSEA_LECE_D3D5, eps = 0)

sumFGSEA_MTLEPTCE <- data.frame(MTPT_D1D1_NES = FGSEA_MTPT_D1D1$NES,
                                MTPT_D5D5_NES = FGSEA_MTPT_D5D5$NES,
                                LECE_D1D1_NES = FGSEA_LECE_D1D1$NES,
                                LECE_D3D5_NES = FGSEA_LECE_D3D5$NES,
                                row.names = FGSEA_LECE_D1D1$pathway)


### Add information on agreement in directional of NES between donors
sumFGSEA_MTLEPTCE$Quadrant <- apply(sumFGSEA_MTLEPTCE[, 1:4] > 0, 1, \(x){
    x <- paste(as.numeric(x),collapse="")
    if(x == "1111"){
        "TR"
    }
    else if(x == "0000"){
        "BL"
    }
    else if(x == "0011"){
        "BR"
    }
    else if(x == "1100"){
        "TL"
    }
    else{
        FALSE
    }
})

sumFGSEA_MTLEPTCE$LECE_MTPT_specific <- sapply(sumFGSEA_MTLEPTCE$Quadrant, \(x){
    if(x=="TR"){
        "MTLE"
    }
    else if(x=="BL"){
        "PTCE"
    }
    else{
        FALSE
    }
})

## Clustering pathways
### Estimating jacc distances
estimateJaccDist <- function(ClusterList, clusterSize = 6){
    stopifnot(is.list(ClusterList),clusterSize >= 1)
    require(parallel)
    cl <- makeCluster(clusterSize)
    tryCatch({
        clusterExport(cl ,"ClusterList", environment())
        jaccMat <- do.call(rbind, parLapply(cl, names(ClusterList), \(i){
            do.call(cbind, lapply(names(ClusterList), \(j){
                length(intersect(ClusterList[[i]],ClusterList[[j]]))/
                    length(union(ClusterList[[i]],ClusterList[[j]]))
            }))
        }))
    }, error = function(e) print(e), finally = stopCluster(cl))
    rownames(jaccMat) <- colnames(jaccMat) <- names(ClusterList)
    return(jaccMat)
}

jaccDist_C2CP <- estimateJaccDist(pathwaysC2CP)

### Clustering pathways higher expressed in PTCE and MTLE
jaccDist_C2CPclust <- jaccDist_C2CP
jaccDist_C2CPclust <- 1-jaccDist_C2CPclust
jaccDist_C2CPclust[jaccDist_C2CPclust == 0] <- 0.001 # convert 0 to 0.001 because algorithm cannot handle 0
PTCECluster <- row.names(sumFGSEA_MTLEPTCE)[sumFGSEA_MTLEPTCE$LECE_PTPT_specific == "PTCE"]
MTLECluster <- row.names(sumFGSEA_MTLEPTCE)[sumFGSEA_MTLEPTCE$LECE_PTPT_specific == "MTLE"]

isoMDS_C2CP_PTCE <- isoMDS(as.dist(jaccDist_C2CPclust[PTCECluster, PTCECluster]))
isoMDS_C2CP_MTLE <- isoMDS(as.dist(jaccDist_C2CPclust[MTLECluster, MTLECluster]))

## Choosing number of clusters & performing k-means clustring
### PTCE pathways
set.seed(42)
wssPTCE <- fviz_nbclust(isoMDS_C2CP_PTCE$points, kmeans, method = "wss", k.max = 15)  # elbow plot
swPTCE <- fviz_nbclust(isoMDS_C2CP_PTCE$points, kmeans, method = "silhouette", k.max = 15) # silhouette width
gapPTCE <- fviz_gap_stat(clusGap(isoMDS_C2CP_PTCE$points, FUN = kmeans, nstart = 25, K.max = 15, B = 50))

set.seed(42)
kmeans_C2CP_PTCE <- kmeans(isoMDS_C2CP_PTCE$points, 8, nstart = 1)

### MTLE pathways
set.seed(42)
wssMTLE <- fviz_nbclust(isoMDS_C2CP_MTLE$points, kmeans, method = "wss", k.max = 15)  # elbow plot
swMTLE <- fviz_nbclust(isoMDS_C2CP_MTLE$points, kmeans, method = "silhouette", k.max = 15) # silhouette width
gapMTLE <- fviz_gap_stat(clusGap(isoMDS_C2CP_MTLE$points, FUN = kmeans, nstart = 25, K.max = 15, B = 50))

set.seed(42)
kmeans_C2CP_MTLE <- kmeans(isoMDS_C2CP_MTLE$points, 5, nstart = 1)


# FGSEA external data
## mArray data
FGSEA_mArray <- setNames(mArrayDE$logFC, row.names(mArrayDE))
FGSEA_mArray <- fgseaMultilevel(pathways = pathwaysC2CP, minSize = 10,
                                stats = FGSEA_mArray, eps = 0)
FGSEA_mArray <- as.data.frame(FGSEA_mArray)
row.names(FGSEA_mArray) <- FGSEA_mArray$pathway


## external scRNAseq data
FGSEA_extScRNA_MTPT <- setNames(sumExtScRNAseq$meanlogFC_ExtSc_MTPT, row.names(sumExtScRNAseq))
FGSEA_extScRNA_MTPT <- fgseaMultilevel(pathways = pathwaysC2CP, minSize = 10,
                                       stats = FGSEA_extScRNA_MTPT, eps = 0)
FGSEA_extScRNA_MTPT <- as.data.frame(FGSEA_extScRNA_MTPT)
row.names(FGSEA_extScRNA_MTPT) <- FGSEA_extScRNA_MTPT$pathway

FGSEA_extSCRNA_MTPT_singleDonors <- apply(sumExtScRNAseq, 2, \(x){
    x <- fgseaMultilevel(pathways = pathwaysC2CP, minSize = 10,
                         stats = x, eps = 0)
    x <- as.data.frame(x)
    row.names(x) <- x$pathway
    return(x[,"NES", drop = F])
})
namesFGSEA_extSCRNA_MTPT_singleDonors <- names(FGSEA_extSCRNA_MTPT_singleDonors)
FGSEA_extSCRNA_MTPT_singleDonors <- do.call(cbind, FGSEA_extSCRNA_MTPT_singleDonors)
colnames(FGSEA_extSCRNA_MTPT_singleDonors) <- namesFGSEA_extSCRNA_MTPT_singleDonors



## Figure 5a: Heatmap all MTLE & PTCE pathways
plotData <- data.frame(Pathways = c(names(kmeans_C2CP_MTLE$cluster),
                                    names(kmeans_C2CP_PTCE$cluster)),
                       Cluster = c(paste0("MTLE", unname(kmeans_C2CP_MTLE$cluster)),
                                   paste0("PTCE", unname(kmeans_C2CP_PTCE$cluster))))
plotData <- cbind(plotData, sumFGSEA_MTLEPTCE[plotData$Pathways,])
plotData$extSc_MTPT <- FGSEA_extScRNA_MTPT[row.names(plotData), "NES"]
plotData$mArray <- FGSEA_mArray[row.names(plotData), "NES"]
plotData2 <- plotData$Cluster[order(plotData$Cluster)]
plotData <- melt(plotData)
plotData$Pathways <- paste0(plotData$Cluster, "_", plotData$Pathways)

p1 <- ggplot(plotData, aes(y = variable, fill = value, x = Pathways)) +
    geom_tile() +
    scale_fill_gradient2(low = "#440154", mid = "#ffffff", high = "#7ad151", 
                         midpoint = 0) +
    xlab("") +
    ylab("") +
    coord_flip() +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", axis.text.y = element_blank(),
          axis.ticks = element_blank())


plotData2 <- as.data.frame(plotData2)
plotData2$Number <- seq(1, nrow(plotData2))

p2 <- ggplot(plotData2, aes(x = Number, color = plotData2, fill = plotData2)) +
    geom_bar() +
    xlab("Cluster") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    scale_color_manual(values = c("#440154", "#3e4b89", "#25848e", "#38b977", "#bbdf24",
                                  "#000004", "#1e1149", "#53137d", "#872781", "#bc3978", "#eb5760", "#fd9166", "#fece91")) +
    scale_fill_manual(values = c("#440154", "#3e4b89", "#25848e", "#38b977", "#bbdf24",
                                 "#000004", "#1e1149", "#53137d", "#872781", "#bc3978", "#eb5760", "#fd9166", "#fece91")) +
    coord_flip() +
    theme_bw(base_size = 11) +
    theme(axis.text = element_blank(), legend.position = "none",
          axis.ticks = element_blank(), axis.title = element_blank())

p2 + p1 + plot_layout(widths = c(0.5, 10))

## Plotting example pathways
### Figure 5b, 5c, 5d
plotPW <- "BIOCARTA_ARENRF2_PATHWAY" ## Figure 5b
plotPW <- "REACTOME_REGULATION_OF_LIPID_METABOLISM_BY_PPARALPHA" ## Figure 5c
plotPW <- "KEGG_OXIDATIVE_PHOSPHORYLATION" ## Figure 5d

 
plotGenes <- pathwaysC2CP[[plotPW]]
plotData <- data.frame(MTPT_D1D1 = sumDE_LECE_MTPT[plotGenes, "MTPT_D1D1"],
                       MTPT_D5D5 = sumDE_LECE_MTPT[plotGenes, "MTPT_D5D5"],
                       LECE_D1D1 = sumDE_LECE_MTPT[plotGenes, "LECE_D1D1"],
                       LECE_D3D5 = sumDE_LECE_MTPT[plotGenes, "LECE_D3D5"],
                       Genes = plotGenes,
                       row.names = plotGenes)

plotData <- plotData[!is.na(rowSums(plotData[, 1:4])),]
geneOrder <- plotData$Genes[order(rowSums(plotData[, 1:4]), decreasing = T)]
plotData$extScRNA <- sumExtScRNAseq[plotData$Genes, "meanlogFC_ExtSc_MTPT"]
plotData$mArray <- mArrayDE[plotData$Genes, "logFC"]

plotData <- melt(plotData)
plotData$value <- sapply(plotData$value, \(x){ifelse(x>1, 1, x)})
plotData$value <- sapply(plotData$value, \(x){ifelse(x<(-1), -1, x)})

ggplot(plotData, aes(y = variable, x = Genes, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "#440154", mid = "#ffffff", high = "#7ad151", 
                         midpoint = 0, breaks = c(-6, -1, 0, 1, 6)) +
    xlab("") +
    ylab("") +
    theme_bw(base_size = 11) +
    #theme(axis.text.x = element_text(angle = 70, vjust = 0.5, size = 5))
    theme(axis.text.x = element_text(angle = 70, vjust = 0.5, size = 8))


### Extended Data Figure 5a, 5b, 5c: Histogram of all genes in example pathways measured in mArray data
plotmArray <- mArrayDE[intersect(row.names(mArrayDE), plotGenes), "logFC"]
max(abs(min(plotmArray)), max(plotmArray))

ggplot(mapping = aes(x = plotmArray)) +
    geom_histogram(color = "#440154", fill = "#440154", breaks = seq(-1, 1, 0.05)) +
    #geom_histogram(color = "#7ad151", fill = "#7ad151", breaks = seq(-1.6, 1.6, 0.05)) +
    geom_vline(xintercept = 0, color = "#a3a3a3", linetype = "longdash", linewidth = 1) +
    geom_vline(xintercept = mean(plotmArray), color = "#ae5bc2", linewidth = 1) +
    #geom_vline(xintercept = mean(plotmArray), color = "#2e780b", linewidth = 1) +
    xlab("log FC - mArray") +
    ylab("Counts") +
    theme_bw(base_size = 11) 


## Extended Data Table 6
extDataTab6 <- sumFGSEA_MTLEPTCE_addPval
extDataTab6$MTPT_specific <- sumFGSEA_MTPT[row.names(extDataTab6), "MTPT_specific"]
extDataTab6 <- extDataTab6[, c(1:4,11,13,6,7:10,12)]
extDataTab6$ExtmArray_MTPT_NES <- FGSEA_mArray[row.names(extDataTab6), "NES"]
extDataTab6$ExtmArray_MTPT_Padj <- FGSEA_mArray[row.names(extDataTab6), "padj"]
extDataTab6$ExtSc_MTPT_NES <- FGSEA_extScRNA_MTPT[row.names(extDataTab6), "NES"]
extDataTab6$ExtSc_MTPT_Padj <- FGSEA_extScRNA_MTPT[row.names(extDataTab6), "padj"]
extDataTab6 <- extDataTab6[order(extDataTab6$MTPT_specific, decreasing = T),]
extDataTab6 <- extDataTab6[order(extDataTab6$LECE_MTPT_specific, decreasing = T),]
