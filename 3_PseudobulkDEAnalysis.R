## Code from Luise Nagel

source("3b_Functions4DEAnalysis.R")
library(corrplot)
library(deming)
library(ggplot2)
library(ggpubr)
library(SummarizedExperiment)
library(patchwork)
library(reshape2)


# Loading data and filtering for DE genes analysis
allSE <- readRDS("../../Data/AllDonors_Filtered/SE_filteredAndPreprocessed.rds")

fSE <- allSE[, !(allSE$Tissue == "lMT" & allSE$Celltype == "LiverEpithelial")]
fSE <- fSE[, !(fSE$Donor == "Donor02" | fSE$Donor == "Donor04")]
fSE <- fSE[, !(fSE$Donor == "Donor03" & fSE$Celltype == "Metastasis")]

# Performing pseudobulk DE analysis 
## Primary tumor vs paired metastasis
### Estimating and centering pseudo-bulk gene expression 
pseudoB_PTMT <- estimateMeanExp(SE = fSE,
                                matName = "counts",
                                logT = F, 
                                ctName = "Celltype",
                                ct2compare = c("PrimaryTumor", "Metastasis"),
                                ct2filter = c("PrimaryTumor", "Metastasis"),
                                donorName = "Donor",
                                minFeaturePerc = 10,
                                normLib = 10000,
                                minCells = 49)

### Estimating donor-wise log2FC & adding classification
sumDE_MTPT <- data.frame(MTPT_D1D1 = pseudoB_PTMT$Donor01_Metastasis - pseudoB_PTMT$Donor01_PrimaryTumor,
                         MTPT_D5D5 = pseudoB_PTMT$Donor05_Metastasis - pseudoB_PTMT$Donor05_PrimaryTumor,
                         row.names = row.names(pseudoB_PTMT))
sumDE_MTPT$meanlogFC_MTPT <- rowMeans(sumDE_MTPT)
sumDE_MTPT$MTPT_Candidate <- apply(sumDE_MTPT, 1, \(x){
    if(as.numeric(x["MTPT_D1D1"])>0.1 & as.numeric(x["MTPT_D5D5"])>0.1){
        "MT"
    } else if(as.numeric(x["MTPT_D1D1"])<(-0.1) & as.numeric(x["MTPT_D5D5"])<(-0.1)){
        "PT"
    } else{
        FALSE
    }
})

### Plotting
### Figure 2b: Density plot PT and MT genes
ggplot(sumDE_MTPT, aes(x = meanlogFC_MTPT, color = MTPT_Candidate, fill = MTPT_Candidate)) +
    geom_density(alpha = 0.5, linewidth = 0.7) +
    scale_color_manual(name = "Candidate", values = c("#b3b3b3", "#7ad151", "#440154")) +
    scale_fill_manual(name = "Candidate", values = c("#b3b3b3", "#7ad151", "#440154")) +
    xlab("mean logFC MT/PT") +
    ylab("Density") +
    theme_bw(base_size = 11)

### Figure 3a: log2FC of PT & MT candidates in colon vs liver epithelium
plotData <- sumDE_LECE_MTPT
plotData$MTPT_Candidate_MTPTonly <- sumDE_MTPT[row.names(plotData), "MTPT_Candidate"]
plotData <- plotData[plotData$MTPT_Candidate_MTPTonly !=F,]

ggplot(plotData, aes(x = meanlogFC_LECE, color = MTPT_Candidate_MTPTonly, 
                     alpha = MTPT_Candidate_MTPTonly, fill = MTPT_Candidate_MTPTonly)) +
    geom_density(linewidth = 0.7) +
    scale_color_manual(name = "", values = c("#7ad151","#440154")) +
    scale_fill_manual(name = "", values = c("#7ad151","#440154")) +
    scale_alpha_manual(name = "", values = c(0.8, 0.5)) +
    geom_vline(xintercept = 0, linetype = "longdash", color = "#000000") +
    xlab("log2FC - liver vs. colon epithelial cells") +
    ylab("Density") +
    theme_bw(base_size = 11)


## Tissue adaptive genes: Primary tumor vs paired metastasis & colon epithelium vs liver epithelium 
### Estimating and centering pseudo-bulk gene expression 
pseudoB_PTCEMTLE <- estimateMeanExp(SE = fSE,
                                    matName = "counts", 
                                    logT = F, 
                                    ctName = "Celltype", 
                                    ct2compare = c("PrimaryTumor", "Metastasis", "ColonEpithelial", "LiverEpithelial"),
                                    ct2filter = c("PrimaryTumor", "Metastasis", "ColonEpithelial", "LiverEpithelial"),
                                    donorName = "Donor",
                                    minFeaturePerc = 10,
                                    normLib = 10000,
                                    minCells = 49)

### Estimating donor-wise log2FC & adding classification
sumDE_LECE_MTPT <- data.frame(MTPT_D1D1 = pseudoB_PTCEMTLE$Donor01_Metastasis - pseudoB_PTCEMTLE$Donor01_PrimaryTumor,
                              LECE_D1D1 = pseudoB_PTCEMTLE$Donor01_LiverEpithelial - pseudoB_PTCEMTLE$Donor01_ColonEpithelial,
                              MTPT_D5D5 = pseudoB_PTCEMTLE$Donor05_Metastasis - pseudoB_PTCEMTLE$Donor05_PrimaryTumor,
                              LECE_D3D5 = pseudoB_PTCEMTLE$Donor03_LiverEpithelial - pseudoB_PTCEMTLE$Donor05_ColonEpithelial,
                              row.names = row.names(pseudoB_PTCEMTLE))
sumDE_LECE_MTPT$meanlogFC_MTPT <- rowMeans(sumDE_LECE_MTPT[, c("MTPT_D1D1", "MTPT_D5D5")])
sumDE_LECE_MTPT$meanlogFC_LECE <- rowMeans(sumDE_LECE_MTPT[, c("LECE_D1D1", "LECE_D3D5")])
sumDE_LECE_MTPT$MTPT_Candidate <- apply(sumDE_LECE_MTPT, 1, \(x){
    if(as.numeric(x["MTPT_D1D1"])>0.1 & as.numeric(x["MTPT_D5D5"])>0.1){
        "MT"
    } else if(as.numeric(x["MTPT_D1D1"])<(-0.1) & as.numeric(x["MTPT_D5D5"])<(-0.1)){
        "PT"
    } else{
        FALSE
    }
})

sumDE_LECE_MTPT$LECE_Candidate <- apply(sumDE_LECE_MTPT, 1, \(x){
    if(as.numeric(x["LECE_D1D1"])>0.1 & as.numeric(x["LECE_D3D5"])>0.1){
        "LE"
    } else if(as.numeric(x["LECE_D1D1"])<(-0.1) & as.numeric(x["LECE_D3D5"])<(-0.1)){
        "CE"
    } else{
        FALSE
    }
})

sumDE_LECE_MTPT$LECE_MTPT_Candidate <- apply(sumDE_LECE_MTPT, 1, \(x){
    if(x["MTPT_Candidate"]=="PT" & x["LECE_Candidate"]=="CE"){
        "PTCE"
    } else if(x["MTPT_Candidate"]=="MT" & x["LECE_Candidate"]=="LE"){
        "MTLE"
    } else{
        FALSE
    }
})

sumDE_LECE_MTPT$Quadrant01  <- apply(sumDE_LECE_MTPT, 1, \(x){
    if(as.numeric(x["MTPT_D1D1"])>0.1 & as.numeric(x["MTPT_D5D5"])>0.1 & as.numeric(x["LECE_D1D1"])>0.1 & as.numeric(x["LECE_D3D5"])>0.1){
        "TR"
    } else if(as.numeric(x["MTPT_D1D1"])<(-0.1) & as.numeric(x["MTPT_D5D5"])<(-0.1) & as.numeric(x["LECE_D1D1"])>0.1 & as.numeric(x["LECE_D3D5"])>0.1){
        "BR"
    } else if(as.numeric(x["MTPT_D1D1"])<(-0.1) & as.numeric(x["MTPT_D5D5"])<(-0.1) & as.numeric(x["LECE_D1D1"])<(-0.1) & as.numeric(x["LECE_D3D5"])<(-0.1)){
        "BL"
    } else if(as.numeric(x["MTPT_D1D1"])>0.1 & as.numeric(x["MTPT_D5D5"])>0.1 & as.numeric(x["LECE_D1D1"])<(-0.1) & as.numeric(x["LECE_D3D5"])<(-0.1)){
        "TL"
    } else{
        FALSE
    }
})


### Figure 3c, Extended Data Figure 2b, 2c:
### log2FC in different tissue comparisons, highlighting MTLE and PTCE genes
plotData <- sumDE_LECE_MTPT
plotData$Quadrant01 <- factor(plotData$Quadrant01, levels = c("FALSE", "TL", "BR", "TR", "BL"))
plotData <- plotData[order(plotData$Quadrant01), ] 

demingCand <- deming(plotData$meanlogFC_MTPT[plotData$Quadrant01 != FALSE]~
                         plotData$meanlogFC_LECE[plotData$Quadrant01 != FALSE])

ggplot(plotData, aes(x = meanlogFC_LECE, y = meanlogFC_MTPT, color = Quadrant01, alpha = Quadrant01)) + # Figure 3c
#ggplot(plotData, aes(x = LECE_D1D1, y = MTPT_D1D1, color = Quadrant01, alpha = Quadrant01)) + # Extended Data Figure 2b
#ggplot(plotData, aes(x = LECE_D3D5, y = MTPT_D5D5, color = Quadrant01, alpha = Quadrant01)) + # Extended Data Figure 2c
    geom_point(size = 2) +
    scale_color_manual(values = c("#b3b3b3", "#000000", "#000000", "#7ad151", "#440154")) +
    scale_alpha_manual(values = c(0.3, 0.5, 0.5, 0.5, 0.5)) +
    xlim(-3, 3) +
    ylim(-3, 3) +
    geom_abline(intercept = demingCand$coefficients[1], 
                slope = demingCand$coefficients[2], linewidth = 1, color = "#595959") +
    xlab("log2FC - Epithelial cells") +
    ylab("log2FC - Tumor cells") +
    theme_bw(base_size = 11) +
    theme(aspect.ratio = 1)


## Comparing MT&PT with MTLE&PTCE candidates
plotData <- sumDE_LECE_MTPT
plotData$meanlogFC_MTPTLECE <- rowMeans(plotData[, c(5,6)])
plotData$MTPT_Candidate_inOnlyMTPT <- sumDE_MTPT[row.names(plotData), "MTPT_Candidate"]
plotData$meanlogFC_MTPT_inOnlyMTPT <- sumDE_MTPT[row.names(plotData), "meanlogFC_MTPT"]

plotData$Colors <- sapply(c(1:nrow(plotData)), \(i){
    if(plotData$LECE_MTPT_Candidate[i]!="FALSE"&plotData$MTPT_Candidate_inOnlyMTPT[i]!="FALSE"){
        "Both"
    }
    else if(plotData$MTPT_Candidate_inOnlyMTPT[i]!="FALSE"){
        "MTPT"
    }
    else if(plotData$LECE_MTPT_Candidate[i]!="FALSE"){
        "MTLEPTCE"
    }
    else{
        "FALSE"
    }
})

plotData$Colors_MTPT <- sapply(c(1:nrow(plotData)), \(i){
    if(plotData$LECE_MTPT_Candidate[i]=="PTCE"&plotData$MTPT_Candidate_inOnlyMTPT[i]=="PT"){
        "PTBoth"
    }
    else if(plotData$LECE_MTPT_Candidate[i]=="MTLE"&plotData$MTPT_Candidate_inOnlyMTPT[i]=="MT"){
        "MTBoth"
    }
    else if(plotData$MTPT_Candidate_inOnlyMTPT[i]!="FALSE"){
        plotData$MTPT_Candidate_inOnlyMTPT[i]
    }
    else{
        "FALSE"
    }
})

plotData$Colors_LECE_MTPT <- sapply(c(1:nrow(plotData)), \(i){
    if(plotData$LECE_MTPT_Candidate[i]=="PTCE"&plotData$MTPT_Candidate_inOnlyMTPT[i]=="PT"){
        "PTBoth"
    }
    else if(plotData$LECE_MTPT_Candidate[i]=="MTLE"&plotData$MTPT_Candidate_inOnlyMTPT[i]=="MT"){
        "MTBoth"
    }
    else if(plotData$LECE_MTPT_Candidate[i]!="FALSE"){
        plotData$LECE_MTPT_Candidate[i]
    }
    else{
        "FALSE"
    }
})

plotData$Colors_MTPT <- factor(plotData$Colors_MTPT, levels = c("FALSE", "PT", "MT", "MTBoth", "PTBoth"))
plotData$Colors_LECE_MTPT <- factor(plotData$Colors_LECE_MTPT, levels = c("FALSE", "PTCE", "MTLE", "MTBoth", "PTBoth"))
plotData$Colors <- factor(plotData$Colors, levels = c( "FALSE", "MTLEPTCE", "MTPT", "Both"))
plotData <- plotData[order(plotData$Colors),]


### Extended Data Figure 3a: log2FC of MTPT in MTvsPT and MTLEvsPTCE analysis
doScatter <- function(xAxis, yAxis, xlab = "log2FC MTPT - MTPT analysis",
                      ylab = "log2FC MTPT - MTLEPTCE analysis"){
    ggplot(mapping = aes(x = xAxis, y = yAxis)) +
        geom_point(size = 1.5, alpha = 0.1) +
        geom_vline(xintercept = 0, linewidth = 0.7, color = "#000000", linetype = "dashed") +
        geom_hline(yintercept = 0, linewidth = 0.7, color = "#000000", linetype = "dashed") +
        xlab(xlab) +
        ylab(ylab) +
        theme_bw(base_size = 11)
}

doScatter(xAxis = sumDE_MTPT[row.names(sumDE_LECE_MTPT), "MTPT_D1D1"],
          yAxis = sumDE_LECE_MTPT[, "MTPT_D1D1"])

doScatter(xAxis = sumDE_MTPT[row.names(sumDE_LECE_MTPT), "MTPT_D5D5"],
          yAxis = sumDE_LECE_MTPT[, "MTPT_D5D5"])

doScatter(xAxis = sumDE_MTPT[row.names(sumDE_LECE_MTPT), "meanlogFC_MTPT"],
          yAxis = sumDE_LECE_MTPT[, "meanlogFC_MTPT"])


## Scatter and density of log2FC MTPT & LECE, colored by candidates in different analysis 
doScatterDensityPlots <- function(xAxis, yAxis, color, 
                                  xlab = "log2FC - Epithelial cells",
                                  ylab = "log2FC - Tumor cells"){
    p1 <- ggplot(mapping = aes(x = xAxis, color = color, 
                               alpha = color, fill = color)) +    
        geom_density(linewidth = 0.5) +
        scale_color_manual(values = c("#b3b3b3", "#b684c2", "#7ad151", "#2d7a09", "#440154")) +
        scale_fill_manual(values = c("#b3b3b3", "#b684c2", "#7ad151", "#2d7a09", "#440154")) +
        scale_alpha_manual(values = c(0.1, 0.3, 0.3, 0.3, 0.3)) +
        xlim(-3, 3) +
        ylab("Density") +
        theme_bw(base_size = 11) +
        theme(legend.position = "none", axis.title.x = element_blank(), 
              axis.text.x = element_blank(), axis.ticks.x = element_blank())
    
    p2 <- ggplot(mapping = aes(x = yAxis, color = color, 
                               alpha = color, fill = color)) +    
        geom_density(linewidth = 0.5) +
        coord_flip() +
        scale_color_manual(values = c("#b3b3b3", "#b684c2", "#7ad151", "#2d7a09", "#440154")) +
        scale_fill_manual(values = c("#b3b3b3", "#b684c2", "#7ad151", "#2d7a09", "#440154")) +
        scale_alpha_manual(values = c(0.1, 0.3, 0.3, 0.3, 0.3)) +
        xlim(-3, 3) +
        ylab("Density") +
        theme_bw(base_size = 11)  +
        theme(legend.position = "none",  axis.title.y = element_blank(), 
              axis.text.y = element_blank(), axis.ticks.y = element_blank())
    
    p3 <- ggplot(mapping = aes(x = xAxis, y = yAxis, 
                               color = color, alpha = color)) +
        geom_point(size = 1.5) +
        scale_color_manual(values = c("#b3b3b3", "#b684c2", "#7ad151", "#2d7a09", "#440154")) +
        scale_alpha_manual(values = c(0.1, 0.3, 0.3, 0.3, 0.3)) +
        xlim(-3, 3) +
        ylim(-3, 3) +
        geom_vline(xintercept = 0, linewidth = 0.7, color = "#000000", linetype = "dashed") +
        geom_hline(yintercept = 0, linewidth = 0.7, color = "#000000", linetype = "dashed") +
        xlab(xlab) +
        ylab(ylab) +
        theme_bw(base_size = 11) +
        theme(legend.position = "none")
    
    plot <- p1 + plot_spacer() + p3 + p2 + 
        plot_layout(nrow = 2, ncol = 2, widths = c(5, 2), heights = c(2, 5))
    
    return(plot)
}

### Extended Data Figure 3b
doScatterDensityPlots(xAxis = plotData$meanlogFC_LECE,
                      yAxis = plotData$meanlogFC_MTPT, 
                      color = plotData$Colors_MTPT)

### Extended Data Figure 3c
doScatterDensityPlots(xAxis = plotData$meanlogFC_LECE,
                      yAxis = plotData$meanlogFC_MTPT, 
                      color = plotData$Colors_LECE_MTPT)




