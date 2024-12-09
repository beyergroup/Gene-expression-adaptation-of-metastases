# Code from Luise Nagel

## Summarizing data.frames from earlier analysis required (see previous R scripts)


library(ggplot2)
library(ggpubr)
library(limma)
library(FactoMineR)
library(patchwork)

# Loading integrated & normalized array data
mArray <- readRDS("../Data/mArrayData/Data150set.rds")

## Adjusting gene names
nameProbeGene <- read.csv("../Data/mArrayData/ToLoadR_HTA-2_0.r3.na36.hg19.a1.transcript.csv", row.names = 1)
nameProbeGene$GeneName <- sapply(nameProbeGene$gene_assignment, function(x){
    strsplit(x, " // ")[[1]][2]
})
row.names(mArray) <- nameProbeGene[row.names(mArray), "GeneName"]
mArray <- mArray[!is.na(row.names(mArray)),]

## Adding study information
mArrayAnnot <- data.frame(row.names = colnames(mArray),
                          Study = substr(colnames(mArray), 1, 6))

mArrayAnnot$TumorStudy <- sapply(mArrayAnnot$Study, function(x){
    if(x == "GSM482"){
        "Liver MT - Sveen 2020"
    }
    else if(x == "GSM253"){
        "Colon PT - Sveen 2017"
    }
    else if(x == "GSM210"|x == "GSM233"){
        "Colon PT - Sveen 2016"
    }
})

mArrayAnnot$Origin <- sapply(mArrayAnnot$Study, function(x){
    if(x == "GSM482"){
        "lMT"
    }
    else{
        "cPT"
    }
})

## PCA of microarray data
mArrayPCA <- PCA(t(mArray), graph = F)

## Extended Data Figure 4b: PCA of microarray data
ggplot(mapping = aes(x = mArrayPCA$ind$coord[, 1], y = mArrayPCA$ind$coord[, 2], 
                     color = mArrayAnnot$TumorStudy, alpha = mArrayAnnot$TumorStudy)) +
    geom_point(size = 2) +
    scale_color_viridis_d(name = "Study", end = 0.9) +
    scale_alpha_manual(values = c(0.8, 0.6, 0.6)) +
    xlab("PC 1") +
    ylab("PC 2") +
    theme_bw(base_size = 11) +
    theme(aspect.ratio = 1)


## DE analysis mArray
do_DEmArray <- function(mArray, origin, lvl1 = "cPT", lvl2 = "lMT"){
    y <- factor(origin, levels = c("cPT", "lMT"))   
    mArray_filter <- data.frame(rowNumber = 1:nrow(mArray),
                                name = row.names(mArray),
                                mean = rowMeans(mArray))
    
    mArray_filter <- do.call(rbind, lapply(split(mArray_filter, mArray_filter$name), function(x){
        if(nrow(x) > 1){
            x$keep <- x$mean == max(x$mean)
        }
        else(
            x$keep <- TRUE
        )
        return(x)
    }))
    mArray_filter <- mArray_filter[order(mArray_filter$rowNumber),]
    mArray_filter <- mArray[mArray_filter$keep,]
    
    mArray_mm <- model.matrix(~y, data = as.data.frame(mArray_filter))
    mArray_fit <- lmFit(mArray_filter, mArray_mm)
    mArray_tmp <- eBayes(mArray_fit)
    mArray_topTable <- topTable(mArray_tmp, sort.by = "P", n = Inf)
    
    mArrayDE <- mArray_topTable
    mArrayDE <- mArrayDE[sapply(mArrayDE$ID, function(x){
        ! x %in% (mArrayDE$ID[duplicated(mArrayDE$ID)])
    }), ]
    row.names(mArrayDE) <- mArrayDE$ID
    return(mArrayDE)
}
mArrayDE <- do_DEmArray(mArray = mArray, origin = mArrayAnnot$Origin)

### Extended Data Figure 3c: Volcano plot of DE analysis of microarray data
p1 <- ggplot(mArrayDE, aes(y = -log10(adj.P.Val), x = logFC, color = adj.P.Val<0.05)) +
    geom_point(size = 1.5, alpha = 0.3) +
    scale_color_manual(values = c("#000000", "#de7600")) +
    scale_fill_manual(values = c("#000000", "#de7600")) +
    xlab("logFC - mArray") +
    ylim(0, 100) +
    ylab("-log10(adj. P value) - mArray") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")

p2 <- ggplot(mArrayDE, mapping = aes(x = -log10(adj.P.Val), color = adj.P.Val<0.05, fill = adj.P.Val<0.05)) +
    geom_histogram(alpha = 0.5, linewidth = 0.7, breaks = seq(0, 100, 2), position = "identity") +
    coord_flip() +
    scale_color_manual(values = c("#000000", "#de7600")) +
    scale_fill_manual(values = c("#000000", "#de7600")) +
    xlab("") +
    ylab("Density") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")

p1 + p2 + plot_layout(widths = c(2, 1))


### Extended Data Table 7
mArrayUsed <- data.frame(SampleID = colnames(mArray),
                         Study = mArrayAnnot$Study,
                         Tissue = mArrayAnnot$Origin)


# Comparing to DE analysis
## Primary tumor vs paired metastasis
### Adding mArray data to previous summarizing data.frames
sumDE_MTPT$mArray_logFC <- mArrayDE[row.names(sumDE_MTPT), "logFC"]
sumDE_MTPT$mArray_AjdPval <- mArrayDE[row.names(sumDE_MTPT), "adj.P.Val"]
sumDE_MTPT$mArray_AjdPval05 <- sumDE_MTPT$mArray_AjdPval<0.05
sumDE_MTPT$mArrayAggreeDirection <- sapply(seq(1, nrow(sumDE_MTPT)), \(i){
    if(sumDE_MTPT$MTPT_Candidate[i] == "PT"){
        sumDE_MTPT$mArray_logFC[i]<0
    }
    else if(sumDE_MTPT$MTPT_Candidate[i] == "MT"){
        sumDE_MTPT$mArray_logFC[i]>0
    }
    else{
        NA
    }
})

sumDE_MTPT$mArrayAggreeDirectionAndAjdPval <- sapply(seq(1, nrow(sumDE_MTPT)), \(i){
    if(is.na(sumDE_MTPT$mArrayAggreeDirection[i])){
        NA
    }
    else if(sumDE_MTPT$mArrayAggreeDirection[i]){
        sumDE_MTPT$mArray_AjdPval05[i]
    }
    else if(!sumDE_MTPT$mArrayAggreeDirection[i]){
        FALSE
    }
})

### Extended Data Figure 4d: log2 FC of microarray data in our PT & MT classified peptides
plotData <- sumDE_MTPT
plotData$MTPT_CandidateSig <- sapply(1:nrow(plotData), \(i){
    ifelse(plotData$mArray_AjdPval05[i], plotData$MTPT_Candidate[i], "not significant")
})
plotData$MTPT_Candidate <- factor(plotData$MTPT_Candidate, levels = c("FALSE", "PT", "MT"))
plotData <- plotData[order(plotData$MTPT_Candidate), ]
plotData$MTPT_CandidateSig <- factor(plotData$MTPT_CandidateSig, levels = c("not significant", "FALSE", "PT", "MT"))

ggplot(plotData, aes(y = mArray_logFC, x = MTPT_Candidate, 
                     fill = MTPT_CandidateSig, color = MTPT_CandidateSig)) +
    geom_violin(alpha = 0.5, linewidth = 0.5) +
    scale_color_manual(values = c("#000000", "#b3b3b3", "#440154", "#7ad151")) +
    scale_fill_manual(values = c("#000000", "#b3b3b3", "#440154", "#7ad151")) +
    scale_x_discrete(labels = c("Not\nspecific", "PT", "MT")) +
    xlab("Classification in scRNAseq") +
    ylab("log2FC(PT/MT) - mArray") +
    theme_bw(base_size = 11) 



## Tissue adaptive genes: Primary tumor vs paired metastasis & colon epithelium vs liver epithelium 
### Adding mArray data to previous summarizing data.frames
sumDE_LECE_MTPT$mArray_logFC <- mArrayDE[row.names(sumDE_LECE_MTPT), "logFC"]
sumDE_LECE_MTPT$mArray_AjdPval <- mArrayDE[row.names(sumDE_LECE_MTPT), "adj.P.Val"]
sumDE_LECE_MTPT$mArray_AjdPval05 <- sumDE_LECE_MTPT$mArray_AjdPval<0.05
sumDE_LECE_MTPT$mArrayAggreeDirection <- sapply(seq(1, nrow(sumDE_LECE_MTPT)), \(i){
    if(sumDE_LECE_MTPT$MTPT_Candidate[i] == "PT"){
        sumDE_LECE_MTPT$mArray_logFC[i]<0
    }
    else if(sumDE_LECE_MTPT$MTPT_Candidate[i] == "MT"){
        sumDE_LECE_MTPT$mArray_logFC[i]>0
    }
    else{
        NA
    }
})
sumDE_LECE_MTPT$mArrayAggreeDirectionAndAjdPval <- sapply(seq(1, nrow(sumDE_LECE_MTPT)), \(i){
    if(is.na(sumDE_LECE_MTPT$mArrayAggreeDirection[i])){
        NA
    }
    else if(sumDE_LECE_MTPT$mArrayAggreeDirection[i]){
        sumDE_LECE_MTPT$mArray_AjdPval05[i]
    }
    else if(!sumDE_LECE_MTPT$mArrayAggreeDirection[i]){
        FALSE
    }
})


### Figure 4d: log2 FC of microarray data in our PT and MT DE classified peptides
plotData <- sumDE_LECE_MTPT
plotData$LECE_MTPT_CandidateSig <- sapply(1:nrow(plotData), \(i){
    ifelse(plotData$mArray_AjdPval05[i], plotData$LECE_MTPT_Candidate[i], "not significant")
})
plotData$LECE_MTPT_Candidate <- factor(plotData$LECE_MTPT_Candidate, levels = c("FALSE", "PTCE", "MTLE"))
plotData$LECE_MTPT_CandidateSig <- factor(plotData$LECE_MTPT_CandidateSig, levels = c("not significant", "FALSE", "PTCE", "MTLE"))

ggplot(plotData, aes(y = mArray_logFC, x = LECE_MTPT_Candidate,
                     fill = LECE_MTPT_CandidateSig, color = LECE_MTPT_CandidateSig)) +
    geom_violin(alpha = 0.5, linewidth = 0.5) +
    scale_color_manual(values = c("#000000", "#b3b3b3", "#440154", "#7ad151")) +
    scale_fill_manual(values = c("#000000", "#b3b3b3", "#440154", "#7ad151")) +
    scale_x_discrete(labels = c("Not\nspecific", "PTCE\ncandidates", "MTLE \ncandidates")) +
    xlab("Classification in scRNAseq") +
    ylab("log2FC(PT/MT) - mArray") +
    theme_bw(base_size = 11)



# Comparing to tissue exclusive analysis
## Add mArray data
sumTE_LECE_MTPT$mArray_logFC <- mArrayDE[row.names(sumTE_LECE_MTPT), "logFC"]
sumTE_LECE_MTPT$mArray_AjdPval <- mArrayDE[row.names(sumTE_LECE_MTPT), "adj.P.Val"]

### Extended Data Figure 4e: log2 FC of microarray data in our PT and MT tissue exclusive classified peptides
plotData <- sumTE_LECE_MTPT[sumTE_LECE_MTPT$MTPT_Candidate != F, ]
ggplot(plotData, aes(x = mArray_logFC, color = MTPT_Candidate, 
                     fill = MTPT_Candidate, alpha = MTPT_Candidate)) +
    geom_histogram(breaks = seq(-1.5, 2, 0.1), position = "identity") +
    geom_vline(xintercept = 0, linewidth = 0.7, color = "black") +
    scale_color_manual(name = "", values = c("#7ad151", "#440154")) +
    scale_fill_manual(name = "", values = c("#7ad151", "#440154")) +
    scale_alpha_manual(values =  c(0.75, 0.4)) +
    xlab("log2FC(MT/PT) - mArray") +
    ylab("Count") +
    theme_bw(base_size = 11)


### Figure 4d: log2 FC of microarray data in our PTCE and MTLE tissue exclusive classified peptides
plotData <- sumTE_LECE_MTPT[sumTE_LECE_MTPT$LECE_MTPT_Candidate != F, ]
ggplot(plotData, aes(x = mArray_logFC, color = LECE_MTPT_Candidate, 
                     fill = LECE_MTPT_Candidate, alpha = LECE_MTPT_Candidate)) +
    geom_histogram(breaks = seq(-1.5, 2, 0.1), position = "identity") +
    geom_vline(xintercept = 0, linewidth = 0.7, color = "black") +
    scale_color_manual(name = "", values = c("#7ad151", "#440154")) +
    scale_fill_manual(name = "", values = c("#7ad151", "#440154")) +
    scale_alpha_manual(values =  c(0.75, 0.4)) +
    xlab("log2FC(MT/PT) - mArray") +
    ylab("Count") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")


## Extended Data Table
extDataTab5 <- sumTE_LECE_MTPT[sumTE_LECE_MTPT$MTPT_Candidate != FALSE,]
extDataTab5 <- extDataTab5[order(extDataTab5$MTPT_Candidate, decreasing = T),]
extDataTab5 <- extDataTab5[order(extDataTab5$LECE_MTPT_Candidate, decreasing = T),]
