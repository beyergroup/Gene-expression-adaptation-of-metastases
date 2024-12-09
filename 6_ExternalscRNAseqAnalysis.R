# Code from Luise Nagel

## Summarizing data.frames from earlier analysis required (see previous R scripts)

source("3b_Functions4DEAnalysis.R")
library(ggplot2)
library(MethComp)
library(SummarizedExperiment)

# Loading external data
matExtScRNAseq <- read.delim("../Data/extscRNAseq_PTMT/non_immune_counts.txt", row.names = 1)
colnames(matExtScRNAseq) <- gsub("[.]", "-", colnames(matExtScRNAseq))
annotExtScRNAseq <- read.delim("../Data/extscRNAseq_PTMT/non_immune_meta.txt", row.names = 1)
annotExtScRNAseq$Celltype <- ifelse(grepl("Tu", annotExtScRNAseq$cluster), "Tumor", "IDC")
annotExtScRNAseq$CelltypeTissue <- paste0(annotExtScRNAseq$organs, "_", annotExtScRNAseq$Celltype)

# Applying my DE approach
## Converting Seurat into summarized Experiment
SEextScRNAseq <- SummarizedExperiment(assays = list("counts" = matExtScRNAseq[, row.names(annotExtScRNAseq)]),
                                      colData = annotExtScRNAseq)
SEextScRNAseq <- SEextScRNAseq[!grepl("MT-", row.names(SEextScRNAseq)),]
SEextScRNAseq <- SEextScRNAseq[apply(assay(SEextScRNAseq, "counts"), 1, \(x) {sum(x[x!=0])>3}),]

rm(matExtScRNAseq); gc() ## clear memory space

## Getting mean exp. per sample
meanExpExtScRNAseq <- EstimateMeanExp(SE = SEextScRNAseq,
                                      matName = "counts", 
                                      logT = F, 
                                      ctName = "CelltypeTissue",
                                      ct2compare = c("CCT_Tumor", "LCT_Tumor"),
                                      ct2filter = c("CCT_Tumor", "LCT_Tumor"),
                                      donorName = "patients",
                                      minFeaturePerc = 5,
                                      normLib = FALSE,
                                      minCells = 120) ## removes s0115, not enough tumor cells
meanExpExtScRNAseq <-  as.data.frame(apply(meanExpExtScRNAseq, 2, function(x) x/sum(x)*10000))

## Estimating logFC
sumExtScRNAseq <- data.frame(meanlogFC_ExtSc_MTPT = rowMeans(meanExpExtScRNAseq[, grepl("LCT", colnames(meanExpExtScRNAseq))])-
                                 rowMeans(meanExpExtScRNAseq[, grepl("CCT", colnames(meanExpExtScRNAseq))]),
                             logFC_ExtSc_MTPT_s0107 = meanExpExtScRNAseq$s0107_LCT_Tumor-meanExpExtScRNAseq$s0107_CCT_Tumor,
                             logFC_ExtSc_MTPT_s0813 = meanExpExtScRNAseq$s0813_LCT_Tumor-meanExpExtScRNAseq$s0813_CCT_Tumor,
                             logFC_ExtSc_MTPT_s0920 = meanExpExtScRNAseq$s0920_LCT_Tumor-meanExpExtScRNAseq$s0920_CCT_Tumor,
                             logFC_ExtSc_MTPT_s1231 = meanExpExtScRNAseq$s1231_LCT_Tumor-meanExpExtScRNAseq$s1231_CCT_Tumor,
                             row.names = row.names(meanExpExtScRNAseq))

# Compare to our results
sumDE_LECE_MTPT$meanlogFC_MTPT_extSc <- sumExtScRNAseq[row.names(sumDE_LECE_MTPT), "meanlogFC_ExtSc_MTPT"]
sumDE_MTPT$meanlogFC_MTPT_extSc <- sumExtScRNAseq[row.names(sumDE_MTPT), "meanlogFC_ExtSc_MTPT"]

## MT & PT candidates
### Extended Data Figure 4a: log2FC external scRNAseq to our MT & PT candidates
plotData <- sumDE_MTPT[sumDE_MTPT$MTPT_Candidate!=F,]

ggplot(plotData, aes(y = meanlogFC_MTPT_extSc, x = MTPT_Candidate, 
                     fill = MTPT_Candidate, color = MTPT_Candidate)) +
    geom_boxplot(alpha = 0.5, linewidth = 0.5) +
    scale_color_manual(values = c("#7ad151", "#440154")) +
    scale_fill_manual(values = c("#7ad151", "#440154")) +
    xlab("Classification in our scRNA-seq data") +
    ylab("log2FC(MT/PT) - external scRNAseq") +
    theme_bw(base_size = 11) 


### Figure 4a: log2FC external scRNAseq to our MTLE & PTCE candidates
plotData <- sumDE_LECE_MTPT[sumDE_LECE_MTPT$LECE_MTPT_Candidate!=F,]

ggplot(plotData, aes(y = meanlogFC_MTPT_extSc, x = LECE_MTPT_Candidate, 
                     fill = LECE_MTPT_Candidate, color = LECE_MTPT_Candidate)) +
    geom_boxplot(alpha = 0.5, linewidth = 0.5) +
    scale_color_manual(values = c("#7ad151", "#440154")) +
    scale_fill_manual(values = c("#7ad151", "#440154")) +
    xlab("Classification in our scRNA-seq data") +
    ylab("log2FC(MT/PT) - external scRNAseq") +
    theme_bw(base_size = 11)


### Figure 4b: log2FC external scRNAseq vs log2FC our scRNAseq data (MTLE & PTCE candidates colored)
plotData <- sumDE_LECE_MTPT[order(sumDE_LECE_MTPT$LECE_MTPT_Candidate),]
ggplot(plotData, aes(x = meanlogFC_MTPT, y = meanlogFC_MTPT_extSc,  color = LECE_MTPT_Candidate,
                     alpha = LECE_MTPT_Candidate)) +
    geom_point(size = 1) +
    scale_color_manual(values = c("#a7a7a7", "#7ad151", "#440154")) + 
    scale_alpha_manual(values = c(0.3, 1, 1)) +
    xlab("log2FC(MT/PT) - our scRNAseq") +
    ylab("log2FC(MT/PT) - external scRNAseq") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")
dev.off()


### Extended Data Table
extDataTab3 <- sumDE_MTPT[, c(1:6,10)]
extDataTab3 <- extDataTab3[order(abs(extDataTab3$meanlogFC_MTPT), decreasing = T),]
extDataTab3 <- extDataTab3[order(extDataTab3$MTPT_Candidate, decreasing = T),]

extDataTab4 <- sumDE_LECE_MTPT[, c(1:9,11,12,16)]
extDataTab4 <- extDataTab4[order(abs(extDataTab4$meanlogFC_MTPT), decreasing = T),]
extDataTab4 <- extDataTab4[order(extDataTab4$LECE_MTPT_Candidate, decreasing = T),]
