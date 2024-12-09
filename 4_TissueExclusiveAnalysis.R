# Retrieving tissue exclusive genes
## Filtering data, removing some cell types, removing groups of cells if not detected often enough in specific donor (min 49 cells)
f2SE <- allSE
f2SE <- f2SE[, !(f2SE$Tissue == "lMT" & f2SE$Celltype == "LiverEpithelial")]
f2SE  <- f2SE[, f2SE$Celltype == "PrimaryTumor"|f2SE$Celltype == "Metastasis"|f2SE$Celltype == "ColonEpithelial"|f2SE$Celltype == "LiverEpithelial"]
f2SE  <- f2SE[, !f2SE$Donor == "Donor02"]
f2SE  <- f2SE[, !(f2SE$Celltype == "LiverEpithelial"&(f2SE$Donor == "Donor04"|f2SE$Donor == "Donor05"))]

## Percent of cells genes are expressed in split by cell type and donor
percGen_perCTDonor <- do.call(cbind, lapply(list("ColonEpithelial", "LiverEpithelial", "PrimaryTumor", "Metastasis"), \(i){
    sapply(c("Donor01", "Donor02", "Donor03", "Donor04", "Donor05"), \(j){
        MyAssay <- as.matrix(assay(f2SE)[, f2SE$Celltype == i & f2SE$Donor == j])
        apply(MyAssay, 1, \(x) sum(x != 0)/length(x)*100)
    })
}))
colnames(percGen_perCTDonor) <- paste0(rep(c("CE", "LE", "PT", "MT"), each = 5), "_", colnames(percGen_perCTDonor))
percGen_perCTDonor <- as.data.frame(percGen_perCTDonor[, !is.nan(colMeans(percGen_perCTDonor))])
percGen_perCTDonor <- as.data.frame(percGen_perCTDonor[rowSums(percGen_perCTDonor)!=0,])


## Keeping only tissue groups with matching donor
#percGen_perCTDonor <- percGen_perCTDonor[, c("CE_Donor01", "CE_Donor05", "LE_Donor01", "LE_Donor03",
#                                             "PT_Donor01", "PT_Donor05", "MT_Donor01", "MT_Donor05")]
## Estimating mean number of times a genes was expressed per tissue over donors
sumTE_LECE_MTPT <- percGen_perCTDonor
sumTE_LECE_MTPT <- as.data.frame(do.call(cbind, lapply(split(as.data.frame(t(sumTE_LECE_MTPT)), 
                                                             gsub("_.*", "", colnames(sumTE_LECE_MTPT))), \(x){
                                                                 rowMeans(t(x))
                                                             })))
sumTE_LECE_MTPT <- sumTE_LECE_MTPT[, c("CE", "LE", "PT", "MT")]

## Defining exclusively expressed genes
sumTE_LECE_MTPT$MTPT_Candidate <- apply(sumTE_LECE_MTPT[, 1:4], 1, \(x){
    if(x["PT"]>20 & x["MT"]<10){
        "PT"
    }
    else if(x["PT"]<10 & x["MT"]>20){
        "MT"
    }
    else{
        FALSE
    }
})

sumTE_LECE_MTPT$LECE_MTPT_Candidate <- apply(sumTE_LECE_MTPT[, 1:4], 1, \(x){
    if(x["PT"]>20 & x["CE"]>20 & x["MT"]<10 & x["LE"]<10){
        "PTCE"
    }
    else if(x["PT"]<10 & x["CE"]<10 & x["MT"]>20 & x["LE"]>20){
        "MTLE"
    }
    else{
        FALSE
    }
})