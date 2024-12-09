# Getting pseudo-bulk expression values 
estimateMeanExp <- function(SE, matName, logT = F, ctName, ct2compare = NULL, ct2filter = NULL, donorName = NULL, minFeature = NULL, minFeaturePerc = NULL, minCells = 1, normLib = 10000){
    if(logT == F){
        SE <- logData(SE, matName)
    }
    if(is.null(ct2compare)){
        ct2compare <- names(table(SE@colData[, ctName]))
        message(ct2compare[1], " and ", ct2compare[2], " will be compared.")
    }
    if(is.null(ct2filter)){
        ct2filter <- ct2compare
    }
    SE <- filterData(SE, matName, ctName, ct2compare, ct2filter, donorName, minFeature, minFeaturePerc, minCells)
    meanExp <- getMeanExp(SE, ctName, ct2compare, donorName, minCells)
    
    if(!is.null(normLib)){
        meanExp <- as.data.frame(apply(meanExp, 2, function(x)
            {x/sum(x)*normLib}
            ))
    }
    return(meanExp)
}


## Function to log transform data
logData <- function(SE, matName){
    message("Log-transforming data.")
    assay(SE, matName) <- log2(as.matrix(assay(SE, matName))+1)
    return(SE)
}

## Filtering cells and genes based on choosen cell types and minimum amount of times a genes is expressed
filterData <- function(SE, matName, ctName, ct2compare, ct2filter, donorName, minFeature, minFeaturePerc, minCells){
    if(!is.null(minFeature)){
        message("Filtering genes based on specific amount of detections.")
        if(is.null(donorName)){
            genesKeep <- Reduce(intersect, lapply(as.list(ct2filter), function(i){
                row.names(SE)[rowSums(assay(SE, matName)[, SE@colData[, ctName] == i, drop = F]>0)>=minFeature]
            }))
        }
        else{
            donors2compare <- unique(SE@colData[, donorName])
            message("Accounting for donor information while filtering genes.")
            genesKeep <- lapply(as.list(donors2compare), function(i){
                list <- lapply(as.list(ct2filter), function(j){
                    MyAssay <- assay(SE, matName)[, SE@colData[, donorName] == i & (SE@colData[, ctName] == j), drop = F]
                    if(ncol(MyAssay)<minCells){
                        return(NULL)
                    }
                    else{
                        
                        row.names(SE)[rowSums(MyAssay>0)>=minFeature]
                    }
                })
                list <- Reduce(intersect, list[unlist(lapply(list, length))!=0])
            })
            genesKeep <- Reduce(intersect, genesKeep[unlist(lapply(genesKeep, length))!=0])
        }
    }
    else if(!is.null(minFeaturePerc)){
        message("Filtering genes based on percent of cells in which they were detected.")
        if(is.null(donorName)){
            genesKeep <- Reduce(intersect, lapply(as.list(ct2filter), function(i){
                MyAssay <- row.names(SE)[rowSums(MyAssay>0) >= minFeature+Perc/100*ncol(MyAssay)]
            }))
        }
        else{
            donors2compare <- unique(SE@colData[, donorName])
            message("Accounting for donor information while filtering genes.")
            genesKeep <- lapply(as.list(donors2compare), function(i){
                list <- lapply(as.list(ct2filter), function(j){
                    MyAssay <- assay(SE, matName)[, SE@colData[, donorName] == i & (SE@colData[, ctName] == j), drop = F]
                    MyAssay <- assay(SE, matName)[, SE@colData[, donorName] == i & (SE@colData[, ctName] == j), drop = F]
                    if(ncol(MyAssay)<minCells){
                        return(NULL)
                    }
                    else{
                        row.names(SE)[rowSums(MyAssay>0) >= minFeaturePerc/100*ncol(MyAssay)]
                    }
                })
                list <- Reduce(intersect, list[unlist(lapply(list, length))!=0])
            })
            genesKeep <- Reduce(intersect, genesKeep[unlist(lapply(genesKeep, length))!=0])
        }
    }
    else{
        genesKeep <- row.names(SE)
    }
    message(length(genesKeep), " genes remaining in analysis.")
    SE <- SummarizedExperiment(assays = list(log2 = assay(SE, matName)[genesKeep, SE@colData[, ctName] %in% ct2compare]),
                               colData = as.data.frame(SE@colData)[SE@colData[, ctName] %in% ct2compare, , drop = F])
    return(SE)
}

## Estimating mean expression values not taking zeros into account
getMeanExp <- function(SE, ctName, ct2compare, donorName, minCells){
    if(is.null(donorName)){
        message("Assuming only one Donor/not differentiating between donors.")
        SEmean <- data.frame(CT1 = apply(as.data.frame(assay(SE, "log2")[, SE@colData[, ctName] == ct2compare[1]]), 1,  function(x) mean(x[x!=0])),
                             CT2 = apply(as.data.frame(assay(SE, "log2")[, SE@colData[, ctName] == ct2compare[2]]), 1,  function(x) mean(x[x!=0])))
        SEmean <- data.frame(CT1 = apply(assay(SE, "log2")[, SE@colData[, ctName] == ct2compare[1]], 1,  function(x) mean(x[x!=0])),
                             CT2 = apply(assay(SE, "log2")[, SE@colData[, ctName] == ct2compare[2]], 1,  function(x) mean(x[x!=0])))
        colnames(SEmean) <- ct2compare
    }
    else{
        donors2compare <- names(table(SE@colData[, donorName]))
        cat(cat(donors2compare, sep = " & "), " will be compared.\n")
        message("Estimating mean per donor befor adjusting and summarizing expression data.")
        SEdonor <- split(as.data.frame(t(assay(SE, "log2"))), paste0(SE@colData[, donorName], "_", SE@colData[, ctName]))
        SEdonor <- SEdonor[lapply(SEdonor, nrow)>= minCells]
        SEmean <- do.call(cbind, lapply(SEdonor, function(x){
            apply(x, 2, function(y) mean(y[y!=0]))
        }))
    }
    return(SEmean)
}