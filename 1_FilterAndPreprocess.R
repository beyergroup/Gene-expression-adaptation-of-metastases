## Code from Marten Wenzel


##### PROCESS RAW DATA #####


  # CRCLM01; batch 1
  # reading in raw data
  sc_raw_data_crclm01_cst <- Seurat::Read10X(data.dir = paste0(raw_data_directory,"/CRCLM01/SID115169/outs/filtered_feature_bc_matrix/")) # CRCLM01 C.sig. Tumor
  sc_raw_data_crclm01_csn <- Read10X(data.dir = paste0(raw_data_directory,"/CRCLM01/SID115170/outs/filtered_feature_bc_matrix/")) # CRCLM01 C.sig. normal
  sc_raw_data_crclm01_can <- Read10X(data.dir = paste0(raw_data_directory,"/CRCLM01//SID115172/outs/filtered_feature_bc_matrix/")) # CRCLM01 C.asc. normal
  sc_raw_data_crclm01_lit <- Read10X(data.dir = paste0(raw_data_directory,"/CRCLM01//SID115173/outs/filtered_feature_bc_matrix/")) # CRCLM01 Liver Metastasis
  sc_raw_data_crclm01_lin <- Read10X(data.dir = paste0(raw_data_directory,"/CRCLM01//SID115174/outs/filtered_feature_bc_matrix/")) # CRCLM01 Liver Normal
  
  # sample annotation
  ident_string_crclm01_cst <- "crclm01_cst"   # CRCLM01 C.sig. Tumor
  ident_string_crclm01_csn <- "crclm01_csn"  # CRCLM01 C.sig. normal
  ident_string_crclm01_can <- "crclm01_can"  # CRCLM01 C.asc. normal
  ident_string_crclm01_lit <- "crclm01_lit" # CRCLM01 Liver Metastasis
  ident_string_crclm01_lin <- "crclm01_lin" # CRCLM01 Liver Normal
  
  # Creating Seurat objects using all genes appearing in at least 3 cells and using all cells with at least 1000 different expressed genes
  seu_crclm01_cst <- CreateSeuratObject(counts = sc_raw_data_crclm01_cst, min.cells = 3, min.features = 1000, project = ident_string_crclm01_cst)
  seu_crclm01_csn <- CreateSeuratObject(counts = sc_raw_data_crclm01_csn, min.cells = 3, min.features = 1000, project = ident_string_crclm01_csn)
  seu_crclm01_can <- CreateSeuratObject(counts = sc_raw_data_crclm01_can, min.cells = 3, min.features = 1000, project = ident_string_crclm01_can)
  seu_crclm01_lit <- CreateSeuratObject(counts = sc_raw_data_crclm01_lit, min.cells = 3, min.features = 1000, project = ident_string_crclm01_lit)
  seu_crclm01_lin <- CreateSeuratObject(counts = sc_raw_data_crclm01_lin, min.cells = 3, min.features = 1000, project = ident_string_crclm01_lin)
  
  # CRCLM01; batch 2
  # reading in raw data
  sc_raw_data_crclm02_cxt <- Read10X(data.dir = paste0(raw_data_directory,"/CRCLM02_CRCLM03//SID127153/outs/filtered_feature_bc_matrix/")) # CRCLM02 Colon Tumor
  sc_raw_data_crclm02_cxn <- Read10X(data.dir = paste0(raw_data_directory,"/CRCLM02_CRCLM03//SID127154/outs/filtered_feature_bc_matrix/")) # CRCLM02 Colon Normal
  sc_raw_data_crclm03_lit <- Read10X(data.dir = paste0(raw_data_directory,"/CRCLM02_CRCLM03//SID127157/outs/filtered_feature_bc_matrix/")) # CRCLM03 Liver Metastasis
  sc_raw_data_crclm03_lin <- Read10X(data.dir = paste0(raw_data_directory,"/CRCLM02_CRCLM03//SID127158/outs/filtered_feature_bc_matrix/")) # CRCLM03 Liver Normal
  
  # sample annotation
  ident_string_crclm02_cxt <- "crclm02_cxt" # CRCLM02 Colon Tumor
  ident_string_crclm02_cxn <- "crclm02_cxn" # CRCLM02 Colon Normal
  ident_string_crclm03_lit <- "crclm03_lit" # CRCLM03 Liver Metastasis
  ident_string_crclm03_lin <- "crclm03_lin" # CRCLM03 Liver Normal
  
  # Creating Seurat objects using all genes appearing in at least 3 cells and using all cells with at least 1000 different expressed genes
  seu_crclm02_cxt <- CreateSeuratObject(counts = sc_raw_data_crclm02_cxt, min.cells = 3, min.features = 1000, project = ident_string_crclm02_cxt)
  seu_crclm02_cxn <- CreateSeuratObject(counts = sc_raw_data_crclm02_cxn, min.cells = 3, min.features = 1000, project = ident_string_crclm02_cxn)
  seu_crclm03_lit <- CreateSeuratObject(counts = sc_raw_data_crclm03_lit, min.cells = 3, min.features = 1000, project = ident_string_crclm03_lit)
  seu_crclm03_lin <- CreateSeuratObject(counts = sc_raw_data_crclm03_lin, min.cells = 3, min.features = 1000, project = ident_string_crclm03_lin)
  
  # CRCLM01; batch 3
  # reading in raw data
  sc_raw_data_crclm04_cxt <- Read10X(data.dir = paste0(raw_data_directory,"/CRCLM04_CRCLM05/SID139399/outs/filtered_feature_bc_matrix/"))  # CRCLM04 Colon Tumor
  sc_raw_data_crclm04_cxn <- Read10X(data.dir = paste0(raw_data_directory,"/CRCLM04_CRCLM05/SID139401/outs/filtered_feature_bc_matrix/"))  # CRCLM04 Colon Normal
  sc_raw_data_crclm04_lin <- Read10X(data.dir = paste0(raw_data_directory,"/CRCLM04_CRCLM05/SID139403/outs/filtered_feature_bc_matrix/"))  # CRCLM04 Liver Normal
  sc_raw_data_crclm05_cxt <- Read10X(data.dir = paste0(raw_data_directory,"/CRCLM04_CRCLM05/SID139405/outs/filtered_feature_bc_matrix/"))  # CRCLM05 Colon Tumor
  sc_raw_data_crclm05_cxn <- Read10X(data.dir = paste0(raw_data_directory,"/CRCLM04_CRCLM05/SID139407/outs/filtered_feature_bc_matrix/"))  # CRCLM05 Colon Normal 
  sc_raw_data_crclm05_lit <- Read10X(data.dir = paste0(raw_data_directory,"/CRCLM04_CRCLM05/SID139409/outs/filtered_feature_bc_matrix/"))  # CRCLM05 Liver Metastasis
  sc_raw_data_crclm05_lin <- Read10X(data.dir = paste0(raw_data_directory,"/CRCLM04_CRCLM05/SID139411/outs/filtered_feature_bc_matrix/"))  # CRCLM05 Liver Normal
  
  # sample annotation
  ident_string_crclm04_cxt <- "crclm04_cxt"  # CRCLM04 Colon Tumor
  ident_string_crclm04_cxn <- "crclm04_cxn"  # CRCLM04 Colon Normal
  ident_string_crclm04_lin <- "crclm04_lin"  # CRCLM04 Liver Normal
  ident_string_crclm05_cxt <- "crclm05_cxt"  # CRCLM05 Colon Tumor
  ident_string_crclm05_cxn <- "crclm05_cxn"  # CRCLM05 Colon Normal #### HIER WAR DER FEHLER
  ident_string_crclm05_lit <- "crclm05_lit"  # CRCLM05 Liver Metastasis
  ident_string_crclm05_lin <- "crclm05_lin"  # CRCLM05 Liver Normal
  
  # Creating Seurat objects using all genes appearing in at least 3 cells and using all cells with at least 1000 different expressed genes
  seu_crclm04_cxt <- CreateSeuratObject(counts = sc_raw_data_crclm04_cxt, min.cells = 3, min.features = 1000, project = ident_string_crclm04_cxt)
  seu_crclm04_cxn <- CreateSeuratObject(counts = sc_raw_data_crclm04_cxn, min.cells = 3, min.features = 1000, project = ident_string_crclm04_cxn)
  seu_crclm04_lin <- CreateSeuratObject(counts = sc_raw_data_crclm04_lin, min.cells = 3, min.features = 1000, project = ident_string_crclm04_lin)
  seu_crclm05_cxt <- CreateSeuratObject(counts = sc_raw_data_crclm05_cxt, min.cells = 3, min.features = 1000, project = ident_string_crclm05_cxt)
  seu_crclm05_cxn <- CreateSeuratObject(counts = sc_raw_data_crclm05_cxn, min.cells = 3, min.features = 1000, project = ident_string_crclm05_cxn)
  seu_crclm05_lit <- CreateSeuratObject(counts = sc_raw_data_crclm05_lit, min.cells = 3, min.features = 1000, project = ident_string_crclm05_lit)
  seu_crclm05_lin <- CreateSeuratObject(counts = sc_raw_data_crclm05_lin, min.cells = 3, min.features = 1000, project = ident_string_crclm05_lin)
  
  # create list for qc
  seurat_list <- list(
    "crclm01_cst" = seu_crclm01_cst,
    "crclm01_csn" = seu_crclm01_csn,
    "crclm01_can" = seu_crclm01_can,
    "crclm01_lit" = seu_crclm01_lit,
    "crclm01_lin" = seu_crclm01_lin,
    
    "crclm02_cxt" = seu_crclm02_cxt,
    "crclm02_cxn" = seu_crclm02_cxn,
    "crclm03_lit" = seu_crclm03_lit,
    "crclm03_lin" = seu_crclm03_lin,
    
    "crclm04_cxt" = seu_crclm04_cxt,
    "crclm04_cxn" = seu_crclm04_cxn,
    "crclm04_lin" = seu_crclm04_lin,
    "crclm05_cxt" = seu_crclm05_cxt,
    "crclm05_cxn" = seu_crclm05_cxn,
    "crclm05_lit" = seu_crclm05_lit,
    "crclm05_lin" = seu_crclm05_lin
  )
  
  # qc annotations
  for(i in 1:length(seurat_list)){
    seurat_list[[i]][["sample"]] <- seurat_list[[i]][["orig.ident"]]
    seurat_list[[i]][["barcode"]] <- str_sub(rownames(seurat_list[[i]][["nCount_RNA"]]), end = -3) # cropping barcodes
    seurat_list[[i]][["percent_mitochondrial_genes"]] <- PercentageFeatureSet(seurat_list[[i]], pattern="^MT-") # calculating the percentage of mitochondrial genes with Seurat function
    seurat_list[[i]][["patient"]] <- paste0("patient_",str_sub(levels(seurat_list[[i]]@active.ident),6,7)) # patient annotation
    switch(str_sub(levels(seurat_list[[i]]@active.ident),6,7),
           "01" = seurat_list[[i]][["batch"]] <- "batch_01",
           "02" = seurat_list[[i]][["batch"]] <- "batch_02",
           "03" = seurat_list[[i]][["batch"]] <- "batch_02",
           "04" = seurat_list[[i]][["batch"]] <- "batch_03",
           "05" = seurat_list[[i]][["batch"]] <- "batch_03",
    )
  }                       
  
  
##### QUALITY CONTROL #####
  
  # QC VlnPlots of merged objects
  seu_merge_batch_01 <- merge(x=seurat_list[["crclm01_cst"]], 
                              y=c(seurat_list[["crclm01_csn"]],
                                  seurat_list[["crclm01_can"]],
                                  seurat_list[["crclm01_lit"]],
                                  seurat_list[["crclm01_lin"]]
                              ),
                              add.cell.ids=c(levels(seurat_list[["crclm01_cst"]]$sample),
                                             levels(seurat_list[["crclm01_csn"]]$sample),
                                             levels(seurat_list[["crclm01_can"]]$sample),
                                             levels(seurat_list[["crclm01_lit"]]$sample),
                                             levels(seurat_list[["crclm01_lin"]]$sample)
                              ),
                              merge.data=TRUE,
                              project="crclm_batch_01")
  
  
  seu_merge_batch_02 <- merge(x=seurat_list[["crclm02_cxt"]], 
                              y=c(seurat_list[["crclm02_cxn"]],
                                  seurat_list[["crclm03_lit"]],
                                  seurat_list[["crclm03_lin"]]
                              ),
                              add.cell.ids=c(levels(seurat_list[["crclm02_cxt"]]$sample),
                                             levels(seurat_list[["crclm02_cxn"]]$sample),
                                             levels(seurat_list[["crclm03_lit"]]$sample),
                                             levels(seurat_list[["crclm03_lin"]]$sample)
                              ), merge.data=TRUE, project="crclm_batch_02")
  
  
  seu_merge_batch_03 <- merge(x=seurat_list[["crclm04_cxt"]], 
                              y=c(seurat_list[["crclm04_cxn"]],
                                  seurat_list[["crclm04_lin"]],
                                  seurat_list[["crclm05_cxt"]],
                                  seurat_list[["crclm05_cxn"]],
                                  seurat_list[["crclm05_lit"]],
                                  seurat_list[["crclm05_lin"]]
                              ),
                              add.cell.ids=c(levels(seurat_list[["crclm04_cxt"]]$sample),
                                             levels(seurat_list[["crclm04_cxn"]]$sample),
                                             levels(seurat_list[["crclm04_lin"]]$sample),
                                             levels(seurat_list[["crclm05_cxt"]]$sample),
                                             levels(seurat_list[["crclm05_cxn"]]$sample),
                                             levels(seurat_list[["crclm05_lit"]]$sample),
                                             levels(seurat_list[["crclm05_lin"]]$sample)
                              ),
                              merge.data=TRUE,
                              project="crclm_batch_03")
  
    
  seurat_list <- list("crclm_batch_01"=seu_merge_batch_01,
                      "crclm_batch_02"=seu_merge_batch_02,
                      "crclm_batch_03"=seu_merge_batch_03)
  

  threshold_list <- list()
  threshold_fraction_mito_gene <- 30 # should range between 10 and 30
  
  for (i in 1:length(seurat_list)) {
    threshold_list[[names(seurat_list)[i]]] <- subset_methodically(seurat_list[[i]], bool_mito = TRUE, bool_nCount = TRUE, bool_nFeature = TRUE, bool_subset = FALSE, bool_toList = TRUE, bool_print = FALSE) # save thresholds to list according to procedure as described in methods with self-developed function subset_methodically
    seurat_list[[i]] <- subset_methodically(seurat_list[[i]], bool_mito = FALSE, bool_nCount = TRUE, bool_nFeature = TRUE, bool_subset = TRUE)
    seurat_list[[i]] <- subset(seurat_list[[i]], subset = percent_mitochondrial_genes < threshold_fraction_mito_gene)
  }
  
 
  
  
##### MERGE SEURAT OBJECTS ####                     


  # simple merge of batch 01
  seu_merge_batch_01 <- merge(x=seurat_list[["crclm01_cst"]], 
                              y=c(seurat_list[["crclm01_csn"]],
                                  seurat_list[["crclm01_can"]],
                                  seurat_list[["crclm01_lit"]],
                                  seurat_list[["crclm01_lin"]]
                              ),
                              add.cell.ids=c(levels(seurat_list[["crclm01_cst"]]$sample),
                                             levels(seurat_list[["crclm01_csn"]]$sample),
                                             levels(seurat_list[["crclm01_can"]]$sample),
                                             levels(seurat_list[["crclm01_lit"]]$sample),
                                             levels(seurat_list[["crclm01_lin"]]$sample)
                              ),
                              merge.data=TRUE,
                              project="crclm_batch_01")
  
  
  # simple merge of batch 02
  seu_merge_batch_02 <- merge(x=seurat_list[["crclm02_cxt"]], 
                              y=c(seurat_list[["crclm02_cxn"]],
                                  seurat_list[["crclm03_lit"]],
                                  seurat_list[["crclm03_lin"]]
                              ),
                              add.cell.ids=c(levels(seurat_list[["crclm02_cxt"]]$sample),
                                             levels(seurat_list[["crclm02_cxn"]]$sample),
                                             levels(seurat_list[["crclm03_lit"]]$sample),
                                             levels(seurat_list[["crclm03_lin"]]$sample)
                              ), merge.data=TRUE, project="crclm_batch_02")
  
  
  # simple merge of batch 03
  seu_merge_batch_03 <- merge(x=seurat_list[["crclm04_cxt"]], 
                              y=c(seurat_list[["crclm04_cxn"]],
                                  seurat_list[["crclm04_lin"]],
                                  seurat_list[["crclm05_cxt"]],
                                  seurat_list[["crclm05_cxn"]],
                                  seurat_list[["crclm05_lit"]],
                                  seurat_list[["crclm05_lin"]]
                              ),
                              add.cell.ids=c(levels(seurat_list[["crclm04_cxt"]]$sample),
                                             levels(seurat_list[["crclm04_cxn"]]$sample),
                                             levels(seurat_list[["crclm04_lin"]]$sample),
                                             levels(seurat_list[["crclm05_cxt"]]$sample),
                                             levels(seurat_list[["crclm05_cxn"]]$sample),
                                             levels(seurat_list[["crclm05_lit"]]$sample),
                                             levels(seurat_list[["crclm05_lin"]]$sample)
                              ),
                              merge.data=TRUE,
                              project="crclm_batch_03")
  
  
  seurat_list <- list("crclm_batch_01"=seu_merge_batch_01,
                      "crclm_batch_02"=seu_merge_batch_02,
                      "crclm_batch_03"=seu_merge_batch_03)
  
  
  pca_met="pca"
  number_of_var_genes <- 5000 # number of most variable genes to be found
  cluster_resolution <- 0.4 # parameter for cluster resolution
  de_test="MAST"
  
  # Selcting Integration Features
  features <- SelectIntegrationFeatures(seurat_list, nfeatures = 6000, verbose = FALSE)
  # Preparing SCT Integration
  seurat_list <- PrepSCTIntegration(seurat_list, anchor.features = features, verbose = FALSE)
  #Finding Integration Anchors
  cca_anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:40, reduction = "cca", anchor.features = features, normalization.method = "LogNormalize", verbose = FALSE) # includes normalization of the merged data
  
  #Integrating data
  seu_merge <- IntegrateData(anchorset = cca_anchors, dims=1:40, normalization.method = "LogNormalize", verbose = FALSE)
  seu_merge@project.name <- "CRCLM-COMPLETE"
  

##### WORKFLOW FOR CRCLM COMPLETE #####

  # normalization, scaling and variable feature calculations on unmerged data
  seu_merge <- NormalizeData(seu_merge, assay = "RNA" )
  seu_merge <- ScaleData(seu_merge, verbose=TRUE, assay = "RNA")
  seu_merge <- FindVariableFeatures(seu_merge, selection.method="vst", nfeatures=6000, assay = "RNA", verbose=TRUE)
  
  # scaling and variable feature calculations of the merged data
  seu_merge <- ScaleData(seu_merge, verbose=TRUE, assay = seu_merge@active.assay)
  seu_merge <- FindVariableFeatures(seu_merge, selection.method="vst", nfeatures=6000, assay = seu_merge@active.assay, verbose=TRUE)
  
  pca_met <- "pca"
  seu_merge <- RunPCA(seu_merge, verbose = TRUE, npcs = 50, reduction.name = "pca",assay = seu_merge@active.assay) # calculating principal components
  number_of_pca <- estimate_inflection_point(seu_merge) # estimate the inflection point in the PCA data
  
  seu_merge <- FindNeighbors(seu_merge, dims = 1:number_of_pca, reduction = pca_met) # knn-algorithm
  seu_merge <- FindClusters(seu_merge, resolution = 0.35)  # resolution controls granularity of clustering
  seu_merge <- RunUMAP(seu_merge,
                       dims = 1:number_of_pca, 
                       reduction = pca_met,n.epochs = 500)
  
  # Calculate cell cycle scores according to Seurat vignette
  s_genes <- cc.genes.updated.2019$s.genes
  g2m_genes <-  cc.genes.updated.2019$g2m.genes
  seu_merge <- CellCycleScoring(seu_merge, s.features=s_genes, g2m.features=g2m_genes, set.ident=FALSE)
  seu_merge$cc_difference <- seu_merge$S.Score - seu_merge$G2M.Score
  tmp <- add_rel_prolif(seu_merge) # calculates the fraction of cycling cells (G2M, M cells) regarding all cells per cluster
  seu_merge <- add_prolif_per_cluster(seu_merge,tmp[[1]])
  
  # numeric batch annotation for MAST
  batch <- seu_merge$batch
  batch[batch=="batch_01"] <- 1
  batch[batch=="batch_02"] <- 2
  batch[batch=="batch_03"] <- 3
  seu_merge$numeric_batch <- batch
  
  tmp <- FindAllMarkers(seu_merge, test.use="MAST", min.pct = 0.40, logfc.threshold = 0.1, only.pos = FALSE, slot="data", features = seu_merge[[seu_merge@active.assay]]@var.features, assay = "RNA", latent.vars = "numeric_batch") # find differentially expressed genes with MAST
  tmp <- add_pct_ratio(tmp) # add to table the ratio of pct.1 and pct.2, e.g. how specific is the rel. expression within a cluster compared to all others

  
  Idents(seu_merge) <- seu_merge$seurat_clusters

  # pos_de_marker table with pos log fold changes
  # list with top marker genes grouped by highest log fold change, highest pct ratio and highest relative expression within a cluster
  top_marker_genes <- subset(seu_merge[[seu_merge@active.assay]]@misc[["pos_de_marker"]], seu_merge[[seu_merge@active.assay]]@misc[["pos_de_marker"]]$avg_log2FC > 0.5) %>% 
    group_by(cluster) %>%
    top_n(n=75, wt=avg_log2FC) %>% 
    dplyr::arrange(desc(avg_log2FC), .by_group=TRUE)
  
  top_marker_genes2 <- subset(seu_merge[[seu_merge@active.assay]]@misc[["pos_de_marker"]], seu_merge[[seu_merge@active.assay]]@misc[["pos_de_marker"]]$avg_log2FC > 0.5) %>% 
    group_by(cluster) %>%
    top_n(n=75, wt=pct.ratio) %>% 
    dplyr::arrange(desc(pct.ratio), .by_group=TRUE)
  
  top_marker_genes3 <- subset(seu_merge[[seu_merge@active.assay]]@misc[["pos_de_marker"]], seu_merge[[seu_merge@active.assay]]@misc[["pos_de_marker"]]$avg_log2FC > 0.5) %>% 
    group_by(cluster) %>%
    top_n(n=75, wt=pct.1) %>% 
    dplyr::arrange(desc(pct.1), .by_group=TRUE)
  
  # find a set of most expressive marker genes per cluster used for cell type annotation
  marker_gene_proposal <- list()
  for (cluster in sort(unique(seu_merge$seurat_clusters))) {
    marker_gene_proposal[[paste0("cluster.",cluster)]] <- propose_gene_marker_set(seu_merge, cluster_id = cluster, assay="RNA", top_n_number = 5)
  }
  
  
  
  
  
  
  
  