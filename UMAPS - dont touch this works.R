#install and load libraries
remotes::install_github("satijalab/seurat", ref = "develop", force = TRUE)
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force=TRUE)
library(Seurat)
library(Matrix)
library(dplyr)
library(data.table)

if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install("rhdf5")
install.packages('hdf5r', type='binary')
library(hdf5r)

install.packages('sctransform')
library(sctransform)
library(ggplot2)

if (!requireNamespace("SingleR", quietly = TRUE)) {
  BiocManager::install("SingleR")
}
library(SingleR)

BiocManager::install("celldex")
library(celldex)

install.packages('harmony')
library(harmony)
library(SeuratDisk)
library(SeuratObject)

#set root folder containing all data
gse128169root <- 'GSE128169_RAW'

#geo sample groups
ctrl_gsms <- c("GSM3666096","GSM3666097","GSM3666098","GSM3666099","GSM3666100")
ssc_gsms  <- c("GSM3666101","GSM3666102","GSM3666103","GSM3666104",
               "GSM3666105","GSM3666106","GSM3666107","GSM3666108")

#function to infer region
infer_region <- function(sc_name) {
  
  if (grepl("LOW", sc_name, ignore.case=TRUE)) return("Lower")
  
  if (grepl("UP", sc_name, ignore.case=TRUE)) return("Upper")
  
  return(NA)
}

#start list
seurat_list <- list()

#load in hpca reference for labelling
reflung <- celldex::HumanPrimaryCellAtlasData()

#process .mtx samples
barcode_files <- list.files(gse128169root, pattern="_barcodes.tsv.gz$", 
                            full.names=TRUE, recursive=TRUE)

for (i in seq_along(barcode_files)) {
  
  prefix <- sub("_barcodes.tsv.gz", "", basename(barcode_files[i]))
  gsm_id <- unlist(strsplit(prefix, "_"))[1]
  
  if (!(gsm_id %in% c(ctrl_gsms, ssc_gsms))) {
    message("skipping sample: ", gsm_id)
    next
  }
  
  message("processing .mtx sample: ", gsm_id)
  sample_dir <- dirname(barcode_files[i])
  
  barcodes_file <- file.path(sample_dir, paste0(prefix, "_barcodes.tsv.gz"))
  genes_file    <- file.path(sample_dir, paste0(prefix, "_genes.tsv.gz"))
  matrix_file   <- file.path(sample_dir, paste0(prefix, "_matrix.mtx.gz"))
  
  barcodes <- fread(barcodes_file, header=FALSE)$V1
  genes    <- fread(genes_file, header=FALSE)$V2
  mat_triplets <- fread(matrix_file, header=FALSE)
  
  counts <- sparseMatrix(
    i = mat_triplets$V1,
    j = mat_triplets$V2,
    x = mat_triplets$V3,
    dims = c(length(genes), length(barcodes))
  )
  
  genes_unique <- make.unique(genes)
  rownames(counts) <- genes_unique
  colnames(counts) <- paste0(gsm_id, "_", barcodes)
  
  so <- CreateSeuratObject(counts=counts, project=gsm_id)
  so$sample_id <- gsm_id
  so$disease <- ifelse(gsm_id %in% ctrl_gsms, "Control", "SSc-ILD")
  so$lung_region <- infer_region(prefix)
  
  so <- subset(so, subset=nFeature_RNA > 200 & nCount_RNA > 500)
  
  so <- NormalizeData(so)
  so <- FindVariableFeatures(so)
  so <- ScaleData(so)
  
  seurat_list[[gsm_id]] <- so
}

#function to get dense expression matrix for SingleR annotation
get_dense_expr_matrix <- function(seurat_obj) {
  
  sample_layers <- grep("^data\\.", names(seurat_obj[["RNA"]]@layers), value = TRUE)
  
  expr_matrix_list <- lapply(sample_layers, function(layer) {
    
    mat <- seurat_obj[["RNA"]]@layers[[layer]]
    
    if (is.null(rownames(mat))) {
      rownames(mat) <- rownames(seurat_obj[["RNA"]])
    }
    
    return(mat)
  })
  
  expr_matrix <- do.call(cbind, expr_matrix_list)
  
  expr_matrix_dense <- as.matrix(expr_matrix)
  
  return(expr_matrix_dense)
}

#split by condition and region
ctrl_list  <- seurat_list[names(seurat_list) %in% ctrl_gsms]
ssc_list   <- seurat_list[names(seurat_list) %in% ssc_gsms]
# upper_list <- Filter(function(x) x$lung_region == "Upper", ssc_list)
# lower_list <- Filter(function(x) x$lung_region == "Lower", ssc_list)
upper_list <- Filter(function(x) unique(x$lung_region) == "Upper", ssc_list)
lower_list <- Filter(function(x) unique(x$lung_region) == "Lower", ssc_list)


#merge control samples
ctrl_combined <- ctrl_list[[1]]
if (length(ctrl_list) > 1) {
  for (i in 2:length(ctrl_list)) {
    ctrl_combined <- merge(ctrl_combined, y = ctrl_list[[i]], merge.data = TRUE)
  }
}
sapply(upper_list, function(x) dim(x)[2])

#scale and pca
ctrl_combined <- ScaleData(ctrl_combined)
ctrl_combined <- RunPCA(ctrl_combined)
#run harmony
ctrl_combined <- RunHarmony(ctrl_combined, group.by.vars = "sample_id")
#dimensional reduction
ctrl_combined <- RunUMAP(ctrl_combined, reduction = "harmony", dims = 1:30)
ctrl_combined <- FindNeighbors(ctrl_combined, reduction = "harmony", dims = 1:30)
ctrl_combined <- FindClusters(ctrl_combined, resolution = 0.5)

#singleR annotation for controls
ctrl_expr_dense <- get_dense_expr_matrix(ctrl_combined)

ctrl_combined$SingleR.labels <- SingleR(
  test = ctrl_expr_dense,
  ref = reflung,
  labels = reflung$label.main
)$labels

#plot UMAP
DimPlot(ctrl_combined, group.by = "SingleR.labels", label = TRUE, repel = TRUE) + 
  ggtitle("Control UMAP (Harmony)")



#merge upper lobe ssc-ild samples
upper_combined <- upper_list[[1]]
if (length(upper_list) > 1) {
  for (i in 2:length(upper_list)) {
    upper_combined <- merge(upper_combined, y = upper_list[[i]], merge.data = TRUE)
  }
}

#scale and pca
upper_combined <- ScaleData(upper_combined)
upper_combined <- RunPCA(upper_combined)
#run harmony
upper_combined <- RunHarmony(upper_combined, group.by.vars = "sample_id")
#dimensional reduction
upper_combined <- RunUMAP(upper_combined, reduction = "harmony", dims = 1:30)
upper_combined <- FindNeighbors(upper_combined, reduction = "harmony", dims = 1:30)
upper_combined <- FindClusters(upper_combined, resolution = 0.5)

#singleR annotation for upper ssc-ild
rm(expr_matrix_dense) #space saving - wipe dense matrix after used each time
upper_expr_dense <- get_dense_expr_matrix(upper_combined)

upper_combined$SingleR.labels <- SingleR(
  test = upper_expr_dense,
  ref = reflung,
  labels = reflung$label.main
)$labels

#plot UMAP
DimPlot(upper_combined, group.by = "SingleR.labels", label = TRUE, repel = TRUE) + 
  ggtitle("Upper SSC-ILD UMAP (Harmony)")

#merge lower lobe ssc-ild samples
lower_combined <- lower_list[[1]]
if (length(lower_list) > 1) {
  for (i in 2:length(lower_list)) {
    lower_combined <- merge(lower_combined, y = lower_list[[i]], merge.data = TRUE)
  }
}

#scale and pca
lower_combined <- ScaleData(lower_combined)
lower_combined <- RunPCA(lower_combined)
#run harmony
lower_combined <- RunHarmony(lower_combined, group.by.vars = "sample_id")
#dimensional reduction
lower_combined <- RunUMAP(lower_combined, reduction = "harmony", dims = 1:30)
lower_combined <- FindNeighbors(lower_combined, reduction = "harmony", dims = 1:30)
lower_combined <- FindClusters(lower_combined, resolution = 0.5)

#singleR annotation for lower ssc-ild
rm(upper_expr_dense) #space saving - wipe dense matrix after used each time
lower_expr_dense <- get_dense_expr_matrix(lower_combined)

lower_combined$SingleR.labels <- SingleR(
  test = lower_expr_dense,
  ref = reflung,
  labels = reflung$label.main
)$labels

#plot UMAP
DimPlot(lower_combined, group.by = "SingleR.labels", label = TRUE, repel = TRUE) + 
  ggtitle("Lower SSC-ILD UMAP (Harmony)")

#save for later use!!
#combined seurats 
saveRDS(ctrl_combined, file = "ctrl_combined.rds")
saveRDS(upper_combined, file = "upper_combined.rds")
saveRDS(lower_combined, file = "lower_combined.rds")
#individual objs for each sample
for (name in names(seurat_list)) {
  saveRDS(seurat_list[[name]], file = paste0("seurat_", name, ".rds"))
}


