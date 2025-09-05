# --- Install packages if not already installed ---

# CRAN packages
cran_packages <- c("data.table", "Matrix", "dplyr", "ggplot2", "tidyr", "devtools", "sctransform", "hdf5r", "harmony","pheatmap")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}
# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
bioc_packages <- c("rhdf5", "SingleR", "celldex")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}
# GitHub packages
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
# Install GitHub packages with force=TRUE to update
remotes::install_github("satijalab/seurat", ref = "develop", force = TRUE)
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", force = TRUE)
devtools::install_github("immunogenomics/presto")
# Additional packages
if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
  install.packages("SeuratDisk")
}
if (!requireNamespace("SeuratObject", quietly = TRUE)) {
  install.packages("SeuratObject")
}

# --- Load libraries ---
library(Seurat)
library(Matrix)
library(dplyr)
library(data.table)
library(rhdf5)
library(hdf5r)
library(sctransform)
library(ggplot2)
library(SingleR)
library(celldex)
library(harmony)
library(SeuratDisk)
library(SeuratObject)
library(tidyr)
library(devtools)
library(presto)
library(pheatmap)


library(data.table)
library(Seurat)
library(dplyr)
library(ggplot2)

# 1. Load Seurat objects
epi_clean <- readRDS("epi_clean.rds")
endo_clean <- readRDS("endo_clean.rds")
fib_clean <- readRDS("fib_clean.rds")

# 2. Load CellTypist prediction CSVs (from Python output)
epi_pred <- fread("epi_predictions.csv", header = TRUE)
endo_pred <- fread("endo_predictions.csv", header = TRUE)
fib_pred <- fread("fib_predictions.csv", header = TRUE)

# 3. Rename the first column from 'V1' (barcode) to 'Cell'
setnames(epi_pred, "V1", "Cell")
setnames(endo_pred, "V1", "Cell")
setnames(fib_pred, "V1", "Cell")

# 4. Set keys for joining
setkey(epi_pred, Cell)
setkey(endo_pred, Cell)
setkey(fib_pred, Cell)

# 5. Extract metadata from Seurat objects with 'Cell' barcode as a column
epi_meta <- as.data.table(epi_clean@meta.data, keep.rownames = "Cell")
endo_meta <- as.data.table(endo_clean@meta.data, keep.rownames = "Cell")
fib_meta <- as.data.table(fib_clean@meta.data, keep.rownames = "Cell")

setkey(epi_meta, Cell)
setkey(endo_meta, Cell)
setkey(fib_meta, Cell)

# 6. Join CellTypist predictions into metadata tables by Cell barcode
epi_annotated <- epi_pred[epi_meta]
endo_annotated <- endo_pred[endo_meta]
fib_annotated <- fib_pred[fib_meta]

# 7. Add predicted labels back to Seurat metadata (excluding 'Cell' column)
epi_clean <- AddMetaData(epi_clean, metadata = epi_annotated[, !c("Cell"), with = FALSE])
endo_clean <- AddMetaData(endo_clean, metadata = endo_annotated[, !c("Cell"), with = FALSE])
fib_clean <- AddMetaData(fib_clean, metadata = fib_annotated[, !c("Cell"), with = FALSE])

# 8. Optional: Normalize, find variable features, scale data, PCA, UMAP
for (obj in list(epi_clean, endo_clean, fib_clean)) {
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, npcs = 30)
  obj <- RunUMAP(obj, dims = 1:20)
}

# Save updated objects if desired
saveRDS(epi_clean, file = "epi_clean_annotated.rds")
saveRDS(endo_clean, file = "endo_clean_annotated.rds")
saveRDS(fib_clean, file = "fib_clean_annotated.rds")

# 9. Plot UMAPs colored by CellTypist predicted labels
DimPlot(epi_clean, reduction = "umap", group.by = "predicted_labels", label = TRUE) + ggtitle("Epithelial Cell Types")
DimPlot(endo_clean, reduction = "umap", group.by = "predicted_labels", label = TRUE) + ggtitle("Endothelial Cell Types")
DimPlot(fib_clean, reduction = "umap", group.by = "predicted_labels", label = TRUE) + ggtitle("Fibroblast Cell Types")

# 10. Functions for filtering and plotting heatmaps by predicted labels

filter_large_clusters <- function(seurat_obj, group_col, min_cells = 100) {
  cluster_sizes <- table(seurat_obj[[group_col]])
  large_clusters <- names(cluster_sizes[cluster_sizes >= min_cells])
  return(large_clusters)
}

plot_markers_heatmap_filtered <- function(seurat_obj, group_col, prefix, min_cells = 100) {
  large_clusters <- filter_large_clusters(seurat_obj, group_col, min_cells)
  cells_to_keep <- WhichCells(seurat_obj, expression = get(group_col) %in% large_clusters)
  seurat_sub <- subset(seurat_obj, cells = cells_to_keep)
  
  # Make sure data is scaled for heatmap
  if (is.null(seurat_sub@assays$RNA@scale.data) || ncol(seurat_sub@assays$RNA@scale.data) == 0) {
    seurat_sub <- ScaleData(seurat_sub)
  }
  
  markers <- FindAllMarkers(seurat_sub, group.by = group_col, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  top_genes <- markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) %>%
    pull(gene) %>%
    unique()
  
  heatmap_plot <- DoHeatmap(seurat_sub, features = top_genes, group.by = group_col) +
    NoLegend() + ggtitle(paste0("Top markers - ", prefix, " (filtered clusters)"))
  
  print(heatmap_plot)
  return(markers)
}

# 11. Plot heatmaps for all three cell types
epi_markers_filtered <- plot_markers_heatmap_filtered(epi_clean, "predicted_labels", "Epi", min_cells = 100)
endo_markers_filtered <- plot_markers_heatmap_filtered(endo_clean, "predicted_labels", "Endo", min_cells = 100)
fib_markers_filtered <- plot_markers_heatmap_filtered(fib_clean, "predicted_labels", "Fib", min_cells = 100)


