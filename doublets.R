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
library(DoubletFinder)

#Load cleaned Seurat objs
epi_clean <- readRDS("epi_clean.RDS")
endo_clean <- readRDS("endo_clean.RDS")
fib_clean <- readRDS("fib_clean.RDS")
# For epithelial cells
epi_clean@meta.data$lung_region[is.na(epi_clean@meta.data$lung_region)] <- "Control"
# For endothelial cells
endo_clean@meta.data$lung_region[is.na(endo_clean@meta.data$lung_region)] <- "Control"
# For fibroblast cells
fib_clean@meta.data$lung_region[is.na(fib_clean@meta.data$lung_region)] <- "Control"


epi_clean@meta.data$lung_region <- ifelse(
  epi_clean@meta.data$disease == "Control",
  "Control",
  epi_clean@meta.data$lung_region
)

endo_clean@meta.data$lung_region <- ifelse(
  endo_clean@meta.data$disease == "Control",
  "Control",
  endo_clean@meta.data$lung_region
)

fib_clean@meta.data$lung_region <- ifelse(
  fib_clean@meta.data$disease == "Control",
  "Control",
  fib_clean@meta.data$lung_region
)

view(endo_clean@meta.data)

process_doublets <- function(seurat_obj) {
  # 1. Preprocessing (if not done yet)
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  
  # 2. Parameter sweep to find optimal pK
  sweep.res.list <- paramSweep(seurat_obj, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  optimal_pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  cat("Optimal pK:", optimal_pK, "\n")
  
  # 3. Estimate expected doublets (5% here, adjust as needed)
  nExp <- round(0.05 * ncol(seurat_obj))
  
  # 4. Run DoubletFinder
  seurat_obj <- doubletFinder(seurat_obj, PCs = 1:10, pN = 0.25, pK = optimal_pK, nExp = nExp)
  
  return(seurat_obj)
}

# Run for each object
epi_clean <- process_doublets(epi_clean)
endo_clean <- process_doublets(endo_clean)
fib_clean <- process_doublets(fib_clean)

# Extract the DF classification column name for each object (should be a single string)
epi_df_col <- grep("DF.classifications", colnames(epi_clean@meta.data), value = TRUE)
endo_df_col <- grep("DF.classifications", colnames(endo_clean@meta.data), value = TRUE)
fib_df_col <- grep("DF.classifications", colnames(fib_clean@meta.data), value = TRUE)

# Subset doublets from each Seurat object
epi_doublets <- subset(epi_clean, subset = !!sym(epi_df_col) == "Doublet")
endo_doublets <- subset(endo_clean, subset = !!sym(endo_df_col) == "Doublet")
fib_doublets <- subset(fib_clean, subset = !!sym(fib_df_col) == "Doublet")

epi_clean
#SVAE!!!
saveRDS(epi_doublets, file = "epi_doublets.rds")
saveRDS(endo_doublets, file = "endo_doublets.rds")
saveRDS(fib_doublets, file = "fib_doublets.rds")
saveRDS(epi_clean, file='epi_clean_df.rds')
saveRDS(endo_clean, file='endo_clean_df.rds')
saveRDS(fib_clean, file='fib_clean_df.rds')

# Epithelial coloured by doublets
DimPlot(epi_clean, group.by = epi_df_col, pt.size = 1) + ggtitle("Epithelial: Singlets vs Doublets")
# Endothelial
DimPlot(endo_clean, group.by = endo_df_col, pt.size = 1) + ggtitle("Endothelial: Singlets vs Doublets")
# Fibroblast
DimPlot(fib_clean, group.by = fib_df_col, pt.size = 1) + ggtitle("Fibroblast: Singlets vs Doublets")
#plot rna metrics qc
VlnPlot(epi_clean, features = c("nFeature_RNA", "nCount_RNA"), group.by = epi_df_col, pt.size = 0) + ggtitle("Epithelial QC metrics")
VlnPlot(endo_clean, features = c("nFeature_RNA", "nCount_RNA"), group.by = endo_df_col, pt.size = 0) + ggtitle("Endothelial QC metrics")
VlnPlot(fib_clean, features = c("nFeature_RNA", "nCount_RNA"), group.by = fib_df_col, pt.size = 0) + ggtitle("Fibroblast QC metrics")

#by region - controls vs ssc-ild upper and lower
# Extract the doublet classification column name for each object
epi_df_col <- grep("DF.classifications", colnames(epi_clean@meta.data), value = TRUE)
endo_df_col <- grep("DF.classifications", colnames(endo_clean@meta.data), value = TRUE)
fib_df_col <- grep("DF.classifications", colnames(fib_clean@meta.data), value = TRUE)

# Summarize doublets by lung_region for epithelial cells
doublet_summary_epi <- epi_clean@meta.data %>%
  group_by(lung_region) %>%
  summarise(
    total_cells = n(),
    doublets = sum(get(epi_df_col) == "Doublet"),
    doublet_rate = doublets / total_cells
  )

# Summarize doublets by lung_region for endothelial cells
doublet_summary_endo <- endo_clean@meta.data %>%
  group_by(lung_region) %>%
  summarise(
    total_cells = n(),
    doublets = sum(get(endo_df_col) == "Doublet"),
    doublet_rate = doublets / total_cells
  )

# Summarize doublets by lung_region for fibroblast cells
doublet_summary_fib <- fib_clean@meta.data %>%
  group_by(lung_region) %>%
  summarise(
    total_cells = n(),
    doublets = sum(get(fib_df_col) == "Doublet"),
    doublet_rate = doublets / total_cells
  )

# Plot doublet rates for epithelial cells
ggplot(doublet_summary_epi, aes(x = lung_region, y = doublet_rate)) +
  geom_col(fill = "steelblue") +
  labs(title = "Doublet Rate by Lung Region (Epithelial)", y = "Doublet Proportion", x = "Lung Region") +
  theme_minimal()

# Plot doublet rates for endothelial cells
ggplot(doublet_summary_endo, aes(x = lung_region, y = doublet_rate)) +
  geom_col(fill = "forestgreen") +
  labs(title = "Doublet Rate by Lung Region (Endothelial)", y = "Doublet Proportion", x = "Lung Region") +
  theme_minimal()

# Plot doublet rates for fibroblast cells
ggplot(doublet_summary_fib, aes(x = lung_region, y = doublet_rate)) +
  geom_col(fill = "darkorange") +
  labs(title = "Doublet Rate by Lung Region (Fibroblast)", y = "Doublet Proportion", x = "Lung Region") +
  theme_minimal()


#set predicted_labels as idents
# For epithelial
epi_clean$celltype <- epi_clean$predicted_labels
Idents(epi_clean) <- epi_clean$celltype
# For endothelial
endo_clean$celltype <- endo_clean$predicted_labels
Idents(endo_clean) <- endo_clean$celltype
# For fibroblasts
fib_clean$celltype <- fib_clean$predicted_labels
Idents(fib_clean) <- fib_clean$celltype

#findallmarkers for each from pred labels
# Epithelial markers
epi_markers <- FindAllMarkers(
  object = epi_clean,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
# Endothelial markers
endo_markers <- FindAllMarkers(
  object = endo_clean,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
# Fibroblast markers
fib_markers <- FindAllMarkers(
  object = fib_clean,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

#save markers for future use
write.csv(epi_markers, file = "epi_predicted_labels_markers.csv")
write.csv(endo_markers, file = "endo_predicted_labels_markers.csv")
write.csv(fib_markers, file = "fib_predicted_labels_markers.csv")

#step 1: top 50 for each pred cell type
top_markers_epi <- epi_markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) %>%
  pull(gene) %>% unique()
top_markers_endo <- endo_markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) %>%
  pull(gene) %>% unique()
top_markers_fib <- fib_markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) %>%
  pull(gene) %>% unique()

# --- Step 2: Extract normalized expression matrix from doublets ---
get_doublets_expression <- function(doublet_obj, markers) {
  mat <- GetAssayData(doublet_obj, assay = "RNA", slot = "data")
  genes_present <- intersect(rownames(mat), markers)
  mat[genes_present, , drop = FALSE]
}

epi_exp <- get_doublets_expression(epi_doublets, top_markers_epi)
endo_exp <- get_doublets_expression(endo_doublets, top_markers_endo)
fib_exp <- get_doublets_expression(fib_doublets, top_markers_fib)

# --- Step 3: Calculate scores per cell (sum expression of marker genes) ---
calc_scores <- function(exp_mat) {
  Matrix::colSums(exp_mat)
}

epi_scores <- calc_scores(epi_exp)
endo_scores <- calc_scores(endo_exp)
fib_scores <- calc_scores(fib_exp)

write.csv(epi_scores, file="epi_scores.csv")
write.csv(endo_scores, file="endo_scores.csv")
write.csv(fib_scores, file="fib_scores.csv")

# --- Step 4: Combine scores for mixed expression detection ---

# Get union of markers to cover for cross-checks (example epi vs endo)
markers_epi_endo <- unique(c(top_markers_epi, top_markers_endo))

# Extract expression for epi_doublets for those combined markers
epi_exp_all <- GetAssayData(epi_doublets, assay = "RNA", slot = "data")
epi_exp_subset <- epi_exp_all[intersect(rownames(epi_exp_all), markers_epi_endo), , drop = FALSE]

# Calculate scores separately
epi_score_in_epi_doublets <- Matrix::colSums(epi_exp_subset[intersect(top_markers_epi, rownames(epi_exp_subset)), , drop = FALSE])
endo_score_in_epi_doublets <- Matrix::colSums(epi_exp_subset[intersect(top_markers_endo, rownames(epi_exp_subset)), , drop = FALSE])

scores_epi_doublets <- data.frame(
  cell = colnames(epi_exp_subset),
  epi_score = epi_score_in_epi_doublets,
  endo_score = endo_score_in_epi_doublets
)

# Identify "mixed" doublets as those above 75th percentile in both
threshold_epi <- quantile(scores_epi_doublets$epi_score, 0.75)
threshold_endo <- quantile(scores_epi_doublets$endo_score, 0.75)

mixed_epi_doublets <- subset(scores_epi_doublets, epi_score > threshold_epi & endo_score > threshold_endo)

# --- Step 5: Plotting example for epi doublets co-expression ---

ggplot(scores_epi_doublets, aes(x = epi_score, y = endo_score)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = threshold_endo, linetype = "dashed", color = "red") +
  geom_vline(xintercept = threshold_epi, linetype = "dashed", color = "red") +
  labs(title = "Co-expression Scores in Epi/Endo Doublets",
       x = "Epi Marker Expression Score",
       y = "Endo Marker Expression Score")

# --- Repeat for endo_doublets (e.g. endo vs fib markers) ---

markers_endo_fib <- unique(c(top_markers_endo, top_markers_fib))
endo_exp_all <- GetAssayData(endo_doublets, assay = "RNA", slot = "data")
endo_exp_subset <- endo_exp_all[intersect(rownames(endo_exp_all), markers_endo_fib), , drop = FALSE]

endo_score_in_endo_doublets <- Matrix::colSums(endo_exp_subset[intersect(top_markers_endo, rownames(endo_exp_subset)), , drop = FALSE])
fib_score_in_endo_doublets <- Matrix::colSums(endo_exp_subset[intersect(top_markers_fib, rownames(endo_exp_subset)), , drop = FALSE])

scores_endo_doublets <- data.frame(
  cell = colnames(endo_exp_subset),
  endo_score = endo_score_in_endo_doublets,
  fib_score = fib_score_in_endo_doublets
)

threshold_endo_2 <- quantile(scores_endo_doublets$endo_score, 0.75)
threshold_fib <- quantile(scores_endo_doublets$fib_score, 0.75)

mixed_endo_doublets <- subset(scores_endo_doublets, endo_score > threshold_endo_2 & fib_score > threshold_fib)

ggplot(scores_endo_doublets, aes(x = endo_score, y = fib_score)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = threshold_fib, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = threshold_endo_2, linetype = "dashed", color = "blue") +
  labs(title = "Co-expression Scores in Endo/Fib Doublets",
       x = "Endo Marker Expression Score",
       y = "Fib Marker Expression Score")

# --- Repeat for fib_doublets (e.g. fib vs epi markers) ---

markers_fib_epi <- unique(c(top_markers_fib, top_markers_epi))
fib_exp_all <- GetAssayData(fib_doublets, assay = "RNA", slot = "data")
fib_exp_subset <- fib_exp_all[intersect(rownames(fib_exp_all), markers_fib_epi), , drop = FALSE]

fib_score_in_fib_doublets <- Matrix::colSums(fib_exp_subset[intersect(top_markers_fib, rownames(fib_exp_subset)), , drop = FALSE])
epi_score_in_fib_doublets <- Matrix::colSums(fib_exp_subset[intersect(top_markers_epi, rownames(fib_exp_subset)), , drop = FALSE])

scores_fib_doublets <- data.frame(
  cell = colnames(fib_exp_subset),
  fib_score = fib_score_in_fib_doublets,
  epi_score = epi_score_in_fib_doublets
)

threshold_fib_2 <- quantile(scores_fib_doublets$fib_score, 0.75)
threshold_epi_2 <- quantile(scores_fib_doublets$epi_score, 0.75)

mixed_fib_doublets <- subset(scores_fib_doublets, fib_score > threshold_fib_2 & epi_score > threshold_epi_2)

ggplot(scores_fib_doublets, aes(x = fib_score, y = epi_score)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = threshold_epi_2, linetype = "dashed", color = "green") +
  geom_vline(xintercept = threshold_fib_2, linetype = "dashed", color = "green") +
  labs(title = "Co-expression Scores in Fib/Epi Doublets",
       x = "Fib Marker Expression Score",
       y = "Epi Marker Expression Score")

#saved mixed doublet df for downstream use
saveRDS(mixed_epi_doublets, "mixed_epi_doublets.rds")
saveRDS(mixed_endo_doublets, "mixed_endo_doublets.rds")
saveRDS(mixed_fib_doublets, "mixed_fib_doublets.rds")

# Extract cell names of mixed doublets
mixed_epi_cells <- mixed_epi_doublets$cell
mixed_endo_cells <- mixed_endo_doublets$cell
mixed_fib_cells <- mixed_fib_doublets$cell
mixed_epi_cells
mixed_endo_cells
mixed_fib_cells

#ADD LUNG REGION****
# For epithelial mixed doublets
epi_region <- epi_clean@meta.data[mixed_epi_cells, "lung_region", drop = FALSE]
epi_region$cell <- rownames(epi_region)
# For endothelial mixed doublets
endo_region <- endo_clean@meta.data[mixed_endo_cells, "lung_region", drop = FALSE]
endo_region$cell <- rownames(endo_region)
# For fibroblast mixed doublets
fib_region <- fib_clean@meta.data[mixed_fib_cells, "lung_region", drop = FALSE]
fib_region$cell <- rownames(fib_region)
# Add lung_region to epithelial mixed doublets
mixed_epi_doublets <- merge(mixed_epi_doublets, epi_region, by = "cell")
# Add lung_region to endothelial mixed doublets
mixed_endo_doublets <- merge(mixed_endo_doublets, endo_region, by = "cell")
# Add lung_region to fibroblast mixed doublets
mixed_fib_doublets <- merge(mixed_fib_doublets, fib_region, by = "cell")
#SAVE
saveRDS(mixed_epi_doublets, "mixed_epi_doublets_with_region.rds")
saveRDS(mixed_endo_doublets, "mixed_endo_doublets_with_region.rds")
saveRDS(mixed_fib_doublets, "mixed_fib_doublets_with_region.rds")
write.csv(mixed_epi_doublets, file='epi_mixed_scored.csv')
write.csv(mixed_endo_doublets, file='endo_mixed_scored.csv')
write.csv(mixed_fib_doublets, file='fib_mixed_scored.csv')
write.csv(epi_mixed_doublets_obj@meta.data, file='epi_mixed_meta.csv')
write.csv(endo_mixed_doublets_obj@meta.data, file='endo_mixed_meta.csv')
write.csv(fib_mixed_doublets_obj@meta.data, file='fib_mixed_meta.csv')
write.csv(epi_clean@meta.data, file='epicleanmeta.csv')
write.csv(endo_clean@meta.data, file='endocleanmeta.csv')
write.csv(fib_clean@meta.data, file='fibcleanmeta.csv')

#colour by region!
ggplot(mixed_epi_doublets, aes(x = epi_score, y = endo_score, color = lung_region)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_hline(yintercept = threshold_endo, linetype = "dashed", color = "red") +
  geom_vline(xintercept = threshold_epi, linetype = "dashed", color = "red") +
  labs(title = "Co-expression Scores in Epi/Endo Doublets by Lung Region",
       x = "Epi Marker Expression Score",
       y = "Endo Marker Expression Score") +
  theme_minimal() +
  scale_color_manual(values = c("Upper" = "#1f78b4", "Lower" = "#33a02c", "Control" = "#fb9a99"))
#endo
ggplot(mixed_endo_doublets, aes(x = endo_score, y = fib_score, color = lung_region)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_hline(yintercept = threshold_fib, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = threshold_endo_2, linetype = "dashed", color = "blue") +
  labs(title = "Co-expression Scores in Endo/Fib Doublets by Lung Region",
       x = "Endo Marker Expression Score",
       y = "Fib Marker Expression Score") +
  theme_minimal() +
  scale_color_manual(values = c("Upper" = "#1f78b4", "Lower" = "#33a02c", "Control" = "#fb9a99"))
#fib
ggplot(mixed_fib_doublets, aes(x = fib_score, y = epi_score, color = lung_region)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_hline(yintercept = threshold_epi_2, linetype = "dashed", color = "green") +
  geom_vline(xintercept = threshold_fib_2, linetype = "dashed", color = "green") +
  labs(title = "Co-expression Scores in Fib/Epi Doublets by Lung Region",
       x = "Fib Marker Expression Score",
       y = "Epi Marker Expression Score") +
  theme_minimal() +
  scale_color_manual(values = c("Upper" = "#1f78b4", "Lower" = "#33a02c", "Control" = "#fb9a99"))

# Extract metadata for doublets using cell names
epi_doublet_meta <- epi_clean@meta.data[mixed_epi_cells, , drop = FALSE]
endo_doublet_meta <- endo_clean@meta.data[mixed_endo_cells, , drop = FALSE]
fib_doublet_meta <- fib_clean@meta.data[mixed_fib_cells, , drop = FALSE]

#COMPARING PREDICTED CELLTYPIST VS DOUBLET PROFILE
# Extract CellTypist predicted labels for mixed epithelial doublets
epi_doublet_cells <- mixed_epi_doublets$cell
epi_doublet_labels <- epi_clean@meta.data[epi_doublet_cells, "predicted_labels", drop = TRUE]
# Add predicted labels back into the mixed doublets dataframe for easier comparison
mixed_epi_doublets$predicted_label <- epi_doublet_labels
# Endothelial
endo_doublet_cells <- mixed_endo_doublets$cell
endo_doublet_labels <- endo_clean@meta.data[endo_doublet_cells, "predicted_labels", drop = TRUE]
mixed_endo_doublets$predicted_label <- endo_doublet_labels
# Fibroblast
fib_doublet_cells <- mixed_fib_doublets$cell
fib_doublet_labels <- fib_clean@meta.data[fib_doublet_cells, "predicted_labels", drop = TRUE]
mixed_fib_doublets$predicted_label <- fib_doublet_labels
# Example summary for epithelial mixed doublets
table(mixed_epi_doublets$predicted_label)
table(mixed_endo_doublets$predicted_label)
table(mixed_fib_doublets$predicted_label)

#subset mixed with clean to add rna data back
epi_mixed_doublets_obj <- subset(epi_clean, cells = mixed_epi_doublets$cell)
endo_mixed_doublets_obj <- subset(endo_clean, cells = mixed_endo_doublets$cell)
fib_mixed_doublets_obj <- subset(fib_clean, cells = mixed_fib_doublets$cell)
# Add a metadata column marking mixed doublets vs others
epi_clean$mixed_doublet_status <- ifelse(colnames(epi_clean) %in% mixed_epi_doublets$cell, "MixedDoublet", "Other")
Idents(epi_clean) <- epi_clean$mixed_doublet_status
# Find markers that distinguish mixed doublets
epi_mixed_markers <- FindMarkers(epi_clean, ident.1 = "MixedDoublet", ident.2 = "Other", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(epi_mixed_markers, 10)
endo_clean$mixed_doublet_status <- ifelse(colnames(endo_clean) %in% mixed_endo_doublets$cell, "MixedDoublet", "Other")
Idents(endo_clean) <- endo_clean$mixed_doublet_status
#endo mixed markers
endo_mixed_markers <- FindMarkers(endo_clean, ident.1 = "MixedDoublet", ident.2 = "Other", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
fib_clean$mixed_doublet_status <- ifelse(colnames(fib_clean) %in% mixed_fib_doublets$cell, "MixedDoublet", "Other")
#fib mixed markers 
Idents(fib_clean) <- fib_clean$mixed_doublet_status
fib_mixed_markers <- FindMarkers(fib_clean, ident.1 = "MixedDoublet", ident.2 = "Other", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

epi_mixed_markers
endo_mixed_markers
fib_mixed_markers

# Take top 50 marker genes from your existing marker results
top_epi_genes <- head(rownames(epi_mixed_markers), 50)
top_endo_genes <- head(rownames(endo_mixed_markers), 50)
top_fib_genes <- head(rownames(fib_mixed_markers), 50)
# Subset cells for epi and endo and fib clean objects
epi_cells <- colnames(epi_clean)
endo_cells <- colnames(endo_clean)
fib_cells <- colnames(fib_clean)

# Extract scaled expression data for top genes in epi and endo
epi_clean_scaled <- ScaleData(epi_clean, features = rownames(epi_clean))
endo_clean_scaled <- ScaleData(endo_clean, features = rownames(endo_clean))
fib_clean_scaled <- ScaleData(fib_clean, features = rownames(fib_clean))
epi_scaled <- GetAssayData(epi_clean_scaled, assay = "RNA", slot = "scale.data")[top_epi_genes, epi_cells]
endo_scaled <- GetAssayData(endo_clean_scaled, assay = "RNA", slot = "scale.data")[top_endo_genes, endo_cells]
fib_scaled <- GetAssayData(fib_clean_scaled, assay = "RNA", slot = "scale.data")[top_fib_genes, fib_cells]
# # Create annotation for columns (cells) using your mixed_doublet_status metadata
# epi_annotation <- data.frame(MixedDoubletStatus = epi_clean$mixed_doublet_status)
# rownames(epi_annotation) <- 'epi/endo_cells'
# endo_annotation <- data.frame(MixedDoubletStatus = endo_clean$mixed_doublet_status)
# rownames(endo_annotation) <- 'endo/fib_cells'
# fib_annotation <- data.frame(MixedDoubletStatus = fib_clean$mixed_doublet_status)
# rownames(fib_annotation) <- 'fib/epi_cells'

# Subset mixed doublet cells only
epi_mixed_cells <- mixed_epi_doublets$cell
endo_mixed_cells <- mixed_endo_doublets$cell
fib_mixed_cells <-mixed_fib_doublets$cell
# Subset scaled expression matrix for top marker genes & mixed doublets only
epi_mixed_scaled <- epi_scaled[top_epi_genes, epi_mixed_cells, drop = FALSE]
endo_mixed_scaled <- endo_scaled[top_endo_genes, endo_mixed_cells, drop = FALSE]
fib_mixed_scaled <- fib_scaled[top_fib_genes, fib_mixed_cells, drop = FALSE]
# Create annotation for cells using predicted labels or mixed doublet status
#and predicted label metadata - DID NOT WORK not fully needed?
# epi_mixed_annotation <- data.frame(PredictedLabel = mixed_epi_doublets$predicted_labels)
# rownames(epi_mixed_annotation) <- mixed_epi_doublets$cell
# 
# endo_mixed_annotation <- data.frame(PredictedLabel = mixed_endo_doublets$predicted_label)
# rownames(endo_mixed_annotation) <- mixed_endo_doublets$cell
# 
# fib_mixed_annotation <- data.frame(PredictedLabel = mixed_fib_doublets$predicted_label)
# rownames(fib_mixed_annotation) <- mixed_fib_doublets$cell
# Plot heatmap for Epi mixed doublets
pheatmap(epi_mixed_scaled,
       annotation_col = sample_id,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = FALSE,
         main = "Epithelial/Endothelial Mixed Doublets")

# Plot heatmap for Endo mixed doublets
pheatmap(endo_mixed_scaled,
        annotation_col = endo_mixed_doublets_obj@meta.data$lung_region,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = FALSE,
         main = "Endothelial/Fibroblast Mixed Doublets")

pheatmap(fib_mixed_scaled,
         annotation_col = sample_id,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = FALSE,
         main = "Fibroblast/Epithelial Mixed Doublets")



lung_regions <- unique(c(epi_mixed_doublets_obj$lung_region,
                         endo_mixed_doublets_obj$lung_region,
                         fib_mixed_doublets_obj$lung_region))

sample_ids <- unique(c(epi_mixed_doublets_obj$sample_id,
                       endo_mixed_doublets_obj$sample_id,
                       fib_mixed_doublets_obj$sample_id))

# Number of unique categories
n_lung <- length(lung_regions)
n_samples <- length(sample_ids)
# Define palettes
colours_lung <- setNames(glasbey.colors(n_lung), lung_regions)
colours_sample <- setNames(firstpal.colors(n_samples), sample_ids)
# ----------------------------
# Epi/Endo Mixed Doublets
# ----------------------------
anno_epi <- epi_mixed_doublets_obj@meta.data[colnames(epi_mixed_scaled), c("sample_id", "lung_region"), drop = FALSE]
pheatmap(epi_mixed_scaled,
         annotation_col = anno_epi,
         annotation_colors = list(lung_region = colours_lung, sample_id = colours_sample),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = FALSE,
         main = "Epithelial/Endothelial Mixed Doublets")


# ----------------------------
# Endo/Fib Mixed Doublets
# ----------------------------
anno_endo <- endo_mixed_doublets_obj@meta.data[colnames(endo_mixed_scaled), c("sample_id", "lung_region"), drop = FALSE]
pheatmap(endo_mixed_scaled,
         annotation_col = anno_endo,
         annotation_colors = list(lung_region = colours_lung, sample_id = colours_sample),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = FALSE,
         main = "Endothelial/Fibroblast Mixed Doublets")

# ----------------------------
# Fib/Epi Mixed Doublets
# ----------------------------
anno_fib <- fib_mixed_doublets_obj@meta.data[colnames(fib_mixed_scaled), c("sample_id", "lung_region"), drop = FALSE]
pheatmap(fib_mixed_scaled,
         annotation_col = anno_fib,
         annotation_colors = list(lung_region = colours_lung, sample_id = colours_sample),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = FALSE,
         main = "Fibroblast/Epithelial Mixed Doublets")

# ----------------------------
# Epi/Endo Mixed Doublets
# ----------------------------
anno_epi <- epi_mixed_doublets_obj@meta.data[colnames(epi_mixed_scaled), c("sample_id", "lung_region"), drop = FALSE]
pheatmap(epi_mixed_scaled,
         annotation_col = anno_epi,
         annotation_colors = list(lung_region = colours_lung, sample_id = colours_sample),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = FALSE,
         main = "Epithelial/Endothelial Mixed Doublets",
         fontsize = 7,       # reduces legend/key text size
         fontsize_row = 6,
         fontsize_col = 6)

# ----------------------------
# Endo/Fib Mixed Doublets
# ----------------------------
anno_endo <- endo_mixed_doublets_obj@meta.data[colnames(endo_mixed_scaled), c("sample_id", "lung_region"), drop = FALSE]
pheatmap(endo_mixed_scaled,
         annotation_col = anno_endo,
         annotation_colors = list(lung_region = colours_lung, sample_id = colours_sample),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = FALSE,
         main = "Endothelial/Fibroblast Mixed Doublets",
         fontsize = 7,
         fontsize_row = 6,
         fontsize_col = 6)

# ----------------------------
# Fib/Epi Mixed Doublets
# ----------------------------
anno_fib <- fib_mixed_doublets_obj@meta.data[colnames(fib_mixed_scaled), c("sample_id", "lung_region"), drop = FALSE]
pheatmap(fib_mixed_scaled,
         annotation_col = anno_fib,
         annotation_colors = list(lung_region = colours_lung, sample_id = colours_sample),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = FALSE,
         main = "Fibroblast/Epithelial Mixed Doublets",
         fontsize = 7,
         fontsize_row = 6,
         fontsize_col = 6)


# Save epithelial mixed doublets object
saveRDS(epi_mixed_doublets_obj, file = "epi_mixed_doublets_obj.rds")
# Save endothelial mixed doublets object
saveRDS(endo_mixed_doublets_obj, file = "endo_mixed_doublets_obj.rds")
# Same for fib
saveRDS(fib_mixed_doublets_obj, file = "fib_mixed_doublets_obj.rds")
# Export metadata for epithelial mixed doublets
write.csv(epi_mixed_doublets_obj@meta.data, file = "epi_mixed_doublets_metadata.csv")
# Export metadata for endothelial mixed doublets
write.csv(endo_mixed_doublets_obj@meta.data, file = "endo_mixed_doublets_metadata.csv")
# Same for fib
write.csv(fib_mixed_doublets_obj@meta.data, file = "fib_mixed_doublets_metadata.csv")

# save the Seurat object as RDS
saveRDS(fib_mixed_doublets_obj, file = "fib_mixed_doublets_obj.rds")

#save top markers for mixed doublets
write.csv(epi_mixed_markers, file='epi_mixed_markers.csv')
write.csv(endo_mixed_markers, file="endo_mixed_markers.csv")
write.csv(fib_mixed_markers, file = 'fib_mixed_markers.csv')
#scores saved and metadata under doubllet score folder- 75% percentile
write.csv(epi_mixed_scaled, file='epi_mixed_scaled_markers.csv')
write.csv(endo_mixed_scaled, file='endo_mixed_scaled_markers.csv')
write.csv(fib_mixed_scaled, file='fib_mixed_scaled_markers.csv')



mixed_epi_doublets<-readRDS('epi_mixed_doublets.rds')
mixed_endo_doublets<-readRDS('endo_mixed_doublets.rds')
mixed_fib_doublets<-readRDS('fib_mixed_doublets.rds')
#ensure df
mixed_epi_doublets <- as.data.frame(mixed_epi_doublets)
mixed_endo_doublets <- as.data.frame(mixed_endo_doublets)
mixed_fib_doublets <- as.data.frame(mixed_fib_doublets)
#checke xist
head(mixed_epi_doublets)
head(mixed_endo_doublets)
head(mixed_fib_doublets)
library(dplyr)

# For mixed_epi_doublets
mixed_epi_doublets_df <- mixed_epi_doublets@meta.data %>%
  select(epi_score, endo_score, lung_region, disease)

# For mixed_endo_doublets
mixed_endo_doublets_df <- mixed_endo_doublets@meta.data %>%
  select(endo_score, fib_score, lung_region, disease)

# For mixed_fib_doublets
mixed_fib_doublets_df <- mixed_fib_doublets@meta.data %>%
  select(fib_score, epi_score, lung_region, disease)

library(dplyr)

# summarise so with type
summarize_doublets <- function(obj, label) {
  obj@meta.data %>%
    dplyr::group_by(sample_id, lung_region) %>%
    dplyr::summarise(count = n(), .groups = "drop") %>%
    dplyr::mutate(doublet_type = label) %>%
    dplyr::select(doublet_type, everything())
}
# summarize each doublet type
epi_summary <- summarize_doublets(epi_mixed_doublets_obj, "epi/endo")
endo_summary <- summarize_doublets(endo_mixed_doublets_obj, "endo/fib")
fib_summary <- summarize_doublets(fib_mixed_doublets_obj, "fib/epi")
# Combine into a single table
all_doublets_summary <- dplyr::bind_rows(epi_summary, endo_summary, fib_summary)
all_doublets_summary

write.csv(all_doublets_summary, file='all_doublets_summary.csv')
library(dplyr)
library(tidyr)

# Helper function to summarize by a column
summarize_by_column <- function(obj, label, column) {
  obj@meta.data %>%
    group_by(across(all_of(column))) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(doublet_type = label)
}

# Summarize each doublet type
epi_sample <- summarize_by_column(epi_mixed_doublets_obj, "epi/endo", "sample_id")
endo_sample <- summarize_by_column(endo_mixed_doublets_obj, "endo/fib", "sample_id")
fib_sample <- summarize_by_column(fib_mixed_doublets_obj, "fib/epi", "sample_id")
epi_region <- summarize_by_column(epi_mixed_doublets_obj, "epi/endo", "lung_region")
endo_region <- summarize_by_column(endo_mixed_doublets_obj, "endo/fib", "lung_region")
fib_region <- summarize_by_column(fib_mixed_doublets_obj, "fib/epi", "lung_region")
# Combine all doublets
all_sample <- bind_rows(epi_sample, endo_sample, fib_sample) %>%
  pivot_wider(names_from = doublet_type, values_from = count, values_fill = 0)
all_region <- bind_rows(epi_region, endo_region, fib_region) %>%
  pivot_wider(names_from = doublet_type, values_from = count, values_fill = 0)
# Print combined tables
all_sample
all_region

write.csv(all_sample, file='all_sample.csv')
write.csv(all_region, file='all_region.csv')

library(ggplot2)

# Convert to long format for ggplot
all_sample_long <- all_sample %>%
  pivot_longer(cols = c(`epi/endo`, `endo/fib`, `fib/epi`),
               names_to = "doublet_type",
               values_to = "count")

ggplot(all_sample_long, aes(x = sample_id, y = count, fill = doublet_type)) +
  geom_col(position = "dodge") +
  theme_minimal() +
  labs(title = "Doublet Distribution by Sample",
       x = "Sample ID",
       y = "Count",
       fill = "Doublet Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

all_region_long <- all_region %>%
  pivot_longer(cols = c(`epi/endo`, `endo/fib`, `fib/epi`),
               names_to = "doublet_type",
               values_to = "count")

ggplot(all_region_long, aes(x = lung_region, y = count, fill = doublet_type)) +
  geom_col(position = "dodge") +
  theme_minimal() +
  labs(title = "Doublet Distribution by Lung Region",
       x = "Lung Region",
       y = "Count",
       fill = "Doublet Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#statistical sig? fishers as some small counts
# Sample ID table
sample_table <- as.matrix(all_sample[,-1])  # drop the first column (sample_id)
rownames(sample_table) <- all_sample$sample_id
# Fisher's exact test for sample_id
fisher_sample<-fisher.test(sample_table, simulate.p.value = TRUE, B = 1e5)
fisher_sample
# Lung region table
region_table <- as.matrix(all_region[,-1])
rownames(region_table) <- all_region$lung_region
# Fisher's exact test for lung_region
fisher_region <- fisher.test(region_table)
fisher_region

#look at scores across doublet types - save for downstream
scores_epi_doublets$doublet_type <- "epi/endo"
scores_endo_doublets$doublet_type <- "endo/fib"
scores_fib_doublets$doublet_type <- "fib/epi"
write.csv(scores_epi_doublets, file='epiendo_scores.csv')
write.csv(scores_endo_doublets, file='endofib_scores.csv')
write.csv(scores_fib_doublets, file='fibepi_scores.csv')

#is proportion of mixed cells bbased on thesholds greater than chance? fishers 2x2
epi_table_mixed <- table(
  epi_high = scores_epi_doublets$epi_score > threshold_epi,
  endo_high = scores_epi_doublets$endo_score > threshold_endo
)
epiendo_fishers<-fisher.test(epi_table_mixed)

endo_table_mixed <- table(
  endo_high = scores_endo_doublets$endo_score > threshold_endo_2,
  fib_high = scores_endo_doublets$fib_score > threshold_fib
)
endofib_fishers<-fisher.test(endo_table_mixed)

fib_table_mixed <- table(
  fib_high = scores_fib_doublets$fib_score > threshold_fib_2,
  epi_high = scores_fib_doublets$epi_score > threshold_epi_2
)
fibepi_fishers<-fisher.test(fib_table_mixed)


# extract fisher.test results into a data frame
extract_fisher <- function(test, name) {
  data.frame(
    comparison = name,
    p_value = test$p.value,
    odds_ratio = unname(test$estimate),
    conf_low = test$conf.int[1],
    conf_high = test$conf.int[2],
    method = test$method
  )
}

# collect results
epi_endo_res  <- extract_fisher(epiendo_fishers, "epi_vs_endo")
endo_fib_res  <- extract_fisher(endofib_fishers, "endo_vs_fib")
fib_epi_res   <- extract_fisher(fibepi_fishers, "fib_vs_epi")
# combine into one table
all_results <- rbind(epi_endo_res, endo_fib_res, fib_epi_res)
# export all together
write.csv(all_results, "doublet_fishers_summary.csv", row.names = FALSE)
# abd seperate
write.csv(epi_endo_res, "fisher_epi_vs_endo.csv", row.names = FALSE)
write.csv(endo_fib_res, "fisher_endo_vs_fib.csv", row.names = FALSE)
write.csv(fib_epi_res, "fisher_fib_vs_epi.csv", row.names = FALSE)

#repeat - Zscore to standardise scales of scores between all 3
zscore <- function(x) {
  (x - mean(x)) / sd(x)
}

#epi vs endo - epi doubelts
epi_exp_all <- GetAssayData(epi_doublets, assay = "RNA", slot = "data")
markers_epi_endo <- unique(c(top_markers_epi, top_markers_endo))
epi_exp_subset <- epi_exp_all[intersect(rownames(epi_exp_all), markers_epi_endo), , drop = FALSE]

epi_score_in_epi_doublets <- Matrix::colSums(epi_exp_subset[intersect(top_markers_epi, rownames(epi_exp_subset)), , drop = FALSE])
endo_score_in_epi_doublets <- Matrix::colSums(epi_exp_subset[intersect(top_markers_endo, rownames(epi_exp_subset)), , drop = FALSE])

# z-score
epi_score_in_epi_doublets <- zscore(epi_score_in_epi_doublets)
endo_score_in_epi_doublets <- zscore(endo_score_in_epi_doublets)

scores_epi_doublets <- data.frame(
  cell = colnames(epi_exp_subset),
  epi_score = epi_score_in_epi_doublets,
  endo_score = endo_score_in_epi_doublets
)

threshold_epi <- quantile(scores_epi_doublets$epi_score, 0.75)
threshold_endo <- quantile(scores_epi_doublets$endo_score, 0.75)

epi_table_mixed <- table(
  epi_high = scores_epi_doublets$epi_score > threshold_epi,
  endo_high = scores_epi_doublets$endo_score > threshold_endo
)
epiendo_fishers <- fisher.test(epi_table_mixed)
#endo vs fib - endo doublets
endo_exp_all <- GetAssayData(endo_doublets, assay = "RNA", slot = "data")
markers_endo_fib <- unique(c(top_markers_endo, top_markers_fib))
endo_exp_subset <- endo_exp_all[intersect(rownames(endo_exp_all), markers_endo_fib), , drop = FALSE]

endo_score_in_endo_doublets <- Matrix::colSums(endo_exp_subset[intersect(top_markers_endo, rownames(endo_exp_subset)), , drop = FALSE])
fib_score_in_endo_doublets <- Matrix::colSums(endo_exp_subset[intersect(top_markers_fib, rownames(endo_exp_subset)), , drop = FALSE])

# z-score
endo_score_in_endo_doublets <- zscore(endo_score_in_endo_doublets)
fib_score_in_endo_doublets <- zscore(fib_score_in_endo_doublets)

scores_endo_doublets <- data.frame(
  cell = colnames(endo_exp_subset),
  endo_score = endo_score_in_endo_doublets,
  fib_score = fib_score_in_endo_doublets
)

threshold_endo_2 <- quantile(scores_endo_doublets$endo_score, 0.75)
threshold_fib <- quantile(scores_endo_doublets$fib_score, 0.75)

endo_table_mixed <- table(
  endo_high = scores_endo_doublets$endo_score > threshold_endo_2,
  fib_high = scores_endo_doublets$fib_score > threshold_fib
)
endofib_fishers <- fisher.test(endo_table_mixed)

#epi vs fib - fib doublets
fib_exp_all <- GetAssayData(fib_doublets, assay = "RNA", slot = "data")
markers_fib_epi <- unique(c(top_markers_fib, top_markers_epi))
fib_exp_subset <- fib_exp_all[intersect(rownames(fib_exp_all), markers_fib_epi), , drop = FALSE]

fib_score_in_fib_doublets <- Matrix::colSums(fib_exp_subset[intersect(top_markers_fib, rownames(fib_exp_subset)), , drop = FALSE])
epi_score_in_fib_doublets <- Matrix::colSums(fib_exp_subset[intersect(top_markers_epi, rownames(fib_exp_subset)), , drop = FALSE])

# z-score
fib_score_in_fib_doublets <- zscore(fib_score_in_fib_doublets)
epi_score_in_fib_doublets <- zscore(epi_score_in_fib_doublets)

scores_fib_doublets <- data.frame(
  cell = colnames(fib_exp_subset),
  fib_score = fib_score_in_fib_doublets,
  epi_score = epi_score_in_fib_doublets
)

threshold_fib_2 <- quantile(scores_fib_doublets$fib_score, 0.75)
threshold_epi_2 <- quantile(scores_fib_doublets$epi_score, 0.75)

fib_table_mixed <- table(
  fib_high = scores_fib_doublets$fib_score > threshold_fib_2,
  epi_high = scores_fib_doublets$epi_score > threshold_epi_2
)
fibepi_fishers <- fisher.test(fib_table_mixed)

#extract results
extract_fisher <- function(test, name) {
  data.frame(
    comparison = name,
    p_value = test$p.value,
    odds_ratio = unname(test$estimate),
    conf_low = test$conf.int[1],
    conf_high = test$conf.int[2],
    method = test$method
  )
}

#save - export csv
epi_endo_res  <- extract_fisher(epiendo_fishers, "epi_vs_endo")
endo_fib_res  <- extract_fisher(endofib_fishers, "endo_vs_fib")
fib_epi_res   <- extract_fisher(fibepi_fishers, "fib_vs_epi")
all_results <- rbind(epi_endo_res, endo_fib_res, fib_epi_res)

write.csv(all_results, "doublet_fishers_summary.csv", row.names = FALSE)
write.csv(epi_endo_res, "fisher_epi_vs_endo.csv", row.names = FALSE)
write.csv(endo_fib_res, "fisher_endo_vs_fib.csv", row.names = FALSE)
write.csv(fib_epi_res, "fisher_fib_vs_epi.csv", row.names = FALSE)

#forest plot ORs
ggplot(all_results, aes(x = comparison, y = odds_ratio)) +
  geom_point(size = 3, color = "darkblue") +
  geom_errorbar(aes(ymin = conf_low, ymax = conf_high), width = 0.2, color = "darkblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_y_log10() +
  coord_flip() +
  labs(
    title = "Odds Ratios of Mixed Doublet Enrichment",
    y = "Odds Ratio (log scale)",
    x = "Comparison"
  ) +
  theme_minimal(base_size = 14)

#ggsave("forest_plot_odds_ratios.png", forest_plot, width = 6, height = 4, dpi = 300)


#p val sig plots
all_results <- all_results %>%
  mutate(log10_p = -log10(p_value))

ggplot(all_results, aes(x = comparison, y = log10_p, fill = comparison)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(
    title = "Significance of Mixed Doublet Enrichment",
    y = "-log10(p-value)",
    x = "Comparison"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

#save coepression scores!
write.csv(scores_epi_doublets, file='all_epi_coexpr_scores.csv')
write.csv(scores_endo_doublets, file='all_endo_coexpr_scores.csv')
write.csv(scores_fib_doublets, file='all_fib_coexpr_scores.csv')
#log reg model - mixed doublet 1 other 0
epi_clean$doublet_status <- ifelse(
  epi_clean$epi_score > threshold_epi & epi_clean$endo_score > threshold_endo,
  1, 0
)

endo_clean$doublet_status <- ifelse(
  endo_clean$endo_score > threshold_endo_2 & endo_clean$fib_score > threshold_fib,
  1, 0
)

fib_clean$doublet_status <- ifelse(
  fib_clean$fib_score > threshold_fib_2 & fib_clean$epi_score > threshold_epi_2,
  1, 0
)

#add coexpr scores to metadata
epi_clean$epi_score <- scores_epi_doublets[Cells(epi_clean), "epi_score"]
epi_clean$endo_score <- scores_epi_doublets[Cells(epi_clean), "endo_score"]
endo_clean$endo_score <- scores_endo_doublets[Cells(endo_clean), "endo_score"]
endo_clean$fib_score <- scores_endo_doublets[Cells(endo_clean), "fib_score"]
fib_clean$fib_score <- scores_fib_doublets[Cells(fib_clean), "fib_score"]
fib_clean$epi_score <- scores_fib_doublets[Cells(fib_clean), "epi_score"]

epi_scores_df <- data.frame(
  cell = names(epi_scores),
  epi_score = as.numeric(epi_scores),
  row.names = NULL
)
head(epi_scores_df)
endo_scores_df <- data.frame(
  cell = names(endo_scores),
  endo_score = as.numeric(endo_scores),
  row.names = NULL
)
head(endo_scores_df)
fib_scores_df <- data.frame(
  cell = names(fib_scores),
  fib_score = as.numeric(fib_scores),
  row.names = NULL
)
head(fib_scores_df)

all(Cells(epi_doublets) %in% names(epi_scores))
all(Cells(endo_doublets) %in% names(endo_scores))
all(Cells(fib_doublets) %in% names(fib_scores))

library(pheatmap)
library(Seurat)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

# ----------------------------
# extract clusters from heatmap
# ----------------------------
extract_clusters <- function(expr_matrix, k_rows = 4, k_cols = 3){
  ph <- pheatmap(expr_matrix, silent = TRUE)
  
  row_clusters <- cutree(ph$tree_row, k = k_rows)
  col_clusters <- cutree(ph$tree_col, k = k_cols)
  
  list(
    row_clusters = row_clusters,
    col_clusters = col_clusters
  )
}

# ----------------------------
# add module scores for each gene cluster
# ----------------------------
add_module_scores <- function(seurat_obj, expr_matrix, row_clusters){
  for(clust in unique(row_clusters)){
    genes <- names(row_clusters[row_clusters == clust])
    mod_name <- paste0("Module_", clust)
    seurat_obj <- AddModuleScore(
      seurat_obj,
      features = list(genes),
      name = mod_name
    )
  }
  return(seurat_obj)
}

# ----------------------------
#plot module scores by lung_region
# ----------------------------
plot_module_scores <- function(seurat_obj, module_prefix = "Module_"){
  module_cols <- grep(module_prefix, colnames(seurat_obj@meta.data), value = TRUE)
  plots <- lapply(module_cols, function(mod){
    VlnPlot(seurat_obj, features = mod, group.by = "lung_region") +
      ggtitle(mod) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  return(plots)
}

# ----------------------------
# 1. Process each mixed doublet type
# ----------------------------
# Example: epi/endo
epi_clusters <- extract_clusters(epi_mixed_scaled, k_rows = 4, k_cols = 3)
epi_mixed_doublets_obj <- add_module_scores(epi_mixed_doublets_obj, epi_mixed_scaled, epi_clusters$row_clusters)
epi_plots <- plot_module_scores(epi_mixed_doublets_obj)

# Example: endo/fib
endo_clusters <- extract_clusters(endo_mixed_scaled, k_rows = 4, k_cols = 3)
endo_mixed_doublets_obj <- add_module_scores(endo_mixed_doublets_obj, endo_mixed_scaled, endo_clusters$row_clusters)
endo_plots <- plot_module_scores(endo_mixed_doublets_obj)

# Example: fib/epi
fib_clusters <- extract_clusters(fib_mixed_scaled, k_rows = 4, k_cols = 3)
fib_mixed_doublets_obj <- add_module_scores(fib_mixed_doublets_obj, fib_mixed_scaled, fib_clusters$row_clusters)
fib_plots <- plot_module_scores(fib_mixed_doublets_obj)

# ----------------------------
# 2. Compare modules across doublet types
# ----------------------------
# Compute average expression per module per doublet type
compute_avg_module <- function(seurat_obj, module_prefix = "Module_"){
  module_cols <- grep(module_prefix, colnames(seurat_obj@meta.data), value = TRUE)
  avg_module <- sapply(module_cols, function(mod) mean(seurat_obj@meta.data[[mod]]))
  return(avg_module)
}

avg_epi <- compute_avg_module(epi_mixed_doublets_obj)
avg_endo <- compute_avg_module(endo_mixed_doublets_obj)
avg_fib <- compute_avg_module(fib_mixed_doublets_obj)

avg_modules_df <- data.frame(
  module = names(avg_epi),
  epi = avg_epi,
  endo = avg_endo,
  fib = avg_fib
)

#correlation heatmap between doublet types
library(corrplot)
corr_matrix <- cor(avg_modules_df[,2:4])
# Rename rows and columns
rownames(corr_matrix) <- colnames(corr_matrix) <- c("Epi/Endo", "Endo/Fib", "Fib/Epi")
corrplot(corr_matrix, method = "color", addCoef.col = "black")
title(main = "Correlation of Module Expression Across Mixed Doublet Types")



# ----------------------------
# 3. annotate modules with GO enrichment
# ----------------------------
enrich_modules <- function(expr_matrix, row_clusters){
  module_results <- list()
  for(clust in unique(row_clusters)){
    genes <- names(row_clusters[row_clusters == clust])
    entrez_ids <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
    go_res <- enrichGO(
      gene = entrez_ids$ENTREZID,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pvalueCutoff = 0.05
    )
    module_results[[paste0("Module_", clust)]] <- go_res
  }
  return(module_results)
}
epi_go <- enrich_modules(epi_mixed_scaled, epi_clusters$row_clusters)
# # Combine into a single dataframe
# combined_go <- do.call(rbind, lapply(names(epi_go), function(mod){
#   df <- as.data.frame(epi_go[[mod]])
#   df$Module <- mod
#   df
# }))
# # Take top 5 GO terms per module
# library(dplyr)
# top_go <- combined_go %>%
#   group_by(Module) %>%
#   slice_max(order_by = -log10(p.adjust), n = 5)
# # Dotplot using ggplot2
# ggplot(top_go, aes(x = Module, y = Description, size = Count, color = -log10(p.adjust))) +
#   geom_point() +
#   scale_color_viridis_c() +
#   theme_minimal() +
#   labs(x = "Module", y = "GO Term", color = "-log10(adj p)", size = "Gene count") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))



