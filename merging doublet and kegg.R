# Load libraries
library(data.table)  # fast fread
library(dplyr)
library(tibble)      # for rownames_to_column
library(uwot)

# -----------------------------
# 1. Load scaled marker CSVs
# -----------------------------
epi_df <- fread("epi_mixed_scaled_markers.csv")
endo_df <- fread("endo_mixed_scaled_markers.csv")
fib_df <- fread("fib_mixed_scaled_markers.csv")

# Make sure first column is "gene"
colnames(epi_df)[1] <- "gene"
colnames(endo_df)[1] <- "gene"
colnames(fib_df)[1] <- "gene"

# -----------------------------
# 2. Extract gene names
# -----------------------------
epi_genes <- epi_df$gene
endo_genes <- endo_df$gene
fib_genes <- fib_df$gene

# -----------------------------
# 3. Transpose so cells are rows
# -----------------------------
epi_mat <- as.data.frame(t(epi_df[,-1]))
colnames(epi_mat) <- epi_genes
rownames(epi_mat) <- colnames(epi_df)[-1]

endo_mat <- as.data.frame(t(endo_df[,-1]))
colnames(endo_mat) <- endo_genes
rownames(endo_mat) <- colnames(endo_df)[-1]

fib_mat <- as.data.frame(t(fib_df[,-1]))
colnames(fib_mat) <- fib_genes
rownames(fib_mat) <- colnames(fib_df)[-1]

# -----------------------------
# 4. Load metadata
# -----------------------------
epi_meta <- read.csv("epi_mixed_doublets_metadata.csv", check.names = FALSE)
endo_meta <- read.csv("endo_mixed_doublets_metadata.csv", check.names = FALSE)
fib_meta <- read.csv("fib_mixed_doublets_metadata.csv", check.names = FALSE)

# Rename first column in metadata to "cell_id" if blank
colnames(epi_meta)[1] <- "cell_id"
colnames(endo_meta)[1] <- "cell_id"
colnames(fib_meta)[1] <- "cell_id"
#***stop here
# -----------------------------
# 5. Add cell_id as a column to matrices
# -----------------------------
epi_mat <- epi_mat %>% rownames_to_column(var = "cell_id")
endo_mat <- endo_mat %>% rownames_to_column(var = "cell_id")
fib_mat <- fib_mat %>% rownames_to_column(var = "cell_id")

# -----------------------------
# 6. Merge metadata with markers
# -----------------------------
epi_full <- epi_meta %>% left_join(epi_mat, by = "cell_id")
endo_full <- endo_meta %>% left_join(endo_mat, by = "cell_id")
fib_full <- fib_meta %>% left_join(fib_mat, by = "cell_id")

# -----------------------------
# 7. Add a cell type column
# -----------------------------
epi_full$cell_type <- "epi/endo"
endo_full$cell_type <- "endo/fib"
fib_full$cell_type <- "fib/epi"

# -----------------------------
# 8. Combine all three into one data frame
# -----------------------------
all_cells <- bind_rows(epi_full, endo_full, fib_full)

# -----------------------------
# 9. Optional: check
# -----------------------------
dim(all_cells)
head(all_cells[,1:10])
write.csv(all_cells, file="all_doublets_with_meta.csv")
write.csv(epi_full, file="epi_with_meta.csv")
write.csv(endo_full, file="endo_with_meta.csv")
write.csv(fib_full, file="fib_with_meta.csv")
write.csv(epi_df, file="epi_df.csv")
write.csv(endo_df, file="endo_df.csv")
write.csv(fib_df, file="fib_df.csv")

library(dplyr)
library(tibble)

# Save gene names
epi_genes <- epi_df$gene
# Remove gene column and transpose
epi_mat <- t(epi_df[,-1])
# Convert to data.frame
epi_mat <- as.data.frame(epi_mat)
# Set column names to gene names
colnames(epi_mat) <- epi_genes
# Set rownames to cell barcodes
rownames(epi_mat) <- colnames(epi_df)[-1]
# Optional: add cell_id column for merging with metadata
epi_mat <- epi_mat %>% rownames_to_column(var = "cell_id")

# Save gene names
endo_genes <- endo_df$gene
# Remove gene column and transpose
endo_mat <- t(endo_df[,-1])
# Convert to data.frame
endo_mat <- as.data.frame(endo_mat)
# Set column names to gene names
colnames(endo_mat) <- endo_genes
# Set rownames to cell barcodes
rownames(endo_mat) <- colnames(endo_df)[-1]
# Optional: add cell_id column for merging with metadata
endo_mat <- endo_mat %>% rownames_to_column(var = "cell_id")

epi_clean<-readRDS('epi_clean.rds')
fib_clean<-readRDS('fib_clean.rds')
endo_clean<-readRDS('endo_clean.rds')
# Run UMAP on PCA embeddings
epi_clean <- RunUMAP(epi_clean, dims = 1:20)
# Plot UMAP
DimPlot(epi_clean, reduction = "umap", group.by = "mixed_doublet_status",
        cols = c("MixedDoublet" = "#F8766D", "Other" = "#00BFC4")) + ggtitle ('Epi Subset UMAP - Mixed Doublets')
DimPlot(endo_clean, reduction = "umap", group.by = "mixed_doublet_status",
        cols = c("MixedDoublet" = "#F8766D", "Other" = "#00BFC4")) + ggtitle ('Endo Subset UMAP - Mixed Doublets')
DimPlot(fib_clean, reduction = "umap", group.by = "mixed_doublet_status",
        cols = c("MixedDoublet" = "#F8766D", "Other" = "#00BFC4")) + ggtitle ('Fib Subset UMAP - Mixed Doublets')


DimPlot(epi_mixed_doublets_obj, reduction = "umap", group.by='sample_id') +ggtitle('Epi Mixed Doublets by Sample ID')
DimPlot(endo_mixed_doublets_obj, reduction = "umap", group.by = 'sample_id') +ggtitle('Endo Mixed Doublets by Sample ID')
DimPlot(fib_mixed_doublets_obj, reduction='umap',group.by='sample_id')+ggtitle('Fib Mixed Doublets by Sample ID')


# Save gene names
fib_genes <- fib_df$gene
# Remove gene column and transpose
fib_mat <- t(fib_df[,-1])
# Convert to data.frame
fib_mat <- as.data.frame(fib_mat)
# Set column names to gene names
colnames(fib_mat) <- fib_genes
# Set rownames to cell barcodes
rownames(fib_mat) <- colnames(fib_df)[-1]
# Optional: add cell_id column for merging with metadata
fib_mat <- fib_mat %>% rownames_to_column(var = "cell_id")

write.csv(epi_mat, file="epi_mat.csv")
write.csv(endo_mat, file="endo_mat.csv")
write.csv(fib_mat, file="fib_mat.csv")

library(dplyr)

epi_mat<-read.csv('epi_mat.csv')
endo_mat<-read.csv('endo_mat.csv')
fib_mat<-read.csv('fib_mat.csv')

(epi_mat)

epi_with_meta<-read.csv('epi_with_meta.csv')
endo_with_meta<-read.csv('endo_with_meta.csv')
fib_with_meta<-read.csv('fib_with_meta.csv')
#new metadata col doublet_type
# epi_mixed
markers_epi <- read.csv("epi_mixed_scaled_markers.csv")
metadata_epi <- read.csv("epi_mixed_doublets_metadata.csv")
metadata_epi$doublet_type <- "epi/endo"
write.csv(metadata_epi, "epi_mixed_doublets_metadata_with_type.csv", row.names = FALSE)
# fib_mixed
metadata_fib <- read.csv("fib_mixed_doublets_metadata.csv")
metadata_fib$doublet_type <- "fib/epi"
write.csv(metadata_fib, "fib_mixed_doublets_metadata_with_type.csv", row.names = FALSE)
# endo_mixed
metadata_endo <- read.csv("endo_mixed_doublets_metadata.csv")
metadata_endo$doublet_type <- "endo/fib"
write.csv(metadata_endo, "endo_mixed_doublets_metadata_with_type.csv", row.names = FALSE)

#must have epi and endo _mixed_doublets_obj loaded in!
avg_epi <- AverageExpression(epi_mixed_doublets_obj, group.by = "lung_region")$RNA
avg_endo <- AverageExpression(endo_mixed_doublets_obj, group.by = "lung_region")$RNA


# =============================
# Load libraries
# =============================
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# =============================
# Function: KEGG enrichment from avg expression
# =============================
run_kegg <- function(avg_expr, lung_region, top_n = 200){
  
  # 1. Select top expressed genes in this region
  top_genes <- rownames(avg_expr)[order(avg_expr[, lung_region], decreasing = TRUE)[1:top_n]]
  
  # 2. Convert gene symbols to Entrez IDs
  entrez_ids <- bitr(top_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  
  # 3. KEGG enrichment
  kegg_res <- enrichKEGG(gene = entrez_ids$ENTREZID,
                         organism = "hsa",
                         pvalueCutoff = 0.05)
  
  # 4. Save gene list
  write.table(top_genes, paste0(lung_region, "_top_genes.txt"), 
              quote=F, row.names=F, col.names=F)
  
  # 5. Return results
  return(kegg_res)
}

# =============================
# Example: Epithelial doublets
# =============================
# avg_epi should be a matrix/data.frame: rows = genes, columns = regions (control, upper, lower)
epi_upper_kegg <- run_kegg(avg_epi, lung_region = "Upper", top_n = 200)
epi_lower_kegg <- run_kegg(avg_epi, lung_region = "Lower", top_n = 200)
epi_control_kegg <- run_kegg(avg_epi, lung_region = "Control", top_n = 200)

# =============================
# Example: Endothelial doublets
# =============================
# avg_endo should be a matrix/data.frame: rows = genes, columns = regions
endo_upper_kegg <- run_kegg(avg_endo, lung_region = "Upper", top_n = 200)
endo_lower_kegg <- run_kegg(avg_endo, lung_region = "Lower", top_n = 200)
endo_control_kegg <- run_kegg(avg_endo, lung_region = "Control", top_n = 200)

# =============================
# Visualize KEGG results
# =============================
# Dotplot for epithelial upper
dotplot(epi_upper_kegg, showCategory=12) + ggtitle("Epi/Endo Upper Lobe KEGG")
# Dotplot for endothelial upper
dotplot(endo_upper_kegg, showCategory=12) + ggtitle("Endo/Fib Upper Lobe KEGG")
# Dotplot for epithelial upper
dotplot(epi_lower_kegg, showCategory=12) + ggtitle("Epi/Endo Lower Lobe KEGG")
# Dotplot for endothelial upper
dotplot(endo_lower_kegg, showCategory=12) + ggtitle("Endo/Fib Lower Lobe KEGG")
# Dotplot for epithelial upper
dotplot(epi_control_kegg, showCategory=12) + ggtitle("Epi/Endo Control KEGG")
# Dotplot for endothelial upper
dotplot(endo_control_kegg, showCategory=12) + ggtitle("Endo/Fib Control KEGG")

# Repeat dotplot() for other conditions as needed

# =============================
# Load libraries
# =============================
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# =============================
# Load libraries
# =============================
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# =============================
# Function: KEGG enrichment from a Seurat object
# =============================
run_kegg_all_cells <- function(seurat_obj, top_n = 200){
  
  # 1. Average expression across all cells
  avg_expr <- AverageExpression(seurat_obj)$RNA
  
  # 2. Select top expressed genes
  top_genes <- rownames(avg_expr)[order(avg_expr[,1], decreasing = TRUE)[1:top_n]]
  
  # 3. Convert gene symbols to Entrez IDs
  entrez_ids <- bitr(top_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  
  # 4. KEGG enrichment
  kegg_res <- enrichKEGG(gene = entrez_ids$ENTREZID,
                         organism = "hsa",
                         pvalueCutoff = 0.05)
  
  # 5. Save top genes
  write.table(top_genes, "top_genes.txt", quote=F, row.names=F, col.names=F)
  
  # 6. Return enrichment results
  return(kegg_res)
}

# =============================
# Run for all epithelial doublets
# =============================
epi_kegg <- run_kegg_all_cells(epi_mixed_doublets_obj, top_n = 200)

# =============================
# Run for all endothelial doublets
# =============================
endo_kegg <- run_kegg_all_cells(endo_mixed_doublets_obj, top_n = 200)

# =============================
# Visualize KEGG results
# =============================
dotplot(epi_kegg, showCategory=12) + ggtitle("Epi/Endo Doublets KEGG")
dotplot(endo_kegg, showCategory=12) + ggtitle("Endo/Fib Doublets KEGG")

# =============================
# Add doublet_type to metadata
# =============================

# Epi/Endo doublets
epi_mixed_doublets_obj$doublet_type <- "epi/endo"
# Endo/Fib doublets
endo_mixed_doublets_obj$doublet_type <- "endo/fib"
# Fib/Epi doublets
fib_mixed_doublets_obj$doublet_type <- "fib/epi"

# =============================
# Check results
# =============================
cat("Epi:\n")
print(table(epi_mixed_doublets_obj$doublet_type))

cat("\nEndo:\n")
print(table(endo_mixed_doublets_obj$doublet_type))

cat("\nFib:\n")
print(table(fib_mixed_doublets_obj$doublet_type))



#dotplots - top genes
# Function: make dotplot of top expressed genes
plot_top_genes <- function(seurat_obj, group.by = "doublet_type", top_n = 20){
  
  # 1. Get average expression
  avg_expr <- AverageExpression(seurat_obj, group.by = group.by)$RNA
  
  # 2. Select top N genes (by overall expression across groups)
  gene_means <- rowMeans(avg_expr)
  top_genes <- names(sort(gene_means, decreasing = TRUE))[1:top_n]
  
  # 3. Make dot plot
  p <- DotPlot(seurat_obj, features = top_genes, group.by = group.by) +
    RotatedAxis() +
    ggtitle(paste("Top", top_n, "genes by", group.by))
  
  return(p)
}

# =============================
# Example for each doublet set
# =============================
# Epi/Endo doublets
plot_top_genes(epi_mixed_doublets_obj, group.by = "doublet_type", top_n = 20)
# Endo/Fib doublets
plot_top_genes(endo_mixed_doublets_obj, group.by = "doublet_type", top_n = 20)
# Fib/Epi doublets
plot_top_genes(fib_mixed_doublets_obj, group.by = "doublet_type", top_n = 20)


# DotPlot(epi_mixed_doublets_obj, features = top_genes, group.by = "doublet_type") + 
#   RotatedAxis() + 
#   ggtitle("Top Expressed Genes - Epi/Endo Doublets") +
#   theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))
# 
# DotPlot(endo_mixed_doublets_obj, features = top_genes, group.by = "doublet_type") + 
#   RotatedAxis() + 
#   ggtitle("Top Expressed Genes - Endo/Fib Doublets") +
#   theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))
# 
# DotPlot(fib_mixed_doublets_obj, features = top_genes, group.by = "doublet_type") + 
#   RotatedAxis() + 
#   ggtitle("Top Expressed Genes - Fib/Epi Doublets") +
#   theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))


# 1. Add doublet_type to metadata for each object
epi_mixed_doublets_obj$doublet_type <- "epi/endo"
endo_mixed_doublets_obj$doublet_type <- "endo/fib"
fib_mixed_doublets_obj$doublet_type <- "fib/epi"

# 2. Merge all three into one object
all_doublets <- merge(
  epi_mixed_doublets_obj, 
  y = list(endo_mixed_doublets_obj, fib_mixed_doublets_obj),
  add.cell.ids = c("EpiEndo", "EndoFib", "FibEpi")
)

# 3. Normalize & scale if not already
all_doublets <- NormalizeData(all_doublets)
all_doublets <- FindVariableFeatures(all_doublets)
all_doublets <- ScaleData(all_doublets)

# 4. Define a plotting function for top expressed genes per group
plot_top_genes <- function(seurat_obj, group.by, top_n = 50) {
  avg_expr <- AverageExpression(seurat_obj, group.by = group.by)$RNA
  top_genes <- apply(avg_expr, 2, function(x) {
    head(names(sort(x, decreasing = TRUE)), top_n)
  })
  top_genes <- unique(unlist(top_genes))
  DoHeatmap(seurat_obj, features = top_genes, group.by = group.by) + 
    ggtitle(paste("Top", top_n, "genes by", group.by))
}

# 5. Plot all together
plot_top_genes(all_doublets, group.by = "doublet_type", top_n = 50)

# --- Function to get top genes per group ---
get_top_genes <- function(seurat_obj, group.by, top_n = 30) {
  avg_expr <- AverageExpression(seurat_obj, group.by = group.by)$RNA
  top_genes <- apply(avg_expr, 2, function(x) {
    head(names(sort(x, decreasing = TRUE)), top_n)
  })
  # flatten + keep unique genes
  top_genes <- unique(as.vector(top_genes))
  return(top_genes)
}

# --- Get top genes across doublet types ---
top_genes <- get_top_genes(all_doublets, group.by = "doublet_type", top_n = 50)
# --- Dotplot ---
DotPlot(all_doublets, features = top_genes, group.by = "doublet_type") + 
  RotatedAxis() + 
  ggtitle("Top Expressed Genes Across Doublet Types") +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))


# GO for epi/endo and endo/fib
run_go_all_cells <- function(seurat_obj, top_n = 200, ont = "BP") {
  
  # 1. Average expression
  avg_expr <- AverageExpression(seurat_obj)$RNA
  
  # 2. Top expressed genes
  top_genes <- rownames(avg_expr)[order(avg_expr[,1], decreasing = TRUE)[1:top_n]]
  
  # 3. Convert to Entrez IDs
  entrez_ids <- bitr(top_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  
  # 4. GO enrichment
  go_res <- enrichGO(
    gene = entrez_ids$ENTREZID,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = ont,            # BP, MF, or CC
    pvalueCutoff = 0.05,
    readable = TRUE       # converts back to gene symbols
  )
  
  # 5. Save top genes
  write.table(top_genes, paste0("top_genes_", ont, ".txt"), quote=F, row.names=F, col.names=F)
  
  return(go_res)
}
# Epi/Endo mixed doublets
epi_go <- run_go_all_cells(epi_mixed_doublets_obj, top_n = 200, ont = "BP")
dotplot(epi_go, showCategory = 15) + ggtitle("Epi/Endo Doublets GO-BP")
# Endo/Fib mixed doublets
endo_go <- run_go_all_cells(endo_mixed_doublets_obj, top_n = 200, ont = "BP")
dotplot(endo_go, showCategory = 15) + ggtitle("Endo/Fib Doublets GO-BP")

# Epi/Endo mixed doublets
epi_go <- run_go_all_cells(epi_mixed_doublets_obj, top_n = 200, ont = "MF")
dotplot(epi_go, showCategory = 12) + ggtitle("Epi/Endo Doublets GO-MF")
# Endo/Fib mixed doublets
endo_go <- run_go_all_cells(endo_mixed_doublets_obj, top_n = 200, ont = "MF")
dotplot(endo_go, showCategory = 12) + ggtitle("Endo/Fib Doublets GO-MF")

# Epi/Endo mixed doublets
epi_go <- run_go_all_cells(epi_mixed_doublets_obj, top_n = 200, ont = "CC")
dotplot(epi_go, showCategory = 12) + ggtitle("Epi/Endo Doublets GO-CC")
# Endo/Fib mixed doublets
endo_go <- run_go_all_cells(endo_mixed_doublets_obj, top_n = 200, ont = "CC")
dotplot(endo_go, showCategory = 12) + ggtitle("Endo/Fib Doublets GO-CC")


Idents(epi_mixed_doublets_obj) <- epi_mixed_doublets_obj$doublet_type  
Idents(endo_mixed_doublets_obj) <- endo_mixed_doublets_obj$doublet_type     
Idents(fib_mixed_doublets_obj) <- fib_mixed_doublets_obj$doublet_type  
# top X variable features for each obj
epi_top15 <- head(VariableFeatures(epi_mixed_doublets_obj), 15)
endo_top15 <- head(VariableFeatures(endo_mixed_doublets_obj), 15)
fib_top15 <- head(VariableFeatures(fib_mixed_doublets_obj), 15)

# violin plots
p1 <- VlnPlot(epi_mixed_doublets_obj, features = epi_top15, ncol = 5, cols='#0CB702') + 
  #ggtitle("Epi/Endo - Top 15 Variable Genes") +
  theme(plot.title = element_text(hjust = 0.5, vjust = 1, size = 14, face = "bold"))

p2 <- VlnPlot(endo_mixed_doublets_obj, features = endo_top15, ncol = 5) + 
  #ggtitle("Endo/Fib - Top 15 Variable Genes") +
  theme(plot.title = element_text(hjust = 0.5, vjust = 1, size = 14, face = "bold"))

p3 <- VlnPlot(fib_mixed_doublets_obj, features = fib_top15, ncol = 5, cols='#00A9FF') + 
  #ggtitle("Fib/Epi - Top 15 Variable Genes") +
  theme(plot.title = element_text(hjust = 0.5, vjust = 1, size = 14, face = "bold"))
p1
p2
p3
# Combine plots
combined_plot <- p1 / p2 / p3
print(combined_plot)

epi_top50 <- head(VariableFeatures(epi_mixed_doublets_obj), 50)
endo_top50 <- head(VariableFeatures(endo_mixed_doublets_obj), 50)
fib_top50 <- head(VariableFeatures(fib_mixed_doublets_obj), 50)


write.csv(epi_top50, file='epitop.csv')
write.csv(endo_top50, file='endotop.csv')
write.csv(fib_top50, file='fibtop.csv')

