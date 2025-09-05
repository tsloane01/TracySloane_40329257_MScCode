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
library(tidyr)
install.packages('devtools')
devtools::install_github('immunogenomics/presto')



#read in ssc-ild and ctrl seurat objs
upper_combined <- readRDS("upper_combined.rds")
lower_combined <- readRDS("lower_combined.rds")
ctrl_combined  <- readRDS("ctrl_combined.rds")

#subset fibroblasts and epithelial cells from all datasets
cellssubset2 <- c("Epithelial_cells", "Fibroblasts")
upper_subset2 <- subset(upper_combined, subset = SingleR.labels %in% cellssubset2)
lower_subset2 <- subset(lower_combined, subset = SingleR.labels %in% cellssubset2)
ctrl_subset2  <- subset(ctrl_combined,  subset = SingleR.labels %in% cellssubset2)
#add condition label
upper_subset2$condition <- "Upper"
lower_subset2$condition <- "Lower"
ctrl_subset2$condition  <- "Control"

#subset epithelial cells only
epi_upper <- subset(upper_subset2, subset = SingleR.labels == "Epithelial_cells")
epi_lower <- subset(lower_subset2, subset = SingleR.labels == "Epithelial_cells")
epi_ctrl  <- subset(ctrl_subset2,  subset = SingleR.labels == "Epithelial_cells")
#merge epithelial cells
epi_subset <- merge(epi_ctrl, y = list(epi_upper, epi_lower))
#normalise and reduce dimensions
epi_subset <- NormalizeData(epi_subset)
epi_subset <- FindVariableFeatures(epi_subset)
epi_subset <- ScaleData(epi_subset)
epi_subset <- RunPCA(epi_subset)

#cluster and umap
epi_subset <- FindNeighbors(epi_subset, dims = 1:20)
epi_subset <- FindClusters(epi_subset, resolution = 0.5)
epi_subset <- RunUMAP(epi_subset, dims = 1:20)

#set identities to clusters and join layers
Idents(epi_subset) <- epi_subset$seurat_clusters
epi_subset <- JoinLayers(epi_subset)
#find markers for epithelial subclusters
markers_epi <- FindAllMarkers(epi_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#get top 5 markers per cluster
top_markers_epi <- markers_epi %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
#heatmap of top epithelial markers
DoHeatmap(epi_subset, features = unique(top_markers_epi$gene)) +
  ggtitle("Top Marker Genes: Epithelial Subclusters") +
  theme(axis.text.y = element_text(size = 8))

#umap plot of epithelial clusters
DimPlot(epi_subset, label = TRUE, repel = TRUE) +
  ggtitle("UMAP of Epithelial Subclusters")

#subset fibroblasts only
fib_upper <- subset(upper_subset2, subset = SingleR.labels == "Fibroblasts")
fib_lower <- subset(lower_subset2, subset = SingleR.labels == "Fibroblasts")
fib_ctrl  <- subset(ctrl_subset2,  subset = SingleR.labels == "Fibroblasts")
#merge fibroblasts
fib_subset <- merge(fib_ctrl, y = list(fib_upper, fib_lower))
#normalise and reduce dimensions
fib_subset <- NormalizeData(fib_subset)
fib_subset <- FindVariableFeatures(fib_subset)
fib_subset <- ScaleData(fib_subset)
fib_subset <- RunPCA(fib_subset)
#cluster and umap
fib_subset <- FindNeighbors(fib_subset, dims = 1:20)
fib_subset <- FindClusters(fib_subset, resolution = 0.5)
fib_subset <- RunUMAP(fib_subset, dims = 1:20)

#set identities to clusters and join layers
Idents(fib_subset) <- fib_subset$seurat_clusters
fib_subset <- JoinLayers(fib_subset)
#find markers for fibroblast subclusters
markers_fib <- FindAllMarkers(fib_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#get top 5 markers per cluster
top_markers_fib <- markers_fib %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
#heatmap of top fibroblast markers
DoHeatmap(fib_subset, features = unique(top_markers_fib$gene)) +
  ggtitle("Top Marker Genes: Fibroblast Subclusters") +
  theme(axis.text.y = element_text(size = 8))

#umap plot of fibroblast clusters
DimPlot(fib_subset, label = TRUE, repel = TRUE) +
  ggtitle("UMAP of Fibroblast Subclusters")

#heatmap grouped by condition - epi
DoHeatmap(epi_subset, group.by = "condition", features = unique(top_markers_epi$gene)) +
  ggtitle("Top Epithelial Markers by Condition")
#fib by condition
DoHeatmap(fib_subset, group.by = "condition", features = unique(top_markers_fib$gene)) +
  ggtitle("Top Fibroblast Markers by Condition")

#epithelial cells heatmap split by both cluster and condition
DoHeatmap(
  epi_subset,
  features = unique(top_markers_epi$gene),
  group.by = c("seurat_clusters", "condition")
) +
  ggtitle("Top Epithelial Markers by Cluster and Condition") +
  theme(axis.text.y = element_text(size = 8))

#fibroblast cells heatmap split by both cluster and condition
DoHeatmap(
  fib_subset,
  features = unique(top_markers_fib$gene),
  group.by = c("seurat_clusters", "condition")
) +
  ggtitle("Top Fibroblast Markers by Cluster and Condition") +
  theme(axis.text.y = element_text(size = 8))

#proportions of clusters per condition?
epi_prop_table <- prop.table(table(epi_subset$seurat_clusters, epi_subset$condition), margin = 2)
fib_prop_table <- prop.table(table(fib_subset$seurat_clusters, fib_subset$condition), margin =2)
#as data frames
epi_df <- as.data.frame(epi_prop_table)
colnames(epi_df) <- c("Cluster", "Condition", "Proportion")
fib_df <- as.data.frame(fib_prop_table)
colnames(fib_df) <- c("Cluster", "Condition", "Proportion")
#plot using ggplot2
plot_prop_bar <- function(df, title) {
  ggplot(df, aes(x = Condition, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(title = title, y = "Proportion", x = "Condition") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
#plotn epi props
plot_prop_bar(epi_df, "Epithelial Subcluster Proportions by Condition")
#plot fib props
plot_prop_bar(fib_df, "Fibroblast Subcluster Proportions by Condition")

#stacked bar
plot_prop_stacked <- function(df, title) {
  ggplot(df, aes(x = Condition, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity", position = "stack") +   # stacked bars
    labs(title = title, y = "Proportion", x = "Condition") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
#plot epi with stacked bars
plot_prop_stacked(epi_df, "Epithelial Subcluster Proportions by Condition")
#fib with stacked bars
plot_prop_stacked(fib_df, "Fibroblast Subcluster Proportions by Condition")


# #statistical significant diff number cells per clusters?
# #for epi
#df with cluster and condition for each cell
epi_meta <- data.frame(
   Cluster = epi_subset$seurat_clusters,
   Condition = epi_subset$condition
 )
# 
# #cells per Cluster per Condition
# epi_counts <- epi_meta %>%
#   group_by(Condition, Cluster) %>%
#   summarise(Count = n()) %>%
#   tidyr::pivot_wider(names_from = Condition, values_from = Count, values_fill = 0)
# 

#for fib
 fib_meta <- data.frame(
   Cluster = fib_subset$seurat_clusters,
   Condition = fib_subset$condition
 )
# 
# fib_counts <- fib_meta %>%
#   group_by(Cluster, Condition) %>%
#   summarise(Count = n()) %>%
#   tidyr::pivot_wider(names_from = Condition, values_from = Count, values_fill = 0)
# 
# #function to run test per cluster
# run_proportion_test <- function(counts_row) {
#   #counts_row is one row with counts per condition
#   mat <- as.matrix(counts_row[ , -1])  #remove cluster name col
#   #chi-sq or fisher depends on counts 
#   if (any(mat < 5)) {
#     test <- fisher.test(mat)
#   } else {
#     test <- chisq.test(mat)
#   }
#   return(c(p.value = test$p.value))
# }
# 
# #epi clusters
# epi_pvals <- apply(epi_counts[ , -1], 1, function(x) {
#   #contingency table: rows = cluster presence (yes/no), cols = conditions
#   #cluster vs rest cells per condition, need total cell counts per condition
#   total_cells <- table(epi_meta$Condition)
#   cluster_cells <- x
#   rest_cells <- total_cells - cluster_cells
#   cont_table <- rbind(cluster_cells, rest_cells)
#   if (any(cont_table < 5)) {
#     fisher.test(cont_table)$p.value
#   } else {
#     chisq.test(cont_table)$p.value
#   }
# })
# 
# #combine with cluster names
# epi_results <- data.frame(
#   Cluster = epi_counts$Cluster,
#   p_value = epi_pvals
# )
# epi_results$adj_pvalue <- p.adjust(epi_results$p_value, method = "BH")
# print(epi_results)
# sig_clusters <- epi_results[epi_results$adj_pvalue < 0.05, ]
# print(sig_clusters)
# 
# 

#filter metadata for upper and lower epithelial cells only
epi_meta_ul <- epi_meta %>% filter(Condition %in% c("Upper", "Lower"))
#filter metadata for upper and lower fibroblast cells only
fib_meta_ul <- fib_meta %>% filter(Condition %in% c("Upper", "Lower"))

#count cells per cluster per condition for epithelial cells
epi_counts_ul <- epi_meta_ul %>%
  group_by(Cluster, Condition) %>%
  summarise(Count = n()) %>%
  pivot_wider(names_from = Condition, values_from = Count, values_fill = 0)
#count cells per cluster per condition for fibroblast cells
fib_counts_ul <- fib_meta_ul %>%
  group_by(Cluster, Condition) %>%
  summarise(Count = n()) %>%
  pivot_wider(names_from = Condition, values_from = Count, values_fill = 0)

#define function to run fisher or chi-square test per cluster for two conditions
test_two_conditions <- function(cluster_counts, total_per_condition) {
  rest_counts <- total_per_condition - cluster_counts
  cont_table <- rbind(cluster_counts, rest_counts)
  if (any(cont_table < 5)) {
    p <- fisher.test(cont_table)$p.value
  } else {
    p <- chisq.test(cont_table)$p.value
  }
  return(p)
}

#get total cells per condition for epithelial subset
total_epi <- table(epi_meta_ul$Condition)
#apply test to epithelial clusters
epi_pvals_ul <- apply(epi_counts_ul[ , c("Upper", "Lower")], 1, function(x) {
  test_two_conditions(cluster_counts = x, total_per_condition = total_epi)
})
#create results df for epithelial clusters
epi_results_ul <- data.frame(
  Cluster = epi_counts_ul$Cluster,
  p_value = epi_pvals_ul
)
#adjust p-values for multiple testing - epithelial
epi_results_ul$adj_pvalue <- p.adjust(epi_results_ul$p_value, method = "BH")

#get total cells per condition for fibroblast subset
total_fib <- table(fib_meta_ul$Condition)
#apply test to fibroblast clusters
fib_pvals_ul <- apply(fib_counts_ul[ , c("Upper", "Lower")], 1, function(x) {
  test_two_conditions(cluster_counts = x, total_per_condition = total_fib)
})
#create results dataframe for fibroblast clusters
fib_results_ul <- data.frame(
  Cluster = fib_counts_ul$Cluster,
  p_value = fib_pvals_ul
)
#adjust p-values for multiple testing - fibroblast
fib_results_ul$adj_pvalue <- p.adjust(fib_results_ul$p_value, method = "BH")

#print significant epithelial clusters between upper and lower
print("significant epithelial clusters (upper vs lower):")
print(epi_results_ul[epi_results_ul$adj_pvalue < 0.05, ])
print(epi_results_ul)
#print significant fibroblast clusters between upper and lower
print("significant fibroblast clusters (upper vs lower):")
print(fib_results_ul[fib_results_ul$adj_pvalue < 0.05, ])
print(fib_results_ul)
#export epithelial cluster results
write.csv(epi_results_ul, file = "epi_cluster_significance_upper_vs_lower.csv", row.names = FALSE)
#export fibroblast cluster results
write.csv(fib_results_ul, file = "fib_cluster_significance_upper_vs_lower.csv", row.names = FALSE)

#save epithelial subset with annotations and clustering
saveRDS(epi_subset, file = "epi_subset.rds")
#save fibroblast subset with annotations and clustering
saveRDS(fib_subset, file = "fib_subset.rds")

#subset endothelial cells only
endo_upper <- subset(upper_combined, subset = SingleR.labels == "Endothelial_cells")
endo_lower <- subset(lower_combined, subset = SingleR.labels == "Endothelial_cells")
endo_ctrl  <- subset(ctrl_combined,  subset = SingleR.labels == "Endothelial_cells")
#add condition label
endo_upper$condition <- "Upper"
endo_lower$condition <- "Lower"
endo_ctrl$condition  <- "Control"
#merge endothelial cells
endo_subset <- merge(endo_ctrl, y = list(endo_upper, endo_lower))
#normalise and reduce dimensions
endo_subset <- NormalizeData(endo_subset)
endo_subset <- FindVariableFeatures(endo_subset)
endo_subset <- ScaleData(endo_subset)
endo_subset <- RunPCA(endo_subset)
#cluster and umap
endo_subset <- FindNeighbors(endo_subset, dims = 1:20)
endo_subset <- FindClusters(endo_subset, resolution = 0.5)
endo_subset <- RunUMAP(endo_subset, dims = 1:20)
#set identities to clusters and join layers
Idents(endo_subset) <- endo_subset$seurat_clusters
endo_subset <- JoinLayers(endo_subset)
#find markers for endothelial subclusters
markers_endo <- FindAllMarkers(endo_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#get top 5 markers per cluster
top_markers_endo <- markers_endo %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
#heatmap of top endothelial markers
DoHeatmap(endo_subset, features = unique(top_markers_endo$gene)) +
  ggtitle("Top Marker Genes: Endothelial Subclusters") +
  theme(axis.text.y = element_text(size = 8))
#umap plot of endothelial clusters
DimPlot(endo_subset, label = TRUE, repel = TRUE) +
  ggtitle("UMAP of Endothelial Subclusters")
#save endothelial subset with annotations and clustering
saveRDS(endo_subset, file = "endo_subset.rds")

#proportions of clusters per condition - endothelial
endo_prop_table <- prop.table(table(endo_subset$seurat_clusters, endo_subset$condition), margin = 2)
#convert to data frame
endo_df <- as.data.frame(endo_prop_table)
colnames(endo_df) <- c("Cluster", "Condition", "Proportion")
#bar plot of endothelial proportions
plot_prop_bar(endo_df, "Endothelial Subcluster Proportions by Condition")
#stacked bar plot of endothelial proportions
plot_prop_stacked(endo_df, "Endothelial Subcluster Proportions by Condition")

#generate metadata for endothelial subset
endo_meta <- data.frame(
  Cluster = endo_subset$seurat_clusters,
  Condition = endo_subset$condition
)
#filter metadata for upper and lower endothelial cells only
endo_meta_ul <- endo_meta %>% filter(Condition %in% c("Upper", "Lower"))
#count cells per cluster per condition for endothelial cells
endo_counts_ul <- endo_meta_ul %>%
  group_by(Cluster, Condition) %>%
  summarise(Count = n()) %>%
  pivot_wider(names_from = Condition, values_from = Count, values_fill = 0)
#get total cells per condition for endothelial subset
total_endo <- table(endo_meta_ul$Condition)
#apply test to endothelial clusters using same function as before
endo_pvals_ul <- apply(endo_counts_ul[ , c("Upper", "Lower")], 1, function(x) {
  test_two_conditions(cluster_counts = x, total_per_condition = total_endo)
})
#create results dataframe for endothelial clusters
endo_results_ul <- data.frame(
  Cluster = endo_counts_ul$Cluster,
  p_value = endo_pvals_ul
)
#adjust p-values for multiple testing - endothelial
endo_results_ul$adj_pvalue <- p.adjust(endo_results_ul$p_value, method = "BH")
#print significant endothelial clusters between upper and lower
print("significant endothelial clusters (upper vs lower):")
print(endo_results_ul[endo_results_ul$adj_pvalue < 0.05, ])
print(endo_results_ul)
#export endothelial cluster results
write.csv(endo_results_ul, file = "endo_cluster_significance_upper_vs_lower.csv", row.names = FALSE)

#make sure the label and condition fields are present
epi_subset$cell_label <- paste0(epi_subset$SingleR.labels, "_", epi_subset$condition)
fib_subset$cell_label <- paste0(fib_subset$SingleR.labels, "_", fib_subset$condition)
endo_subset$cell_label <- paste0(endo_subset$SingleR.labels, "_", endo_subset$condition)
#merge all subsets and normalise (scale)
combined_isg15 <- merge(epi_subset, y = list(fib_subset, endo_subset))

# fetch raw or normalized expression for ISG15
expr <- FetchData(combined_isg15, vars = "ISG15")
# scale (z-score) manually
expr_scaled <- scale(expr$ISG15)
# add back scaled values as metadata for plotting
combined_isg15$ISG15_scaled <- expr_scaled
df <- FetchData(combined_isg15, vars = c("ISG15_scaled", "SingleR.labels", "condition"))
df$group <- factor(paste(df$SingleR.labels, df$condition, sep = "_"),
                   levels = c("Epithelial_cells_Control", "Epithelial_cells_Upper", "Epithelial_cells_Lower",
                              "Fibroblasts_Control", "Fibroblasts_Upper", "Fibroblasts_Lower",
                              "Endothelial_cells_Control", "Endothelial_cells_Upper", "Endothelial_cells_Lower"))

ggplot(df, aes(x = group, y = ISG15_scaled, fill = sub("_[^_]+$", "", group))) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.shape = NA) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cell type and Condition", y = "Scaled ISG15 expression", fill = "Cell Type")


#violin plot grouped by SingleR.labels + condition
VlnPlot(
  combined_isg15,
  features = "ISG15",
  group.by = "cell_label",
  pt.size = 0
) +
  ggtitle("ISG15 Expression by Cell Type and Condition (SingleR)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#fetch ISG15 expression from all three subsets and combine into one vector
epi_isg15  <- FetchData(epi_subset, vars = "ISG15")$ISG15
fib_isg15  <- FetchData(fib_subset, vars = "ISG15")$ISG15
endo_isg15 <- FetchData(endo_subset, vars = "ISG15")$ISG15
#combine all values and get max (ignoring NA values)
all_isg_values <- c(epi_isg15, fib_isg15, endo_isg15)
ymax <- max(all_isg_values, na.rm = TRUE)

#epithelial
VlnPlot(
  epi_subset,
  features = "ISG15",
  group.by = "condition",
  pt.size = 0
) +
  ggtitle("ISG15 in Epithelial Cells") +
  ylim(0, ymax)

#fibroblasts
VlnPlot(
  fib_subset,
  features = "ISG15",
  group.by = "condition",
  pt.size = 0
) +
  ggtitle("ISG15 in Fibroblast Cells") +
  ylim(0, ymax)

#endothelial
VlnPlot(
  endo_subset,
  features = "ISG15",
  group.by = "condition",
  pt.size = 0
) +
  ggtitle("ISG15 in Endothelial Cells") +
  ylim(0, ymax)

#97 cells in fib lower ie small sample group can skew
table(combined_isg15$SingleR.labels, combined_isg15$condition)["Fibroblasts", "Lower"]



#labelling clusters 
library(panglaoDBdata)
library(dplyr)

data("panglao_markers")

label_cluster <- function(cluster_genes, species = "Human", organ = "lung") {
  panglao_sub <- panglao_markers %>%
    filter(species == species, organ == organ)
  
  scores <- panglao_sub %>%
    group_by(cell_type) %>%
    summarise(overlap = sum(gene %in% cluster_genes)) %>%
    arrange(desc(overlap))
  
  return(scores)
}

# function to get top label per cluster from your markers df
get_top_labels <- function(top_markers_df) {
  clusters <- unique(top_markers_df$cluster)
  labels_list <- lapply(clusters, function(cl) {
    genes <- top_markers_df %>% filter(cluster == cl) %>% pull(gene)
    res <- label_cluster(genes)
    res$cluster <- cl
    res
  })
  labels_df <- do.call(rbind, labels_list)
  top_labels <- labels_df %>%
    group_by(cluster) %>%
    slice_max(order_by = overlap, n = 1) %>%
    ungroup()
  return(top_labels)
}

# run for epithelial
top_labels_epi <- get_top_labels(top_markers_epi)
print(top_labels_epi)
# run for fibroblast
top_labels_fib <- get_top_labels(top_markers_fib)
print(top_labels_fib)
# run for endothelial
top_labels_endo <- get_top_labels(top_markers_endo)
print(top_labels_endo)

#save for use in python
library(Seurat)
library(SeuratDisk)

library(Seurat)
library(SeuratDisk)
library(dplyr)

# Save each subset as an .rds file
saveRDS(epi_subset, file = "epi_clusters_saved.rds")
saveRDS(endo_subset, file = "endo_clusters_saved.rds")
saveRDS(fib_subset, file = "fib_clusters_saved.rds")

library(Matrix)

#export to make own clean new seurat obj with cells genes and metadata w no complications
export_clean_data <- function(seurat_obj, prefix) {
  # Extract counts matrix (usually sparse dgCMatrix)
  counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
  
  # Extract cell metadata
  meta <- seurat_obj@meta.data
  
  # Save counts as Matrix Market file
  Matrix::writeMM(counts, file = paste0(prefix, "_counts.mtx"))
  
  # Save gene names (features)
  write.csv(data.frame(Gene=rownames(counts)), file = paste0(prefix, "_genes.csv"), row.names = FALSE)
  
  # Save cell barcodes
  write.csv(data.frame(Cell=colnames(counts)), file = paste0(prefix, "_cells.csv"), row.names = FALSE)
  
  # Save metadata
  write.csv(meta, file = paste0(prefix, "_metadata.csv"))
  
  message("Export complete for ", prefix)
}

# Example: export for your three subsets
export_clean_data(epi_subset, "epi")
export_clean_data(endo_subset, "endo")
export_clean_data(fib_subset, "fib")

library(Seurat)
library(Matrix)

build_seurat_from_export <- function(prefix, reference_obj = NULL) {
  counts <- readMM(paste0(prefix, "_counts.mtx"))
  genes <- read.csv(paste0(prefix, "_genes.csv"))
  rownames(counts) <- genes$Gene
  cells <- read.csv(paste0(prefix, "_cells.csv"))
  colnames(counts) <- cells$Cell
  
  seurat_obj <- CreateSeuratObject(counts = counts, project = prefix)
  
  meta <- read.csv(paste0(prefix, "_metadata.csv"), row.names = 1)
  meta <- meta[colnames(seurat_obj), , drop = FALSE]
  seurat_obj <- AddMetaData(seurat_obj, metadata = meta)
  
  # Normalize data, this stores normalized counts in the "data" layer
  seurat_obj <- NormalizeData(seurat_obj)
  
  # If reference object provided, copy counts and normalized data from it
  if (!is.null(reference_obj)) {
    # Use GetAssayData to fetch counts and normalized data
    counts_ref <- GetAssayData(reference_obj, assay = "RNA", slot = "counts")
    data_ref <- GetAssayData(reference_obj, assay = "RNA", slot = "data")
    
    # Set these data in the layers of new object via SetAssayData
    seurat_obj <- SetAssayData(seurat_obj, assay = "RNA", slot = "counts", new.data = counts_ref)
    seurat_obj <- SetAssayData(seurat_obj, assay = "RNA", slot = "data", new.data = data_ref)
  }
  
  return(seurat_obj)
}

#fix assays before saving as rds
fix_assay_names <- function(seurat_obj, reference_obj) {
  # Get gene and cell names from reference (original) object
  gene_names <- rownames(reference_obj)
  cell_names <- colnames(reference_obj)
  
  # Extract counts and data matrices
  counts <- seurat_obj@assays$RNA@layers[["counts"]]
  data <- seurat_obj@assays$RNA@layers[["data"]]
  
  # Assign row and column names
  rownames(counts) <- gene_names
  colnames(counts) <- cell_names
  rownames(data) <- gene_names
  colnames(data) <- cell_names
  
  # Put matrices back in assay layers
  seurat_obj@assays$RNA@layers[["counts"]] <- counts
  seurat_obj@assays$RNA@layers[["data"]] <- data
  
  return(seurat_obj)
}


# Apply for all
epi_clean <- build_seurat_from_export("epi")
endo_clean <- build_seurat_from_export("endo")
fib_clean <- build_seurat_from_export("fib")
epi_clean <- fix_assay_names(epi_clean, epi_subset)
endo_clean <- fix_assay_names(endo_clean, endo_subset)
fib_clean <- fix_assay_names(fib_clean, fib_subset)

saveRDS(epi_clean, file='epi_clean.rds')
saveRDS(endo_clean, file='endo_clean.rds')
saveRDS(fib_clean, file='fib_clean.rds')



#export for use in celltypist
# clean_and_export <- function(seurat_obj, prefix) {
#   message("Starting cleanup for ", prefix)
#   
#   # Remove unwanted assay layers, keep only counts and data
#   layers <- names(seurat_obj@assays$RNA@layers)
#   layers_to_remove <- setdiff(layers, c("counts", "data"))
#   if(length(layers_to_remove) > 0) {
#     for(layer in layers_to_remove) {
#       seurat_obj@assays$RNA@layers[[layer]] <- NULL
#       message("Removed layer: ", layer)
#     }
#   }
#   
#   # Make sure counts and data are dgCMatrix
#   seurat_obj <- SetAssayData(seurat_obj, assay = "RNA", slot = "counts",
#                              new.data = as(GetAssayData(seurat_obj, assay = "RNA", slot = "counts"), "dgCMatrix"))
#   seurat_obj <- SetAssayData(seurat_obj, assay = "RNA", slot = "data",
#                              new.data = as(GetAssayData(seurat_obj, assay = "RNA", slot = "data"), "dgCMatrix"))
#   
#   # Slim object to essentials
#   seurat_obj <- DietSeurat(seurat_obj, assays = "RNA", dimreducs = NULL, graphs = NULL)
#   
#   # Try saving with error catching
#   tryCatch({
#     SaveH5Seurat(seurat_obj, filename = paste0(prefix, ".h5Seurat"), overwrite = TRUE)
#     Convert(paste0(prefix, ".h5Seurat"), dest = "h5ad", overwrite = TRUE)
#     message("Exported ", prefix, " successfully.")
#   }, error = function(e) {
#     message("Error exporting ", prefix, ": ", e$message)
#   })
#   
#   return(seurat_obj)
# }
# 
# # Then run:
# epi_clean <- clean_and_export(epi_clean, "epi_clean")
# endo_clean <- clean_and_export(endo_clean, "endo_clean")
# fib_clean <- clean_and_export(fib_clean, "fib_clean")

filter_and_export <- function(seurat_obj, prefix) {
  counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
  
  # Filter out zero-count genes
  keep_genes <- Matrix::rowSums(counts) > 0
  filtered_genes <- rownames(counts)[keep_genes]
  
  seurat_filtered <- subset(seurat_obj, features = filtered_genes)
  
  # Optional: print filtering summary
  message(prefix, ": Removed ", sum(!keep_genes), " zero-count genes; ", length(filtered_genes), " genes remain.")
  
  # Save h5Seurat and convert to h5ad
  SaveH5Seurat(seurat_filtered, paste0(prefix, "_filtered.h5Seurat"), overwrite = TRUE)
  Convert(paste0(prefix, "_filtered.h5Seurat"), dest = "h5ad", overwrite = TRUE)
  
  return(seurat_filtered)
}

# Then run for all three:
epi_filtered <- filter_and_export(epi_clean, "epi")
endo_filtered <- filter_and_export(endo_clean, "endo")
fib_filtered <- filter_and_export(fib_clean, "fib")

#convert step keeps failing but do have h5seurat objs for each
