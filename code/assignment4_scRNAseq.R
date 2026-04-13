# BINF*6110 Assignment 4
# scRNA-seq Analysis of IAV Nasal Mucosa
# Dataset: Kazer et al. 2024
# Comparing Naive vs D05 in Respiratory Mucosa (RM)


# ---- Load libraries ----

library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)
library(plyr)
library(DESeq2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(patchwork)


# ---- Load and subset data ----

# Load the Seurat object provided for the assignment
seurat_obj <- readRDS("seurat_ass4.rds")

# Check what the object looks like and what metadata is available
seurat_obj
head(seurat_obj@meta.data)

# Check what timepoints and tissues are in the data
table(seurat_obj@meta.data$time)
table(seurat_obj@meta.data$organ_custom)
table(seurat_obj@meta.data$time, seurat_obj@meta.data$organ_custom)

# Subset to Respiratory Mucosa only, keeping Naive and D05
# RM is the primary site of IAV infection
# D05 = peak myeloid recruitment, so we expect the strongest immune response here
seurat_sub <- subset(seurat_obj,
                     subset = organ_custom == "RM" & time %in% c("Naive", "D05"))

# Remove the full object to free up memory
rm(seurat_obj)
gc()

# Confirm the subset looks right
seurat_sub
table(seurat_sub@meta.data$time)


# ---- Quality control ----

# Calculate mitochondrial gene percentage per cell
# High mitochondrial content indicates stressed or dying cells
# Mouse mitochondrial genes start with "mt-"
seurat_sub[["percent.mt"]] <- PercentageFeatureSet(seurat_sub, pattern = "^mt-")

# Visualize QC metrics before filtering
VlnPlot(seurat_sub,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3,
        pt.size = 0)

# Filter out low quality cells based on the QC plots above
# Keeping cells with more than 200 genes detected (removes empty droplets)
# Removing cells with more than 4000 genes detected (potential doublets)
# Removing cells with more than 15% mitochondrial reads (stressed/dying cells)
seurat_sub <- subset(seurat_sub,
                     subset = nFeature_RNA > 200 &
                       nFeature_RNA < 4000 &
                       percent.mt < 15)

# Check how many cells remain after filtering
seurat_sub


# ---- Normalization, variable features, scaling, and PCA ----

# Log-normalize the data to account for differences in sequencing depth
seurat_sub <- NormalizeData(seurat_sub)

# Find the top 2000 most variable genes
# These genes carry the most biological signal across cells
seurat_sub <- FindVariableFeatures(seurat_sub, nfeatures = 2000)

# Scale the data so that each gene has mean 0 and variance 1
# This prevents highly expressed genes from dominating the analysis
seurat_sub <- ScaleData(seurat_sub)

# Run PCA for dimensionality reduction
seurat_sub <- RunPCA(seurat_sub, npcs = 30)

# Elbow plot to decide how many PCs to use
# We use 30 PCs because the standard deviation decreases gradually without a sharp drop
ElbowPlot(seurat_sub, ndims = 30)


# ---- Harmony integration ----

# Run Harmony to correct for batch effects across individual mice
# This ensures that cells cluster by cell type rather than by sample identity
# Harmony was chosen over CCA integration because it is lighter on memory
seurat_sub <- RunHarmony(seurat_sub,
                         group.by.vars = "mouse_id",
                         reduction = "pca",
                         reduction.save = "harmony")


# ---- Clustering and UMAP ----

# Build a shared nearest neighbour graph using the Harmony embedding
seurat_sub <- FindNeighbors(seurat_sub, reduction = "harmony", dims = 1:30)

# Cluster cells using the Louvain algorithm at resolution 0.5
seurat_sub <- FindClusters(seurat_sub, resolution = 0.5)

# Run UMAP for 2D visualization using the Harmony embedding
seurat_sub <- RunUMAP(seurat_sub, reduction = "harmony", dims = 1:30)

# Plot UMAP coloured by cluster number and by timepoint side by side
p1 <- DimPlot(seurat_sub, reduction = "umap", label = TRUE) + ggtitle("Clusters")
p2 <- DimPlot(seurat_sub, reduction = "umap", group.by = "time") + ggtitle("Naive vs D05")
p1 + p2

# Save progress before moving to annotation
saveRDS(seurat_sub, "seurat_sub_clustered.rds")


# ---- Manual cluster annotation ----

# Find the top marker genes for each cluster to identify cell types
# only.pos = TRUE keeps only upregulated markers
# min.pct = 0.25 means the gene must be expressed in at least 25% of cells in the cluster
# logfc.threshold = 0.25 sets a minimum fold change cutoff
all_markers <- FindAllMarkers(seurat_sub,
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25)

# View the top 3 markers per cluster to guide annotation
top_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 3)

print(top_markers, n = 87)

# Assign cell type labels to each cluster based on canonical marker genes
# Annotation was validated using PanglaoDB and Durante et al. 2020
new_labels <- c(
  "0"  = "Olfactory Sensory Neurons",    # Slc1a2, Cnga4
  "1"  = "Respiratory Epithelial Cells", # Defb1, Serpinb5
  "2"  = "Macrophages",                  # Adgre4, Gpr141
  "3"  = "M2 Macrophages",               # Fcrls, Mrc1, Pf4
  "4"  = "NK Cells",                     # Prf1, Gzma, Ncr1
  "5"  = "Immature Neurons",             # Gap43, Ly6h
  "6"  = "Fibroblasts",                  # Dpt, Scara5, Hhip
  "7"  = "Secretory Epithelial Cells",   # Bpifa1, Tff2, Muc4
  "8"  = "Endothelial Cells",            # Sele, Ptprb, Gpihbp1
  "9"  = "Olfactory Sensory Neurons",    # Abca4, Nrcam
  "10" = "Rare/Unknown",                 # Thsd7b, Tex15
  "11" = "Dendritic Cells",              # Clec9a, Flt3, Cd209a
  "12" = "Neutrophils",                  # S100a9, Cxcr2, Retnlg
  "13" = "Rare/Unknown",                 # Mup4
  "14" = "Pericytes",                    # Rgs4, Higd1b
  "15" = "Ionocytes",                    # Ascl3, Kl, Galnt13
  "16" = "Club Cells",                   # Scgb2b27, Car6
  "17" = "B Cells",                      # Ms4a1, Iglc1, Fcmr
  "18" = "Goblet Cells",                 # Muc2, Sec14l3
  "19" = "Sustentacular Cells",          # Calb1, Slc22a20
  "20" = "Neuronal Progenitors",         # Neurog1, Hist1h1b
  "21" = "Ciliated Cells",               # Pih1h3b
  "22" = "Proliferating Cells",          # Cdc20, Plk1
  "23" = "Sensory Support Cells",        # Otoa, Svopl
  "24" = "Chondrocytes",                 # Matn3, Ucma, Snorc
  "25" = "Osteoblasts",                  # Bglap, Bglap2
  "26" = "Pvalb+ Neurons",               # Pvalb, Ankrd63
  "27" = "Schwann Cells",                # Mpz, Cldn19, Foxd3
  "28" = "Glandular Cells"               # Ptgds, Bpifb6, Ces1a
)

seurat_sub$cell_type <- plyr::mapvalues(
  seurat_sub$seurat_clusters,
  from = names(new_labels),
  to = new_labels
)

# Plot the annotated UMAP
DimPlot(seurat_sub,
        reduction = "umap",
        group.by = "cell_type",
        label = TRUE,
        repel = TRUE,
        label.size = 3) +
  ggtitle("Cell Type Annotation") +
  theme(legend.position = "bottom")


# ---- Feature plots ----

# Visualize expression of known marker genes across the UMAP
# This validates that our cluster annotations are biologically correct
FeaturePlot(seurat_sub,
            features = c("C1qa",    # Macrophages
                         "S100a9",  # Neutrophils
                         "Ms4a1",   # B Cells
                         "Ncr1",    # NK Cells
                         "Ptprb",   # Endothelial Cells
                         "Defb1"),  # Respiratory Epithelial Cells
            reduction = "umap",
            ncol = 3)


# ---- Pseudobulk differential expression using DESeq2 ----

# Subset to macrophage clusters for DE analysis
# Macrophages are a key innate immune cell type during IAV infection
macro <- subset(seurat_sub,
                subset = cell_type %in% c("Macrophages", "M2 Macrophages"))

# Pseudobulk: aggregate counts by mouse_id and timepoint
# This approach is recommended over single-cell DE because it reduces false positives
# by treating each biological replicate (mouse) as a single observation
pseudo <- AggregateExpression(macro,
                              assays = "RNA",
                              return.seurat = TRUE,
                              group.by = c("mouse_id", "time"))

# Add condition labels and set as identity
pseudo$time <- ifelse(grepl("Naive", colnames(pseudo)), "Naive", "D05")
Idents(pseudo) <- "time"

# Run DE analysis using DESeq2
de_results <- FindMarkers(pseudo,
                          ident.1 = "D05",
                          ident.2 = "Naive",
                          test.use = "DESeq2",
                          min.pct = 0.1,
                          logfc.threshold = 0.25)

# Check top DE genes
head(de_results, 20)

# Label genes by significance category for plotting
de_results$gene <- rownames(de_results)
de_results$category <- "Not Significant"
de_results$category[de_results$avg_log2FC > 0.5 & de_results$p_val_adj < 0.05] <- "Upregulated"
de_results$category[de_results$avg_log2FC < -0.5 & de_results$p_val_adj < 0.05] <- "Downregulated"

# Get top 15 most significant genes to label on the volcano plot
top_genes <- de_results %>%
  filter(category != "Not Significant") %>%
  arrange(p_val_adj) %>%
  head(15)

# Volcano plot of DE results
ggplot(de_results, aes(x = avg_log2FC, y = -log10(p_val_adj), color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "Not Significant" = "grey")) +
  geom_text(data = top_genes, aes(label = gene),
            size = 3, vjust = -0.5, color = "black") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  theme_classic() +
  labs(title = "Macrophage DE: D05 vs Naive",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value",
       color = "")

# Violin plots of top interferon-stimulated genes across conditions
VlnPlot(macro,
        features = c("Rsad2", "Ifit1", "Isg15", "Ifit3"),
        group.by = "time",
        ncol = 2,
        pt.size = 0.5) &
  theme(plot.title = element_text(size = 11, face = "bold.italic"),
        axis.title.x = element_blank())


# ---- GSEA on GO Biological Processes ----

# Create a ranked gene list sorted by log2 fold change
# Genes most upregulated in D05 are at the top, most downregulated at the bottom
gene_list <- de_results$avg_log2FC
names(gene_list) <- rownames(de_results)
gene_list <- sort(gene_list, decreasing = TRUE)

# Convert mouse gene symbols to Entrez IDs (required by clusterProfiler)
gene_df <- bitr(names(gene_list),
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Mm.eg.db)

# Match Entrez IDs back to the ranked gene list
gene_list_entrez <- gene_list[gene_df$SYMBOL]
names(gene_list_entrez) <- gene_df$ENTREZID
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)

# Run GSEA using GO Biological Process gene sets
gsea_results <- gseGO(geneList = gene_list_entrez,
                      OrgDb = org.Mm.eg.db,
                      ont = "BP",
                      minGSSize = 15,
                      maxGSSize = 500,
                      pvalueCutoff = 0.05,
                      verbose = FALSE)

# Dotplot showing activated and suppressed pathways
dotplot(gsea_results,
        showCategory = 15,
        split = ".sign") +
  facet_grid(. ~ .sign) +
  ggtitle("GSEA: Macrophages D05 vs Naive")


# ---- Save publication-quality figures ----

# Figure 1: QC violin plot
png("results/Fig1_QC_violin.png", width = 6000, height = 2000, res = 250)
p <- VlnPlot(seurat_sub,
             features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
             ncol = 3, pt.size = 0) &
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 13, face = "bold"),
        plot.margin = margin(10, 20, 10, 20))
p + plot_annotation(
  title = "Figure 1: Quality Control Metrics Across Clusters",
  theme = theme(plot.title = element_text(size = 15, face = "bold")))
dev.off()

# Figure 2: Elbow plot
png("results/Fig2_elbow.png", width = 1800, height = 1200, res = 250)
ElbowPlot(seurat_sub, ndims = 30) +
  labs(title = "Figure 2: Elbow Plot of Principal Components",
       x = "Principal Component",
       y = "Standard Deviation") +
  theme_classic() +
  theme(plot.title = element_text(size = 12, face = "bold"))
dev.off()

# Figure 3: UMAP clusters and timepoint
png("results/Fig3_umap_clusters_time.png", width = 3600, height = 1400, res = 250)
p1 <- DimPlot(seurat_sub, reduction = "umap", label = TRUE, repel = TRUE) +
  labs(title = "A: Unsupervised Clusters") +
  theme(plot.title = element_text(size = 12, face = "bold"))
p2 <- DimPlot(seurat_sub, reduction = "umap", group.by = "time") +
  labs(title = "B: Naive vs D05") +
  theme(plot.title = element_text(size = 12, face = "bold"))
print(p1 + p2 + plot_annotation(
  title = "Figure 3: UMAP Visualization of Clusters and Timepoints",
  theme = theme(plot.title = element_text(size = 14, face = "bold"))))
dev.off()

# Figure 4: Annotated UMAP
png("results/Fig4_umap_annotated.png", width = 3000, height = 2200, res = 250)
DimPlot(seurat_sub,
        reduction = "umap",
        group.by = "cell_type",
        label = TRUE,
        repel = TRUE,
        label.size = 3) +
  labs(title = "Figure 4: Cell Type Annotation of Respiratory Mucosa") +
  theme(plot.title = element_text(size = 13, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 8)) +
  guides(color = guide_legend(ncol = 4, override.aes = list(size = 3)))
dev.off()

# Figure 5: Feature plots
png("results/Fig5_feature_plots.png", width = 3600, height = 2400, res = 250)
p2 <- FeaturePlot(seurat_sub,
                  features = c("C1qa", "S100a9", "Ms4a1",
                               "Ncr1", "Ptprb", "Defb1"),
                  reduction = "umap", ncol = 3) &
  theme(plot.title = element_text(size = 11, face = "bold.italic"))
p2 + plot_annotation(
  title = "Figure 5: Feature Plots of Canonical Cell Type Markers",
  theme = theme(plot.title = element_text(size = 14, face = "bold")))
dev.off()

# Figure 6: Volcano plot
png("results/Fig6_volcano.png", width = 2400, height = 2000, res = 250)
ggplot(de_results, aes(x = avg_log2FC, y = -log10(p_val_adj), color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "#D62728",
                                "Downregulated" = "#1F77B4",
                                "Not Significant" = "grey70")) +
  geom_text(data = top_genes, aes(label = gene),
            size = 3, vjust = -0.5, color = "black") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  theme_classic() +
  labs(title = "Figure 6: Differential Gene Expression in Macrophages (D05 vs Naive)",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value",
       color = "") +
  theme(plot.title = element_text(size = 11, face = "bold"),
        legend.position = "top")
dev.off()

# Figure 7: Violin plots of top ISGs
png("results/Fig7_violin_ISGs.png", width = 2400, height = 2000, res = 250)
VlnPlot(macro,
        features = c("Rsad2", "Ifit1", "Isg15", "Ifit3"),
        group.by = "time",
        ncol = 2,
        pt.size = 0.5) &
  theme(plot.title = element_text(size = 11, face = "bold.italic"),
        axis.title.x = element_blank())
dev.off()

# Figure 8: GSEA dotplot
png("results/Fig8_gsea.png", width = 3000, height = 3600, res = 250)
dotplot(gsea_results,
        showCategory = 15,
        split = ".sign") +
  facet_grid(. ~ .sign) +
  labs(title = "Figure 8: GSEA of GO Biological Processes in Macrophages (D05 vs Naive)") +
  theme(plot.title = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 9))
dev.off()


# ---- Save all outputs ----

saveRDS(seurat_sub, "seurat_sub_final.rds")
write.csv(de_results, "code/macrophage_DE_results.csv")
write.csv(as.data.frame(gsea_results), "code/gsea_results.csv")

print("All done!")