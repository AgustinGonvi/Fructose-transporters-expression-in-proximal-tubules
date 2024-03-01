library(Seurat)
library(ggplot2)
library(stringr)
library(writexl)
library(gridExtra)
library(tidyverse)

print(getwd())
dir.create("./Results")
list.dirs()

print(getwd())
dataPath <- "G:\\Shared drives\\Gonzalez-Vicente Lab\\03_Data\\2023_Rat scRNAseq\\01_Resolution 1.6\\Results\\ALL.combined.rds"
print(dataPath)

all.combined <- readRDS(dataPath)

DimPlot(object = all.combined,
              reduction = "umap",
              label = T)
ggsave(filename = paste("./Results/01_ALL clusters_UMAP", ".tiff", sep = ""),
       device = "tiff",
       dpi = 600)
DimPlot(object = all.combined,
              reduction = "tsne",
              label = T)
ggsave(filename = paste("./Results/01_ALL clusters_TSNE", ".tiff", sep = ""),
       device = "tiff",
       dpi = 600)
DimPlot(object = all.combined,
              reduction = "pca",
              label = T)
ggsave(filename = paste("./Results/01_ALL clusters_PCA", ".tiff", sep = ""),
       device = "tiff",
       dpi = 600)

all.combined$Cell.Class <- case_when(
  all.combined$seurat_clusters %in% c(26,29,19) ~ "Endothelial",
  all.combined$seurat_clusters %in% c(1,18,3,0,4,13,27,23,30,12,11,9,22,7,35,25,21,32,15,14,16,6,31,17) ~ "Epithelial",
  all.combined$seurat_clusters %in% c(24,34,10,28,5,2,8,20) ~ "Immune",
  all.combined$seurat_clusters %in% c(33) ~ "Stromal")
table(all.combined$Cell.Class)

DimPlot(object = all.combined,
        reduction = "umap",
        group.by = "Cell.Class",
        label = T)
ggsave(filename = paste("./Results/02_Cell.Class_UMAP", ".tiff", sep = ""),
       device = "tiff",
       dpi = 600)
DimPlot(object = all.combined,
        reduction = "tsne",
        group.by = "Cell.Class",
        label = T)
ggsave(filename = paste("./Results/02_Cell.Class_TSNE", ".tiff", sep = ""),
       device = "tiff",
       dpi = 600)
DimPlot(object = all.combined,
        reduction = "pca",
        group.by = "Cell.Class",
        repel = T,
        label = T)
ggsave(filename = paste("./Results/02_Cell.Class_PCA", ".tiff", sep = ""),
       device = "tiff",
       dpi = 600)

all.combined <- SetIdent(all.combined, value = all.combined$Cell.Class)

Epithelial <- subset(all.combined,
                     idents = "Epithelial")

DimPlot(object = Epithelial,
        reduction = "umap",
        group.by = "Cell.Class",
        repel = T,
        label = T)
ggsave(filename = paste("./Results/03_Cell.Class_UMAP", ".tiff", sep = ""),
       device = "tiff",
       dpi = 600)
DimPlot(object = Epithelial,
        reduction = "tsne",
        group.by = "Cell.Class",
        repel = T,
        label = T)
ggsave(filename = paste("./Results/03_Cell.Class_TSNE", ".tiff", sep = ""),
       device = "tiff",
       dpi = 600)
DimPlot(object = Epithelial,
        reduction = "pca",
        group.by = "Cell.Class",
        repel = T,
        label = T)
ggsave(filename = paste("./Results/03_Cell.Class_PCA", ".tiff", sep = ""),
       device = "tiff",
       dpi = 600)
Epithelial <- SetIdent(Epithelial, value = Epithelial$seurat_clusters)

HubMAP_Curated <- c(                 
                    "Slc5a12",                       # PT_S1
                    "Slc22a8","Slc22a6",             # PT_S1/S2
                    "Slc22a7", "Slc7a13",            # PT_S3
                    "Satb2",                         # PT_S3/DTL
                    "Spp1",                          # DTL
                    "Cldn3",                         # ATL
                    "Prox1",                         # M_TAL
                    "Slc12a1",                       # TAL
                    "Cldn16","Enox1",                # C_TAL
                    "Slc12a3",                       # DCT
                    "Scnn1g",                        # CNT/CD_PC
                    "Slc4a1",                        # CNT/CD_IC.A
                    "Slc26a4",                       # CD_IC.B
                    "Slc14a1"                       # Papillary Tip
                    #"Slc14a2"                        # IMCD/DTL
                    )


DotPlot(Epithelial,
        features = HubMAP_Curated,
        assay = "RNA",
        cols = c("white", "blue"),
        scale = T,
        scale.by = "radius",
        dot.scale = 6) +
    RotatedAxis()
ggsave(filename = paste("./Results/04_Cluster Markers_Epithelial", ".tiff", sep = ""),
       bg = "white",
       device = "tiff",
       dpi = 600)

table(all.combined$seurat_clusters)

all.combined$Cell.Type1 <- case_when(
  all.combined$seurat_clusters %in% c(26,29,19,24,34,10,28,5,2,8,20,33) ~ NA,
  all.combined$seurat_clusters %in% c(1,18)         ~ "PT.S1",
  all.combined$seurat_clusters %in% c(0,3)          ~ "PT.S2",
  all.combined$seurat_clusters %in% c(4,13)         ~ "PT.S3",
  all.combined$seurat_clusters %in% c(12)           ~ "DTL",
  all.combined$seurat_clusters %in% c(30)           ~ "DTL",
  all.combined$seurat_clusters %in% c(23)           ~ "ATL",
  all.combined$seurat_clusters %in% c(22,9)         ~ "MTAL",
  all.combined$seurat_clusters %in% c(7,35)         ~ "CTAL",
  all.combined$seurat_clusters %in% c(25,11)        ~ "CTAL",
  all.combined$seurat_clusters %in% c(21)           ~ "DCT",
  all.combined$seurat_clusters %in% c(16,14)        ~ "IC.A",
  all.combined$seurat_clusters %in% c(6,31)         ~ "IC.B",
  all.combined$seurat_clusters %in% c(15)           ~ "PC",
  all.combined$seurat_clusters %in% c(32)           ~ "IMCD",
  all.combined$seurat_clusters %in% c(17,27)        ~ NA)
table(all.combined$Cell.Type1)

all.combined$Cell.Type2 <- case_when(
  all.combined$seurat_clusters %in% c(26,29,19,24,34,10,28,5,2,8,20,33) ~ NA,
  all.combined$seurat_clusters %in% c(1,18)         ~ "PT",
  all.combined$seurat_clusters %in% c(0,3)          ~ "PT",
  all.combined$seurat_clusters %in% c(4,13)         ~ "PT",
  all.combined$seurat_clusters %in% c(12)           ~ "DTL",
  all.combined$seurat_clusters %in% c(30)           ~ "DTL",
  all.combined$seurat_clusters %in% c(23)           ~ "ATL",
  all.combined$seurat_clusters %in% c(22,9)         ~ "TAL",
  all.combined$seurat_clusters %in% c(7,35)         ~ "TAL",
  all.combined$seurat_clusters %in% c(25,11)        ~ "TAL",
  all.combined$seurat_clusters %in% c(21)           ~ "DCT",
  all.combined$seurat_clusters %in% c(16,14)        ~ "IC",
  all.combined$seurat_clusters %in% c(6,31)         ~ "IC",
  all.combined$seurat_clusters %in% c(15)           ~ "PC",
  all.combined$seurat_clusters %in% c(32)           ~ "IMCD",
  all.combined$seurat_clusters %in% c(17,27)        ~ NA)
table(all.combined$Cell.Type2)

all.combined$Region <- case_when(
  all.combined$seurat_clusters %in% c(26,29,19,24,34,10,28,5,2,8,20,33) ~ NA,
  all.combined$seurat_clusters %in% c(1,18)         ~ "PT",
  all.combined$seurat_clusters %in% c(0,3)          ~ "PT",
  all.combined$seurat_clusters %in% c(4,13)         ~ "PT",
  all.combined$seurat_clusters %in% c(12)           ~ "Thin.limbs",
  all.combined$seurat_clusters %in% c(30)           ~ "Thin.limbs",
  all.combined$seurat_clusters %in% c(23)           ~ "Thin.limbs",
  all.combined$seurat_clusters %in% c(22,9)         ~ "TAL",
  all.combined$seurat_clusters %in% c(7,35)         ~ "TAL",
  all.combined$seurat_clusters %in% c(25,11)        ~ "TAL",
  all.combined$seurat_clusters %in% c(21)           ~ "DCT",
  all.combined$seurat_clusters %in% c(16,14)        ~ "CD",
  all.combined$seurat_clusters %in% c(6,31)         ~ "CD",
  all.combined$seurat_clusters %in% c(15)           ~ "CD",
  all.combined$seurat_clusters %in% c(32)           ~ "CD",
  all.combined$seurat_clusters %in% c(17,27)        ~ NA)
table(all.combined$Region)


rm(Epithelial)

all.combined <- SetIdent(all.combined, value = all.combined$Cell.Class)

Epithelial <- subset(all.combined,
                     idents = "Epithelial")

Epithelial <- SetIdent(Epithelial, value = Epithelial$Cell.Type1)

Epithelial <- subset(Epithelial,
                     idents = c("PT.S1",
                                "PT.S2",
                                "PT.S3",
                                "DTL",
                                "ATL",
                                "MTAL",
                                "CTAL",
                                "DCT",
                                "PC",
                                "IC.A",
                                "IC.B",
                                "IMCD"))

Epithelial@active.ident <- factor(Epithelial@active.ident,
                                  levels=c("PT.S1",
                                           "PT.S2",
                                           "PT.S3",
                                           "DTL",
                                           "ATL",
                                           "MTAL",
                                           "CTAL",
                                           "DCT",
                                           "PC",
                                           "IC.A",
                                           "IC.B",
                                           "IMCD"
                                           ))

DimPlot(object = Epithelial,
        reduction = "umap",
        repel = T,
        label = T)
ggsave(filename = paste("./Results/05_Cell.Type1_UMAP", ".tiff", sep = ""),
       device = "tiff",
       dpi = 600)
DimPlot(object = Epithelial,
        reduction = "tsne",
        repel = T,
        label = T)
ggsave(filename = paste("./Results/05_Cell.Type1_TSNE", ".tiff", sep = ""),
       device = "tiff",
       dpi = 600)
DimPlot(object = Epithelial,
        reduction = "pca",
        repel = T,
        label = T)
ggsave(filename = paste("./Results/05_Cell.Type1_PCA", ".tiff", sep = ""),
       device = "tiff",
       dpi = 600)

table(Epithelial$Cell.Type1)

DotPlot(Epithelial,
        features = HubMAP_Curated,
        assay = "RNA",
        cols = c("white", "blue"),
        scale = T,
        scale.by = "radius",
        dot.scale = 6) +
    RotatedAxis()
ggsave(filename = paste("./Results/06_DotPlot_Cell.Type1", ".tiff", sep = ""),
       bg = "white",
       device = "tiff",
       dpi = 600)

dir.create(file.path("./Results", "Cluster Markers_Cell.Type1"), showWarnings = T)

DefaultAssay(object = Epithelial) <- "RNA"

markers_list <- list()
# Loop over each cell identity (cluster)
for (clus in unique(Epithelial@active.ident)) {
  # Find markers for the current cluster
  markers <- FindMarkers(Epithelial,
                         ident.1 = clus,
                         logfc.threshold = 1,
                         min.pct = 0.25,
                         test.use = "wilcox",
                         only.pos = FALSE,
                         max.cells.per.ident = Inf,
                         random.seed = 1,
                         verbose = TRUE)
  
  # Add the markers to the list, with the cluster identity as the list name
  markers_list[[paste("Cluster", clus, sep = " ")]] <- markers
  
  # Run garbage collection and remove objects from memory
  gc(verbose = TRUE)
}

# Write the markers to separate CSV files
lapply(names(markers_list), function(x) {
  write.csv(markers_list[[x]],
            file = paste0("./Results/Cluster Markers_Cell.Type1/", x, "_Markers.csv"),
            row.names = T)
})

dir.create(file.path("./Results", "Cluster Averages_Cell.Type1"), showWarnings = T)

cluster.averages <- AverageExpression(Epithelial)
head(cluster.averages[["RNA"]][1:10,])
write.csv(cluster.averages$RNA, "./Results/Cluster Averages_Cell.Type1/Cell.Type1 Clusters mean RNA expression.csv")

cluster.averages <- AverageExpression(Epithelial)
head(cluster.averages[["SCT"]][1:10,])
write.csv(cluster.averages$SCT, "./Results/Cluster Averages_Cell.Type1/Cell.Type1 Clusters mean SCT expression.csv")

# Write capitalized genes for comparing to humans
df <- cluster.averages$RNA
rownames(df) <-   str_to_upper(rownames(df))
write.csv(df, "./Results/Cluster Averages_Cell.Type1/Cell.Type1 Clusters mean RNA expression_CAPS.csv")

df <- cluster.averages$SCT
rownames(df) <-   str_to_upper(rownames(df))
write.csv(df, "./Results/Cluster Averages_Cell.Type1/Cell.Type1 Clusters mean SCT expression_CAPS.csv")
rm(df)

Epithelial <- SetIdent(Epithelial, value = Epithelial$Cell.Type2)

Epithelial <- subset(Epithelial,
                     idents = c("PT",
                                "DTL",
                                "ATL",
                                "TAL",
                                "DCT",
                                "PC",
                                "IC"))

Epithelial@active.ident <- factor(Epithelial@active.ident,
                                  levels=c("PT",
                                           "DTL",
                                           "ATL",
                                           "TAL",
                                           "DCT",
                                           "PC",
                                           "IC"))
DimPlot(object = Epithelial,
        reduction = "umap",
        #group.by = "Cell.Type2",
        repel = T,
        label = T)
ggsave(filename = paste("./Results/07_Cell.Type2_UMAP", ".tiff", sep = ""),
       device = "tiff",
       dpi = 600)
DimPlot(object = Epithelial,
        reduction = "tsne",
        #group.by = "Cell.Type2",
        repel = T,
        label = T)
ggsave(filename = paste("./Results/07_Cell.Type2_TSNE", ".tiff", sep = ""),
       device = "tiff",
       dpi = 600)
DimPlot(object = Epithelial,
        reduction = "pca",
        #group.by = "Cell.Type",
        repel = T,
        label = T)
ggsave(filename = paste("./Results/07_Cell.Type2_PCA", ".tiff", sep = ""),
       device = "tiff",
       dpi = 600)

table(Epithelial$Cell.Type2)

DotPlot(Epithelial,
        features = HubMAP_Curated,
        assay = "RNA",
        cols = c("white", "blue"),
        scale = T,
        scale.by = "radius",
        dot.scale = 6) +
    RotatedAxis()
ggsave(filename = paste("./Results/08_DotPlot_Cell.Type2", ".tiff", sep = ""),
       bg = "white",
       device = "tiff",
       dpi = 600)

dir.create(file.path("./Results", "Cluster Markers_Cell.Type2"), showWarnings = T)

DefaultAssay(object = Epithelial) <- "RNA"

markers_list <- list()
# Loop over each cell identity (cluster)
for (clus in unique(Epithelial@active.ident)) {
  # Find markers for the current cluster
  markers <- FindMarkers(Epithelial,
                         ident.1 = clus,
                         logfc.threshold = 1,
                         min.pct = 0.25,
                         test.use = "wilcox",
                         only.pos = FALSE,
                         max.cells.per.ident = Inf,
                         random.seed = 1,
                         verbose = TRUE)
  
  # Add the markers to the list, with the cluster identity as the list name
  markers_list[[paste("Cluster", clus, sep = " ")]] <- markers
  
  # Run garbage collection and remove objects from memory
  gc(verbose = TRUE)
}

# Write the markers to separate CSV files
lapply(names(markers_list), function(x) {
  write.csv(markers_list[[x]],
            file = paste0("./Results/Cluster Markers_Cell.Type2/", x, "_Markers.csv"),
            row.names = T)
})

dir.create(file.path("./Results", "Cluster Averages_Cell.Type2"), showWarnings = T)

cluster.averages <- AverageExpression(Epithelial)
head(cluster.averages[["RNA"]][1:10,])
write.csv(cluster.averages$RNA, "./Results/Cluster Averages_Cell.Type2/Cell.Type2 Clusters mean RNA expression.csv")

cluster.averages <- AverageExpression(Epithelial)
head(cluster.averages[["SCT"]][1:10,])
write.csv(cluster.averages$SCT, "./Results/Cluster Averages_Cell.Type2/Cell.Type2 Clusters mean SCT expression.csv")

# Write capitalized genes for comparing to humans
df <- cluster.averages$RNA
rownames(df) <-   str_to_upper(rownames(df))
write.csv(df, "./Results/Cluster Averages_Cell.Type2/Cell.Type2 Clusters mean RNA expression_CAPS.csv")

df <- cluster.averages$SCT
rownames(df) <-   str_to_upper(rownames(df))
write.csv(df, "./Results/Cluster Averages_Cell.Type2/Cell.Type2 Clusters mean SCT expression_CAPS.csv")
rm(df)

Epithelial <- SetIdent(Epithelial, value = Epithelial$Region)

Epithelial <- subset(Epithelial,
                     idents = c("PT",
                                "Thin.limbs",
                                "TAL",
                                "DCT",
                                "CD"))

Epithelial@active.ident <- factor(Epithelial@active.ident,
                                  levels=c("PT",
                                           "Thin.limbs",
                                           "TAL",
                                           "DCT",
                                           "CD"))
DimPlot(object = Epithelial,
        reduction = "umap",
        repel = T,
        label = T)
ggsave(filename = paste("./Results/09_Region_UMAP", ".tiff", sep = ""),
       device = "tiff",
       dpi = 600)
DimPlot(object = Epithelial,
        reduction = "tsne",
        repel = T,
        label = T)
ggsave(filename = paste("./Results/09_Region_TSNE", ".tiff", sep = ""),
       device = "tiff",
       dpi = 600)
DimPlot(object = Epithelial,
        reduction = "pca",
        #group.by = "Cell.Type",
        repel = T,
        label = T)
ggsave(filename = paste("./Results/09_Region_PCA", ".tiff", sep = ""),
       device = "tiff",
       dpi = 600)

table(Epithelial$Region)

DotPlot(Epithelial,
        features = HubMAP_Curated,
        assay = "RNA",
        cols = c("white", "blue"),
        scale = T,
        scale.by = "radius",
        dot.scale = 6) +
    RotatedAxis()
ggsave(filename = paste("./Results/10_DotPlot_Region", ".tiff", sep = ""),
       bg = "white",
       device = "tiff",
       dpi = 600)

dir.create(file.path("./Results", "Cluster Markers_Region"), showWarnings = T)

DefaultAssay(object = Epithelial) <- "RNA"

markers_list <- list()
# Loop over each cell identity (cluster)
for (clus in unique(Epithelial@active.ident)) {
  # Find markers for the current cluster
  markers <- FindMarkers(Epithelial,
                         ident.1 = clus,
                         logfc.threshold = 1,
                         min.pct = 0.25,
                         test.use = "wilcox",
                         only.pos = FALSE,
                         max.cells.per.ident = Inf,
                         random.seed = 1,
                         verbose = TRUE)
  
  # Add the markers to the list, with the cluster identity as the list name
  markers_list[[paste("Cluster", clus, sep = " ")]] <- markers
  
  # Run garbage collection and remove objects from memory
  gc(verbose = TRUE)
}

# Write the markers to separate CSV files
lapply(names(markers_list), function(x) {
  write.csv(markers_list[[x]],
            file = paste0("./Results/Cluster Markers_Region/", x, "_Markers.csv"),
            row.names = T)
})

dir.create(file.path("./Results", "Cluster Averages_Region"), showWarnings = T)

cluster.averages <- AverageExpression(Epithelial)
head(cluster.averages[["RNA"]][1:10,])
write.csv(cluster.averages$RNA, "./Results/Cluster Averages_Region/Region Clusters mean RNA expression.csv")

cluster.averages <- AverageExpression(Epithelial)
head(cluster.averages[["SCT"]][1:10,])
write.csv(cluster.averages$SCT, "./Results/Cluster Averages_Region/Region Clusters mean SCT expression.csv")

# Write capitalized genes for comparing to humans
df <- cluster.averages$RNA
rownames(df) <-   str_to_upper(rownames(df))
write.csv(df, "./Results/Cluster Averages_Region/Region Clusters mean RNA expression_CAPS.csv")

df <- cluster.averages$SCT
rownames(df) <-   str_to_upper(rownames(df))
write.csv(df, "./Results/Cluster Averages_Region/Region Clusters mean SCT expression_CAPS.csv")
rm(df)

saveRDS(all.combined, file = "../ALL.assigned_r1.rds")

sessionInfo()

