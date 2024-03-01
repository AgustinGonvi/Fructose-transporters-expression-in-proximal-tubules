# Cleaned R Code from R Markdown

library(Seurat)
library(cowplot)
library(stringr)
library(ggplot2)

saved.seed <- 5678

path_WD <- "G:/Shared drives/Gonzalez-Vicente Lab/03_Data/2023_Rat scRNAseq/01_Resolution 1.6"

path_Data <- "G:\\Shared drives\\Gonzalez-Vicente public data\\Paper0501_Caloric Restriction\\00_Data and QC\\ALL animals QC.rds"

print(paste("Working directory :", path_WD, sep = " "))
print(paste("Data directory :", path_Data, sep = " "))

ALL.animals <- readRDS(file = path_Data)


object.list <- SplitObject(ALL.animals,
                           split.by = "Sample")

object.list <- lapply(X = object.list,
                      FUN = SCTransform)

features <- SelectIntegrationFeatures(object.list = object.list,
                                      nfeatures = 3000)


object.list <- PrepSCTIntegration(object.list = object.list,
                                  anchor.features = features)

all.anchors <- FindIntegrationAnchors(object.list = object.list,
                                      normalization.method = "SCT",
                                      anchor.features = features)

all.combined <- IntegrateData(anchorset = all.anchors,
                              normalization.method = "SCT")

Resolution <- 1.6

all.combined <- RunPCA(all.combined, npcs = 30, verbose = FALSE)
all.combined <- RunTSNE(all.combined, dims = 1:30, seed.use = saved.seed, perplexity = 30)
all.combined <- RunUMAP(all.combined, reduction = "pca", dims = 1:30)
all.combined <- FindNeighbors(all.combined, reduction = "pca", dims = 1:30)
all.combined <- FindClusters(all.combined, resolution = Resolution)

p1 <- DimPlot(object = all.combined,
              reduction = "umap",
              group.by = "Sample",
              pt.size = 0.01)
p1
ggsave(filename = paste("./Results/01_UMAP_bySample_Res", Resolution, ".tiff", sep = ""),
       device = "tiff",
       dpi = 600)

p2 <- DimPlot(object = all.combined,
              reduction = "umap",
              split.by = "Sample")
p2
ggsave(filename = paste("./Results/02_UMAP_Sample_Res", Resolution, ".tiff", sep = ""),
       device = "tiff",
       dpi = 600,
       width = 50, height = 20, units = "cm")
       
p3 <- DimPlot(object = all.combined,
              reduction = "umap",
              label = T)
p3
ggsave(filename = paste("./Results/03_UMAP_Labels_Res", Resolution, ".tiff", sep = ""),
       device = "tiff",
       dpi = 600)    

rm(p1, p2, p3)

p1 <- DimPlot(object = all.combined,
              reduction = "pca",
              group.by = "Sample",
              pt.size = 0.01)
p1
ggsave(filename = paste("./Results/04_PCA_bySample_Res", Resolution, ".tiff", sep = ""),
       device = "tiff",
       dpi = 600)

p2 <- DimPlot(object = all.combined,
              reduction = "pca",
              split.by = "Sample")
p2
ggsave(filename = paste("./Results/05_PCA_Sample_Res", Resolution, ".tiff", sep = ""),
       device = "tiff",
       dpi = 600,
       width = 50, height = 20, units = "cm")
       
p3 <- DimPlot(object = all.combined,
              reduction = "pca",
              label = T)
p3
ggsave(filename = paste("./Results/06_PCA_Labels_Res", Resolution, ".tiff", sep = ""),
       device = "tiff",
       dpi = 600)    

rm(p1, p2, p3)

p1 <- DimPlot(object = all.combined,
              reduction = "tsne",
              group.by = "Sample",
              pt.size = 0.01)
p1
ggsave(filename = paste("./Results/07_tsne_bySample_Res", Resolution, ".tiff", sep = ""),
       device = "tiff",
       dpi = 600)

p2 <- DimPlot(object = all.combined,
              reduction = "tsne",
              split.by = "Sample")
p2
ggsave(filename = paste("./Results/08_tsne_Sample_Res", Resolution, ".tiff", sep = ""),
       device = "tiff",
       dpi = 600,
       width = 50, height = 20, units = "cm")
       
p3 <- DimPlot(object = all.combined,
              reduction = "tsne",
              label = T)
p3
ggsave(filename = paste("./Results/09_tsne_Labels_Res", Resolution, ".tiff", sep = ""),
       device = "tiff",
       dpi = 600)    

rm(p1, p2, p3)

dir.create(file.path("./Results", "Cluster Markers"), showWarnings = T)

DefaultAssay(object = all.combined) <- "RNA"
clus = 0
while (clus < length(levels(all.combined@meta.data$seurat_clusters)))
{
    #print(paste("Cluster", clus, sep=" "))
    markers <- FindMarkers(all.combined,
                           ident.1 = clus,
                           logfc.threshold = 0.5,
                           min.pct = 0.10,
                           test.use = "wilcox",
                           only.pos = FALSE,
                           max.cells.per.ident = Inf,
                           random.seed = 1,
                           verbose = TRUE)
    write.csv(markers, file = paste("./Results/Cluster Markers/Cluster ", clus, "_Markers.csv", sep = ""))
    clus = clus + 1
    gc(verbose = T)
    rm(markers)
}

dir.create(file.path("./Results", "Cluster Averages"), showWarnings = T)
           
cluster.averages <- AverageExpression(all.combined)
head(cluster.averages[["RNA"]][1:10,])
write.csv(cluster.averages$RNA, "./Results/Cluster Averages/Numbers mean RNA expression.csv")

cluster.averages <- AverageExpression(all.combined)
head(cluster.averages[["SCT"]][1:10,])
write.csv(cluster.averages$SCT, "./Results/Cluster Averages/Numbers mean SCT expression.csv")

# Write capitalized genes for comparing to humans
df <- cluster.averages$SCT
rownames(df) <-   str_to_upper(rownames(df))
write.csv(df, "./Results/Cluster Averages/Numbers mean SCT expression_CAPS.csv")
rm(df, cluster.averages)

saveRDS(all.combined, file = "./Results/ALL.combined.rds")

sessionInfo()