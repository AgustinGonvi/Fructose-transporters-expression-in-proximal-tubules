# Load required libraries
library(Seurat)
library(cowplot)
library(stringr)
library(ggplot2)
library(tidyverse)
library(DoubletFinder)
library(gridExtra)

# Set seed for reproducibility
saved.seed <- 6789

# Define paths to data directories
path_General <- "G:\\Shared drives\\Gonzalez-Vicente public data\\Paper0501_Caloric Restriction\\00_Data and QC"
path_KMY  <- paste(path_General, "1_GSM4331828_Kidney-M-Y\\00_Matrices", sep = "\\") %>% print()
path_KMO  <- paste(path_General, "2_GSM4331829_Kidney-M-O\\00_Matrices", sep = "\\") %>% print()
path_KMCR <- paste(path_General, "3_GSM4331830_Kidney-M-CR\\00_Matrices", sep = "\\") %>% print()
path_KFY  <- paste(path_General, "4_GSM4331831_Kidney-F-Y\\00_Matrices", sep = "\\") %>% print()
path_KFO  <- paste(path_General, "5_GSM4331832_Kidney-F-O\\00_Matrices", sep = "\\") %>% print()
path_KFCR <- paste(path_General, "6_GSM4331833_Kidney-F-CR\\00_Matrices", sep = "\\") %>% print()

# Create Seurat objects for each dataset
KMY <- CreateSeuratObject(Read10X(data.dir = path_KMY), min.cells = 3, min.features = 200)
KMY$Sample <- "Young.M.AL"
KMY$Sex <- "Male"
KMY$Age <- "Young"
KMY$Diet <- "Ad.Libitum"

KMO <- CreateSeuratObject(Read10X(data.dir = path_KMO), min.cells = 3, min.features = 200)
KMO$Sample <- "Old.M.AL"
KMO$Sex <- "Male"
KMO$Age <- "Old"
KMO$Diet <- "Ad.Libitum"

KMCR <- CreateSeuratObject(Read10X(data.dir = path_KMCR), min.cells = 3, min.features = 200)
KMCR$Sample <- "Old.M.CR"
KMCR$Sex <- "Male"
KMCR$Age <- "Old"
KMCR$Diet <- "Restricted"

KFY <- CreateSeuratObject(Read10X(data.dir = path_KFY), min.cells = 3, min.features = 200)
KFY$Sample <- "Young.F.AL"
KFY$Sex <- "Female"
KFY$Age <- "Young"
KFY$Diet <- "Ad.Libitum"

KFO <- CreateSeuratObject(Read10X(data.dir = path_KFO), min.cells = 3, min.features = 200)
KFO$Sample <- "Old.F.AL"
KFO$Sex <- "Female"
KFO$Age <- "Old"
KFO$Diet <- "Ad.Libitum"

KFCR <- CreateSeuratObject(Read10X(data.dir = path_KFCR), min.cells = 3, min.features = 200)
KFCR$Sample <- "Old.F.CR"
KFCR$Sex <- "Female"
KFCR$Age <- "Old"
KFCR$Diet <- "Restricted"

# Merge all datasets into a single Seurat object
all.merged <- merge(KMY, c(KMO, KMCR, KFY, KFO, KFCR))

# Display dimensions and sample table of the merged dataset
dim(all.merged)
table(all.merged$Sample, all.merged$orig.ident)

# Clean up environment by removing path and Seurat object variables
rm(path_KFCR, path_KFO, path_KFY, path_KMCR, path_KMO, path_KMY)
rm(KMY, KMO, KMCR, KFY, KFO, KFCR)

# Perform garbage collection
gc(verbose = T)

# Identify mitochondrial genes
grep("^Mt-", rownames(all.merged@assays$RNA@counts), value = T)

# Calculate percentage of mitochondrial genes
all.merged <- PercentageFeatureSet(all.merged, "^Mt-", col.name = "percent.MT")
head(all.merged$percent.MT)

# Identify ribosomal proteins
grep("^Rp[sl]", rownames(all.merged@assays$RNA@counts), value = T)

# Calculate percentage of ribosomal proteins
all.merged <- PercentageFeatureSet(all.merged, "^Rp[sl]", col.name = "percent.Ribosomal")
head(all.merged$percent.Ribosomal)

# Exclude "Malat1" from analysis
all.merged[rownames(all.merged) != "Malat1",] -> SObject.F

# Find the maximum count and its index for each cell
apply(SObject.F@assays$RNA@counts, 2, max) -> SObject.F$largest_count
apply(SObject.F@assays$RNA@counts, 2, which.max) -> SObject.F$largest_index
rownames(SObject.F)[SObject.F$largest_index] -> SObject.F$largest_gene

# Calculate percentage of the largest gene
100 * SObject.F$largest_count / SObject.F$nCount_RNA -> SObject.F$percent.Largest.Gene

# Update all.merged with the largest gene and its percentage
SObject.F$largest_gene -> all.merged$largest_gene
SObject.F$percent.Largest.Gene -> all.merged$percent.Largest.Gene

# Remove SObject.F to free up memory
rm(SObject.F)

# Plot violin plots for various metrics
VlnPlot(all.merged, features=c("nFeature_RNA", "nCount_RNA", "percent.MT", "percent.Ribosomal", "percent.Largest.Gene"), ncol = 3)
VlnPlot(all.merged, features=c("percent.Largest.Gene"), ncol = 1) + scale_y_log10()
VlnPlot(all.merged, features=c("nFeature_RNA", "nCount_RNA", "percent.MT", "percent.Ribosomal", "percent.Largest.Gene"), group.by = "Sample", ncol = 3)
VlnPlot(all.merged, features=c("percent.Largest.Gene"), group.by = "Sample", ncol = 1) + scale_y_log10()

# Generate feature scatter plots
grid.arrange(
  FeatureScatter(all.merged, feature1 = "nCount_RNA", feature2 = "percent.Largest.Gene", group.by = "orig.ident", pt.size = 0.5),
  FeatureScatter(all.merged, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5),
  ncol = 2
)

# Generate QC metrics and visualize them
qc.metrics <- as_tibble(all.merged[[]], rownames="Cell.Barcode")
head(qc.metrics)
qc.metrics %>%
  arrange(percent.MT) %>%
  ggplot(aes(nCount_RNA, nFeature_RNA, colour=percent.MT)) +
  geom_point() +
  scale_color_gradientn(colors=c("black", "blue", "green2", "red", "yellow")) +
  ggtitle("QC metrics") +
  geom_hline(yintercept = 500) +
  geom_hline(yintercept = 5000)

# Generate more detailed QC plots
qc.metrics %>%
  arrange(percent.MT) %>%
  ggplot(aes(nCount_RNA, nFeature_RNA, colour=percent.MT)) +
  geom_point(size=0.7) +
  scale_color_gradientn(colors=c("black", "blue", "green2", "red", "yellow")) +
  ggtitle("QC metrics") +
  geom_hline(yintercept = 500) +
  geom_hline(yintercept = 5000) +
  scale_y_continuous() +
  scale_x_log10()

# Analyze gene expression complexity
qc.metrics %>%
  mutate(complexity=log10(nFeature_RNA) / log10(nCount_RNA)) -> qc.metrics
lm(log10(qc.metrics$nFeature_RNA)~log10(qc.metrics$nCount_RNA)) -> complexity.lm
complexity.lm

# Calculate complexity difference
qc.metrics %>%
  mutate(
    complexity_diff = log10(nFeature_RNA) - ((log10(qc.metrics$nCount_RNA)*complexity.lm$coefficients[2])+complexity.lm$coefficients[1])
  ) -> qc.metrics
qc.metrics %>%
  ggplot(aes(x=complexity_diff)) +
  geom_density(fill="yellow")

# Adjust complexity scale and visualize
min(c(max(qc.metrics$complexity_diff),0-min(qc.metrics$complexity_diff))) -> complexity_scale
qc.metrics %>%
  mutate(complexity_diff=replace(complexity_diff, complexity_diff< -0.1, -0.1)) %>%
  ggplot(aes(x=log10(nCount_RNA), y=log10(nFeature_RNA), colour=complexity_diff)) +
  geom_point(size=0.5) +
  geom_abline(slope=complexity.lm$coefficients[2], intercept = complexity.lm$coefficients[1]) +
  scale_colour_gradient2(low="blue2", mid="grey", high="red2")

# Further analysis of gene expression
qc.metrics %>%
  ggplot(aes(x=complexity_diff, y=percent.Largest.Gene)) +
  geom_point()
qc.metrics %>%
  group_by(largest_gene) %>%
  count() %>%
  arrange(desc(n)) -> largest_gene_list
largest_gene_list

# Plot genes with highest expression
largest_gene_list %>%
  filter(n>140) %>%
  pull(largest_gene) -> largest_genes_to_plot
qc.metrics %>%
  filter(largest_gene %in% largest_genes_to_plot) %>%
  mutate(largest_gene=factor(largest_gene, levels=largest_genes_to_plot)) %>%
  arrange(largest_gene) %>%
  ggplot(aes(x=log10(nCount_RNA), y=log10(nFeature_RNA), colour=largest_gene)) +
  geom_point(size=1) +
  scale_colour_manual(values=c("grey", RColorBrewer::brewer.pal(9,"Set1")))

# Visualize QC metrics by mitochondrial content
qc.metrics %>%
  arrange(percent.MT) %>%
  ggplot(aes(x=complexity_diff, y=percent.Largest.Gene, colour=percent.MT)) +
  geom_point() +
  scale_colour_gradient(low="grey", high="red2")

# Visualize QC metrics by ribosomal content
qc.metrics %>%
  arrange(percent.Ribosomal) %>%
  ggplot(aes(x=complexity_diff, y=percent.Largest.Gene, colour=percent.Ribosomal)) +
  geom_point() +
  scale_colour_gradient(low="grey", high="red2")

# Histogram of mitochondrial content distribution
qc.metrics %>%
  ggplot(aes(percent.MT)) +
  geom_histogram(binwidth=0.5, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage Mitochondrion") +
  geom_vline(xintercept=10)

# Histogram of largest gene percentage distribution
qc.metrics %>%
  ggplot(aes(percent.Largest.Gene)) +
  geom_histogram(binwidth=0.7, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage Largest Gene") +
  geom_vline(xintercept=10)

# Filter all.merged based on specified criteria
all.filtered <- subset(all.merged,
                       nFeature_RNA > 560 & nFeature_RNA < 4500 &
                       nCount_RNA < 30000 &
                       percent.MT < 40 &
                       percent.Ribosomal < 30 &
                       percent.Largest.Gene < 25)

# Display the filtered dataset
all.filtered

# Violin plots for filtered data
VlnPlot(all.filtered,
        features=c("nFeature_RNA", "nCount_RNA", "percent.MT", "percent.Ribosomal", "percent.Largest.Gene"),
        ncol=3)
VlnPlot(all.filtered,
        features=c("percent.Largest.Gene"),
        ncol=1) + scale_y_log10()
VlnPlot(all.filtered,
        features=c("nFeature_RNA", "nCount_RNA", "percent.MT", "percent.Ribosomal", "percent.Largest.Gene"),
        group.by="Sample",
        ncol=3)
VlnPlot(all.filtered,
        features=c("percent.Largest.Gene"),
        group.by="Sample",
        ncol=1) + scale_y_log10()

# Feature scatter plots for filtered data
grid.arrange(
  FeatureScatter(all.filtered, feature1="nCount_RNA", feature2="percent.Largest.Gene", group.by="Sample", pt.size=1),
  FeatureScatter(all.filtered, "nCount_RNA", "nFeature_RNA", group.by="Sample", pt.size=1),
  ncol=2
)

# Compare violin plots before and after filtering
VlnPlot(all.merged,
        features=c("nFeature_RNA", "nCount_RNA", "percent.MT", "percent.Ribosomal", "percent.Largest.Gene"),
        ncol=5)
VlnPlot(all.filtered,
        features=c("nFeature_RNA", "nCount_RNA", "percent.MT", "percent.Ribosomal", "percent.Largest.Gene"),
        ncol=5)

# Feature scatter plots comparison before and after filtering
grid.arrange(
  FeatureScatter(all.merged, feature1="nCount_RNA", feature2="percent.Largest.Gene", group.by="orig.ident", pt.size=0.5),
  FeatureScatter(all.filtered, feature1="nCount_RNA", feature2="percent.Largest.Gene", group.by="orig.ident", pt.size=0.5),
  FeatureScatter(all.merged, "nCount_RNA", "nFeature_RNA", group.by="orig.ident", pt.size=0.5),
  FeatureScatter(all.filtered, "nCount_RNA", "nFeature_RNA", group.by="orig.ident", pt.size=0.5),
  ncol=2
)

# Normalize filtered data
all.filtered <- NormalizeData(all.filtered, normalization.method="LogNormalize", verbose=T)

# Analyze gene expression
gene.expression <- apply(all.filtered@assays$RNA@data, 1, mean)
gene.expression <- sort(gene.expression, decreasing=T)
head(gene.expression, n=50)

# Histogram of Gapdh expression
ggplot(mapping=aes(all.filtered@assays$RNA@data["Gapdh",])) +
  geom_histogram(binwidth=0.05, fill="yellow", colour="black") +
  ggtitle("Gapdh expression")

# Seurat cell cycle stage marker genes
cc.genes.updated.2019
# Non-human format
str_to_title(cc.genes.updated.2019)

# Cell cycle scoring
all.filtered <- CellCycleScoring(all.filtered,
                                 s.features=str_to_title(cc.genes.updated.2019$s.genes),
                                 g2m.features=str_to_title(cc.genes.updated.2019$g2m.genes),
                                 set.ident=T)

# Display cell cycle scores
all.filtered[[]]
as_tibble(all.filtered[[]]) %>%
  ggplot(aes(Phase, fill=Sample)) +
  geom_bar()
as_tibble(all.filtered[[]]) %>%
  ggplot(aes(Sample, fill=Phase)) +
  geom_bar()
as_tibble(all.filtered[[]]) %>%
  ggplot(aes(x=S.Score, y=G2M.Score, color=Phase)) +
  geom_point(size=1)

# Continue analysis with further data filtering and exploration...
all.filtered

# Filter based on cell cycle scores
all.filtered <- subset(all.filtered, S.Score < 0.2 & G2M.Score < 0.2)
all.filtered

# Find variable features
all.filtered <- FindVariableFeatures(all.filtered, selection.method="vst", nfeatures=500)
head(VariableFeatures(all.filtered), 15)
top15 <- head(VariableFeatures(all.filtered), 15)

# Plot variable features with and without labels
plot1 <- VariableFeaturePlot(all.filtered)
plot2 <- LabelPoints(plot=plot1, points=top15, repel=T)
plot2

# Clean up environment
rm(top15, plot1, plot2)

# Analyze gene variance
variance.data <- as_tibble(HVFInfo(all.filtered), rownames="Gene")
variance.data <- variance.data %>%
  mutate(hypervariable=Gene %in% VariableFeatures(all.filtered))
head(variance.data, n=10)
variance.data %>%
  ggplot(aes(log(mean), log(variance), color=hypervariable)) +
  geom_point(size=0.5) +
  scale_color_manual(values=c("black", "red"))

# Scale data
all.filtered <- ScaleData(all.filtered)
# Optionally, regress out variables
# all.filtered <- ScaleData(all.filtered, vars.to.regress=c("S.Score", "G2M.Score", "Age", "Diet", "Sex"))

# Run PCA
all.filtered <- RunPCA(all.filtered, nfeatures.print=10)
ElbowPlot(all.filtered, ndims=50, reduction="pca")

# Dimensional reduction plots
grid.arrange(DimPlot(all.filtered, reduction="pca"),
             DimPlot(all.filtered, reduction="pca", group.by="Sex"),
             DimPlot(all.filtered, reduction="pca", group.by="Diet"),
             DimPlot(all.filtered, reduction="pca", group.by="Age"),
             ncol=2, nrow=2)

# Dimensional heatmap
DimHeatmap(all.filtered, dims=1:3, cells=10000)

# t-SNE analysis
all.filtered <- RunTSNE(all.filtered, dims=1:40, seed.use=saved.seed, perplexity=30)
grid.arrange(DimPlot(all.filtered, reduction="tsne"),
             DimPlot(all.filtered, reduction="tsne", group.by="Sex"),
             DimPlot(all.filtered, reduction="tsne", group.by="Diet"),
             DimPlot(all.filtered, reduction="tsne", group.by="Age"),
             ncol=2)

# UMAP analysis
all.filtered = RunUMAP(all.filtered, dims=1:40, verbose=T)
grid.arrange(DimPlot(all.filtered, reduction="umap"),
             DimPlot(all.filtered, reduction="umap", group.by="Sex"),
             DimPlot(all.filtered, reduction="umap", group.by="Diet"),
             DimPlot(all.filtered, reduction="umap", group.by="Age"),
             ncol=2)

# Doublet detection
sweep.list <- paramSweep_v3(all.filtered, PCs=1:40, sct=F, num.cores=1)
sweep.stats <- summarizeSweep(sweep.list, GT=F)
Find.pK.stats <- find.pK(sweep.stats)
barplot(Find.pK.stats$BCmetric, names.arg=Find.pK.stats$pK, las=2)

# Set parameters for doublet detection
nPCs <- 1:30                                 # Relevant principal components
nExp <- round(ncol(all.filtered) * 0.05)     # Expect 5% doublets
pK <- 0.27                                   # Optimal pK from previous plot

# Run DoubletFinder
all.filtered <- doubletFinder_v3(all.filtered, pN=0.25, PCs=nPCs, pK=pK, nExp=nExp)

# Extract the name of the DoubletFinder classification column
DF.name = colnames(all.filtered@meta.data)[grepl("DF.classification", colnames(all.filtered@meta.data))]

# UMAP plots highlighting doublets
grid.arrange(DimPlot(all.filtered, reduction="umap", group.by="orig.ident"),
             DimPlot(all.filtered, reduction="umap", group.by=DF.name),
             ncol=2)

# Violin plot for feature RNA by DoubletFinder classification
VlnPlot(all.filtered, features="nFeature_RNA", group.by=DF.name, pt.size=0.1)

# Display dimensions of the filtered dataset
dim(all.filtered)

# Filter out doublets
all.filtered = all.filtered[, all.filtered@meta.data[, DF.name] == "Singlet"]

# Display dimensions of the singlet dataset
dim(all.filtered)

# Save the filtered dataset
saveRDS(all.filtered, file=paste("ALL animals QC", ".rds", sep=""))

# Display session information
sessionInfo()

# Additional UMAP plots
grid.arrange(DimPlot(all.filtered, reduction="umap", group.by="orig.ident"),
             DimPlot(all.filtered, reduction="umap", group.by=DF.name),
             ncol=2)

# Additional UMAP plots for detailed analysis
grid.arrange(DimPlot(all.filtered, reduction="umap", group.by="orig.ident"),
             DimPlot(all.filtered, reduction="umap", group.by=DF.name),
             ncol=2)

# Final UMAP plot for presentation
grid.arrange(DimPlot(all.filtered, reduction="umap", group.by="orig.ident"),
             DimPlot(all.filtered, reduction="umap", group.by=DF.name),
             ncol=2)
