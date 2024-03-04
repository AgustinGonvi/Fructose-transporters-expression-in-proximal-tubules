# Transporter Location - Differential Gene Expression (DESeq)
# Author: Ronghao Zhang
# Date: Current Date (Sys.Date() is used in Rmd to generate this dynamically)

# Preliminaries
# Packages for Preliminaries 
library(knitr)
library(rmdformats)
# Global options
options(max.print="75")
opts_chunk$set(comment=NA)
opts_knit$set(width=75)

# R Packages
# This code chunk loads up the necessary R packages to manage the data and finish Homework 01. 
library(DESeq2)
library(annotables)
library(tidyverse)

# Themes
# This code chunk sets the theme of plots. 
theme_set(theme_bw())
plot_theme <- theme(text=element_text(size=12, family="serif")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

# Initialization
nsample <- 21

# Add Columns to Metadata
# Load Metadata
metadata <- read.table("./../00_Inputs/00_Meta_Data.csv", header = TRUE, sep = ",")

# Differential Gene Expression
# Load normalized counts
# factorize metadata
metadata$PTSegment <- factor(metadata$PTSegment)

# load saved matrix from SubRead
raw_count <- read.table("./../00_Inputs/raw_count_matrix.csv", header = TRUE, sep = ",") 
colnames(raw_count)[1] <- "ratGene"

# get the gene names from there
gene_name <- raw_count[,1]

# assign the gene names back to dataframe
raw_count <- raw_count %>% select(-ratGene)
rownames(raw_count) <- gene_name

# Import Matrix to DeseqDataSet Object
# check row names and column names are matched
rownames(metadata) <- metadata$ID
all(colnames(raw_count)==rownames(metadata))

# Create the DESeq Object for comparing s1 v.s. s2 and s3
dds_s1_base <- DESeqDataSetFromMatrix(countData = raw_count,
                                      colData = metadata,
                                      design = ~ PTSegment)

# Create the DESeq Object for comparing s2 v.s. s3
dds_s2_base <- dds_s1_base 
dds_s2_base$PTSegment <- relevel(dds_s2_base$PTSegment,"S2")

# Normalize and Transformation
# select the protein coding gene for rats
prot_code_gene <- rnor6 %>% filter(biotype == "protein_coding")

NormalizeDDS <- function(dds) {
  # normalize using size factors
  dds <- DESeq(dds)
  
  # normalize using variance
  dataExp_N <- counts(dds, normalized = TRUE)
  row_var <- apply(dataExp_N,1,var)
  dataExp_V <- as.data.frame(dataExp_N[row_var > 0.395,])
  
  # find the protein coding gene in the count matrix
  index <- match(prot_code_gene$symbol, rownames(dataExp_V))
  index <- index[!is.na(index)]
  dataExp_P <- dataExp_V[index,]
  
  # Obtain the indices of only desired genes
  genesNotToRemove <- which(rownames(dds) %in% rownames(dataExp_P))

  # Cut your desired genes in the DESeq object
  dds <- dds[genesNotToRemove, ]

  return(dds)
}

dds_s1_base <- NormalizeDDS(dds_s1_base)
dds_s2_base <- NormalizeDDS(dds_s2_base)
resultsNames(dds_s1_base)
resultsNames(dds_s2_base)

# Differential Gene Expression Results
res_s1_s2 <- results(dds_s1_base, 
                     contrast = c("PTSegment", "S2","S1"),
                     alpha = 0.05, lfcThreshold = 0)

res_s1_s3 <- results(dds_s1_base, 
                     contrast = c("PTSegment", "S3","S1"),
                     alpha = 0.05, lfcThreshold = 0)

res_s2_s3 <- results(dds_s2_base, 
                     contrast = c("PTSegment", "S3","S2"),
                     alpha = 0.05, lfcThreshold = 0)
summary(res_s1_s2)
summary(res_s1_s3)
summary(res_s2_s3)

# Export the Results to Matrix
S1_S2 <- data.frame(res_s1_s2)
S1_S3 <- data.frame(res_s1_s3)
S2_S3 <- data.frame(res_s2_s3)

write.table(S1_S2,"./../02_Results/DEG_S1_S2.csv",sep = ",",row.names = TRUE)
write.table(S1_S3,"./../02_Results/DEG_S1_S3.csv",sep = ",",row.names = TRUE)
write.table(S2_S3,"./../02_Results/DEG_S2_S3.csv",sep = ",",row.names = TRUE)

# Export the Results for Log Transformation
# Export the metadata with Size Factors
metadata$SizeFactor <- sizeFactors(dds_s1_base)
write.csv(metadata,"./../02_Results/MetaData_With_Size_Factor.csv",row.names = FALSE)
# export the normalized matrix for log transformation
norm_count <- as.data.frame(counts(dds_s1_base, normalized=TRUE)) 
save(norm_count,file = "./../02_Results/Normalized_Filtered_Counts.RData")

# Session Info
sessionInfo()
