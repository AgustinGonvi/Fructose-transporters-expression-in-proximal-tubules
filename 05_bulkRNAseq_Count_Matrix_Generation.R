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
library(Rsubread)
library(ggplot2)
library(tidyverse)

# Themes
# This code chunk sets the theme of plots. 
theme_set(theme_bw())
plot_theme <- theme(text=element_text(size=12, family="serif")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

# Generate Count Matrix
# This is the pipeline that is used to count the numbers read for each genome.

# Step 01: to index all the `bam` files in the directory. 
sorted_align_files <- list.files(path = "./../../../../mRNA_data_processed/bam/", 
                                 pattern = "align_sorted.bam", 
                                 full.names = TRUE)
length(sorted_align_files)

# Step 02: to count reads for all the `bam` files.
bam_counts <- featureCounts(files = sorted_align_files, 
                            annot.inbuilt = NULL, 
                            annot.ext = "./../00_Inputs/ref_genome_ann.gtf", 
                            isGTFAnnotationFile = TRUE, 
                            isPairedEnd = TRUE, 
                            nthreads = 16)

# Step 03: Export the count matrix 
count_mat <- as.data.frame(bam_counts$counts)
colnames(count_mat) <- sub("_align_sorted.bam", "", colnames(count_mat))
write.csv(count_mat,"./../02_Results/raw_count_matrix.csv",row.names = TRUE)

# Session Info
sessionInfo()
