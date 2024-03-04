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
# Loading the normalized and filtered counts data.
load(file = "./../00_Inputs/Normalized_Filtered_Counts.RData")

# Log2 Transformation
# Applying log2 transformation to the normalized counts and saving the result.
norm_dist_count <- log2(norm_count+1)
write.csv(norm_dist_count,"./../02_Results/norm_dist_counts.csv",row.names = TRUE)

# Session Info
# Displaying session info for reproducibility.
sessionInfo()
