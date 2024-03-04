# 04 Renal Segment Preprocessing: Mean and Median per Segment
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
# This code chunk loads up the necessary R packages.
library(tidyverse)

# Themes
# This code chunk sets the theme of plots. 
theme_set(theme_bw())
plot_theme <- theme(text=element_text(size=12, family="serif")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

# Import and Tidy the Data

## Import the Data
# Import the count and mata data matrix. The `raw_count` variable is generated using `rSubread`, 
# while the `meta_data` variable stores the index of each files, this is downloaded directly from ENA website (Project Number: PRJNA244440).
raw_count <- read.table("./../00_Inputs/count_matrix.csv", header = TRUE, sep = ",")
meta_data <- read.table("./../00_Inputs/PRJNA244440_metadata.txt", header = TRUE, sep = "\t")

## Pre-Processing 
# Order the metadata and update the rownames
meta_data <- dplyr::arrange(meta_data, run_accession) %>% 
  select(-sample_accession)

# Update the rownames of the count matrix
row.names(raw_count) <- raw_count$X
raw_count <- raw_count %>% select(-X)

# Log2 transform the data
log_count <- log2(raw_count+1)

# Check whether the run accession number match the metadata number
FALSE %in% (check1 <- colnames(raw_count) %in% meta_data$run_accession)

# Process the Data
# This section generates 2 separate tables: [1] mean by segment [2] median by segment.
# Transpose the log count matrix
log_count_t <- t(log_count)

# Calculate the mean by segment
mean_by_segment <- aggregate(log_count_t[,], list(meta_data$sample_title), mean) %>% 
  t() %>% as.data.frame()

# Calculate the median by segment
median_by_segment <- aggregate(log_count_t[,], list(meta_data$sample_title), median) %>%
  t() %>% as.data.frame()

# Export CSV Files
# Assign the column names and delete the first row (duplicated column name)
colnames(mean_by_segment) <- mean_by_segment[1,]
colnames(median_by_segment) <- median_by_segment[1,]

write.csv(median_by_segment[2:dim(median_by_segment)[1],], "./../02_Results/median_by_segment.csv", row.names = TRUE)
write.csv(mean_by_segment[2:dim(mean_by_segment)[1],], "./../02_Results/mean_by_segment.csv", row.names = TRUE)

# Session Info
sessionInfo()
