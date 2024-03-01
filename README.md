# Profiling Cell Heterogeneity and Fructose Transporter Expression in the Rat Nephron

## Authors
- Ronghao Zhang
- Darshan Aatmaram Jadhav
- Benjamin Kramer
- Agustin Gonzalez-Vicente* (Corresponding Author)
- and the Kidney Precision Medicine Project

## Citation
If you use any of the code or workflows in this repository, please cite our manuscript in xxx (will update).

## Welcome
Welcome to our GitHub repository! Here you will find analysis scripts for our manuscript. Please contact any corresponding authors Dr. Agustin Gonzalez-Vicente.

## Lab Website
[Visit the Gonzalez-Vicente lab website](https://physiology.case.edu/research/labs/gonzalez-vicente-lab/)

## How to Use Our scRNA Script
1. **Preprocess the filtered matrices** using the `01_scRNAseq_QC.R` file. Remember to change the path for the matrices and save the `.Rds` file properly for further analysis.
2. **Perform Comprehensive analysis** of the rat scRNA data. In the `02_scRNAseq_Integration.R` code, load the path for the preprocessed Seurat object from step 1 and then continue with the analysis.
3. **Use the `03_scRNAseq_ClusterAssignment.R` code** to categorize cells into groups like Endothelial, Epithelial, Immune, and Stromal based on their cluster IDs and further analysis.