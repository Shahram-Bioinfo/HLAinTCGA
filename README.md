
# HLAinTCGA

The R package is specifically designed to analyze and visualize **HLA allele and genotype frequencies**, with a focus on **class I HLA alleles**, across various projects in *The Cancer Genome Atlas (TCGA)*. This package not only computes the **frequency of HLA alleles and genotypes**, but also incorporates the calculation of the **mean expression levels** for each specific allele. By leveraging these tools, researchers can explore the distribution of HLA alleles and genotypes along with their expression patterns, providing insights into the role of these alleles in cancer biology.

Additionally, the package offers advanced functionality to **stratify data by key clinical and demographic subgroups**, including **gender, race, ethnicity, clinical stage, and other relevant variables**. This stratification enables users to examine patterns and correlations within specific subgroups, facilitating a deeper understanding of the interplay between HLA diversity and clinical outcomes in cancer.

To enhance usability, the package includes features for **comprehensive and interactive visualizations**, allowing researchers to easily interpret and present their findings. With minimal manual operations, users can efficiently query and visualize allele distributions, genotype frequencies, and expression levels across TCGA projects. The intuitive workflows and visual outputs are designed to help researchers quickly achieve a clear understanding of HLA allele distributions, making this package a valuable tool for exploring the immunogenetic landscape of cancer.

# Installation

You can install the latest version of **HLAinTCGA** directly from GitHub:
# Install the devtools package if you haven't already
```r
install.packages("devtools")
```

# Install HLAinTCGA from GitHub
```r
devtools::install_github("Shahram-Bioinfo/HLAinTCGA")
```

# Usage
**Loading the Package**
```r
library(HLAinTCGA)
```

# Functions Overview:
**genotype_freq**: The genotype_freq function computes and visualizes the frequency of all HLA genotypes associated with a specific allele (e.g. **A*02:01**)across TCGA projects or within a specific project. It can also compute the frequency of a specific genotype (e.g. **A*02:01** - **A*11:01**)if requested. This function helps identify allele-genotype relationships and their distribution patterns across different cancer types or within a particular project.

# Example Usage
```r
# Load required libraries and dataset
library(HLAinTCGA)

# Load the dataset
data("data", package = "HLAinTCGA")

# Example 1: Analyze all genotypes containing the allele "A*01:01" across all TCGA projects
result_all_projects <- genotype_freq(allele1 =  "A*01:01", data = data)

# Example 2: Analyze all genotypes containing the allele "A*01:01" within the "TCGA-BRCA" project
result_filtered_project <- genotype_freq(allele1 = "A*01:01", data = data, project_filter = "TCGA-BRCA")

# Example 3: Analyze the genotype frequency for a specific pair of alleles ("A*01:01" and "A*03:03") across all TCGA projects
result_specific_genotype_all <- genotype_freq(allele1 = "A*01:01", allele2 = "A*03:03", data = data)

# Example 4: Analyze the genotype frequency for the same pair of alleles ("A*01:01" and "A*03:03") within the "TCGA-BRCA" project
result_specific_genotype_filtered <- genotype_freq(allele1 = "A*01:01", allele2 = "A*03:03", data = data, project_filter = "TCGA-BRCA")

# View results: Display the results tables for analysis
cat("\nResults Across All Projects:\n")
print(result_all_projects$results_table)

cat("\nResults for TCGA-BRCA Project:\n")
print(result_filtered_project$results_table)

# Plot results: Visualize the analysis
cat("\nPlot for Results Across All Projects:\n")
print(result_all_projects$plot)

cat("\nPlot for Results in TCGA-BRCA Project:\n")
print(result_filtered_project$plot)
```

**allele_freq**: Analyzes allele frequencies and allele-specific expression (ASE) levels  across TCGA projects.

# Example Usage
```r
# Load the required dataset from the package
data("data", package = "HLAinTCGA")

# Step 1: Analyze allele frequency for the HLA allele "A*01:01"
allele_of_interest <- "A*01:01"
result <- allele_freq(allele = allele_of_interest, data = data)

# Step 2: Display the results table
cat("\n--- Results Table for Allele:", allele_of_interest, "---\n")
print(result$results_table)

# Step 3: Visualize the allele frequency results
cat("\n--- Visualization of Allele Frequency for:", allele_of_interest, "---\n")
print(result$plot)
```

**allele_freq_gender**: Analyzes allele frequencies and expression levels stratified by gender across TCGA projects and in specific project(e.g. TCGA-BRCA)

**allele_freq_race**: Analyzes allele frequencies and expression levels stratified by race across TCGA projects and in specific project(e.g. TCGA-BRCA)

**allele_freq_stage**: Analyzes allele frequencies and expression levels stratified by AJCC stages across TCGA projects and in specific project(e.g. TCGA-BRCA)

**allele_freq_ethnicity**: Analyzes allele frequencies and expression levels stratified by ethnicity across TCGA projects and in specific project(e.g. TCGA-BRCA)

# Example Usage: Allele Frequency Analysis based Gender
```r
# Load the required dataset
data("data", package = "HLAinTCGA")

# Example 1: Analyze allele frequency for "A*01:01" stratified by gender across all TCGA projects
result_all_projects <- allele_freq_gender(allele = "A*01:01", data = data)

# Example 2: Analyze allele frequency for "A*01:01" stratified by gender within the "TCGA-BRCA" project
result_filtered_project <- allele_freq_gender(allele = "A*01:01", data = data, project_filter = "TCGA-BRCA")

# View results: Display the results tables for each analysis
cat("\nResults Across All Projects:\n")
print(result_all_projects$results_table)

cat("\nResults for TCGA-BRCA Project:\n")
print(result_filtered_project$results_table)

# Plot results: Visualize the analysis
cat("\nPlot for Results Across All Projects:\n")
print(result_all_projects$plot)

cat("\nPlot for Results in TCGA-BRCA Project:\n")
print(result_filtered_project$plot)
```

# Data
The package includes preprocessed TCGA data:

-Expression Levels: TPM values for HLA alleles.
-Demographic Data: Information such as ethnicity, gender, and race.
-Clinical Data: AJCC stage 

# Testing
The package includes unit tests for core functions. To run tests, use:
```r
devtools::test()
```

# Contribution
Contributions, suggestions, and bug reports are welcome! Please create a pull request or submit an issue on the GitHub repository.

# License
This project is licensed under the GPL-3 License - see the LICENSE file for details.





