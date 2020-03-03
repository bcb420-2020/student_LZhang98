# A script that performs the data processing and cleaning in A1, but without R
# Markdown documentation. Intended to be sourced in A2_LukeZhang.Rmd.

# The result of this script is a dataframe with unique HGNC symbols as rownames
# representing our cleaned, processed, and normalized data.

# Download packages ===========================================================

if (! requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
if (! requireNamespace("GEOquery", quietly = TRUE)) {
    BiocManager::install("GEOquery")
}
if (! requireNamespace("edgeR", quietly = TRUE)) {
    BiocManager::install("edgeR")
}
if (! requireNamespace("limma", quietly = TRUE)) {
    BiocManager::install("limma")
}
if (! requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    BiocManager::install("org.Hs.eg.db")
}
if (! requireNamespace("annotate", quietly = TRUE)) {
    BiocManager::install("annotate")
}
library(org.Hs.eg.db)
# =============================================================================

# Download dataset GSE110021 from GEO
datafile <- "GSE110021_counts.Aug2015.txt.gz"
if (!file.exists(datafile)) {
    files <- GEOquery::getGEOSuppFiles("GSE110021", makeDirectory = FALSE)
}
data <- read.delim(datafile,header=TRUE,check.names=FALSE)
dim(data)
data[1:10,1:4]


# RNAseq protocol has 4 lanes per sample. Merge these lanes by adding counts
# together for each gene
new_data <- list()
for (i in 0:11) {
    new_data <- c(new_data, list(rowSums(data[,(i+1):(4*(i+1))])))
}
new_data <- as.data.frame(new_data)
rownames(new_data) <- rownames(data)

# data colnames format:
# <dayNum>.<treatment>.<replicateNum>_<sampleNum>_<laneNum>
for (i in 1:12) {
    colnames(new_data)[i] <- unlist(strsplit(colnames(data[i*4]), "_"))[1]
}


# Define 4 groups within our dataset model:
# 1. no TGFb, day 1
# 2. no TGFb, day 20
# 3. with TGFb, day 1
# 4. with TGFb, day 20
samples <- data.frame(lapply(colnames(new_data), FUN=function(x){unlist(strsplit(x,"\\."))[c(1,2)]}))
rownames(samples) <- c("Day", "Treatment")
colnames(samples) <- colnames(new_data)
samples <- data.frame(t(samples))


# Filter out low counts
cpm_data <- edgeR::cpm(new_data)
rownames(cpm_data) <- rownames(new_data)
filtered_data <- new_data[rowSums(cpm_data >= 1) >= 3, ]


# Normalize the data with edgeR TMM with respect to treatment type
filtered_matrix <- as.matrix(filtered_data)
d <- edgeR::DGEList(counts=filtered_matrix, group=samples$Treatment)
d <- edgeR::calcNormFactors(d)
norm_data <- edgeR::cpm(d)


# Map Ensemble Gene IDs to HGNC symbols
symbols <- annotate::lookUp(rownames(filtered_data), "org.Hs.eg", "SYMBOL")

# Remove all unmapped rows
mapped_indices <- which(!is.na(symbols))
mapped_data <- filtered_data[mapped_indices, ]

# Set HUGO symbols as rownames of new dataframe
mapped_symbols <- symbols[!is.na(symbols)]
rownames(mapped_data) <- mapped_symbols

# mapped_data is our result data frame. Use for A2.