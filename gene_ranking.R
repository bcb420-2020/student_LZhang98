# A script that performs the QLF testing and gene ranking in A2
# but without R Markdown documentation. Intended to be sourced in A2_LukeZhang.Rmd.

# The result of this script is a sorted list of differentially expressed genes, with
# upregulated genes at the top and downregulated genes at the bottom.

# Download packages ===========================================================

if (! requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

if (! requireNamespace("edgeR", quietly = TRUE)) {
    BiocManager::install("edgeR")
}
library(edgeR)

# PROCEDURE ===================================================================

# Construct a model and run differential analysis.

d <- DGEList(counts=mapped_data, group=samples$Treatment)
# construct model
design <- model.matrix(~ samples$Day + samples$Treatment)

# estimate dispersion
d <- estimateDisp(d, design)
# calculate normalization factors
d <- calcNormFactors(d)
# fit model to data
fit <- glmQLFit(d, design)
# calculate differential expression
qlf_test <- glmQLFTest(fit, coef="samples$TreatmentTGFb")

# Get the results
qlf_output <- topTags(qlf_test, sort.by="p.value", n=nrow(mapped_data))

ranked_genes <- qlf_output