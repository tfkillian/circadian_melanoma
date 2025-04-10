################################################################################
# 1) INSTALL AND LOAD REQUIRED PACKAGES
#    (Uncomment the install commands if packages are not installed already)
################################################################################

# install.packages("biomaRt")
# install.packages("sva")
# install.packages("limma")
# install.packages("RUVSeq")
# install.packages("edgeR")

library("biomaRt")
library("sva")
library("limma")
library("RUVSeq")
library("edgeR")
library("dplyr")
library("ExperimentHub")
library("depmap")
library("ggplot2")
library("reshape2")
library("pheatmap")

################################################################################
# 2) RETRIEVE CIRCADIAN GENES (GO:0007623) FROM ENSEMBL VIA biomaRt
################################################################################

# Connect to Ensembl (GRCh38) human genes
# mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# 
# # Query biomaRt for circadian genes (Gene Ontology GO:0007623)
# circadian_info <- getBM(
#   attributes = c("hgnc_symbol", "ensembl_gene_id"),
#   filters    = "go",
#   values     = "GO:0007623",
#   mart       = mart
# )
# 
# saveRDS(object = circadian_info, file = "~/tmp/circadian_melanoma/MAE/circadian_info.rds")

# Extract unique HGNC symbols for circadian genes
circadian_info <- readRDS(file = "~/tmp/circadian_melanoma/MAE/circadian_info.rds")
circadian_gene_symbols <- unique(circadian_info$hgnc_symbol)
cat("Number of circadian genes retrieved from GO:0007623:",
    length(circadian_gene_symbols), "\n")

################################################################################
# 3) LOAD GENE TRANSCRIPT COUNT DATA
#    - Assumes row 1 is a header with sample IDs
#    - Assumes first column is 'gene_name'
#    - Adjust 'sep' as needed if the file is truly CSV (sep=",")
################################################################################
## access depmap
eh <- ExperimentHub()
query(eh, "depmap")

# metadata
eh[["EH7558"]] %>%
  dplyr::filter(grepl("melanoma", lineage_subtype) |
                  grepl("elanoma", Cellosaurus_NCIt_disease) |
                  grepl("elanoma", subtype_disease) | 
                  grepl("elanoma", cell_line),
                !is.na(primary_or_metastasis)) %>%
  dplyr::select(-contains("issues")) %>%
  dplyr::mutate(primary_or_metastasis = case_when(
    primary_or_metastasis == "Primary" ~ 0,
    primary_or_metastasis == "Metastasis" ~ 1)) %>%
  as.data.frame() -> melanoma_data
rownames(melanoma_data) <- melanoma_data$depmap_id
# View(melanoma_data)

# tpm expression
eh[["EH7556"]] %>%
  dplyr::filter(depmap_id %in% melanoma_data$depmap_id) %>% ## only melanoma cell lines!
  dplyr::select(depmap_id, rna_expression, gene_name) %>%
  tidyr::pivot_wider(names_from = depmap_id,
                     values_from = rna_expression) %>%
  dplyr::filter(!duplicated(gene_name)) %>%
  tibble::column_to_rownames(var = "gene_name") %>%
  as.data.frame() -> count_data

################################################################################
# 4) LOAD METADATA
#    - Replace with your actual metadata file/path if different
#    - Make sure 'depmap_id' in metadata matches column names of melanoma_data
################################################################################

# Example (adjust path, delimiter, etc. to match your metadata file)
# meta_file <- "~/tmp/circadian_melanoma/MAE/metadata.csv"
# meta_data <- read.csv(meta_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# For illustration: We'll assume you have already loaded 'meta_data' into your R session.
# Ensure that 'depmap_id' in meta_data matches colnames(melanoma_data).
# Let's align the order of columns/samples in melanoma_data with meta_data:

# common_samples <- intersect(colnames(melanoma_data), meta_data$depmap_id)
# melanoma_data  <- melanoma_data[, common_samples]
# meta_data      <- meta_data[match(common_samples, meta_data$depmap_id), ]
# cat("After matching, dimensions of melanoma_data:", dim(melanoma_data), "\n")
# cat("Metadata rows after matching:", nrow(meta_data), "\n")

################################################################################
# 5) FILTER EXPRESSION MATRIX TO CIRCADIAN GENES (FOR COMPARISON PURPOSES)
################################################################################

circadian_data <- count_data[rownames(count_data) %in% circadian_gene_symbols, ]
cat("Dimensions of circadian_data:", dim(circadian_data), "\n")

################################################################################
# 6) IDENTIFY CIRCADIAN PATTERNS USING SVA (Surrogate Variable Analysis)
#    - Build a simple design matrix (e.g., ~ primary_or_metastasis)
#    - Then fit SVA to detect hidden variation (possibly circadian).
################################################################################

# Check that meta_data has a column "primary_or_metastasis" or similar:
# For example:
# design_mod <- model.matrix(~ primary_or_metastasis, data=meta_data)

# For demonstration, we use a placeholder design since code context may vary:
# If you only have "Primary" vs "Metastasis," ensure factors are set correctly.
# Suppose your metadata factor is meta_data$primary_or_metastasis:

# design_mod <- model.matrix(~ 0 + meta_data$primary_or_metastasis)
# We'll do a simpler placeholder below if the variable is not in environment:
# design_mod <- model.matrix(~ 1, data = count_data)   # Minimal design
meta_data <- melanoma_data[match(colnames(circadian_data), melanoma_data$depmap_id), ]
design_mod <- model.matrix(~ primary_or_metastasis, data = meta_data)
design_null <- model.matrix(~ 1, data = meta_data)

# Convert to matrix just in case it’s not:
expr_matrix <- as.matrix(circadian_data)

# Run sva
sva_results <- sva(dat = expr_matrix,
                   mod = design_mod,
                   mod0 = design_null)

# str(sva_results)

# Extract SVs from sva_results
sv_df <- as.data.frame(sva_results$sv)
colnames(sv_df) <- paste0("SV", 1:ncol(sv_df))  # Rename for clarity

# Add metadata for grouping
sv_df$group <- meta_data$primary_or_metastasis  # Assuming coded as 0 (Primary) and 1 (Metastasis)

# Scatter plot of SV1 vs SV2
ggplot(sv_df, aes(x = SV1, y = SV2, color = factor(group))) +
  geom_point(size = 3, alpha = 0.7) +
  labs(x = "Surrogate Variable 1", y = "Surrogate Variable 2", color = "Group") +
  theme_minimal() +
  ggtitle("Scatter Plot of SVA Surrogate Variables")

  # For reshaping data for ggplot

# Reshape for ggplot
sv_melt <- melt(sv_df, id.vars = "group")

# Boxplot by group
ggplot(sv_melt, aes(x = factor(group), y = value, fill = factor(group))) +
  geom_boxplot() +
  facet_wrap(~variable, scales = "free_y") +
  labs(x = "Group", y = "Surrogate Variable Value") +
  theme_minimal() +
  ggtitle("Boxplot of Surrogate Variables by Group")

# Heatmap of SVs
pheatmap(sva_results$sv,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         # annotation_col = meta_data["primary_or_metastasis"],
         main = "Heatmap of Surrogate Variables")

## Correlation Matrix Between SVs and Known Variables
# Correlation with known groups
cor_sv_meta <- cor(sva_results$sv, as.numeric(meta_data$primary_or_metastasis), use = "complete.obs")

# Bar plot of correlation strength
# barplot(cor_sv_meta,
#         main = "Correlation Between Surrogate Variables and Sample Group",
#         ylab = "Correlation Coefficient",
#         names.arg = paste0("SV", 1:length(cor_sv_meta)))

# PCA on Surrogate Variables
pca <- prcomp(sva_results$sv, scale. = TRUE)
pca_df <- as.data.frame(pca$x)

# Add metadata for coloring
pca_df$group <- meta_data$primary_or_metastasis

# PCA Plot
ggplot(pca_df, aes(x = PC1, y = PC2, color = factor(group))) +
  geom_point(size = 3, alpha = 0.7) +
  labs(x = "PC1", y = "PC2", color = "Group") +
  theme_minimal() +
  ggtitle("PCA of Surrogate Variables")

# Interpreting Surrogate Variables (SVs) from sva()
# 
# The sva() function generates surrogate variables (SVs) that capture hidden sources of variation in your data — often unmeasured confounders like batch effects, circadian variation, or other systematic biases.
# Key Concepts for Interpretation
# 
# SVs are data-driven latent variables. They are extracted from the expression data without prior assumptions about what they represent.
# 
# Each SV is a linear combination of the original gene expression data that best explains the largest sources of variation not accounted for by the known variables in your design model.
# 
# SV1 generally captures the strongest unmeasured effect, SV2 the next strongest, and so on.
# 
# SVs may or may not correlate with known biological variables (e.g., sample group, batch, etc.) — part of the analysis is to determine what they represent.

################################################################################
# 7) REMOVE BATCH EFFECTS USING LIMMA WITH SVA SURROGATE VARIABLES
################################################################################

# For instance, using the surrogate variables (sva_results$sv):
sva_corrected_data <- removeBatchEffect(
  # melanoma_data,
  expr_matrix,
  covariates = sva_results$sv,
  design = design_mod
)

# Quick check of dimensions
cat("Dimensions of sva_corrected_data:", dim(sva_corrected_data), "\n")

################################################################################
# 8) ALTERNATIVE APPROACH: RUVSeq TO DETECT HIDDEN FACTORS
#    - RUVg or RUVs typically requires negative control genes or replicate info.
#    - In practice, pick a reliable negative control gene set to represent
#      non-differentially expressed genes (not circadian).
################################################################################

# Convert raw counts to SeqExpressionSet
set_ruv <- newSeqExpressionSet(
  counts = as.matrix(count_data),
  # phenoData = data.frame(row.names = colnames(count_data)) 
  phenoData = data.frame(meta_data) 
  # optionally include meta_data columns here as well
)

# Between-lane normalization (e.g., upper quartile or TMM)
set_ruv <- betweenLaneNormalization(set_ruv, which = "upper")

# Identify negative control genes. 
# In a circadian context, you might exclude known circadian genes, 
# or you might have a known stable housekeeping set.
negative_control_genes <- setdiff(rownames(count_data), circadian_gene_symbols)

# For demonstration, we assume negative_control_genes is large enough
# Use RUVg with k=1 (1 factor of unwanted variation) for simplicity
ruv_results <- RUVg(
  set_ruv,
  cIdx = negative_control_genes,
  k = 1
)

# Extract the RUV factor(s) from ruv_results
W_ruv <- pData(ruv_results)$W_1  # This is the first RUV factor

# Remove batch effect using the RUV factor(s)
# If you have multiple factors (k>1), then you'd pass them as a matrix.
ruv_corrected_data <- removeBatchEffect(
  as.matrix(count_data),
  covariates = W_ruv
  # optionally: design = design_mod
)

################################################################################
# 9) ALTERNATIVE STRATEGY: RANK-BASED CORRECTION USING edgeR
################################################################################

dge_obj <- DGEList(counts = count_data)
dge_obj <- calcNormFactors(dge_obj) 
edgeR_corrected_data <- cpm(dge_obj, normalized.lib.sizes = TRUE)

cat("Dimensions of edgeR_corrected_data:", dim(edgeR_corrected_data), "\n")

################################################################################
# 10) VALIDATE NORMALIZATION:
#     - Compare correlation structures among circadian genes before vs after
################################################################################

# Correlation of circadian genes before normalization
cor_before <- cor(t(circadian_data))

# Correlation after SVA correction
cor_after_sva <- cor(
  t(sva_corrected_data[rownames(sva_corrected_data) %in% circadian_gene_symbols, ])
)

# Correlation after RUV correction
cor_after_ruv <- cor(
  t(ruv_corrected_data[rownames(ruv_corrected_data) %in% circadian_gene_symbols, ])
)

# Correlation after edgeR rank-based normalization
cor_after_edgeR <- cor(
  t(edgeR_corrected_data[rownames(edgeR_corrected_data) %in% circadian_gene_symbols, ])
)

cat("\nCorrelations before/after (example stats):\n")
cat("Mean correlation (before):", mean(cor_before[upper.tri(cor_before)]), "\n")
cat("Mean correlation (SVA):", mean(cor_after_sva[upper.tri(cor_after_sva)]), "\n")
cat("Mean correlation (RUV):", mean(cor_after_ruv[upper.tri(cor_after_ruv)]), "\n")
cat("Mean correlation (edgeR):", mean(cor_after_edgeR[upper.tri(cor_after_edgeR)]), "\n")

################################################################################
# 11) OPTIONAL: COMPARE CIRCADIAN GENE EXPRESSION IN PRIMARY VS METASTASIS
#     - Example: a quick differential expression (DE) analysis with limma
#     - Modify as needed if your metadata has other variable names
################################################################################

# Example if your metadata has 'primary_or_metastasis' with levels "Primary","Metastasis"
# design_limma <- model.matrix(~ 0 + primary_or_metastasis, data=meta_data)
# colnames(design_limma) <- c("Metastasis","Primary")
# 
# # SVA-corrected data for limma:
# fit <- lmFit(sva_corrected_data, design = design_limma)
# contrast_matrix <- makeContrasts(Met_vs_Prim = Metastasis - Primary, levels = design_limma)
# fit2 <- contrasts.fit(fit, contrast_matrix)
# fit2 <- eBayes(fit2)
# top_de <- topTable(fit2, coef = "Met_vs_Prim", number = Inf)
# 
# cat("\nTop differentially expressed circadian genes between Metastasis and Primary:\n")
# head(top_de[rownames(top_de) %in% circadian_gene_symbols, ], 20)

################################################################################
# END OF WORKFLOW
################################################################################
