library("GenomicRanges")
library("IRanges")
library("MultiAssayExperiment")
library("SummarizedExperiment")
library("S4Vectors")
library("dplyr")
library("Matrix") 
library("MOFA2")
library("ggplot2")
library("parallel")  # For multi-core detection
library("reticulate")

# Set Python environment
use_python("/usr/bin/python3")  # Specify your Python path if needed

# Step 1: Load Data
file_path <- "~/tmp/circadian_melanoma/depmap_data/"
# depmap_mae <- readRDS(file = paste0(file_path, "depmap_MAE_2025-03-19.rds"))
depmap_mae <- readRDS(file = paste0(file_path, "depmap_MAE_2025-03-18.rds"))

# Step 2: Prepare Data for MOFA2, Align samples across assays
depmap_mae <- intersectColumns(depmap_mae)

# Convert assay data into sparse matrices for efficiency
convert_to_matrix <- function(se_object) {
  assay_data <- assay(se_object)
  if (is(assay_data, "dgCMatrix")) {
    return(as.matrix(assay_data))  # Convert sparse matrices to dense matrices for MOFA2
  } else if (is.data.frame(assay_data)) {
    return(as.matrix(assay_data))  # Convert data.frame to matrix
  } else {
    return(assay_data)  # Already in matrix format
  }
}

# Extract and format the assays
mofa_data <- lapply(experiments(depmap_mae), convert_to_matrix)

# Confirm all elements are now sparse matrices
str(mofa_data)

# Step 3: Create and Train the MOFA2 Model, Create MOFA object
mofa_model <- create_mofa(mofa_data)

# Define model and training options
data_opts <- get_default_data_options(mofa_model)
model_opts <- get_default_model_options(mofa_model)
train_opts <- get_default_training_options(mofa_model)

# Performance Improvements
train_opts$maxiter <- 1000                   # Increase iterations for convergence

# Train the MOFA model
mofa_model <- prepare_mofa(mofa_model,
                           data_options = data_opts, 
                           model_options = model_opts, 
                           training_options = train_opts)

mofa_model <- run_mofa(mofa_model, use_basilisk = TRUE)

# Step 4: Visualize Results, Visualize variance explained by each factor
plot_variance_explained(mofa_model)

# Plot the first two latent factors
plot_factors(mofa_model, factors = 1:2)

# Heatmap visualization of feature weights in each latent factor
plot_weights(mofa_model, view = "tpm_exp", factor = 1)

# Step 5: Cluster and Analyze the Latent Space, Extract factor scores
factor_scores <- get_factors(mofa_model)[[1]]

# Perform k-means clustering
set.seed(42)
kmeans_clusters <- kmeans(factor_scores, centers = 3)

# Visualize clusters
factor_df <- data.frame(factor_scores, Cluster = factor(kmeans_clusters$cluster))
ggplot(factor_df, aes(x = Factor1, y = Factor2, color = Cluster)) +
  geom_point() +
  theme_minimal() +
  labs(title = "MOFA Latent Space Clustering")
