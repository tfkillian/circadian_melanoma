## Create MultiAssayExperiment object from current MELANOMA Depmap datasets

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install(c("MultiAssayExperiment", 
#                        "SummarizedExperiment",
#                        "ExperimentHub", 
#                        "depmap", 
#                        "Matrix", 
#                        "MOFA2",
#                        "ggplot2",
#                        "dplyr",
#                        "tidyr",
#                        "tibble",
#                        "readr",
#                        "GenomicRanges",
#                        "IRanges",
#                        "S4Vectors",
#                        "reticulate",
#                        "parallel"), force = FALSE)

library("S4Vectors")
library("caret")
library("lattice")
library("dplyr")
library("tidyr")
library("tibble")
library("readr")
library("GenomicRanges")
library("IRanges")
library("MultiAssayExperiment")
library("SummarizedExperiment")
library("ExperimentHub")
library("Matrix")
library("depmap")

## increase buffer size for downloada
Sys.setenv(VROOM_CONNECTION_SIZE = 2000000L)

## output file path
# file_path <- "media/seq-srv-05/vrc/Project/Project_Theo" ## server
file_path <- "~/tmp/circadian_melanoma/MAE/" ## local

#################### depmap `metadata_22Q2` dataset ############################
# NOTE: the most current depmap release doesn't appear to have a "metadata" file
# as with previous releases, therefore, we will rely on this file for now.
eh <- ExperimentHub()
query(eh, "depmap")
eh[["EH7558"]] %>%
  dplyr::filter(grepl("melanoma", lineage_subtype) |
                  grepl("elanoma", Cellosaurus_NCIt_disease) |
                  grepl("elanoma", subtype_disease) | 
                  grepl("elanoma", cell_line)) %>%
  dplyr::select(-contains("issues")) %>%
  as.data.frame() -> melanoma_metadata_22Q2
rownames(melanoma_metadata_22Q2) <- melanoma_metadata_22Q2$depmap_id

# melanoma_metadata_22Q2 %>% 
#   as_tibble() %>% 
#     dplyr::slice(1:10) %>%
#     write_csv(file = paste0("~/tmp/circadian_melanoma/MAE/meta_test.csv"))

######################## depmap dep_2_name  ####################################
### `dep_2_name` to add `depmap_id` or `cell_line` to other datasets
melanoma_metadata_22Q2 %>%
  dplyr::select(depmap_id, cell_line, lineage_subtype) %>% 
  dplyr::filter(grepl("melanoma", lineage_subtype)) %>% ## only melanoma cell lines!
  dplyr::select(-lineage_subtype) -> dep_2_name_22Q2

######################## depmap `crispr_24Q1` dataset ##########################
readr::read_csv("https://plus.figshare.com/ndownloader/files/51064667") %>%
  dplyr::rename(depmap_id = names(.)[1]) %>%
  dplyr::rename_with(~ gsub("&", ";", sub(" \\(.+\\)$", "", .))) %>%
  dplyr::filter(depmap_id %in% melanoma_metadata_22Q2$depmap_id) %>% ## only melanoma cell lines!
  tibble::column_to_rownames(var = "depmap_id") %>%
  t() %>% scale(center = TRUE, scale = TRUE) %>%
  as.data.frame() -> crispr_24Q1_mat

# crispr_names <- names(crispr_24Q1_mat) ## store names

## build SE object
data.frame(Samples = names(crispr_24Q1_mat)) -> crispr_colData
rownames(crispr_colData) <- crispr_colData$Samples
SummarizedExperiment(assays = crispr_24Q1_mat,
                     colData = crispr_colData) -> crispr_se
rm(crispr_24Q1_mat)

######################## depmap `copyNumber_24Q1` dataset ######################
readr::read_csv("https://plus.figshare.com/ndownloader/files/51065324") %>%
  dplyr::rename(depmap_id = names(.)[1]) %>%
  dplyr::rename_with(~ gsub("&", ";", sub(" \\(.+\\)$", "", .))) %>%
  dplyr::filter(depmap_id %in% melanoma_metadata_22Q2$depmap_id) %>% ## only melanoma cell lines!
  tibble::column_to_rownames(var = "depmap_id") %>%
  t() %>% scale(center = TRUE, scale = TRUE) %>%
  as.data.frame() -> copyNumber_24Q1_mat

# copyNumber_names <- names(copyNumber_24Q1_mat) ## store names

## build SE object
data.frame(Samples = names(copyNumber_24Q1_mat)) -> copyNumber_colData
rownames(copyNumber_colData) <- copyNumber_colData$Samples
SummarizedExperiment(assays = copyNumber_24Q1_mat,
                     colData = copyNumber_colData) -> copyNumber_se
rm(copyNumber_24Q1_mat)

######################## depmap `TPM_22Q2` dataset #############################
## NOTE: the new expression dataset has weird Depmap IDs and it's unclear what they are
# readr::read_csv("https://plus.figshare.com/ndownloader/files/51065360") %>%
#   dplyr::rename(depmap_id = names(.)[1]) %>%
#   # dplyr::rename_with(~ gsub("&", ";", sub(" \\(.+\\)$", "", .))) %>%
#   # dplyr::filter(depmap_id %in% melanoma_metadata_22Q2$depmap_id) %>% ## only melanoma cell lines!
#   # tibble::column_to_rownames(var = "depmap_id") %>%
#   # t() %>%
#   as.data.frame() -> TPM_22Q2_mat
# View(TPM_22Q2_mat)

# eh <- ExperimentHub()
# query(eh, "depmap")
eh[["EH7556"]] %>%
  dplyr::filter(depmap_id %in% melanoma_metadata_22Q2$depmap_id) %>% ## only melanoma cell lines!
  dplyr::select(depmap_id, rna_expression, gene_name) %>%
  tidyr::pivot_wider(names_from = depmap_id,
                     values_from = rna_expression) %>%
  dplyr::filter(!duplicated(gene_name)) %>%
  tibble::column_to_rownames(var = "gene_name") %>%
  scale(center = TRUE, scale = TRUE) %>%
  as.data.frame() -> TPM_22Q2_mat

# TPM_22Q2_mat %>% 
#   tibble::rownames_to_column(var = "gene_name") %>% 
#   dplyr::slice(1:10) %>% 
#   write_csv(file = paste0("~/tmp/circadian_melanoma/MAE/tmp_test.csv"))

# TPM_names <- names(TPM_22Q2_mat) ## store names

## build SE object
data.frame(Samples = names(TPM_22Q2_mat)) -> TPM_colData
rownames(TPM_colData) <- TPM_colData$Samples
SummarizedExperiment(assays = TPM_22Q2_mat,
                     colData = TPM_colData) -> TPM_se
rm(TPM_22Q2_mat)

##################### depmap `RPPA_19Q3` dataset ###############################
## download the 19Q3 metadata
# eh <- ExperimentHub()
# query(eh, "depmap")
eh[["EH3086"]] -> metadata_19Q3
metadata_19Q3 %>% dplyr::select(depmap_id, cell_line) -> dep_2_name_19Q3

### loading data (downloading .csv file from online source)
readr::read_csv(
  paste0("https://depmap.org/portal/download/api/download?file_name=ccle%2Fccle",
         "_2019%2FCCLE_RPPA_20181003.csv&bucket=depmap-external-downloads")) %>%
  dplyr::rename(cell_line = names(.)[1]) %>%
  dplyr::left_join(dep_2_name_19Q3, by = "cell_line") %>%
  dplyr::mutate(
    depmap_id = coalesce(
      depmap_id,
      dplyr::case_when(
        cell_line == "NCIH684_LIVER" ~ "ACH-000089",
        cell_line == "KE97_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE" ~ "ACH-000167",
        TRUE ~ NA_character_))) %>%
  dplyr::select(-cell_line) %>%
  dplyr::filter(depmap_id %in% melanoma_metadata_22Q2$depmap_id) %>% ## only melanoma cell lines!
  tibble::column_to_rownames(var = "depmap_id") %>%
  t() %>% as.data.frame() -> RPPA_19Q3_mat

## https://web.expasy.org/cellosaurus/CVCL_9980 # NCIH684_LIVER
## https://web.expasy.org/cellosaurus/CVCL_3386 # KE97_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE

## build SE object
data.frame(Samples = names(RPPA_19Q3_mat)) -> RPPA_colData
rownames(RPPA_colData) <- RPPA_colData$Samples
SummarizedExperiment(assays = RPPA_19Q3_mat,
                     colData = RPPA_colData) -> RPPA_se
rm(RPPA_19Q3_mat)

##################### depmap `proteomic_20Q2` dataset ##########################
readr::read_csv(paste(
  "https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz")) %>% 
  dplyr::select(Gene_Symbol, Uniprot_Acc, MDAMB468_BREAST_TenPx01:last_col()) %>%
  dplyr::mutate(
    Gene_Symbol = case_when(
      Uniprot_Acc == "H7BZ55" ~ "CROCC2",
      Uniprot_Acc == "A6NL28" ~ "LOC101929943",
      Uniprot_Acc == "Q9H552" ~ "KRT8P11",
      Uniprot_Acc == "S4R362" ~ "SEPTIN1",
      Uniprot_Acc == "Q96FF7" ~ "MISP3",
      Uniprot_Acc == "Q2M2H8" ~ "MGAM2",
      Uniprot_Acc == "P01614" ~ "IGKV2D-40",
      Uniprot_Acc == "P01619" ~ "IGKV3-20",
      Uniprot_Acc == "Q69YL0" ~ "NCBP2AS2",
      Uniprot_Acc == "I3L1I5" ~ "ARHGEF18",
      Uniprot_Acc == "E9PS84" ~ "AP000783.1",
      Uniprot_Acc == "P60507" ~ "ERVFC1",
      Uniprot_Acc == "H3BQF6" ~ "RP11-12J10.3",
      Uniprot_Acc == "C9J7I0" ~ "UMAD1",
      Uniprot_Acc == "Q902F9" ~ "HERVK_113",
      Uniprot_Acc == "P63135" ~ "ERVK-7",
      Uniprot_Acc == "H3BMH7" ~ "RP11-178L8.4",
      Uniprot_Acc == "P01780" ~ "IGHV3-7",
      TRUE ~ Gene_Symbol)) %>%
  dplyr::filter(!is.na(Gene_Symbol)) %>%
  dplyr::filter(!duplicated(Gene_Symbol)) %>%
  tibble::column_to_rownames(var = "Gene_Symbol") %>% 
  dplyr::select(-Uniprot_Acc) %>%
  as.data.frame() -> proteomic_20Q2
names(proteomic_20Q2) <- gsub('(.*)_\\w+', '\\1', names(proteomic_20Q2))
proteomic_20Q2 %>%
  t() %>% as.data.frame() %>%
  tibble::rownames_to_column() %>%
  dplyr::rename(cell_line = rowname) %>%
  dplyr::left_join(dep_2_name_22Q2, by = "cell_line") %>%
  dplyr::select(cell_line, depmap_id, everything()) %>%
  dplyr::mutate(
    depmap_id = case_when(
      cell_line == "SW948_LARGE_INTESTINE.1" ~ "ACH-000680",
      cell_line == "CAL120_BREAST.1" ~ "ACH-000212",
      cell_line == "X8505C_THYROID" ~ "ACH-001307",
      cell_line == "HCT15_LARGE_INTESTINE.1" ~ "ACH-000997",
      cell_line == "SW948_LARGE_INTESTINE" ~ "ACH-000680",
      TRUE ~ depmap_id)) %>%
  dplyr::filter(!is.na(depmap_id)) %>%
  dplyr::filter(!duplicated(depmap_id)) %>%
  dplyr::select(-cell_line) %>%
  dplyr::filter(depmap_id %in% melanoma_metadata_22Q2$depmap_id) %>% ## only melanoma cell lines!
  tibble::column_to_rownames(var = "depmap_id") %>%
  t() %>% scale(center = TRUE, scale = TRUE) %>%
  as.data.frame() -> proteomic_20Q2_mat

# proteomic_names <- names(proteomic_20Q2_mat) ## store names

data.frame(Samples = names(proteomic_20Q2_mat)) -> proteomic_colData
rownames(proteomic_colData) <- proteomic_colData$Samples
SummarizedExperiment(assays = proteomic_20Q2_mat,
                     colData = proteomic_colData) -> proteomic_se

rm(proteomic_20Q2_mat)

############################ methylation dataset ###############################
readr::read_csv(
  paste0("~/tmp/circadian_melanoma/MAE/Methylation_(1kb_upstream_TSS)_subsetted_NAsdropped.csv")) %>%
  dplyr::rename(depmap_id = names(.)[1]) %>%
  dplyr::filter(depmap_id %in% melanoma_metadata_22Q2$depmap_id) %>% ## only melanoma cell lines!
  tibble::column_to_rownames(var = "depmap_id") %>%
  t() %>% scale(center = TRUE, scale = TRUE) %>%
  as.data.frame() -> methylation_mat

# methylation_names <- names(methylation_mat) ## store names

data.frame(Samples = names(methylation_mat),
           row.names = methylation_mat$Samples) -> methylation_colData

## build SE object
SummarizedExperiment(assays = methylation_mat,
                     colData = methylation_colData) -> methylation_se
rm(methylation_mat)

###################### depmap `mutationCalls_21Q4` dataset #####################
readr::read_csv("https://ndownloader.figshare.com/files/31315930") %>%
  dplyr::rename_with(~ c(
    "gene_name", "entrez_id", "ncbi_build", "chromosome", "start_pos",
    "end_pos", "strand", "var_class", "var_type", "ref_allele",
    "tumor_seq_allele1", "dbSNP_RS", "dbSNP_val_status", "genome_change",
    "annotation_trans", "depmap_id", "cDNA_change", "codon_change",
    "protein_change", "is_deleterious", "is_tcga_hotspot", "tcga_hsCnt",
    "is_cosmic_hotspot", "cosmic_hsCnt", "ExAC_AF", "var_annotation",
    "CGA_WES_AC", "HC_AC", "RD_AC", "RNAseq_AC", "sanger_WES_AC", "WGS_AC"
    ), .cols = 1:32) %>%
  dplyr::filter(depmap_id %in% melanoma_metadata_22Q2$depmap_id) %>% ## only melanoma cell lines!
  dplyr::select(depmap_id, everything()) -> mutationCalls_21Q4

mutationCalls_names <- mutationCalls_21Q4$depmap_id ## store names

## 1) Create a GRanges object
gr_mutations <- GRanges(
  seqnames = mutationCalls_21Q4[["chromosome"]],
  ranges   = IRanges(
    start = mutationCalls_21Q4[["start_pos"]],
    end   = mutationCalls_21Q4[["end_pos"]]),
  strand   = mutationCalls_21Q4[["strand"]])

################ create encoded mutation call matrices #########################

# Encode mutation data for MOFA2
mutationCalls_21Q4 %>%
  # Create unique mutation feature (gene + genomic position)
  mutate(mutation_feature = paste0(gene_name, "_chr", chromosome, "_", start_pos, "_", end_pos)) %>%
  # Encode binary features
  mutate(is_deleterious = as.integer(is_deleterious == "TRUE"),
         is_tcga_hotspot = as.integer(is_tcga_hotspot == "TRUE"),
         is_cosmic_hotspot = as.integer(is_cosmic_hotspot == "TRUE")) %>%
  # Encode allele changes
  # mutate(allele_change = paste0(ref_allele, ">", tumor_seq_allele1)) %>%
  # 
  # # Normalize continuous features
  # # mutate(ExAC_AF = scale(ExAC_AF, center = TRUE, scale = TRUE)) %>%
  # 
  # # One-hot encode categorical variables (var_class, var_type)
  # # mutate(across(c(var_class, var_type, allele_change), as.factor)) %>%s
  # # pivot_wider(names_from = depmap_id, 
  # #             values_from = c(var_class, var_type, is_deleterious, ExAC_AF), 
  # #             values_fill = 0) %>%
  # dplyr::select(-c(gene_name:strand), -c(dbSNP_RS:protein_change),
  #               -c(CGA_WES_AC:WGS_AC), -tcga_hsCnt) %>%
  as.data.frame() -> mutation_encoded

# mutation_encoded %>% 
#   dplyr::slice(1:10) %>% 
#   readr::write_csv(file = "~/tmp/circadian_melanoma/depmap_data/mutation_encoded.csv")

##
## 1) var_class (categorical: one-hot encode)
##
var_class_mat <- mutation_encoded %>%
  # Keep only the relevant columns
  select(depmap_id, mutation_feature, var_class) %>%
  # Create a dummy value of 1 for each row
  mutate(value = 1) %>%
  # Pivot so that each distinct var_class becomes a separate column
  pivot_wider(
    id_cols      = c(depmap_id, mutation_feature),   # rows identified by (sample, feature)
    names_from   = var_class,                        # each unique var_class → new column
    values_from  = value,
    values_fill  = 0                                 # fill missing combinations with 0
  )

# var_class_mat %>% 
#   select(depmap_id, mutation_feature, Missense_Mutation) %>%
#   pivot_wider(
#     id_cols      = mutation_feature,   # rows identified by (sample, feature)
#     names_from   = depmap_id,          # each unique var_class → new column
#     values_from = Missense_Mutation) %>% 
#   tibble::column_to_rownames(Missense_Mutation) -> var_class_Missense_Mutation

var_cols <- names(var_class_mat)[3:19]  # adjust the range as needed

# 2) Define a helper function to pivot one column at a time into a sparse matrix
build_sparse_matrix <- function(df, var_col) {
df %>%
    select(depmap_id, mutation_feature, all_of(var_col)) %>%
    pivot_wider(
      id_cols = mutation_feature,    # rows = mutation features
      names_from = depmap_id,        # columns = sample IDs
      values_from = all_of(var_col), # values = the chosen column
      values_fill = 0                # fill missing combos with 0
    ) %>%
    column_to_rownames("mutation_feature") %>%
    as.matrix() %>%
    Matrix(sparse = TRUE)  # convert to a dgCMatrix
}

# 3) Apply this function across all desired columns to build a list of sparse matrices
sparse_list1 <- lapply(var_cols, function(vcol) {
  build_sparse_matrix(var_class_mat, vcol)
})

# Optionally, name each list element by the mutation column
names(sparse_list1) <- paste0("var_class_", gsub("'", "_", var_cols))

SE_list1 <- lapply(names(sparse_list1), function(view_name) {
  SummarizedExperiment(
    assays = list(counts = sparse_list1[[view_name]])
  )
})

names(SE_list1) <- names(sparse_list1)

SE_list1 <- lapply(SE_list1, function(se) {
  
  # Create colData based on sample names in each SummarizedExperiment
  colData <- data.frame(Samples = colnames(assay(se)))
  
  # Assign rownames to colData
  rownames(colData) <- colData$Samples
  
  # Add colData to the SummarizedExperiment
  SummarizedExperiment(
    assays = assay(se),  # Preserve the assay data
    colData = colData    # Add the new colData
  )
})

## 2) var_type (categorical: one-hot encode)
var_type_mat <- mutation_encoded %>%
  select(depmap_id, mutation_feature, var_type) %>%
  mutate(value = 1) %>%
  pivot_wider(
    id_cols     = c(depmap_id, mutation_feature),
    names_from  = var_type,
    values_from = value,
    values_fill = 0
  )

var_cols <- names(var_type_mat)[3:6]  # adjust the range as needed
sparse_list2 <- lapply(var_cols, function(vcol) {
  build_sparse_matrix(var_type_mat, vcol)
})

names(sparse_list2) <- paste0("var_type_", gsub("-", "_", gsub(" ", "_", var_cols)))

SE_list2 <- lapply(names(sparse_list2), function(view_name) {
  SummarizedExperiment(
    assays = list(counts = sparse_list2[[view_name]])
  )
})

names(SE_list2) <- names(sparse_list2)

SE_list2 <- lapply(SE_list2, function(se) {
  
  # Create colData based on sample names in each SummarizedExperiment
  colData <- data.frame(Samples = colnames(assay(se)))
  
  # Assign rownames to colData
  rownames(colData) <- colData$Samples
  
  # Add colData to the SummarizedExperiment
  SummarizedExperiment(
    assays = assay(se),  # Preserve the assay data
    colData = colData    # Add the new colData
  )
})

## 3) is_deleterious (binary: already 0/1, but shown as one-hot for completeness)
is_deleterious_mat <- mutation_encoded %>%
  select(depmap_id, mutation_feature, is_deleterious) %>%
  pivot_wider(
    id_cols     = mutation_feature,
    names_from  = depmap_id,
    values_from = is_deleterious,
    values_fill = 0) %>% 
  column_to_rownames("mutation_feature") %>%
  as.matrix() %>%
  Matrix(sparse = TRUE)

data.frame(Samples = colnames(is_deleterious_mat)) -> is_deleterious_colData
rownames(is_deleterious_colData) <- is_deleterious_colData$Samples
SummarizedExperiment(assays = is_deleterious_mat,
                     colData = is_deleterious_colData) -> is_deleterious_se

## 4) is_tcga_hotspot (binary)
is_tcga_hotspot_mat <- mutation_encoded %>%
  select(depmap_id, mutation_feature, is_tcga_hotspot) %>%
  pivot_wider(
    id_cols     = mutation_feature,
    names_from  = depmap_id,
    values_from = is_tcga_hotspot,
    values_fill = 0) %>% 
  column_to_rownames("mutation_feature") %>%
  as.matrix() %>%
  Matrix(sparse = TRUE)

data.frame(Samples = colnames(is_tcga_hotspot_mat)) -> is_tcga_hotspot_colData
rownames(is_tcga_hotspot_colData) <- is_tcga_hotspot_colData$Samples
SummarizedExperiment(assays = is_tcga_hotspot_mat,
                     colData = is_tcga_hotspot_colData) -> is_tcga_hotspot_se

## 5) is_cosmic_hotspot (binary)
is_cosmic_hotspot_mat <- mutation_encoded %>%
  select(depmap_id, mutation_feature, is_cosmic_hotspot) %>%
  pivot_wider(
    id_cols     = mutation_feature,
    names_from  = depmap_id,
    values_from = is_cosmic_hotspot,
    values_fill = 0) %>% 
  column_to_rownames("mutation_feature") %>%
  as.matrix() %>%
  Matrix(sparse = TRUE)

data.frame(Samples = colnames(is_cosmic_hotspot_mat)) -> is_cosmic_hotspot_colData
rownames(is_cosmic_hotspot_colData) <- is_cosmic_hotspot_colData$Samples
SummarizedExperiment(assays = is_cosmic_hotspot_mat,
                     colData = is_cosmic_hotspot_colData) -> is_cosmic_hotspot_se

## 6) cosmic_hsCnt (continuous/numeric: no one-hot, just pivot)
cosmic_hsCnt_mat <- mutation_encoded %>%
  select(depmap_id, mutation_feature, cosmic_hsCnt) %>%
  pivot_wider(
    id_cols     = mutation_feature,
    names_from  = depmap_id,
    values_from = cosmic_hsCnt,
    values_fill = 0) %>% 
  column_to_rownames("mutation_feature") %>%
  as.matrix() %>%
  Matrix(sparse = TRUE)

data.frame(Samples = colnames(cosmic_hsCnt_mat)) -> cosmic_hsCnt_colData
rownames(cosmic_hsCnt_colData) <- cosmic_hsCnt_colData$Samples
SummarizedExperiment(assays = cosmic_hsCnt_mat,
                     colData = cosmic_hsCnt_colData) -> cosmic_hsCnt_se

## 7) var_annotation (categorical: one-hot encode)
var_annotation_mat <- mutation_encoded %>%
  select(depmap_id, mutation_feature, var_annotation) %>%
  mutate(value = 1) %>%
  pivot_wider(
    id_cols     = c(depmap_id, mutation_feature),
    names_from  = var_annotation,
    values_from = value,
    values_fill = 0
  )

var_cols <- names(var_annotation_mat)[3:6]  # adjust the range as needed
sparse_list3 <- lapply(var_cols, function(vcol) {
  build_sparse_matrix(var_annotation_mat, vcol)
})

names(sparse_list3) <- paste0("var_annotation_", gsub("-", "_", gsub(" ", "_", var_cols)))

SE_list3 <- lapply(names(sparse_list3), function(view_name) {
  SummarizedExperiment(
    assays = list(counts = sparse_list3[[view_name]])
  )
})

names(SE_list3) <- names(sparse_list3)

SE_list3 <- lapply(SE_list3, function(se) {
  
  # Create colData based on sample names in each SummarizedExperiment
  colData <- data.frame(Samples = colnames(assay(se)))
  
  # Assign rownames to colData
  rownames(colData) <- colData$Samples
  
  # Add colData to the SummarizedExperiment
  SummarizedExperiment(
    assays = assay(se),  # Preserve the assay data
    colData = colData    # Add the new colData
  )
})

## 8) ExAC_AF (continuous/numeric)
ExAC_AF_mat <- mutation_encoded %>%
  select(depmap_id, mutation_feature, ExAC_AF) %>%
  pivot_wider(
    id_cols     = mutation_feature,
    names_from  = depmap_id,
    values_from = ExAC_AF) %>%
  column_to_rownames("mutation_feature") %>%
  as.matrix() %>%
  scale(center = TRUE, scale = TRUE) %>%
  Matrix(sparse = TRUE)

data.frame(Samples = colnames(ExAC_AF_mat)) -> ExAC_AF_colData
rownames(ExAC_AF_colData) <- ExAC_AF_colData$Samples
SummarizedExperiment(assays = ExAC_AF_mat,
                     colData = ExAC_AF_colData) -> var_effect_predictor_se

## 2) Store the remaining columns as metadata columns in the GRanges
##    First, pick out which columns you want to exclude from metadata:
exclude_cols <- c("chromosome", "start_pos", "end_pos", "strand")
metadata_cols <- setdiff(colnames(mutationCalls_21Q4), exclude_cols)
mcols(gr_mutations) <- mutationCalls_21Q4[, metadata_cols, drop = FALSE]

rm(mutationCalls_21Q4)

## new Mutation calls.... but there are no Depmap IDs... 
# readr::read_csv("https://plus.figshare.com/ndownloader/files/51065732") %>%
#   as.data.frame() -> b1
# View(b1)

###################### inner join of Depmap IDs ###############################

# length(methylation_names)
# length(proteomic_names)
# length(TPM_names)
# length(copyNumber_names)
# length(crispr_names)
# length(mutationCalls_names)
# 
# # Find common sample IDs (column names) using Reduce and intersect
# common_samples <- Reduce(intersect, list(
#   # methylation_names,
#   # proteomic_names,
#   mutationCalls_names,
#   TPM_names,
#   copyNumber_names,
#   crispr_names
# ))
# 
# length(common_samples)

#################### depmap `drug_dependency_21Q4` dataset #####################
read_csv("https://ndownloader.figshare.com/files/17741420") %>%
  as.data.frame() %>% dplyr::rename(depmap_id = names(.)[1]) %>%
  dplyr::mutate(depmap_id = gsub("_FAILED_STR", "", depmap_id)) %>%
  dplyr::filter(depmap_id %in% melanoma_metadata_22Q2$depmap_id) %>% ## only melanoma cell lines!
  tibble::column_to_rownames(var = "depmap_id") %>%
  as.matrix() %>% t() %>%
  scale(center = TRUE, scale = TRUE) -> drug_dep_19Q3_w

data.frame(Samples = colnames(drug_dep_19Q3_w)) -> drug_dep_colData
rownames(drug_dep_colData) <- drug_dep_colData$Samples
SummarizedExperiment(assays = drug_dep_19Q3_w,
                     colData = drug_dep_colData) -> drug_dep_se

rm(drug_dep_19Q3_w)

#################### depmap `drug_dependency_21Q4` metadata ####################
## data cleaning of `drug_sensativity` dataset and screen info
# read_csv("https://ndownloader.figshare.com/files/20237715") %>%
#   rename(compound = column_name) -> screen_info
# eh[["EH3087"]] -> drug_dep_19Q3

############################# build depmap MAE #################################
## Combine to a named list and call the ExperimentList constructor function
c(
  list(
     crispr      = crispr_se,
     copy_number = copyNumber_se,
     tpm_exp     = TPM_se,
     # proteomic   = proteomic_se,
     # methylation = methylation_se,
     # mutations   = gr_mutations,
     # RPPA        = RPPA_se,
     # drug_dep    = drug_dep_se,
     is_deleterious = is_deleterious_se,
     is_tcga_hotspot = is_tcga_hotspot_se,
     is_cosmic_hotspot = is_cosmic_hotspot_se,
     cosmic_hsCnt = cosmic_hsCnt_se,
     var_effect_predictor = var_effect_predictor_se
     ),
  setNames(SE_list1, names(SE_list1)),
  setNames(SE_list2, names(SE_list2)),
  setNames(SE_list3, names(SE_list3))
  ) -> assayList

## Use the ExperimentList constructor
exp_list <- ExperimentList(assayList)

## build colData, a DataFrame describing the characteristics of biological units
mae_coldata <- DataFrame(melanoma_metadata_22Q2)

# Step 1: Extract sample mappings for SummarizedExperiment objects
sampleMap_se <- lapply(names(assayList), function(assay_name) {
  if (assay_name == "mutations") return(NULL)  # Skip mutations here
  data.frame(
    assay = assay_name,
    primary = colnames(assayList[[assay_name]]),  # Primary sample IDs
    colname = colnames(assayList[[assay_name]])   # Same as colnames for SE objects
  )
}) %>% dplyr::bind_rows()

# Step 2: Extract sample mapping for the GRanges object (mutations)
# sampleMap_gr <- data.frame(
#   assay = "mutations",
#   primary = gr_mutations$depmap_id,  # Primary sample IDs
#   colname = gr_mutations$depmap_id   # Sample IDs directly from 'depmap_id'
# )

# Step 3: Combine sample maps
sampleMap <- dplyr::bind_rows(sampleMap_se#, sampleMap_gr
                              )

## build MAE
MultiAssayExperiment::MultiAssayExperiment(
  experiments = exp_list,
  colData = mae_coldata,
  sampleMap = sampleMap #,
  # checkDims = FALSE  # <-- turns off automatic dimension checking
) -> depmap_mae

## save depmap MAE
saveRDS(depmap_mae, file = paste0(file_path, "depmap_MAE_", Sys.Date(), ".rds"))

## TODO:
# 0) MAE https://www.bioconductor.org/packages/devel/bioc/vignettes/MultiAssayExperiment/inst/doc/QuickStartMultiAssay.html#transformation-reshaping
# 1) you need to download the original "square" omics matrices from Depmap: 
# https://depmap.org/portal/data_page/?tab=allData
# https://plus.figshare.com/articles/dataset/DepMap_24Q1_Public/27993248/1
# 2) are Regulatory or inhibitory compounds that effect of m6A modifications in the drug screen data?