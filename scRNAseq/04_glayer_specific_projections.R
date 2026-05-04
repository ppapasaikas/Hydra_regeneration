## Generate ectoderm/endoderm positional PCA/UMAP projections w. segmentation info
#
# Outputs:
#   sce_ect_with_proj.rds
#   sce_end_with_proj.rds
#
## Usage:
#   Rscript scripts/04_make_ect_end_positional_projections.R \
#     sce_FullAnnot.rds \
#     Segments_SegmentCounts_NewAssembly.rds \
#     /positional_projections

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(scater)
    library(Matrix)
    library(matrixStats)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
    stop("Expected 3 arguments: input_sce segment_counts_rds output_dir", call. = FALSE)
}

input_sce <- args[[1]]
segment_counts_rds <- args[[2]]
output_dir <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

n_threads <- 16

sce <- readRDS(input_sce)
segment_counts <- readRDS(segment_counts_rds)

required_cols <- c("majority_class", "EXP_TIME")
missing_cols <- setdiff(required_cols, colnames(colData(sce)))
if (length(missing_cols)) {
    stop("Missing required colData columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}
if (!"logcounts" %in% assayNames(sce)) {
    stop("The input SCE must contain a logcounts assay.", call. = FALSE)
}

## majority_class coding:
##   1 = interstitial, 2 = ectoderm, 3 = endoderm

common_genes <- intersect(rownames(sce), rownames(segment_counts))
sce <- sce[common_genes, ]
segment_counts <- segment_counts[common_genes, , drop = FALSE]


## Segmentation genes used for positional projections 
segment_cpm <- 1e6 * apply(segment_counts, 2, function(x) {
    x / sum(x[x < stats::quantile(x, 0.95)])
})

body_cols <- grep("body", colnames(segment_cpm), ignore.case = TRUE)
if (!length(body_cols)) {
    stop("No body columns found in segment_counts column names.", call. = FALSE)
}

expressed_genes <- apply(segment_cpm, 1, function(x) stats::quantile(x, 0.05) > 3)
variable_genes <- apply(segment_cpm, 1, function(x) max(x + 1) / min(x + 1) > 5)
body_variable_genes <- apply(segment_cpm[, body_cols, drop = FALSE], 1, function(x) {
    max(x + 1) / min(x + 1) > 2
})

segmentation_genes <- rownames(segment_cpm)[expressed_genes & variable_genes & body_variable_genes]

## Germ-layer-specific expressed genes 
pseudobulks <- cbind(
    PB_int = Matrix::rowSums(counts(sce[, sce$majority_class == 1])),
    PB_ect = Matrix::rowSums(counts(sce[, sce$majority_class == 2])),
    PB_end = Matrix::rowSums(counts(sce[, sce$majority_class == 3]))
)

projection_genes_ect <- intersect(
    segmentation_genes,
    rownames(pseudobulks)[pseudobulks[, "PB_ect"] > stats::quantile(pseudobulks[, "PB_ect"], 0.5)]
)
projection_genes_end <- intersect(
    segmentation_genes,
    rownames(pseudobulks)[pseudobulks[, "PB_end"] > stats::quantile(pseudobulks[, "PB_end"], 0.5)]
)

expressed_in_ect <- rownames(sce)[
    Matrix::rowSums(counts(sce[, sce$majority_class == 2]) > 0) > 0.05 * sum(sce$majority_class == 2)
]
expressed_in_end <- rownames(sce)[
    Matrix::rowSums(counts(sce[, sce$majority_class == 3]) > 0) > 0.05 * sum(sce$majority_class == 3)
]

# Add reference genes for later checks:
marker_genes_to_keep <- c(
    "100203050", "100214868", "101237920", "100204831", "100201806",
    "100209110", "100204611", "100214371", "100208376", "100209580",
    "100199630", "105846072", "101237470", "100199257", "100215038",
    "100214773", "100207084", "100205558", "100213185"
)

sce_ect <- sce[unique(c(expressed_in_ect, marker_genes_to_keep, projection_genes_ect)), sce$majority_class == 2]
sce_end <- sce[unique(c(expressed_in_end, marker_genes_to_keep, projection_genes_end)), sce$majority_class == 3]

## Projection helper
add_positional_projection <- function(sce_layer, projection_genes, params, invert_umap = FALSE) {
    projection_genes <- intersect(projection_genes, rownames(sce_layer))
    if (length(projection_genes) < 10L) {
        stop("Too few projection genes available for this germ layer.", call. = FALSE)
    }
    
    sce_sub <- sce_layer[projection_genes, ]
    
    keep <- colSums(as.matrix(counts(sce_sub)) == 0) > 0.5 * nrow(sce_sub)
    sce_sub <- sce_sub[, keep]
    sce_layer <- sce_layer[, keep]
    
    pca <- stats::prcomp(t(logcounts(sce_sub)), retx = TRUE, center = TRUE, scale. = FALSE)$x
    
    lib_size <- log(sce_sub$total)
    pc_cor <- abs(cor(lib_size, pca, use = "pairwise.complete.obs"))
    keep_pcs <- which(pc_cor <= params$cor_threshold)
    pca_filtered <- pca[, keep_pcs, drop = FALSE]
    
    n_pcs <- min(params$n_pcs, ncol(pca_filtered))
    reducedDim(sce_layer, "PCAsegm") <- pca_filtered[, seq_len(n_pcs), drop = FALSE]
    
    reducedDim(sce_layer, "UMAP_PCAsegm") <- scater::calculateUMAP(
        sce_layer,
        n_threads = n_threads,
        dimred = "PCAsegm",
        min_dist = 0.66,
        spread = params$spread,
        seed = params$seed,
        n_neighbors = params$n_neighbors
    )
    
    if (invert_umap) {
        reducedDim(sce_layer, "UMAP_PCAsegm")[, 1] <- -reducedDim(sce_layer, "UMAP_PCAsegm")[, 1]
        reducedDim(sce_layer, "UMAP_PCAsegm")[, 2] <- -reducedDim(sce_layer, "UMAP_PCAsegm")[, 2]
    }
    
    sce_layer
}

## Parameters
set.seed(1)
seeds <- sample(1:100, 16)
ncomp <- sample(50:90, 16, replace = TRUE)
thr <- sample(c(0.4, 0.5), 16, replace = TRUE)
neighbors <- sample(40:80, 16)
spread <- runif(16, 0.241, 0.28)

ecto_params <- list(
    seed = seeds[4],
    n_pcs = ncomp[4],
    cor_threshold = thr[4],
    n_neighbors = neighbors[4],
    spread = spread[4]
)

endo_params <- list(
    seed = seeds[3],
    n_pcs = ncomp[3],
    cor_threshold = thr[3],
    n_neighbors = neighbors[3],
    spread = spread[3]
)

sce_ect <- add_positional_projection(
    sce_ect,
    projection_genes = projection_genes_ect,
    params = ecto_params,
    invert_umap = TRUE
)

sce_end <- add_positional_projection(
    sce_end,
    projection_genes = projection_genes_end,
    params = endo_params,
    invert_umap = FALSE
)

saveRDS(sce_ect, file.path(output_dir, "sce_ect_with_proj.rds"))
saveRDS(sce_end, file.path(output_dir, "sce_end_with_proj.rds"))

