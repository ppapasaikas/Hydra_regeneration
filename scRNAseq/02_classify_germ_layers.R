### Batch-correct a SingleCellExperiment, assign germ layers from marker scores,
# collapse cell-level assignments to cluster-level majority labels, and save
# full/endoderm/ectoderm/interstitial SCEs.
##
# Usage:
#   Rscript scripts/02_batch_correct_and_split_germ_layers.R \
#     alevin_spliced_introns_separate_sce.rds \
#     AEP_105Hydra2_105ncbi_genemapping.rds \
#     aav9314-Table-S9.txt \
#     annotated_sce
#
# Required colData columns in the input SCE:
#   experiment   sample/batch identifier


suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(S4Vectors)
    library(scater)
    library(scran)
    library(batchelor)
    library(BiocParallel)
    library(BiocNeighbors)
    library(BiocSingular)
    library(Rtsne)
    library(rsvd)
    library(igraph)
    library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
    stop(
        "Need 4 arguments: input_sce gene_mapping_rds juliano_marker_table output_dir",
        call. = FALSE
    )
}

input_sce <- args[[1]]
gene_mapping_rds <- args[[2]]
juliano_marker_table <- args[[3]]
output_dir <- args[[4]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

n_workers <- 4
min_fraction_detected <- .01
n_pcs <- 32
snn_k <- 64

bpparam <- MulticoreParam(n_workers)

sce <- readRDS(input_sce)
if (!"experiment" %in% colnames(colData(sce))) {
    stop("The input SCE must contain colData(sce)$experiment.", call. = FALSE)
}

ID2COND <- c(
    "2869F1" = "REG_HEAD_96H", "2869F2" = "REG_HEAD_96H",
    "3113F1" = "REG_HEAD_0H",  "3113F2" = "REG_HEAD_6H",
    "3113F3" = "REG_HEAD_12H", "3113F4" = "REG_HEAD_24H",
    "3271F1" = "REG_FOOT_96H", "3271F2" = "REG_FOOT_24H",
    "3271F3" = "REG_FOOT_12H", "3271F4" = "REG_FOOT_6H",
    "3271F5" = "REG_FOOT_0H",
    "3357F1" = "REG_FOOT_48H", "3357F2" = "REG_HEAD_48H",
    "3746F1" = "REG_HEAD_72H", "3746F2" = "REG_FOOT_72H",
    "3746F3" = "REG_HEAD_3H",  "3746F4" = "REG_FOOT_3H"
)

missing_conditions <- setdiff(unique(as.character(sce$experiment)), names(ID2COND))
if (length(missing_conditions)) {
    stop(
        "No condition mapping found for sample(s): ",
        paste(missing_conditions, collapse = ", "),
        call. = FALSE
    )
}

sce$EXP_TIME <- unname(ID2COND[as.character(sce$experiment)])
sce$EXP <- sub("_[0-9]+H$", "", sce$EXP_TIME)
sce$TIME <- as.numeric(sub(".*_([0-9]+)H$", "\1", sce$EXP_TIME))


# Expressed-gene selection
is_expressed <- Matrix::rowSums(counts(sce) > 0) > min_fraction_detected * ncol(sce)
expressed_genes <- rownames(sce)[is_expressed]

if ("product" %in% colnames(rowData(sce))) {
    ribosomal_genes <- rownames(sce)[grep("ribosomal protein", rowData(sce)$product, ignore.case = TRUE)]
    expressed_genes <- setdiff(expressed_genes, ribosomal_genes)
}

high_threshold <- stats::quantile(Matrix::rowSums(counts(sce)), 0.995)
high_genes <- rownames(sce)[Matrix::rowSums(counts(sce)) > 5 * high_threshold]
expressed_genes <- setdiff(expressed_genes, high_genes)

sce_gfiltered <- sce[expressed_genes, ]

# Normalization and batch rescaling
set.seed(42)
clusters <- scran::quickCluster(
    sce_gfiltered,
    method = "igraph",
    min.mean = 0.1,
    use.ranks = FALSE,
    BPPARAM = bpparam
)

sce_gfiltered <- scran::computeSumFactors(
    sce_gfiltered,
    min.mean = 0.1,
    cluster = clusters,
    BPPARAM = bpparam
)

sce_gfiltered <- scater::logNormCounts(sce_gfiltered, BPPARAM = bpparam)
sce_rescaled <- batchelor::multiBatchNorm(sce_gfiltered, batch = sce_gfiltered$experiment)

var_by_batch <- lapply(unique(sce_rescaled$experiment), function(batch) {
    scran::modelGeneVar(sce_rescaled[, sce_rescaled$experiment == batch])
})
names(var_by_batch) <- unique(sce_rescaled$experiment)

combined_var <- scran::combineVar(var_by_batch)
hvgs <- scran::getTopHVGs(combined_var, prop = 0.75)
sce_proj <- sce_rescaled[hvgs, ]



## Dimensional reduction and batch correction
reducedDim(sce, "PCA") <- rsvd::rpca(
    t(logcounts(sce_proj)),
    k = n_pcs,
    retx = TRUE,
    center = TRUE,
    scale = FALSE
)$x

set.seed(1)
reducedDim(sce, "TSNE") <- Rtsne::Rtsne(
    reducedDim(sce, "PCA"),
    perplexity = 40,
    initial_dims = min(24, n_pcs),
    pca = FALSE,
    num_threads = n_workers,
    theta = 0.25
)$Y

reducedDim(sce, "UMAP") <- scater::calculateUMAP(
    sce,
    n_threads = n_workers,
    dimred = "PCA",
    n_dimred = n_pcs,
    min_dist = 0.3,
    spread = 0.38,
    seed = 42,
    n_neighbors = 60
)

mnn <- suppressWarnings(batchelor::fastMNN(
    logcounts(sce_proj),
    batch = sce_proj$experiment,
    d = n_pcs,
    BSPARAM = IrlbaParam()
))

reducedDim(sce, "batchcor") <- reducedDim(mnn, "corrected")

set.seed(1)
reducedDim(sce, "TSNE.batchcor") <- Rtsne::Rtsne(
    reducedDim(sce, "batchcor"),
    perplexity = 40,
    initial_dims = min(24, n_pcs),
    pca = FALSE,
    num_threads = n_workers,
    theta = 0.25
)$Y

reducedDim(sce, "UMAP.batchcor") <- scater::calculateUMAP(
    sce,
    n_threads = n_workers,
    dimred = "batchcor",
    n_dimred = n_pcs,
    min_dist = 0.2,
    spread = 0.2,
    seed = 42,
    n_neighbors = 50
)

reducedDims(sce_rescaled) <- reducedDims(sce)

lib_size <- Matrix::colSums(counts(sce))
logcounts(sce) <- log2(1e4 * sweep(counts(sce), 2, lib_size, FUN = "/") + 1)

## Cell-level germ-layer assignment
read_juliano_germ_markers <- function(marker_table, mapping_rds) {
    juliano <- read.delim(marker_table)
    gene_map <- readRDS(mapping_rds)
    
    required_cols <- c(
        "Cluster.ID",
        "Average.Log.Fold.Change",
        "Percent.Positive.Cells.in.Cluster",
        "Percent.Positive.Cells.Outside.of.Cluster",
        "P.value",
        "Gene.ID.Swissprot.Annotation"
    )
    missing_cols <- setdiff(required_cols, colnames(juliano))
    if (length(missing_cols)) {
        stop("Missing required Juliano marker-table columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    if (!"MATCHED_NCBI_gIDs" %in% colnames(gene_map)) {
        stop("The gene mapping RDS must contain a MATCHED_NCBI_gIDs column.", call. = FALSE)
    }
    
    map_ids <- function(ids) {
        ids <- intersect(ids, rownames(gene_map))
        unique(na.omit(gene_map[ids, "MATCHED_NCBI_gIDs"]))
    }
    
    select_markers <- function(cluster_pattern, lfc, pct_in, pct_out) {
        in_cluster <- grepl(cluster_pattern, juliano$Cluster.ID)
        keep <- in_cluster &
            juliano$Average.Log.Fold.Change > log2(lfc) &
            juliano$Percent.Positive.Cells.in.Cluster > pct_in &
            juliano$Percent.Positive.Cells.Outside.of.Cluster < pct_out &
            juliano$P.value < 1e-10
        
        map_ids(juliano$Gene.ID.Swissprot.Annotation[keep])
    }
    
    list(
        inter = select_markers("^i_", lfc = 1.5, pct_in = 0.85, pct_out = 0.6),
        ecto = select_markers("^ec", lfc = 2.0, pct_in = 0.75, pct_out = 0.4),
        endo = select_markers("^en", lfc = 2.0, pct_in = 0.75, pct_out = 0.4)
    )
}

score_marker_sets <- function(sce, markers) {
    lib_size <- Matrix::colSums(counts(sce))
    norm_counts <- stats::median(lib_size) * sweep(counts(sce), 2, lib_size, FUN = "/")
    
    marker_sets <- list(
        inter = setdiff(intersect(markers$inter, rownames(sce)), c(markers$ecto, markers$endo)),
        ecto = setdiff(intersect(markers$ecto, rownames(sce)), c(markers$inter, markers$endo)),
        endo = setdiff(intersect(markers$endo, rownames(sce)), c(markers$inter, markers$ecto))
    )
    
    empty_sets <- names(marker_sets)[vapply(marker_sets, length, integer(1)) == 0L]
    if (length(empty_sets)) {
        stop("No marker genes found for: ", paste(empty_sets, collapse = ", "), call. = FALSE)
    }
    
    scores <- sapply(marker_sets, function(genes) {
        x <- Matrix::colSums(norm_counts[genes, , drop = FALSE])
        log2((x / mean(x)) + 1e-12)
    })
    
    scaled_scores <- t(scale(t(scores)))
    scaled_scores[is.na(scaled_scores)] <- 0
    scaled_scores
}

smooth_scores <- function(scores, projection, k) {
    nn <- BiocNeighbors::findKmknn(projection, k = k, BPPARAM = bpparam)
    smoothed <- scores
    
    for (j in seq_len(ncol(scores))) {
        vals <- scores[, j][as.vector(nn$index)]
        vals <- matrix(vals, nrow = nrow(nn$index), byrow = FALSE)
        smoothed[, j] <- rowMeans(vals)
    }
    
    smoothed
}

assign_germ_layer <- function(scores, projection, k, min_margin = NULL) {
    smoothed <- smooth_scores(scores, projection, k = k)
    layer <- colnames(smoothed)[max.col(smoothed, ties.method = "first")]
    
    if (!is.null(min_margin)) {
        margin <- apply(smoothed, 1, max) - apply(smoothed, 1, median)
        layer[margin < min_margin] <- "undetermined"
    }
    
    unname(c(inter = 1L, ecto = 2L, endo = 3L, undetermined = 4L)[layer])
}

refine_germ_markers <- function(sce, class_col = "class_initial") {
    cls <- colData(sce)[[class_col]]
    count_mat <- counts(sce)
    
    summed <- cbind(
        inter = Matrix::rowSums(count_mat[, cls == 1L, drop = FALSE]),
        ecto = Matrix::rowSums(count_mat[, cls == 2L, drop = FALSE]),
        endo = Matrix::rowSums(count_mat[, cls == 3L, drop = FALSE])
    )
    summed <- apply(summed, 2, function(x) 1e7 * x / sum(x))
    expr <- log2(summed + 8)
    
    frac_pos <- cbind(
        inter = Matrix::rowSums(count_mat[, cls == 1L, drop = FALSE] != 0) / sum(cls == 1L),
        ecto = Matrix::rowSums(count_mat[, cls == 2L, drop = FALSE] != 0) / sum(cls == 2L),
        endo = Matrix::rowSums(count_mat[, cls == 3L, drop = FALSE] != 0) / sum(cls == 3L)
    )
    
    list(
        inter = rownames(sce)[
            expr[, "inter"] - expr[, "ecto"] > 1 &
                expr[, "inter"] - expr[, "endo"] > 1 &
                frac_pos[, "inter"] > 0.05 &
                frac_pos[, "ecto"] < 0.01 &
                frac_pos[, "endo"] < 0.01
        ],
        ecto = rownames(sce)[
            expr[, "ecto"] - expr[, "inter"] > 1 &
                expr[, "ecto"] - expr[, "endo"] > 1 &
                frac_pos[, "ecto"] > 0.25 &
                frac_pos[, "inter"] < 0.02 &
                frac_pos[, "endo"] < 0.02
        ],
        endo = rownames(sce)[
            expr[, "endo"] - expr[, "inter"] > 1 &
                expr[, "endo"] - expr[, "ecto"] > 1 &
                frac_pos[, "endo"] > 0.15 &
                frac_pos[, "inter"] < 0.02 &
                frac_pos[, "ecto"] < 0.02
        ]
    )
}

juliano_markers <- read_juliano_germ_markers(juliano_marker_table, gene_mapping_rds)
initial_scores <- score_marker_sets(sce, juliano_markers)
sce$class_initial <- assign_germ_layer(
    initial_scores,
    projection = reducedDim(sce, "UMAP.batchcor"),
    k = 50,
    min_margin = 1
)

refined_markers <- refine_germ_markers(sce, class_col = "class_initial")
refined_scores <- score_marker_sets(sce, refined_markers)
sce$class <- assign_germ_layer(
    refined_scores,
    projection = reducedDim(sce, "UMAP.batchcor"),
    k = 5,
    min_margin = NULL
)

## Numeric germ-layer coding follows the original analysis:
#   1 = interstitial, 2 = ectoderm, 3 = endoderm, 4 = undetermined
sce$class_label <- factor(
    sce$class,
    levels = c(1L, 2L, 3L, 4L),
    labels = c("inter", "ecto", "endo", "undetermined")
)

saveRDS(juliano_markers, file.path(output_dir, "germ_layer_markers_juliano.rds"))
saveRDS(refined_markers, file.path(output_dir, "germ_layer_markers_refined.rds"))

# Cluster-level majority germ-layer assignment 

graph <- scran::buildSNNGraph(sce, k = snn_k, use.dimred = "PCA")
sce$cluster <- factor(igraph::cluster_louvain(graph)$membership)

## Numeric majority_class coding:
#   1 = interstitial, 2 = ectoderm, 3 = endoderm
sce$majority_class <- sce$class
for (cluster_id in levels(sce$cluster)) {
    in_cluster <- sce$cluster == cluster_id
    sce$majority_class[in_cluster] <- as.integer(names(which.max(table(sce$class[in_cluster]))))
}

sce$majority_class_label <- factor(
    sce$majority_class,
    levels = c(1L, 2L, 3L),
    labels = c("inter", "ecto", "endo")
)

# Germ-layer-specific batch-corrected projections 
add_subset_projection <- function(sce_subset, genes) {
    genes <- intersect(genes, rownames(sce_subset))
    
    mnn <- suppressWarnings(batchelor::fastMNN(
        logcounts(sce_subset[genes, ]),
        batch = sce_subset$experiment,
        d = n_pcs,
        BSPARAM = IrlbaParam()
    ))
    
    reducedDim(sce_subset, "batchcor") <- reducedDim(mnn, "corrected")
    
    set.seed(1)
    reducedDim(sce_subset, "TSNE.batchcor") <- Rtsne::Rtsne(
        reducedDim(sce_subset, "batchcor"),
        perplexity = min(40, floor((ncol(sce_subset) - 1) / 3)),
        initial_dims = min(24, n_pcs),
        pca = FALSE,
        num_threads = n_workers,
        theta = 0.25
    )$Y
    
    reducedDim(sce_subset, "UMAP.batchcor") <- scater::calculateUMAP(
        sce_subset,
        n_threads = n_workers,
        dimred = "batchcor",
        n_dimred = n_pcs,
        min_dist = 0.2,
        spread = 0.2,
        seed = 42,
        n_neighbors = 30
    )
    
    sce_subset
}

projection_genes <- rownames(sce_proj)

sce_inter <- add_subset_projection(sce[, sce$majority_class == 1L], projection_genes)
sce_ecto <- add_subset_projection(sce[, sce$majority_class == 2L], projection_genes)
sce_endo <- add_subset_projection(sce[, sce$majority_class == 3L], projection_genes)

# Save outputs
saveRDS(sce, file.path(output_dir, "sce_FullAnnot.rds"))
saveRDS(sce_endo, file.path(output_dir, "sceEndo_FullAnnot.rds"))
saveRDS(sce_ecto, file.path(output_dir, "sceEcto_FullAnnot.rds"))
saveRDS(sce_inter, file.path(output_dir, "sceInter_FullAnnot.rds"))

