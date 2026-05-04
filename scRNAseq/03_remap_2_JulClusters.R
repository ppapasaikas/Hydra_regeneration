# Generate:
#   Figures/Full_projection_JAnnot.pdf
#
# Usage:
#   Rscript scripts/03_remap_2_JulClusters.R \
#     annotated_sce/sce_FullAnnot.rds \
#     nonDubMarkers.csv \
#     HVAEP_105_ConversionTable.rds \
#     Full_projection_JAnnot.pdf \
#     JulianoClusterAssignment.rds

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(scater)
    library(ggplot2)
    library(AUCell)
    library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
    stop(
        "Expected 5 arguments:\n",
        "  input_sce juliano_markers_csv conversion_table_rds output_pdf assignment_cache_rds",
        call. = FALSE
    )
}

input_sce <- args[[1]]
juliano_markers_csv <- args[[2]]
conversion_table_rds <- args[[3]]
output_pdf <- args[[4]]
assignment_cache_rds <- args[[5]]

dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(assignment_cache_rds), recursive = TRUE, showWarnings = FALSE)

sce <- readRDS(input_sce)

if (!"majority_class" %in% colnames(colData(sce))) {
    stop("The input SCE must contain colData(sce)$majority_class.", call. = FALSE)
}
if (!"UMAP.batchcor" %in% reducedDimNames(sce)) {
    stop("The input SCE must contain reducedDim(sce, 'UMAP.batchcor').", call. = FALSE)
}

## majority_class coding from the original analysis:
# 1 = interstitial, 2 = ectoderm, 3 = endoderm

JulMarkers <- read.csv(juliano_markers_csv, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
ConvTable <- readRDS(conversion_table_rds)

if (!"gene" %in% colnames(JulMarkers)) {
    stop("Juliano marker table must contain a 'gene' column.", call. = FALSE)
}
if (!all(c("HvAEP_GID", "Hv105_GID") %in% colnames(ConvTable))) {
    stop("Conversion table must contain columns 'HvAEP_GID' and 'Hv105_GID'.", call. = FALSE)
}

JulMarkers$gene <- gsub("\\-", "_", JulMarkers$gene)
JulMarkers$gene_105 <- ConvTable$Hv105_GID[match(JulMarkers$gene, ConvTable$HvAEP_GID)]
JulMarkers <- JulMarkers[!is.na(JulMarkers$gene_105), , drop = FALSE]


## Keep only cluster and mapped gene columns.
JulMarkers <- unique(JulMarkers[, c(6, 8)])
colnames(JulMarkers) <- c("cluster", "gene_105")

if (file.exists(assignment_cache_rds)) {
    JulClust <- readRDS(assignment_cache_rds)
} else {
    rankings <- AUCell_buildRankings(
        counts(sce),
        plotStats = FALSE,
        verbose = FALSE,
        splitByBlocks = TRUE,
        BPPARAM = BiocParallel::MulticoreParam(workers = 8)
    )
    
    gene_sets <- split(JulMarkers$gene_105, JulMarkers$cluster)
    
    I_gene_sets  <- gene_sets[grep("^I_", names(gene_sets))]
    Ec_gene_sets <- gene_sets[grep("^Ec", names(gene_sets))][1:5]
    En_gene_sets <- gene_sets[grep("^En", names(gene_sets))][1:4]
    
    assign_juliano <- function(gene_sets, cell_subset) {
        auc <- AUCell_calcAUC(gene_sets, rankings[, cell_subset, drop = FALSE])
        auc_mat <- t(assay(auc))
        colnames(auc_mat)[max.col(auc_mat)]
    }
    
    JulClust <- rep(NA_character_, ncol(sce))
    
    JulClust[sce$majority_class == 1] <- assign_juliano(I_gene_sets,  sce$majority_class == 1)
    JulClust[sce$majority_class == 2] <- assign_juliano(Ec_gene_sets, sce$majority_class == 2)
    JulClust[sce$majority_class == 3] <- assign_juliano(En_gene_sets, sce$majority_class == 3)
    
    names(JulClust) <- colnames(sce)
    saveRDS(JulClust, assignment_cache_rds)
}

sce$JulClust <- JulClust

cluster_colors <- c(
    "Ec_BasalDisk"   = "#D1E9FA",
    "Ec_BodyCol/SC"  = "#31C8FF",
    "Ec_Head"        = "#1E58DF",
    "Ec_Peduncle"    = "#9BAFAE",
    "Ec_Tentacle"    = "#2A1F8F",
    "En_BodyCol/SC"  = "#D8F8D4",
    "En_Foot"        = "#44D04A",
    "En_Head"        = "#A3C61A",
    "En_Tentacle"    = "#0B5A3A",
    "I_DesmoNB"      = "#599E41",
    "I_DesmoNC"      = "#A09C67",
    "I_EarlyNem"     = "#E79A8E",
    "I_Ec1/5N"       = "#F47878",
    "I_Ec1N"         = "#EC4B4C",
    "I_Ec2N"         = "#E31E20",
    "I_Ec3N"         = "#EB4F36",
    "I_Ec4N"         = "#F48954",
    "I_En1N"         = "#FDBC6B",
    "I_En2N"         = "#FDA644",
    "I_En3N"         = "#FE8F1C",
    "I_FemGC"        = "#F98314",
    "I_GlProgen"     = "#E79660",
    "I_GranGl"       = "#D4A8AC",
    "I_ISC"          = "#BA9FCC",
    "I_IsoNB"        = "#9875B7",
    "I_IsoNC"        = "#764CA1",
    "I_MaleGC"       = "#8B6899",
    "I_Neuro"        = "#C0AD99",
    "I_SpumMucGl"    = "#F5F299",
    "I_StenoNB"      = "#E8CE78",
    "I_StenoNC"      = "#CC9350",
    "I_ZymoGl"       = "#B15928"
)

pJ <- scater::plotReducedDim(
    sce,
    dimred = "UMAP.batchcor",
    colour_by = "JulClust",
    point_size = 0.75
) +
    scale_color_manual(
        name = "Juliano Cluster",
        values = cluster_colors
    ) +
    xlab("UMAP1") +
    ylab("UMAP2")

ggsave(
    filename = output_pdf,
    plot = pJ,
    device = cairo_pdf,
    width = 8.125,
    height = 6
)

