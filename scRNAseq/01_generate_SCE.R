## Build SingleCellExperiment object from alevin quantifications

# Usage:
#   Rscript scripts/01_build_alevin_sce.R \
#     data/alevin_samples.tsv \
#     reference/Hydra_vulgaris_105_v3_processed.gff \
#     reference/salmon_linkedtxome.json \
#     reference/alevin_features.tsv \
#     results/alevin_spliced_introns_separate_sce.rds \
#     results/qc
#
# Required sample sheet columns:
#   sample_id    sample name, e.g. 2869F1
#   quant_file   path to alevin/quants_mat.gz

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(S4Vectors)
    library(scater)
    library(tximeta)
    library(rtracklayer)
    library(BiocParallel)
    library(dplyr)
    library(ggplot2)
    library(cowplot)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
    stop(
        "Need 6 arguments:\n",
        "  sample_sheet annotation_gff linked_txome_json feature_table output_rds qc_dir\n",
        call. = FALSE
    )
}

sample_sheet <- args[[1]]
annotation_gff <- args[[2]]
linked_txome_json <- args[[3]]
feature_table <- args[[4]]
output_rds <- args[[5]]
qc_dir <- args[[6]]

output_dir <- dirname(output_rds)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

n_workers <- 16
bpparam <- MulticoreParam(n_workers)

resolve_paths <- function(x, base_dir) {
    is_abs <- grepl("^/", x)
    normalizePath(ifelse(is_abs, x, file.path(base_dir, x)), mustWork = TRUE)
}

save_plot <- function(plot, filename, width = 8, height = 5) {
    ggplot2::ggsave(
        filename = file.path(qc_dir, filename),
        plot = plot,
        width = width,
        height = height,
        units = "in"
    )
}

sample_base <- dirname(normalizePath(sample_sheet, mustWork = TRUE))
samples <- read.delim(sample_sheet, stringsAsFactors = FALSE, check.names = FALSE)

required_cols <- c("sample_id", "quant_file")
missing_cols <- setdiff(required_cols, colnames(samples))
if (length(missing_cols)) {
    stop("Missing required sample sheet columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

samples <- samples %>%
    transmute(
        sample_id = as.character(sample_id),
        quant_file = resolve_paths(as.character(quant_file), sample_base)
    )

stopifnot(!anyDuplicated(samples$sample_id))

# Gene annotation 

gff <- rtracklayer::import(annotation_gff)

gene_annot <- as.data.frame(gff) %>%
    select(any_of(c("gene_id", "product", "gene", "seqnames"))) %>%
    distinct()

if (!"gene_id" %in% colnames(gene_annot)) {
    stop("The annotation file must contain a 'gene_id' field.", call. = FALSE)
}

if ("seqnames" %in% colnames(gene_annot)) {
    gene_annot <- rename(gene_annot, chromosome_name = seqnames)
}

# Read alevin quantifications

message("Loading linked transcriptome metadata")
tximeta::loadLinkedTxome(linked_txome_json)

feature_groups <- read.delim(feature_table, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
colnames(feature_groups)[colnames(feature_groups) == "intron"] <- "unspliced"

read_alevin_sample <- function(i) {
    sample_id <- samples$sample_id[[i]]
    message("Reading sample ", i, "/", nrow(samples), ": ", sample_id)
    
    coldata <- data.frame(
        files = samples$quant_file[[i]],
        names = sample_id,
        stringsAsFactors = FALSE
    )
    
    txi <- tximeta::tximeta(
        coldata = coldata,
        type = "alevin",
        skipMeta = TRUE,
        skipSeqinfo = TRUE,
        dropInfReps = TRUE,
        alevinArgs = list(dropMeanVar = TRUE)
    )
    
    sce_i <- tximeta::splitSE(txi, feature_groups, assayName = "counts")
    sce_i <- as(sce_i, "SingleCellExperiment")
    
    colnames(sce_i) <- paste0(colnames(sce_i), "__", sample_id)
    colData(sce_i) <- DataFrame(
        experiment = factor(rep(sample_id, ncol(sce_i)), levels = samples$sample_id),
        row.names = colnames(sce_i)
    )
    
    sce_i
}

sce <- do.call(cbind, lapply(seq_len(nrow(samples)), read_alevin_sample))
stopifnot(!anyDuplicated(colnames(sce)))

assays(sce) <- list(
    counts = assay(sce, "spliced"),
    spliced = assay(sce, "spliced"),
    unspliced = assay(sce, "unspliced")
)

row_data <- DataFrame(
    data.frame(gene_id = rownames(sce), stringsAsFactors = FALSE) %>%
        left_join(gene_annot, by = "gene_id")
)
rownames(row_data) <- row_data$gene_id
rowData(sce) <- row_data[rownames(sce), , drop = FALSE]

# QC and filtering

message("Calculating QC metrics")
sce <- scater::addPerCellQC(sce, BPPARAM = bpparam)
sce <- scater::addPerFeatureQC(sce, BPPARAM = bpparam)

min_detected <- as.integer(Sys.getenv("MIN_DETECTED_GENES", "850"))
max_detected <- as.integer(Sys.getenv("MAX_DETECTED_GENES", "7000"))

sce$retain <- sce$detected > min_detected & sce$detected < max_detected

qc_summary <- as.data.frame(table(sce$experiment, sce$retain))
colnames(qc_summary) <- c("sample_id", "retain", "n_cells")
write.table(
    qc_summary,
    file = file.path(qc_dir, "cell_filtering_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

qc_df <- as.data.frame(colData(sce))

p <- ggplot(qc_df, aes(x = experiment, fill = retain)) +
    geom_bar() +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    xlab(NULL) +
    ylab("Number of cells")
save_plot(p, "cells_retained_per_sample.pdf")

p <- ggplot(qc_df, aes(x = experiment, fill = retain)) +
    geom_bar(position = "fill") +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    xlab(NULL) +
    ylab("Fraction of cells")
save_plot(p, "fraction_cells_retained_per_sample.pdf")

p <- ggplot(qc_df, aes(x = sum, fill = retain)) +
    geom_histogram(bins = 100, alpha = 0.5) +
    theme_cowplot() +
    xlab("Total UMI count per cell") +
    ylab("Number of cells")
save_plot(p, "umi_count_histogram.pdf")

p <- ggplot(qc_df, aes(x = experiment, y = sum)) +
    geom_violin(scale = "width") +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    xlab(NULL) +
    ylab("Total UMI count per cell")
save_plot(p, "umi_count_per_sample.pdf")

p <- ggplot(qc_df, aes(x = experiment, y = detected)) +
    geom_violin(scale = "width") +
    geom_hline(yintercept = c(min_detected, max_detected), linetype = 3) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    xlab(NULL) +
    ylab("Number of detected genes per cell")
save_plot(p, "detected_genes_per_sample.pdf")

p <- ggplot(qc_df, aes(x = sum, y = detected, color = retain)) +
    geom_point(size = 1, alpha = 0.15) +
    geom_hline(yintercept = c(min_detected, max_detected), linetype = 3) +
    theme_cowplot() +
    xlab("Total UMI count per cell") +
    ylab("Number of detected genes per cell") +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))
save_plot(p, "umi_count_vs_detected_genes.pdf")

p <- qc_df %>%
    arrange(desc(detected)) %>%
    mutate(rank = row_number()) %>%
    ggplot(aes(x = rank, y = detected, color = retain)) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = c(min_detected, max_detected), linetype = 3) +
    theme_cowplot() +
    xlab("Rank") +
    ylab("Number of detected genes per cell")
save_plot(p, "detected_genes_rank_plot.pdf")

sce <- sce[, sce$retain]

message("Saving filtered SCE: ", output_rds)
saveRDS(sce, output_rds)


