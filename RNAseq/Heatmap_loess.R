## Produce Heatmap of Tip Regeneration data using:
# selected variable genes
# loess smoothing
# Ordering of the data according to correlation to Wnt3 (100203050)

library(ComplexHeatmap)
library(circlize)


## Paths and parameters
input_file <- "ref_data/TipReg_var_genes.tsv"
output_pdf <- "Tips_Heatmap_WntOrdering_LoessSmooth.pdf"

target_gene <- "100203050"
ngrid <- 100L
loess_span <- 0.2

ID_annot <- c(
    "101237052", "100203327", "100211959", "100201451", "100202979",
    "100200795", "105844391", "100192284", "100199257", "100206389", "100197018"
)

Gsymb <- c(
    "FAT4", "CELSR2", "JADE1", "Daam2", "ITPR1",
    "Wnt8", "SFRP3", "Ctnnb1", "TCF", "Axin", "Notum"
)
names(Gsymb) <- ID_annot


## row-wise LOESS smoothing
loessSmoothRows <- function(mat, x, newx,
                            span = 0.4,
                            degree = 2,
                            robust = TRUE,
                            robust.iter = 4L) {
    stopifnot(ncol(mat) == length(x))
    
    res <- matrix(NA_real_, nrow(mat), length(newx))
    rownames(res) <- rownames(mat)
    colnames(res) <- paste0("t", signif(newx, 4))
    
    if (robust) {
        family <- "symmetric"
        iterations <- robust.iter
    } else {
        family <- "gaussian"
        iterations <- 1L
    }
    
    ctrl <- loess.control(surface = "direct", iterations = iterations)
    
    for (i in seq_len(nrow(mat))) {
        y <- mat[i, ]
        ok <- is.finite(y) & is.finite(x)
        
        if (sum(ok) < 3L || length(unique(x[ok])) < 3L || sd(y[ok]) == 0) {
            next
        }
        
        fit <- try(
            loess(
                y ~ x,
                data = data.frame(x = x[ok], y = y[ok]),
                span = span,
                degree = degree,
                family = family,
                control = ctrl
            ),
            silent = TRUE
        )
        
        if (!inherits(fit, "try-error")) {
            res[i, ] <- predict(fit, newdata = data.frame(x = newx))
        }
    }
    
    res
}



## Read expression matrix and smooth over real time
PD <- as.matrix(read.delim(input_file, row.names = 1, check.names = TRUE))

## Extract numeric time from column names, e.g. X28A -> 28.
time_vec <- as.numeric(sub("^X([0-9]+).*", "\\1", colnames(PD)))

if (anyNA(time_vec)) {
    stop("Could not extract numeric time from one or more column names.")
}

ord_cols <- order(time_vec)
PD <- PD[, ord_cols, drop = FALSE]
time_vec <- time_vec[ord_cols]

time_grid <- seq(
    min(time_vec, na.rm = TRUE),
    max(time_vec, na.rm = TRUE),
    length.out = ngrid
)

PD_loess <- loessSmoothRows(
    PD,
    x = time_vec,
    newx = time_grid,
    span = loess_span
)

PD_scaled_sm <- t(scale(t(PD_loess)))


## Order rows by correlation to the target gene
if (!target_gene %in% rownames(PD_scaled_sm)) {
    stop(sprintf("Target gene '%s' not found in rownames(PD_scaled_sm).", target_gene))
}

target_vec <- as.numeric(PD_scaled_sm[target_gene, ])

gene_cor <- apply(PD_scaled_sm, 1, function(x) {
    cor(x, target_vec, use = "pairwise.complete.obs", method = "pearson")
})

gene_cor[is.na(gene_cor)] <- -Inf
RO <- order(gene_cor, decreasing = TRUE)

Mplot <- PD_scaled_sm
Mplot[Mplot > 3] <- 3
Mplot[Mplot < -3] <- -3
Mplot_ord <- Mplot[RO, , drop = FALSE]

## Heatmap annotations
time.axis <- pretty(time_grid, n = 10)
time.axis <- time.axis[time.axis >= min(time_grid) & time.axis <= max(time_grid)]

ta <- HeatmapAnnotation(
    time = anno_mark(
        at = round(scales::rescale(time.axis, to = c(1, ncol(Mplot_ord)))),
        labels = paste0(time.axis, "h"),
        labels_rot = 0,
        padding = grid::unit(1, "mm"),
        link_width = grid::unit(1, "mm"),
        extend = grid::unit(0.5, "mm")
    ),
    gap = grid::unit(1, "points")
)

rows_to_mark <- which(rownames(Mplot_ord) %in% ID_annot)
labels_for_rows <- Gsymb[rownames(Mplot_ord)[rows_to_mark]]

row_ha <- rowAnnotation(
    genes = anno_mark(
        at = rows_to_mark,
        labels = labels_for_rows,
        which = "row",
        side = "right",
        link_width = grid::unit(5, "mm"),
        link_gp = grid::gpar(col = "black", lwd = 0.7),
        labels_gp = grid::gpar(fontsize = 10),
        padding = grid::unit(1, "mm"),
        extend = grid::unit(1, "mm")
    )
)


## Final heatmap
H1 <- Heatmap(
    Mplot_ord,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    col = circlize::colorRamp2(
        c(-2, 0, 2),
        c("royalblue", "grey98", "firebrick")
    ),
    width = grid::unit(9, "inches"),
    show_row_dend = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    bottom_annotation = ta,
    right_annotation = row_ha,
    heatmap_legend_param = list(title = "Regen. GExpr.")
)

pdf(output_pdf, height = 10, width = 11, useDingbats = FALSE)
draw(H1)
dev.off()
