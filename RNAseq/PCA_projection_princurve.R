## Preprocessing of BulkRNAseq Tip Regeraration data and generation of:
##   1. Simple Projections (PCA, tSNE)
##   2. Principal curve analysis

library(rsvd)
library(matrixStats)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(tibble)
library(princurve)



#### helper functions: 
# 
select_variable_genes <- function(m, f) {
    zeroes <- which(rowSums(m) <= max(1, min(rowSums(m))))
    all.nz.genes <- rownames(m[-zeroes, ])
    m <- m[all.nz.genes, ]
    
    df <- data.frame(
        mean = rowMeans(m + 1 / ncol(m)),
        cv   = apply(m, 1, sd) / rowMeans(m + 1 / ncol(m)),
        var  = apply(m, 1, var)
    )
    
    df$dispersion <- with(df, var / mean)
    df$mean_bin <- with(
        df,
        cut(mean, breaks = c(-Inf, unique(quantile(mean, seq(0.1, 1, 0.02), na.rm = TRUE)), Inf))
    )
    
    var_by_bin <- data.frame(
        mean_bin   = factor(levels(df$mean_bin), levels = levels(df$mean_bin)),
        bin_median = as.numeric(tapply(df$dispersion, df$mean_bin, stats::median)),
        bin_mad    = as.numeric(tapply(df$dispersion, df$mean_bin, stats::mad))
    )[table(df$mean_bin) > 0, ]
    
    df$bin_disp_median <- var_by_bin$bin_median[match(df$mean_bin, var_by_bin$mean_bin)]
    df$bin_disp_mad    <- var_by_bin$bin_mad[match(df$mean_bin, var_by_bin$mean_bin)]
    df$dispersion_norm <- with(df, (dispersion - bin_disp_median) / (bin_disp_mad + 0.01))
    
    n_genes_keep <- ceiling(f * nrow(m))
    disp_cut_off <- sort(df$dispersion_norm, decreasing = TRUE)[n_genes_keep]
    genes_keep <- which(df$dispersion_norm >= disp_cut_off)
    
    rownames(m)[genes_keep]
}

map2color <- function(x, pal, limits = NULL) {
    if (is.null(limits)) limits <- range(x, na.rm = TRUE)
    pal[findInterval(x, seq(limits[1], limits[2], length.out = length(pal) + 1), all.inside = TRUE)]
}

project_to_curve_arclength <- function(points, curve_points) {
    points <- as.matrix(points)
    curve_points <- as.matrix(curve_points)
    
    seg_start <- curve_points[-nrow(curve_points), , drop = FALSE]
    seg_end   <- curve_points[-1, , drop = FALSE]
    seg_vec   <- seg_end - seg_start
    seg_len   <- sqrt(rowSums(seg_vec^2))
    seg_len2  <- seg_len^2
    cum_len   <- c(0, cumsum(seg_len))
    
    out <- numeric(nrow(points))
    
    for (i in seq_len(nrow(points))) {
        p <- matrix(points[i, ], nrow = nrow(seg_start), ncol = ncol(seg_start), byrow = TRUE)
        t_raw <- rowSums((p - seg_start) * seg_vec) / seg_len2
        t_clamped <- pmax(0, pmin(1, t_raw))
        proj <- seg_start + seg_vec * t_clamped
        d2 <- rowSums((proj - p)^2)
        j <- which.min(d2)
        out[i] <- cum_len[j] + t_clamped[j] * seg_len[j]
    }
    
    out
}

p2star <- function(p) {
    if (p < 0.001) "***"
    else if (p < 0.01) "**"
    else if (p < 0.05) "*"
    else "ns"
}



### sample-count matrices and annotation
SampleInfoTips <- readRDS("ref_data/TipRNAseq_SampleInfo.rds")
SampleCountsM <- readRDS("ref_data/TipRNAseq_SampleCounts.rds")

orderSmpl <- SampleInfoTips$SampleName[order(as.numeric(SampleInfoTips$time), SampleInfoTips$SampleName)]
SampleCountsM <- SampleCountsM[, orderSmpl]

## Remove spurious samples and normalize to CPMs.
rem <- grepl("^1B$", colnames(SampleCountsM)) |
    grepl("^27A$", colnames(SampleCountsM)) |
    grepl("^31A$", colnames(SampleCountsM)) |
    grepl("^33B$", colnames(SampleCountsM)) |
    grepl("^24B$", colnames(SampleCountsM))

NSampleCounts <- 1e6 * apply(
    SampleCountsM[-c(1:4), !rem],
    2,
    function(x) x / sum(x[x < quantile(x, 0.98)])
)

LNCounts <- log2(NSampleCounts + 2)
Log_SampleCounts <- LNCounts
times <- as.numeric(gsub("[a-zA-Z]", "", colnames(Log_SampleCounts)))



## Stepwise batch correction of time-series samples by aligning neighboring time-point gene-wise medians/means.
## Log_SampleCounts.corMd is the median-corrected version.
Delta <- rowMedians(LNCounts[, 20:25]) - rowMedians(LNCounts[, 26:31])
Log_SampleCounts.cor <- LNCounts
Log_SampleCounts.cor[, 1:25] <- sweep(Log_SampleCounts.cor[, 1:25], 1, Delta, FUN = "-")

Delta2 <- rowMedians(Log_SampleCounts.cor[, 42:47]) - rowMedians(Log_SampleCounts.cor[, 48:53])
Log_SampleCounts.corMd <- Log_SampleCounts.cor
Log_SampleCounts.corMd[, 1:47] <- sweep(Log_SampleCounts.cor[, 1:47], 1, Delta2, FUN = "-")


## Log_SampleCounts.corMn is the mean-corrected version.
Delta <- rowMeans(LNCounts[, 20:25]) - rowMeans(LNCounts[, 26:31])
Log_SampleCounts.cor <- LNCounts
Log_SampleCounts.cor[, 1:25] <- sweep(Log_SampleCounts.cor[, 1:25], 1, Delta, FUN = "-")

Log_SampleCounts.corMn <- Log_SampleCounts.cor
Delta2 <- rowMeans(Log_SampleCounts.cor[, 42:47]) - rowMeans(Log_SampleCounts.cor[, 48:53])
Log_SampleCounts.corMn[, 1:47] <- sweep(Log_SampleCounts.cor[, 1:47], 1, Delta, FUN = "-")

## Gene set for tip PCA/pseudotime.
nzeroes <- which(rowSums(NSampleCounts > 0) > 4)
expr <- which(rowSums(NSampleCounts) > 8)

nz.genes <- rownames(NSampleCounts[nzeroes, ])
exp.genes <- intersect(rownames(NSampleCounts[expr, ]), nz.genes)

var.genes.strict <- select_variable_genes(NSampleCounts[exp.genes, ], 0.5919632)

Data <- Log_SampleCounts.corMn[var.genes.strict, ]

## Tip PCA and principal-curve pseudotime.
PC <- rsvd::rpca(
    t(Data),
    k = 15,
    retx = TRUE,
    center = TRUE,
    scale = FALSE,
    rand = 0
)

fitl <- summary(lm(times ~ PC$x[, 1] + PC$x[, 2]))
coef <- fitl$coefficients[2:3, 1]
coef <- coef / sum(abs(coef))

Emb <- cbind(
    coef[1] * PC$x[, 1] + coef[2] * PC$x[, 2],
    PC$x[, 4]
)
rownames(Emb) <- colnames(Data)

prinCRV.tip <- princurve::principal_curve(Emb, smoother = "smooth_spline")

time.assigned <- scales::rescale(prinCRV.tip$lambda, c(range(times)))
names(time.assigned) <- colnames(Data)

dlt.time <- abs(times - time.assigned)
names(dlt.time) <- colnames(Data)

rem2 <- which(dlt.time > 6.5)

pt.ord <- rank(prinCRV.tip$lambda[-rem2])
names(pt.ord) <- colnames(Data)[-rem2]

pt.ord.use <- pt.ord
time.filt <- times[-rem2]
names(time.filt) <- colnames(Data)[-rem2]

TipData <- Log_SampleCounts.corMd[, names(pt.ord.use)]



## Segment counts 
SampleInfo <- read.delim("Segmentation_Data/SampleSheet.txt", stringsAsFactors = FALSE)
Segments <- as.vector(unique(SampleInfo$Description))

SampleCountsS <- readRDS("Segmentation_Data/Analysis/RDS_files/Segments_SampleCounts.rds")

remS <- grep("3_8", colnames(SampleCountsS))
if (length(remS) > 0) {
    SampleCountsS <- SampleCountsS[, -remS]
}

NSampleCountsS <- 1e6 * apply(
    SampleCountsS[-c(1:4), ],
    2,
    function(x) x / sum(x[x < quantile(x, 0.95)])
)

expr.genesS <- which(apply(NSampleCountsS, 1, function(x) quantile(x, 0.05)) > 10)

Ntop <- 500
Top_PC1 <- names(sort(abs(PC$rotation[, 1]), decreasing = TRUE)[1:Ntop])
Top_PC2 <- names(sort(abs(PC$rotation[, 4]), decreasing = TRUE)[1:Ntop])
Top_PC <- unique(c(Top_PC1, Top_PC2))

use.genesALL <- intersect(Top_PC, rownames(NSampleCountsS)[expr.genesS])
use.genesALL <- intersect(use.genesALL, rownames(TipData))

ColorSegm <- colorRampPalette(brewer.pal(8, "Dark2"))(length(Segments))
ColorSegm <- paste0(ColorSegm, "AA")
names(ColorSegm) <- Segments



####################### SIMPLE PROJECTIONS (PCA, tSNE) #########################

PCsegm <- prcomp(
    t(log2(NSampleCountsS[use.genesALL, ] + 8)),
    retx = TRUE,
    center = TRUE,
    scale. = FALSE
)

ColorSegmN <- ColorSegm
ColorSegmN[1:8] <- rev(colorRampPalette(c("#053061", "#E6E6FA", "#67001F"))(8))

PC_x.new <- predict(
    PCsegm,
    t(TipData[use.genesALL, names(pt.ord.use)])
)

mypal <- colorRampPalette(
    c("#439722", "#7ABE81", "#D9F0D3", "#E7D4E8", "#A990BB")
)(16)

ColorTime2 <- map2color(pt.ord.use, mypal)
names(ColorTime2) <- names(pt.ord.use)

ColorSegm2N <- paste0(substr(ColorSegmN, 1, 7), "22")
names(ColorSegm2N) <- names(ColorSegmN)

par(mfrow = c(1, 2))

plot(
    PCsegm$x[, 1], PCsegm$x[, 2],
    pch = 20,
    col = ColorSegmN[SampleInfo$Description[match(colnames(NSampleCountsS), SampleInfo$SampleName)]],
    xlab = "PC1",
    ylab = "PC2",
    cex = 3.0
)

plot(
    PCsegm$x[, 1], PCsegm$x[, 2],
    pch = 20,
    col = ColorSegm2N[SampleInfo$Description[match(colnames(NSampleCountsS), SampleInfo$SampleName)]],
    xlab = "PC1",
    ylab = "PC2",
    cex = 3.0
)

points(
    PC_x.new[, 1], PC_x.new[, 2],
    pch = 21,
    cex = 1.5,
    bg = ColorTime2[rownames(PC_x.new)],
    col = "black"
)

## Fit principal curve on segment PCA.
EMBS <- PCsegm$x[!grepl("[FT]_", rownames(PCsegm$x)), ]

prinCRV <- princurve::principal_curve(
    EMBS[, 1:2],
    smoother = "smooth_spline",
    stretch = 10000,
    thresh = 0.1
)

colSegm <- ColorSegmN[SampleInfo$Description[match(rownames(EMBS), SampleInfo$SampleName)]]

OrdS <- prinCRV$s[prinCRV$ord, , drop = FALSE]

lambda_curve <- c(0, cumsum(sqrt(rowSums(diff(OrdS)^2))))
lambda_tip <- project_to_curve_arclength(PC_x.new[, 1:2], OrdS)
lambda_seg <- project_to_curve_arclength(EMBS[, 1:2], OrdS)

new_positions <- 1 - (lambda_tip - min(lambda_curve)) / diff(range(lambda_curve))
new_positions <- pmax(0, pmin(1, new_positions))
names(new_positions) <- rownames(PC_x.new)

new_positionsS <- 1 - (lambda_seg - min(lambda_curve)) / diff(range(lambda_curve))
new_positionsS <- pmax(0, pmin(1, new_positionsS))
names(new_positionsS) <- rownames(EMBS)

xlim <- range(c(EMBS[, 1], PC_x.new[, 1]))
ylim <- range(c(EMBS[, 2], PC_x.new[, 2]))

layout(matrix(c(1, 2), nrow = 1), widths = c(1, 1))

plot(
    EMBS,
    pch = 20,
    col = colSegm,
    xlab = "PC1",
    ylab = "PC2",
    cex = 3.5,
    xlim = xlim,
    ylim = ylim
)

for (i in 1:(nrow(OrdS) - 1)) {
    segments(
        OrdS[i, 1], OrdS[i, 2],
        OrdS[i + 1, 1], OrdS[i + 1, 2],
        col = "red",
        lwd = 4
    )
}

legend(
    -15, 0,
    unique(names(colSegm)),
    col = ColorSegmN[unique(names(colSegm))],
    pch = 20,
    bty = "n",
    pt.cex = 2,
    title.adj = 0.2,
    x.intersp = 0.2,
    y.intersp = 1
)

points(PC_x.new[, 1:2], pch = 16, col = ColorTime2[rownames(PC_x.new)])

plot(
    time.filt[names(new_positions)],
    new_positions,
    pch = 20,
    col = "darkorange",
    xlab = "Experimental time",
    ylab = "Position on princ. curve"
)

fit <- lm(new_positions ~ time.filt[names(new_positions)])
abline(fit, col = "darkgrey", lwd = 2)

fit_summary <- summary(fit)
r_squared <- round(fit_summary$r.squared, 3)
p_value <- signif(fit_summary$coefficients[2, 4], 2)

legend(
    "topleft",
    legend = c(paste0("R²: ", r_squared), paste0("p-value: ", p_value)),
    bty = "n",
    text.col = "black"
)



## Alternative version with loess fit:
layout(matrix(c(1, 2), nrow = 1), widths = c(1, 1))

plot(
    EMBS,
    pch = 20,
    col = colSegm,
    xlab = "PC1",
    ylab = "PC2",
    cex = 3.5,
    xlim = xlim,
    ylim = ylim
)

for (i in 1:(nrow(OrdS) - 1)) {
    segments(
        OrdS[i, 1], OrdS[i, 2],
        OrdS[i + 1, 1], OrdS[i + 1, 2],
        col = "red",
        lwd = 4
    )
}

legend(
    -15, 0,
    unique(names(colSegm)),
    col = ColorSegmN[unique(names(colSegm))],
    pch = 20,
    bty = "n",
    pt.cex = 2,
    title.adj = 0.2,
    x.intersp = 0.2,
    y.intersp = 1
)

points(PC_x.new[, 1:2], pch = 16, col = ColorTime2[rownames(PC_x.new)])

plot(
    time.filt[names(new_positions)],
    new_positions,
    pch = 20,
    col = "darkorange",
    xlab = "Experimental time",
    ylab = "Position on princ. curve"
)

lo <- loess(new_positions ~ time.filt[names(new_positions)], span = 0.4, degree = 2, family = "symmetric")
pred <- predict(lo)
ord <- order(time.filt[names(new_positions)])
lines(time.filt[names(new_positions)][ord], pred[ord], col = "steelblue", lwd = 2)



#Plot v2: boxplot with sign.
SEGMu <- c("B_", "1_", "2_", "3_", "4_", "H_")
names(SEGMu) <- c("Bz", "B1", "B2", "B3", "B4", "H")

MeanSegmP <- rep(NA, length(SEGMu))
names(MeanSegmP) <- names(SEGMu)

for (segm in names(SEGMu)) {
    pat <- SEGMu[segm]
    MeanSegmP[segm] <- mean(new_positionsS[grep(pat, names(new_positionsS))], na.rm = TRUE)
}

ref_vals <- MeanSegmP[is.finite(MeanSegmP)]

phase_cols <- c(
    "Phase A" = "#FDC863",
    "Phase B" = "#80CDC1",
    "Phase C" = "#35978F"
)

BOX_W <- 0.55
POINT_W <- 0.52

df_plot <- tibble(
    time = as.numeric(time.filt[names(new_positions)]),
    coord = as.numeric(new_positions)
) %>%
    mutate(
        group = case_when(
            time <= 15 ~ "Phase A",
            time > 15 & time <= 35 ~ "Phase B",
            TRUE ~ "Phase C"
        ),
        group = factor(group, levels = c("Phase A", "Phase B", "Phase C")),
        x_phase = as.numeric(group)
    ) %>%
    group_by(group) %>%
    mutate(
        time_scaled = ifelse(
            max(time) > min(time),
            (time - min(time)) / (max(time) - min(time)),
            0.5
        ),
        x_point = x_phase + POINT_W * (time_scaled - 0.5)
    ) %>%
    ungroup()

p_ab <- t.test(
    coord ~ group,
    data = filter(df_plot, group %in% c("Phase A", "Phase B"))
)$p.value

p_bc <- t.test(
    coord ~ group,
    data = filter(df_plot, group %in% c("Phase B", "Phase C"))
)$p.value

y_top <- max(df_plot$coord, na.rm = TRUE)

stat_df <- tibble(
    x = c(1, 2),
    xend = c(2, 3),
    y = c(y_top + 0.04, y_top + 0.08),
    label = c(p2star(p_ab), p2star(p_bc))
)

p <- ggplot(df_plot) +
    geom_boxplot(
        aes(x = x_phase, y = coord, group = group, fill = group),
        width = BOX_W,
        outlier.shape = NA,
        alpha = 0.35,
        colour = "grey30"
    ) +
    geom_point(
        aes(x = x_point, y = coord, colour = group),
        size = 1.6,
        alpha = 0.8
    ) +
    geom_segment(
        data = stat_df,
        aes(x = x, xend = xend, y = y, yend = y),
        inherit.aes = FALSE,
        linewidth = 0.6
    ) +
    geom_text(
        data = stat_df,
        aes(x = (x + xend) / 2, y = y + 0.01, label = label),
        inherit.aes = FALSE,
        vjust = 0
    ) +
    geom_hline(
        data = enframe(ref_vals, value = "y"),
        aes(yintercept = y),
        linetype = "dashed",
        colour = "firebrick"
    ) +
    scale_fill_manual(values = phase_cols, guide = "none") +
    scale_colour_manual(values = phase_cols, guide = "none") +
    scale_x_continuous(
        breaks = 1:3,
        labels = levels(df_plot$group),
        expand = expansion(add = 0.3)
    ) +
    scale_y_continuous(
        limits = c(0, max(stat_df$y) + 0.05),
        breaks = sort(ref_vals),
        labels = names(sort(ref_vals))
    ) +
    labs(
        x = NULL,
        y = "Position"
    ) +
    theme_bw(base_size = 12)

print(p)
