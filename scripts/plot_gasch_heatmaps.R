#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(pheatmap)
    library(grid)
    library(gtable)
})

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

infile <- file.path("data", "gasch2000.txt")

if (!file.exists(infile)) {
    stop("cannot find ", infile,
         "\nplease download the example data into the data/ directory")
}

raw <- fread(infile, sep = "\t", header = TRUE, data.table = FALSE)

rownames(raw) <- raw$UID

# Remove descriptor columns
expr <- raw[, -c(1:3)]

# Force numeric
expr <- apply(expr, 2, as.numeric)
rownames(expr) <- raw$UID

# Remove rows/cols that are all NA
expr <- expr[rowSums(is.na(expr)) < ncol(expr), , drop = FALSE]
expr <- expr[, colSums(is.na(expr)) < nrow(expr), drop = FALSE]

# ------------------------------------------------------------
# Compute CV and take top 10 genes
# ------------------------------------------------------------

cv <- apply(expr, 1, function(x) {
    if (all(is.na(x))) return(NA_real_)
    sd(x, na.rm = TRUE) / abs(mean(x, na.rm = TRUE))
})

top10 <- names(sort(cv, decreasing = TRUE, na.last = NA))[1:10]

expr <- expr[top10, , drop = FALSE]


# ------------------------------------------------------------
# Heatmap function
# ------------------------------------------------------------

draw_heatmap <- function(mat,
                         out_prefix,
                         fontsize_row = 12,
                         fontsize_col = 12,
                         fontsize = 16,
                         main_title = NULL,
                         xlab = "Gene ID",
                         ylab = "Experimental condition",
                         ...) {

    # Remove infinite values
    mat[is.infinite(mat)] <- NA

    # Impute NA with row means
    row_means <- rowMeans(mat, na.rm = TRUE)

    na_idx <- which(is.na(mat), arr.ind = TRUE)

    if (nrow(na_idx) > 0) {
        mat[na_idx] <- row_means[na_idx[, 1]]
    }

    # Remove zero-variance rows/cols
    var_r <- apply(mat, 1, var)
    var_c <- apply(mat, 2, var)

    mat <- mat[var_r > 0, var_c > 0, drop = FALSE]

    if (nrow(mat) < 2 || ncol(mat) < 2) {
        warning("Too few valid entries: ", out_prefix)
        return(invisible(NULL))
    }

    # Transpose: genes -> X, conditions -> Y
    mat <- t(mat)

    # Color palette
    pal <- colorRampPalette(c("white", "blue"))(50)


    # --------------------------------------------------------
    # Plot function
    # --------------------------------------------------------

    open_and_draw <- function(filename) {

        # PNG only
        png(filename,
            width = 1800,
            height = 1800,
            res = 150)

        # Create heatmap object
        ht <- pheatmap(
            mat,
            fontsize = fontsize,
            scale = "row",
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            show_rownames = TRUE,
            show_colnames = TRUE,
            fontsize_row = fontsize_row,
            fontsize_col = fontsize_col,
            angle_col = 45,
            color = pal,
            main = main_title,
            silent = TRUE,
            ...
        )


        # Add space for X label
        ht$gtable <- gtable_add_rows(
            ht$gtable,
            unit(2, "lines"),
            pos = nrow(ht$gtable)
        )

        ht$gtable <- gtable_add_grob(
            ht$gtable,
            textGrob(
                xlab,
                gp = gpar(cex = fontsize / 10)
            ),
            t = nrow(ht$gtable),
            l = 1,
            r = ncol(ht$gtable)
        )


        # Add space for Y label
        ht$gtable <- gtable_add_cols(
            ht$gtable,
            unit(2, "lines"),
            pos = 0
        )

        ht$gtable <- gtable_add_grob(
            ht$gtable,
            textGrob(
                ylab,
                rot = 90,
                gp = gpar(cex = fontsize / 10)
            ),
            t = 1,
            b = nrow(ht$gtable),
            l = 1
        )


        # --------------------------------------------------
        # CRITICAL PART: Scale inside viewport (no clipping)
        # --------------------------------------------------

        grid.newpage()

        pushViewport(
            viewport(
                width  = unit(0.96, "npc"),
                height = unit(0.96, "npc"),
                x = 0.5,
                y = 0.5
            )
        )

        grid.draw(ht$gtable)

        popViewport()

        dev.off()
    }
    open_and_draw(paste0(out_prefix, ".png"))
}


# ------------------------------------------------------------
# Main heatmaps
# ------------------------------------------------------------

out1 <- file.path("results", "heatmap_top10_genes")

draw_heatmap(
    expr,
    out1,
    main_title = "Heat map of top 10 variable genes across all conditions",
    xlab = "Gene identifier",
    ylab = "Experimental condition"
)


# First 30 conditions

if (ncol(expr) >= 30) {

    out_first30 <- file.path("results", "heatmap_first30conds")

    draw_heatmap(
        expr[, 1:30],
        out_first30,
        main_title = "Heat map of first 30 conditions (top 10 genes)",
        xlab = "Gene identifier",
        ylab = "Experimental condition"
    )
}


# Heat shock subset

hs_cols <- grep("^Heat Shock", colnames(expr), value = TRUE)

if (length(hs_cols) > 0) {

    out2 <- file.path("results", "heatmap_heatshock")

    draw_heatmap(
        expr[, hs_cols, drop = FALSE],
        out2,
        main_title = "Heat-induced expression changes (top 10 genes)",
        xlab = "Gene identifier",
        ylab = "Heat-shock condition"
    )
}


# ------------------------------------------------------------
# Category heatmaps
# ------------------------------------------------------------

categories <- c(
    HeatShock = "^Heat Shock",
    Temperature = "shock|deg|Heat",
    Osmotic = "sorbitol",
    Oxidative = "H2O2|Menadione|DTT|diamide",
    Starvation = "starv|Nitrogen|Diauxic|YPD"
)

for (cat in names(categories)) {

    pat <- categories[cat]

    cols <- grep(
        pat,
        colnames(expr),
        ignore.case = TRUE,
        value = TRUE
    )

    if (length(cols) > 1) {

        outc <- file.path("results", paste0("heatmap_", cat))

        draw_heatmap(
            expr[, cols, drop = FALSE],
            outc,
            main_title = paste("Top 10 genes under", cat, "conditions"),
            xlab = "Gene identifier",
            ylab = paste(cat, "condition"),
            # increase font sizes for these smaller subset plots
            fontsize_row = 10,
            fontsize_col = 10
        )
    }
}
# ------------------------------------------------------------

cat("Heat maps written to results/\n")