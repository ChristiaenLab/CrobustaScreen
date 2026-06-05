aovs <- function(dat, grp_vec) {
    # 1. Ensure dat is a data frame and has columns
    dat <- as.data.frame(dat)
    factors <- colnames(dat)
    if (is.null(factors)) stop("Data must have column names.")

    # 2. If grp_vec is a matrix, pick the first column
    if (!is.null(dim(grp_vec))) grp_vec <- grp_vec[, 1]

    # 3. Create a temporary data frame for ANOVA
    # Using a name like '.grp' to avoid collisions with data columns
    tmp_dat <- dat
    tmp_dat$.grp <- as.factor(grp_vec)
    
    # 4. Map over factors
    stats <- factors |>
        map(~ {
            # Use backticks in formula to handle non-standard column names
            form <- as.formula(paste0("`", .x, "` ~ .grp"))
            as.data.frame(summary(aov(form, data = tmp_dat))[[1]])
        })
    
    # 5. Extract results using index-based row selection to be safe
    # Row 1 is the group (.grp), Row 2 is Residuals
    res <- map_dfr(stats, ~ .x[1, ])
    resid <- map_dfr(stats, ~ .x[2, ])
    
    # 6. Final assembly
    colnames(resid) <- paste0(colnames(resid), "_residual")
    res <- cbind(res, resid[, 2:3])
    rownames(res) <- factors
    res$FDR <- p.adjust(res[["Pr(>F)"]], method = "fdr")
    
    return(res)
}

write.aovs <- function(params, z, encoded, groups, clusts, out, ...) {
    # Ensure we use the dataframe part of embeddings and ONE column of clusts
    dat_list <- list(params = params, z = z, embeddings = encoded)
    
    # Pick the best clustering result
    best_clust <- if(!is.null(dim(clusts))) clusts[,1] else clusts
    
    grp_list <- list(
        cluster   = best_clust, 
        condition = groups[, "Condition"], 
        pheno     = groups[, "phenotype.label"]
    )
    
    names(dat_list) |>
        walk(\(x) {
            names(grp_list) |>
                walk(\(y) {
                    # Run and save
                    result <- aovs(dat_list[[x]], grp_list[[y]])
                    filename <- paste0(x, "_", y)
                    # Use a subfolder 'anova' inside 'out'
                    dir.csv(result, filename, path = paste0(out, "/anova"), append.date = F, row.names = T, ...)
                })
        })
}
