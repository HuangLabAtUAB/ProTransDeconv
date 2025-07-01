#' Run EDec Stage 2 for each gene
#'
#' This function iteratively performs deconvolution on each gene across samples using
#' a specified stage 2 deconvolution method (e.g., `EDec::run_edec_stage_2` or a custom one).
#'
#' @param expression_matrix Matrix or data frame of gene expression (rows = genes, columns = samples).
#' @param cell_type_proportions Matrix of cell type proportions (rows = samples, columns = cell types).
#' @param min_samples Minimum number of non-NA samples required to perform deconvolution for a gene.
#' @param method A function that performs stage 2 deconvolution (must accept gene_exp_bulk_samples and cell_type_props).
#'
#' @return A list containing:
#'   \describe{
#'     \item{means}{Matrix of estimated means for each gene and cell type}
#'     \item{std_errors}{Matrix of standard errors}
#'     \item{residuals}{Matrix of residuals for each gene}
#'     \item{explained_variances}{Named numeric vector of R-squared per gene}
#'     \item{degrees_of_freedom}{Named numeric vector of degrees of freedom per gene}
#'   }
#' @export
#' @importFrom stats na.omit
run_EDec_stage2_for_genes <- function(expression_matrix, cell_type_proportions, min_samples = 30, method) {
  if (is.null(expression_matrix) || nrow(expression_matrix) == 0 || ncol(expression_matrix) == 0) {
    message("Expression matrix is empty or NULL. Skipping EDec stage 2 for genes.")
    return(list(means = NULL, std_errors = NULL, residuals = NULL,
                explained_variances = NULL, degrees_of_freedom = NULL))
  }
  
  all_means <- vector("list", length = nrow(expression_matrix))
  all_std_errors <- vector("list", length = nrow(expression_matrix))
  all_explained_variances <- vector("list", length = nrow(expression_matrix))
  all_residuals <- vector("list", length = nrow(expression_matrix))
  all_df <- vector("list", length = nrow(expression_matrix))
  names(all_means) <- rownames(expression_matrix)
  names(all_std_errors) <- rownames(expression_matrix)
  names(all_explained_variances) <- rownames(expression_matrix)
  names(all_residuals) <- rownames(expression_matrix)
  names(all_df) <- rownames(expression_matrix)
  
  for (i in seq_along(rownames(expression_matrix))) {
    gene <- rownames(expression_matrix)[i]
    # Ensure gene_expr is a named numeric vector
    gene_expr <- expression_matrix[gene, , drop = TRUE]
    gene_expr_valid <- gene_expr[!is.na(gene_expr)]
    
    if (length(gene_expr_valid) < min_samples) next
    
    valid_samples <- names(gene_expr_valid)
    # Ensure matched_props is a matrix and not just a vector if only one sample
    matched_props <- cell_type_proportions[valid_samples, , drop = FALSE]
    
    if (nrow(matched_props) == 0 || nrow(matched_props) <= ncol(matched_props)) {
      next # Not enough valid samples or samples <= cell types for deconvolution
    }
    
    stage2_result <- tryCatch({
      method(
        gene_exp_bulk_samples = matrix(gene_expr_valid, nrow = 1, dimnames = list(gene, names(gene_expr_valid))),
        cell_type_props = matched_props
      )
    }, error = function(e) {
      warning(paste("EDec stage 2 failed for gene", gene, ":", e$message))
      NULL
    })
    
    if (!is.null(stage2_result)) {
      all_means[[gene]] <- stage2_result$means
      all_std_errors[[gene]] <- stage2_result$std.errors
      all_explained_variances[[gene]] <- stage2_result$explained.variances
      all_residuals[[gene]] <- stage2_result$residuals
      all_df[[gene]] <- stage2_result$degrees.of.freedom
    }
  }
  # Remove NULL entries
  all_means <- all_means[!sapply(all_means, is.null)]
  all_std_errors <- all_std_errors[!sapply(all_std_errors, is.null)]
  all_explained_variances <- all_explained_variances[!sapply(all_explained_variances, is.null)]
  all_residuals <- all_residuals[!sapply(all_residuals, is.null)]
  all_df <- all_df[!sapply(all_df, is.null)]
  
  if (length(all_means) == 0) {
    message("No genes successfully deconvolved. Returning empty results.")
    return(list(means = NULL, std_errors = NULL, residuals = NULL,
                explained_variances = NULL, degrees_of_freedom = NULL))
  }
  
  # Collect all unique column names for residuals for standardization
  all_residual_cols <- unique(unlist(lapply(all_residuals, colnames)))
  if (is.null(all_residual_cols)) { # If no residuals generated (e.g., all failed)
    all_residual_cols <- character(0)
  }
  
  # Helper for standardizing residuals to a common set of columns (internal)
  standardize_residual_matrix <- function(mat, all_cols) {
    if (is.null(mat) || nrow(mat) == 0) {
      return(matrix(NA, nrow = 1, ncol = length(all_cols), dimnames = list(NULL, all_cols)))
    }
    missing_cols <- setdiff(all_cols, colnames(mat))
    if (length(missing_cols) > 0) {
      temp_mat <- matrix(NA, nrow = nrow(mat), ncol = length(missing_cols), dimnames = list(rownames(mat), missing_cols))
      mat <- cbind(mat, temp_mat)
    }
    mat <- mat[, all_cols, drop = FALSE]
    return(mat)
  }
  # Apply standardization only if there are columns to standardize to
  standardized_residuals <- if (length(all_residual_cols) > 0) {
    lapply(all_residuals, standardize_residual_matrix, all_cols = all_residual_cols)
  } else {
    list()
  }
  
  # Combine results, handling potential NULL or empty matrices
  combined_means <- if(length(all_means) > 0) do.call(rbind, all_means) else NULL
  combined_std_errors <- if(length(all_std_errors) > 0) do.call(rbind, all_std_errors) else NULL
  combined_residuals <- if(length(standardized_residuals) > 0) do.call(rbind, standardized_residuals) else NULL
  combined_explained_variances <- unlist(all_explained_variances)
  combined_df <- unlist(all_df)
  
  list(
    means = combined_means,
    std_errors = combined_std_errors,
    residuals = combined_residuals,
    explained_variances = combined_explained_variances,
    degrees_of_freedom = combined_df
  )
}

#' Batch EDec deconvolution on transformed data
#'
#' @param transformed_list List of transformed expression matrices.
#' @param cell_proportion Cell type proportion matrix.
#' @return Named list of deconvolution results (specifically, the 'means' matrix).
#' @export
#' @importFrom EDec run_edec_stage_2
deconvolve_transformed_list <- function(transformed_list, cell_proportion) {
  if (is.null(cell_proportion) || nrow(cell_proportion) == 0 || ncol(cell_proportion) == 0) {
    message("No cell proportion data provided or it is empty. Deconvolution skipped.")
    return(NULL)
  }
  
  if (!requireNamespace("EDec", quietly = TRUE)) {
    stop("Package 'EDec' is required but not installed. Please install it using install.packages('BiocManager'); BiocManager::install('EDec')")
  }
  
  EDec_results <- list()
  # Define which methods use EDec's default stage 2 vs. the custom one
  use_EDec_ns <- c("Min_score", "Inverse", "Original", "Quantile", "Ratio_shift")
  use_alt_func <- c("Logistic_like", "Tanh", "Min_max")
  
  for (method in names(transformed_list)) {
    mat <- transformed_list[[method]] # Already converted to matrix by ProTransDeconv
    if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
      warning(paste("Skipping", method, "- matrix is NULL or empty."))
      next
    }
    
    # Determine the effective number of samples for min_samples calculation
    common_samples <- intersect(colnames(mat), rownames(cell_proportion))
    effective_num_samples <- length(common_samples)
    
    if (effective_num_samples < 2) {
      warning(paste0("Skipping deconvolution for '", method, "': Not enough common samples (", effective_num_samples, ") with cell proportions."))
      next
    }
    
    # min_samples_for_gene should be at least (num_cell_types + 1) for a unique solution
    min_samples_for_gene <- max(3, ncol(cell_proportion) + 1, floor(effective_num_samples * 0.3))
    
    selected_edec_method <- if (method %in% use_EDec_ns) {
      EDec::run_edec_stage_2
    } else if (method %in% use_alt_func) {
      run_edec_stage_2 # This refers to the custom function defined below
    } else {
      warning(paste("Unknown method:", method, "- using default EDec::run_edec_stage_2"))
      EDec::run_edec_stage_2
    }
    
    result <- tryCatch({
      run_EDec_stage2_for_genes(
        expression_matrix = mat,
        cell_type_proportions = cell_proportion,
        min_samples = min_samples_for_gene,
        method = selected_edec_method
      )
    }, error = function(e) {
      warning(paste("EDec failed for", method, ":", e$message))
      NULL
    })
    
    if (!is.null(result) && !is.null(result$means)) {
      result$means[result$means < 0] <- 0 # Ensure non-negative expression values
      EDec_results[[method]] <- result$means
    } else {
      EDec_results[[method]] <- NULL # Store NULL if deconvolution failed for this method
    }
  }
  
  message("EDec deconvolution finished for all valid transformations.")
  return(EDec_results)
}

#' Custom EDec Stage 2 Implementation (Non-negative least squares with bounds 0-1)
#'
#' This function estimates average expression profiles of constituent cell types
#' by solving constrained least squares problems using quadratic programming.
#' This is a custom implementation inspired by EDec's principles, specifically
#' enforcing solutions between 0 and 1 (if appropriate for the data type).
#'
#' @param gene_exp_bulk_samples Numeric matrix of bulk gene expression (genes x samples).
#' @param cell_type_props Numeric matrix of cell type proportions (samples x cell types).
#' @return A list containing estimated means, standard errors, residuals,
#'   explained variances, and degrees of freedom.
#' @export
#' @importFrom stats na.omit
#' @importFrom quadprog solve.QP
run_edec_stage_2 <- function(gene_exp_bulk_samples, cell_type_props) {
  
  if (!is.matrix(gene_exp_bulk_samples) && !is.data.frame(gene_exp_bulk_samples)) {
    gene_exp_bulk_samples <- as.matrix(gene_exp_bulk_samples)
  }
  
  if (sum(is.na(gene_exp_bulk_samples)) > 0) {
    warning("Your input expression profiles contain NA values.
             Loci with NA values in any samples will not be included
             in the analysis, and will not be present in cell type
             specific methylation profiles.")
    gene_exp_bulk_samples <- as.matrix(stats::na.omit(gene_exp_bulk_samples))
    if (nrow(gene_exp_bulk_samples) == 0) {
      stop("No valid gene expression data remaining after NA omission.")
    }
  }
  
  num_cell_types <- ncol(cell_type_props)
  num_genes <- nrow(gene_exp_bulk_samples)
  num_samples <- nrow(cell_type_props)
  
  if (num_samples <= num_cell_types) {
    stop("Number of samples must be greater than the number of cell types for deconvolution.")
  }
  
  # Constraints for quadratic programming: 0 <= x <= 1
  a_matrix <- cbind(diag(rep(1, num_cell_types)), diag(rep(-1, num_cell_types)))
  b_vector <- c(rep(0, num_cell_types), rep(-1, num_cell_types))
  
  # Pre-calculate Dmat, which is constant for all genes
  d_matrix <- t(cell_type_props) %*% cell_type_props
  
  estimate_cell_type_exp_single_gene <- function(x) {
    d_vector <- as.numeric(x %*% cell_type_props) # Ensure d_vector is numeric
    result <- tryCatch(
      quadprog::solve.QP(Dmat = d_matrix,
                         dvec = d_vector,
                         Amat = a_matrix,
                         bvec = b_vector,
                         meq = 0),
      error = function(e) {
        warning(paste("quadprog::solve.QP failed:", e$message))
        list(solution = rep(NA, num_cell_types))
      }
    )
    result$solution
  }
  
  estimated_cell_type_gene_exp <- t(apply(gene_exp_bulk_samples, 1, estimate_cell_type_exp_single_gene))
  
  rownames(estimated_cell_type_gene_exp) <- rownames(gene_exp_bulk_samples)
  colnames(estimated_cell_type_gene_exp) <- colnames(cell_type_props)
  
  # Compute goodness of fit metrics
  predicted_expression <- estimated_cell_type_gene_exp %*% t(cell_type_props)
  
  common_cols <- intersect(colnames(gene_exp_bulk_samples), colnames(predicted_expression))
  if (length(common_cols) == 0) {
    warning("No common samples between gene expression and predicted expression for residuals calculation. Residuals will be NA.")
    residuals <- matrix(NA, nrow = nrow(gene_exp_bulk_samples), ncol = ncol(gene_exp_bulk_samples),
                        dimnames = dimnames(gene_exp_bulk_samples))
  } else {
    residuals <- gene_exp_bulk_samples[, common_cols, drop = FALSE] -
      predicted_expression[, common_cols, drop = FALSE]
  }
  
  df_val <- max(0, num_samples - num_cell_types)
  
  mean_squared_residuals <- apply(residuals^2, 1, sum, na.rm = TRUE) / df_val
  mean_squared_residuals[!is.finite(mean_squared_residuals)] <- NA
  
  ss_total <- apply(gene_exp_bulk_samples, 1, function(x) sum((x - mean(x, na.rm = TRUE))^2, na.rm = TRUE))
  ss_residual <- apply(residuals^2, 1, sum, na.rm = TRUE)
  
  explained_variances <- 1 - (ss_residual / ss_total)
  explained_variances[!is.finite(explained_variances)] <- NA
  
  # Estimate standard errors
  m <- tryCatch(solve(d_matrix), error = function(e) {
    warning(paste("Failed to solve Dmat for standard error calculation:", e$message))
    matrix(NA, nrow = num_cell_types, ncol = num_cell_types) # Return NA matrix on error
  })
  diag_m <- diag(m)
  
  compute_variances <- function(x){
    variance <- x * diag_m
    return(variance)
  }
  variances <- t(sapply(mean_squared_residuals, compute_variances))
  
  std_devs <- sqrt(variances)
  rownames(std_devs) <- rownames(estimated_cell_type_gene_exp)
  colnames(std_devs) <- colnames(estimated_cell_type_gene_exp)
  
  result <- list(means = estimated_cell_type_gene_exp,
                 std.errors = std_devs,
                 degrees.of.freedom = df_val,
                 explained.variances = explained_variances,
                 residuals = residuals)
  
  return(result)
}

#' Run csSAMfit for each gene
#'
#' This function iteratively applies the `csSAMfit` function from the `csSAM` package
#' to individual genes, handling data preparation and error catching.
#'
#' @param expression_matrix Expression matrix (genes x samples).
#' @param cell_type_proportions Cell type proportion matrix (samples x cell types).
#' @param min_samples Minimum number of non-NA samples required for a gene to be processed.
#' @return A list containing combined basis matrix and a list of basisfit results.
#' @export
#' @importFrom csSAM csSAMfit
run_csSAMfit_for_genes <- function(expression_matrix, cell_type_proportions, min_samples = 30) {
  
  all_basis <- list()
  all_basisfit <- list()
  all_genes <- rownames(expression_matrix)
  
  for (gene in all_genes) {
    
    gene_expr <- expression_matrix[gene, ]
    valid_samples <- !is.na(gene_expr)
    
    if (sum(valid_samples) < min_samples) {
      warning(sprintf("Gene %s skipped: not enough non-NA samples (%d < %d).", gene, sum(valid_samples), min_samples))
      next
    }
    
    gene_expr_clean <- gene_expr[valid_samples]
    cell_props_clean <- cell_type_proportions[valid_samples, , drop = FALSE] # ensure samples x celltypes
    
    if (nrow(cell_props_clean) <= ncol(cell_props_clean)) {
      warning(sprintf("Gene %s skipped: not enough samples (%d) after NA filtering for deconvolution relative to cell types (%d).", gene, nrow(cell_props_clean), ncol(cell_props_clean)))
      next
    }
    
    tryCatch({
      # csSAMfit expects x as genes x samples and cc as celltypes x samples
      # Here, x is a single gene (1 x N_samples) and cc is N_celltypes x N_samples
      cs_result <- csSAM::csSAMfit(
        x = matrix(gene_expr_clean, nrow = 1, dimnames = list(gene, names(gene_expr_clean))),
        cc = t(cell_props_clean), # Transpose cell proportions for csSAMfit
        nperms = 0, # Not performing permutations for differential expression
        fit = "block", # Use block-wise fitting
        verbose = FALSE
      )
      
      this_basis <- cs_result$basis
      this_basis[this_basis < 0] <- 0  # Enforce non-negative values
      # rownames(this_basis) <- gene  # Ensure gene name is preserved
      all_basis[[gene]] <- this_basis
      all_basisfit[[gene]] <- cs_result$basisfit # Store basisfit if needed
      
    }, error = function(e) {
      warning(sprintf("csSAMfit failed for gene %s: %s", gene, e$message))
    })
  }
  
  combined_basis <- if(length(all_basis) > 0) do.call(rbind, all_basis) else NULL
  
  return(list(
    basis = combined_basis,
    basisfit = all_basisfit
  ))
}

#' Run bMIND for each gene
#'
#' This function iteratively applies the `bMIND` function from the `MIND` package
#' to individual genes, handling data preparation and error catching.
#'
#' @param expression_matrix Expression matrix (genes x samples).
#' @param cell_type_proportions Cell type proportion matrix (samples x cell types).
#' @param min_samples Minimum number of valid samples required for a gene to be processed.
#' @return A list containing the estimated mean expression matrix per cell type.
#' @export
#' @importFrom MIND bMIND
run_bMIND_for_genes <- function(expression_matrix, cell_type_proportions,
                                min_samples = 30) {
  all_mu <- lapply(rownames(expression_matrix), function(gene) {
    gene_expr <- expression_matrix[gene, ]
    
    # Find samples with non-NA gene expression
    valid_expr_idx <- !is.na(gene_expr)
    
    # Subset both expression and cell proportion for those samples
    valid_samples_names <- colnames(expression_matrix)[valid_expr_idx]
    matched_props <- cell_type_proportions[valid_samples_names, , drop = FALSE]
    
    # Now also filter out rows with NA or Inf in cell proportions
    keep_props_finite <- apply(matched_props, 1, function(row) all(is.finite(row)))
    final_samples_names <- valid_samples_names[keep_props_finite]
    
    # Final sample counts for this gene (non-NA expression & finite proportions)
    if (length(final_samples_names) < min_samples) {
      # message(sprintf("Gene '%s' skipped for bMIND: Not enough valid samples (%d < %d).", gene, length(final_samples_names), min_samples))
      return(NULL)
    }
    
    gene_expr_filtered <- gene_expr[final_samples_names]
    props_filtered <- cell_type_proportions[final_samples_names, , drop = FALSE]
    
    # Try bMIND
    stage2_result <- tryCatch({
      MIND::bMIND(gene_expr_filtered, props_filtered)
    }, error = function(e) {
      warning(sprintf("bMIND failed for gene '%s': %s", gene, e$message))
      return(NULL)
    })
    
    if (is.null(stage2_result) || is.null(stage2_result$mu)) return(NULL)
    if (any(!is.finite(stage2_result$mu))) return(NULL)
    
    return(stage2_result$mu)
  })
  
  names(all_mu) <- rownames(expression_matrix)
  valid_all_mu <- Filter(Negate(is.null), all_mu)
  
  if (length(valid_all_mu) == 0) return(NULL)
  
  mean_matrix <- do.call(rbind, valid_all_mu)
  rownames(mean_matrix) <- names(valid_all_mu)
  
  # Ensure non-negativity as per typical expression data
  mean_matrix[mean_matrix < 0] <- 1e-6 # Small positive value to avoid zero for logs later
  
  return(list(
    mean = mean_matrix
  ))
}