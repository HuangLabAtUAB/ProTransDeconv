#' Compute Coefficient of Variation (CV)
#'
#' Computes the Coefficient of Variation (CV) for each row of a matrix.
#' CV is calculated as standard deviation divided by the mean.
#' Rows with mean <= 0 will return NA.
#'
#' @param mat A numeric matrix (genes x samples).
#' @return A numeric vector of CV values for each row.
#' @export
#' @importFrom stats sd
compute_cv <- function(mat) {
  apply(mat, 1, function(x) {
    m <- mean(x, na.rm = TRUE)
    if (m <= 0) return(NA)
    stats::sd(x, na.rm = TRUE) / m
  })
}

#' Min-Score Normalization
#'
#' Normalizes data by subtracting the minimum value and dividing by the standard deviation for each row.
#'
#' @param data A numeric matrix or data frame.
#' @return A normalized matrix.
#' @export
#' @importFrom stats sd
Min_score_normalize <- function(data) {
  data <- as.matrix(data)
  t(apply(data, 1, function(x) {
    min_val <- min(x, na.rm = TRUE)
    sd_val <- stats::sd(x, na.rm = TRUE)
    if (sd_val == 0 || is.na(sd_val)) { # Handle cases with zero standard deviation
      normed <- rep(0, length(x))
    } else {
      normed <- (x - min_val) / sd_val
    }
    normed[is.na(x)] <- NA
    return(normed)
  }))
}

#' Min-Max Normalization
#'
#' Scales data to a specified range (default 0 to 1) column-wise.
#'
#' @param data A numeric matrix or data frame.
#' @param a Lower bound of the target range.
#' @param b Upper bound of the target range.
#' @return A normalized matrix.
#' @export
Min_max_normalize <- function(data, a = 0, b = 1) {
  data <- as.matrix(data)
  apply(data, 2, function(x) {
    rng <- range(x, na.rm = TRUE)
    if (rng[2] == rng[1]) { # Handle case where range is zero
      scaled <- rep(a, length(x))
    } else {
      scaled <- (x - rng[1]) / (rng[2] - rng[1]) * (b - a) + a
    }
    scaled[is.na(x)] <- NA
    return(scaled)
  })
}

#' Logistic-like Transformation
#'
#' Applies a logistic-like transformation (1 - exp(-a*x)).
#'
#' @param x A numeric vector.
#' @return A transformed numeric vector.
#' @export
Logistic_like <- function(x) {
  result <- rep(NA, length(x))
  if (all(is.na(x))) return(result)
  max_val_x <- max(x, na.rm = TRUE)
  if (max_val_x <= 0 || !is.finite(max_val_x)) { # Avoid division by zero or Inf
    a <- 0 # If max_val_x is problematic, 'a' becomes 0, so 1 - exp(0) = 0.
  } else {
    a <- 1 / max_val_x
  }
  result[!is.na(x)] <- 1 - exp(-a * x[!is.na(x)])
  return(result)
}

#' Scaled Tanh Transformation
#'
#' Applies a scaled hyperbolic tangent transformation ((tanh(x) + 1) / 2).
#'
#' @param x A numeric vector.
#' @return A transformed numeric vector.
#' @export
Tanh_scaled <- function(x) {
  result <- rep(NA, length(x))
  result[!is.na(x)] <- (tanh(x[!is.na(x)]) + 1) / 2
  return(result)
}

#' Quantile Normalization
#'
#' Performs quantile normalization on a matrix.
#' Requires the 'preprocessCore' package.
#'
#' @param mat A numeric matrix.
#' @return A quantile normalized matrix.
#' @export
#' @importFrom preprocessCore normalize.quantiles
Quantile_normalize <- function(mat) {
  if (!requireNamespace("preprocessCore", quietly = TRUE)) {
    stop("Package 'preprocessCore' is required for Quantile normalization. Please install it using install.packages('BiocManager'); BiocManager::install('preprocessCore')")
  }
  mat <- as.matrix(mat)
  normed <- preprocessCore::normalize.quantiles(mat)
  dimnames(normed) <- dimnames(mat)
  return(normed)
}

#' Ratio Shift Transformation
#'
#' Shifts values by subtracting the minimum value for each row.
#'
#' @param data A numeric matrix or data frame.
#' @return A shifted matrix.
#' @export
Ratio_shift <- function(data) {
  data <- as.matrix(data)
  t(apply(data, 1, function(x) {
    min_val <- min(x, na.rm = TRUE)
    shifted <- x - min_val
    shifted[is.na(x)] <- NA
    return(shifted)
  }))
}

# Helper function for logistic/spectra count normalization (internal to the package)
# Not exported, so no @export tag needed.
process_logistic_spectra <- function(data) {
  dat <- as.matrix(data)
  original_colnames <- colnames(dat)
  original_rownames <- rownames(dat)
  
  dat <- t(apply(t(dat), 2, Logistic_like))
  
  # Handle rows with all NAs or zeros after transformation
  valid_rows <- rowSums(dat, na.rm = TRUE) > 0
  if (sum(valid_rows) == 0) {
    warning("No valid rows remaining after logistic/spectra count transformation. Returning empty matrix.")
    return(matrix(NA, nrow = 0, ncol = ncol(dat), dimnames = list(NULL, original_colnames)))
  }
  dat <- dat[valid_rows, , drop = FALSE]
  
  # Avoid division by zero if max is 0 or NA
  max_val <- max(dat, na.rm = TRUE)
  if (is.finite(max_val) && max_val > 0) {
    dat[!is.na(dat)] <- dat[!is.na(dat)] / max_val
  } else {
    warning("Maximum value after logistic/spectra count transformation is zero or non-finite. Skipping final scaling.")
  }
  
  colnames(dat) <- original_colnames
  rownames(dat) <- original_rownames[valid_rows]
  
  return(dat)
}

# Helper function for compute_cv_summary (internal to the package)
# Not exported.
compute_cv_summary <- function(mat, use_markers_only, marker_genes) {
  if (use_markers_only && !is.null(marker_genes)) {
    mat <- mat[rownames(mat) %in% marker_genes, , drop = FALSE]
  }
  # Handle case where mat might be empty after filtering
  if (nrow(mat) == 0) return(numeric(0))
  vals <- compute_cv(mat)
  vals[is.finite(vals)]
}