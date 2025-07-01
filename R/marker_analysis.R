#' Compute Specificity from Deconvolution Means
#'
#' Calculates the specificity of each gene for each cell type based on estimated
#' cell-type-specific gene expression means. Specificity is calculated by
#' dividing the gene's expression in a given cell type by its total expression across all cell types.
#'
#' @param mat A matrix of estimated gene expression means (genes x cell types).
#' @return A matrix of specificity values (genes x cell types).
#' @export
compute_specificity <- function(mat) {
  mat <- as.matrix(mat)
  row_sums <- rowSums(mat, na.rm = TRUE)
  # Handle genes with zero sum to avoid NaNs (specificity will be 0)
  row_sums[row_sums == 0] <- 1 # Set to 1 to avoid division by zero
  spec_mat <- sweep(mat, 1, row_sums, "/")
  spec_mat[!is.finite(spec_mat)] <- 0 # Convert Inf/NaN to 0
  return(spec_mat)
}

#' Calculate Recall Matrix for Deconvolution Markers
#'
#' Compares identified markers from deconvolution methods against a gold standard
#' marker list and calculates recall rates.
#'
#' @param marker_tables_list A named list of data frames, where each data frame
#'   contains marker information for a deconvolution method (e.g., from EDec, Rodeo).
#'   Each data frame should have a 'Gene' column and columns named `cell_type_marker`
#'   indicating if a gene is a marker for that cell type ("yes"/"no").
#' @param gold_standard_markers A named list, where names are cell types and
#'   values are character vectors of gold standard marker gene names for each cell type.
#' @param cell_types A character vector of all cell types to consider.
#' @param methods A character vector of all deconvolution methods to consider.
#' @return A list containing:
#'   \describe{
#'     \item{recall}{Matrix of recall values (cell types x methods).}
#'     \item{tp}{Matrix of true positives counts (cell types x methods).}
#'     \item{total}{Matrix of total gold standard markers (cell types x methods).}
#'   }
#'   Returns NULL if no gold standard markers are found.
#' @export
calculate_recall_matrix <- function(marker_tables_list, gold_standard_markers, cell_types, methods) {
  if (is.null(gold_standard_markers) || length(gold_standard_markers) == 0 || all(sapply(gold_standard_markers, length) == 0)) {
    message("No gold standard markers found for recall calculation.")
    return(NULL)
  }
  
  recall_mat <- matrix(NA, nrow = length(cell_types), ncol = length(methods),
                       dimnames = list(cell_types, methods))
  tp_mat <- matrix(NA, nrow = length(cell_types), ncol = length(methods),
                   dimnames = list(cell_types, methods))
  total_mat <- matrix(NA, nrow = length(cell_types), ncol = length(methods),
                      dimnames = list(cell_types, methods))
  
  for (method in methods) {
    current_marker_table <- marker_tables_list[[method]]
    if (is.null(current_marker_table) || nrow(current_marker_table) == 0 || !("Gene" %in% colnames(current_marker_table))) next
    
    for (cell_type in cell_types) {
      gs_markers <- gold_standard_markers[[cell_type]]
      if (is.null(gs_markers) || length(gs_markers) == 0) {
        recall_mat[cell_type, method] <- NA # No gold standard to compare against
        tp_mat[cell_type, method] <- NA
        total_mat[cell_type, method] <- NA
        next
      }
      
      marker_col_name <- paste0(cell_type, "_marker")
      if (marker_col_name %in% colnames(current_marker_table)) {
        predicted_markers <- current_marker_table$Gene[current_marker_table[[marker_col_name]] == "yes"]
        
        true_positives <- intersect(predicted_markers, gs_markers)
        recall <- length(true_positives) / length(gs_markers)
        
        recall_mat[cell_type, method] <- recall
        tp_mat[cell_type, method] <- length(true_positives)
        total_mat[cell_type, method] <- length(gs_markers)
      } else {
        # If the marker column doesn't exist, it means no markers were identified for this cell type/method
        recall_mat[cell_type, method] <- 0 # Recall is 0 if no predictions for this cell type
        tp_mat[cell_type, method] <- 0
        total_mat[cell_type, method] <- length(gs_markers)
      }
    }
  }
  # Remove rows/cols that are all NA if no gold standard or no method results
  # If a row/column has any non-NA value, keep it.
  rows_to_keep <- rowSums(!is.na(recall_mat)) > 0
  cols_to_keep <- colSums(!is.na(recall_mat)) > 0
  
  recall_mat_filtered <- recall_mat[rows_to_keep, cols_to_keep, drop = FALSE]
  tp_mat_filtered <- tp_mat[rows_to_keep, cols_to_keep, drop = FALSE]
  total_mat_filtered <- total_mat[rows_to_keep, cols_to_keep, drop = FALSE]
  
  if(nrow(recall_mat_filtered) == 0 || ncol(recall_mat_filtered) == 0) return(NULL)
  return(list(recall = recall_mat_filtered, tp = tp_mat_filtered, total = total_mat_filtered))
}

# Helper function to process marker tables (internal to the package)
# Not exported.
process_marker_table <- function(res_means) {
  if (is.null(res_means) || nrow(res_means) == 0) return(NULL)
  spec <- compute_specificity(res_means) # Calls the exported compute_specificity
  genes <- rownames(spec)
  cells <- colnames(spec)
  bin_mat <- matrix("no", nrow = length(genes), ncol = length(cells), dimnames = list(genes, cells))
  for (cell in cells) bin_mat[spec[, cell] > 0.5, cell] <- "yes"
  colnames(spec) <- paste0(colnames(spec), "_specificity")
  colnames(bin_mat) <- paste0(colnames(bin_mat), "_marker")
  combined <- cbind(spec[rownames(bin_mat), , drop = FALSE], bin_mat)
  data.frame(Gene = rownames(combined), combined, row.names = NULL)
}