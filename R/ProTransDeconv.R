#' Process and Report CV with Optional Deconvolution
#'
#' Computes CV across multiple normalization methods, generates a report,
#' and performs deconvolution if cell type proportions are provided.
#'
#' @param data Expression matrix (genes x samples)
#' @param type One of "intensity", "ratio", or "spectra count"
#' @param run_bmind Logical, whether to run bMIND deconvolution. Default is FALSE.
#' @param cv_threshold Threshold for flagging high CV. Default is 0.25.
#' @param output_html Path to save the generated HTML report. Default is "summary_report.html" in current working directory.
#' @param cor_threshold_gold_standard Correlation threshold for gold standard marker identification. Default is 0.5.
#' @param pval_threshold_gold_standard P-value threshold for gold standard marker identification. Default is 0.05.
#' @param cell_proportion Optional cell type proportion matrix for deconvolution (samples x cell types).
#'
#' @return A list containing transformed data, CV summary, ridge plot data,
#'         deconvolution results (EDec, Rodeo, csSAM, bMIND),
#'         gene-cell correlation, and identified gold standard markers.
#' @export
#' @import ggplot2
#' @import ggridges
#' @import dplyr
#' @import knitr
#' @import rmarkdown
#' @import Rodeo
#' @import csSAM
#' @import MASS
#' @import MIND
#' @importFrom purrr map_dfr
#' @importFrom reshape2 melt
#' @importFrom parallel detectCores makeCluster clusterExport clusterEvalQ parLapply stopCluster
#' @importFrom grDevices png dev.off
#' @importFrom graphics plot
#' @importFrom stats sd median complete.cases cor.test quantile na.omit cmdscale
#' @importFrom utils head tail
ProTransDeconv <- function(data,
                           type = c("intensity", "ratio", "spectra count"),
                           run_bmind = FALSE,
                           cv_threshold = 0.25, cor_threshold_gold_standard = 0.5,
                           output_html = NULL, pval_threshold_gold_standard = 0.05,
                           cell_proportion = NULL) {
  # Package startup messages are suppressed globally by the R package system
  # when the package is loaded, so `suppressPackageStartupMessages`
  # is not strictly needed here for the `library` calls.
  # However, it's good practice to ensure all used packages are in DESCRIPTION Imports.
  
  type <- match.arg(type)
  
  if (is.null(output_html)) {
    output_html <- file.path(getwd(), "summary_report.html")
  }
  
  dir_path <- dirname(output_html)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
  
  output_html <- tryCatch(
    normalizePath(output_html, mustWork = FALSE),
    error = function(e) output_html
  )
  
  # Check if input data is empty
  if (is.null(data) || nrow(data) == 0 || ncol(data) == 0) {
    warning("Input data is empty or NULL. Skipping transformations and downstream analysis.")
    return(invisible(list(
      transformed_list = list(),
      cv_summary = data.frame(),
      cv_plot_data = data.frame(),
      EDec = list(EDec_results = NULL, EDec_marker_table = NULL),
      rodeo = list(rodeo_results = NULL, rodeo_marker_table = NULL),
      csSAM = list(csSAM_result = NULL, csSAM_marker_table = NULL),
      bMIND = list(bMIND_result = NULL, bMIND_marker_table = NULL),
      gene_cell_correlation = NULL,
      gold_standard_markers = list()
    )))
  }
  
  # --- Data Transformation Methods ---
  transform_methods <- list(
    Min_score = Min_score_normalize,
    Min_max = Min_max_normalize
  )
  
  if (type == "intensity") {
    transform_methods$Logistic_like <- process_logistic_spectra
    transform_methods$Quantile <- Quantile_normalize
    transform_methods$Inverse <- function(data) 2 ^ as.matrix(data)
    transform_methods$Original <- as.matrix
  }
  
  if (type == "ratio") {
    transform_methods$Ratio_shift <- Ratio_shift
    transform_methods$Inverse <- function(data) 2 ^ as.matrix(data)
    transform_methods$Tanh <- function(data) {
      dat <- as.matrix(data)
      t(apply(dat, 1, Tanh_scaled))
    }
  }
  
  if (type == "spectra count") {
    transform_methods$Logistic_like <- process_logistic_spectra
    transform_methods$Quantile <- Quantile_normalize
    transform_methods$Original <- as.matrix
  }
  
  transformed_list <- lapply(transform_methods, function(f) {
    res <- tryCatch(f(data), error = function(e) {
      warning(paste("Transformation failed for method:", deparse(substitute(f)), "-", e$message))
      return(NULL)
    })
    # Ensure it's a matrix with dimnames, even if empty or failed
    if (is.null(res) || !is.matrix(res)) {
      res <- matrix(NA, nrow = 0, ncol = ncol(data), dimnames = list(NULL, colnames(data)))
    }
    res
  })
  # Remove any failed or empty transformations
  transformed_list <- transformed_list[!sapply(transformed_list, is.null)]
  transformed_list <- transformed_list[sapply(transformed_list, function(x) nrow(x) > 0)]
  
  if (length(transformed_list) == 0) {
    warning("No valid transformations were performed. Skipping CV calculation and reporting.")
    return(invisible(list(
      transformed_list = list(),
      cv_summary = data.frame(),
      cv_plot_data = data.frame(),
      EDec = list(EDec_results = NULL, EDec_marker_table = NULL),
      rodeo = list(rodeo_results = NULL, rodeo_marker_table = NULL),
      csSAM = list(csSAM_result = NULL, csSAM_marker_table = NULL),
      bMIND = list(bMIND_result = NULL, bMIND_marker_table = NULL),
      gene_cell_correlation = NULL,
      gold_standard_markers = list()
    )))
  }
  
  
  # --- CV Calculation and Summary ---
  prop_col_name <- paste0("Prop_GT_", cv_threshold)
  
  cv_summary <- lapply(names(transformed_list), function(method) {
    vals <- compute_cv_summary(transformed_list[[method]])
    if (length(vals) == 0) return(NULL) # Skip if no CV values computed
    qtiles <- quantile(vals, probs = c(0.25, 0.75), na.rm = TRUE)
    row <- data.frame(
      Method = method,
      Group = type,
      Mean = round(mean(vals), 3),
      Median = round(median(vals), 3),
      Min = round(min(vals), 3),
      Max = round(max(vals), 3),
      Q1 = round(qtiles[1], 3),
      Q3 = round(qtiles[2], 3),
      stringsAsFactors = FALSE
    )
    row[[prop_col_name]] <- round(mean(vals > cv_threshold), 3)
    row
  }) %>% bind_rows()
  
  if (nrow(cv_summary) == 0) {
    warning("CV summary is empty. Skipping CV plot and report generation.")
    return(invisible(list(
      transformed_list = transformed_list,
      cv_summary = data.frame(),
      cv_plot_data = data.frame(),
      EDec = list(EDec_results = NULL, EDec_marker_table = NULL),
      rodeo = list(rodeo_results = NULL, rodeo_marker_table = NULL),
      csSAM = list(csSAM_result = NULL, csSAM_marker_table = NULL),
      bMIND = list(bMIND_result = NULL, bMIND_marker_table = NULL),
      gene_cell_correlation = NULL,
      gold_standard_markers = list()
    )))
  }
  rownames(cv_summary) <- seq_len(nrow(cv_summary))
  
  cv_df <- lapply(names(transformed_list), function(method) {
    vals <- compute_cv_summary(transformed_list[[method]] )
    if (length(vals) == 0) return(NULL)
    data.frame(Method = method, CV = vals, Group = type)
  }) %>% bind_rows() %>% na.omit()
  
  # Ensure cv_df is not empty before arranging levels
  if (nrow(cv_df) > 0) {
    method_order <- cv_summary %>%
      arrange(desc(Median)) %>%
      pull(Method)
    cv_df$Method <- factor(cv_df$Method, levels = method_order)
    cv_summary$Method <- factor(cv_summary$Method, levels = method_order)
  }
  
  
  # --- Deconvolution Results Initialization ---
  EDec_results <- NULL
  EDec_marker_table <- NULL
  rodeo_results <- list()
  csSAM_results <- list()
  bMIND_results <- list()
  rodeo_marker_table <- list()
  csSAM_marker_table <- list()
  bMIND_marker_table <- list()
  
  combined_mat <- NULL # Initialize combined_mat for the correlation results
  gold_standard_markers <- list() # Initialize gold_standard_markers
  
  
  # --- EDec Deconvolution ---
  if (!is.null(cell_proportion) && nrow(cell_proportion) > 0) {
    EDec_results_full <- deconvolve_transformed_list(transformed_list, cell_proportion)
    EDec_marker_table <- lapply(names(EDec_results_full), function(method) {
      res <- EDec_results_full[[method]]
      process_marker_table(res) # res here is the 'means' matrix from EDec_results_full
    })
    names(EDec_marker_table) <- names(EDec_results_full)
    # Store just the means for overall results
    # EDec_results <- lapply(EDec_results_full, `[[`, "means")
    EDec_results <- EDec_results_full
    
  } else {
    message("No cell proportion data provided or it is empty. Skipping EDec deconvolution.")
  }

  # --- Rodeo Deconvolution ---
  if (requireNamespace("Rodeo", quietly = TRUE)) {
    for (method in names(transformed_list)) {
      mat <- transformed_list[[method]]
      if (is.null(mat) || nrow(mat) == 0 || is.null(cell_proportion) || nrow(cell_proportion) == 0) {
        next
      }
      library(MASS)
      tryCatch({

        result <- Rodeo::Rodeo(mat, t(cell_proportion))
        if (!is.null(result)) {
          spec <- compute_specificity(result)
          genes <- rownames(spec)
          cells <- colnames(spec)
          bin_mat <- matrix("no", nrow = length(genes), ncol = length(cells), dimnames = list(genes, cells))
          for (cell in cells) bin_mat[spec[, cell] > 0.5, cell] <- "yes"
          colnames(spec) <- paste0(colnames(spec), "_specificity")
          colnames(bin_mat) <- paste0(colnames(bin_mat), "_marker")
          combined <- cbind(spec[rownames(bin_mat), , drop = FALSE], bin_mat)
          marker_table <- data.frame(Gene = rownames(combined), combined, row.names = NULL)
          rodeo_results[[method]] <- result
          rodeo_marker_table[[method]] <- marker_table
        }
        
        else {
          warning("Rodeo result is NULL for method: ", method)
        }
        
      }, error = function(e) {
        message("Rodeo failed for method '", method, "': ", conditionMessage(e))
      }
      )
      
      
    }
    message("Rodeo deconvolution finished for all valid transformations.")
  } else {
    warning("Package 'Rodeo' is not installed. Skipping Rodeo deconvolution.")
  }

  
  # --- csSAM Deconvolution ---
  if (requireNamespace("csSAM", quietly = TRUE)) {
    for (method in names(transformed_list)) {
     
      mat <- transformed_list[[method]]
      
      if (is.null(mat) || nrow(mat) == 0 || is.null(cell_proportion) || nrow(cell_proportion) == 0) next
  
      tryCatch({
        # csSAMfit expects x as genes x samples and cc as celltypes x samples
        common_samples_csSAM <- intersect(colnames(mat), rownames(cell_proportion))
        
        if (length(common_samples_csSAM) < 2) {
          message(paste0("Skipping csSAM for method '", method, "': Not enough common samples (min 2 required)."))
          next
        }
        
        aligned_mat <- mat[, common_samples_csSAM, drop = FALSE]
        aligned_cell_props <- cell_proportion[common_samples_csSAM, , drop = FALSE]  # samples x celltypes
        
        # Run csSAM fit
        result_list <- run_csSAMfit_for_genes(
          expression_matrix = aligned_mat,
          cell_type_proportions = aligned_cell_props,
          min_samples = 30
        )
        
        if (!is.null(result_list) && !is.null(result_list$basis)) {
          marker_table <- process_marker_table(result_list$basis)
          csSAM_results[[method]] <- result_list$basis
          csSAM_marker_table[[method]] <- marker_table
        }
      }, error = function(e) {
        warning(paste0("csSAM failed for method '", method, "': ", conditionMessage(e)))
      })
      
    }
    message("csSAM deconvolution finished for all valid transformations.")
  } else {
    warning("Package 'csSAM' is not installed. Skipping csSAM deconvolution.")
  }
  
  
  # --- bMIND Deconvolution ---
  if (run_bmind && requireNamespace("MIND", quietly = TRUE)) {
    num_cores_to_use <- max(1, parallel::detectCores() - 1) # Leave one core free
    cl <- parallel::makeCluster(num_cores_to_use)
    on.exit(parallel::stopCluster(cl), add = TRUE) # Ensure cluster stops on exit
    
    # Export necessary objects and functions to the cluster
    parallel::clusterExport(cl, varlist = c(
      "transformed_list", "cell_proportion", "run_bMIND_for_genes",
      "bMIND" # Explicitly export MIND::bMIND if needed by worker nodes
    ), envir = environment())
    parallel::clusterEvalQ(cl, library(MIND)) # Load MIND on workers
    
    # Filter out "Inverse" transformation for bMIND as suggested in the original code
    transformed_list_bMIND <- transformed_list
    if ("Inverse" %in% names(transformed_list_bMIND)) {
      transformed_list_bMIND[["Inverse"]] <- NULL
    }
    
    bmind_raw_results <- parallel::parLapply(cl, names(transformed_list_bMIND), function(method) {
      mat <- transformed_list_bMIND[[method]]
      if (is.null(mat) || nrow(mat) == 0 || is.null(cell_proportion) || nrow(cell_proportion) == 0) return(NULL)
      
      # Ensure common samples and sufficient samples for bMIND
      common_samples <- intersect(colnames(mat), rownames(cell_proportion))
      if (length(common_samples) < 30) { # bMIND needs sufficient samples
        message(sprintf("Skipping bMIND for method %s: Not enough common samples (%d < 30).", method, length(common_samples)))
        return(NULL)
      }
      
      mat_use <- as.data.frame(mat[, common_samples, drop = FALSE])
      props_use <- cell_proportion[common_samples, , drop = FALSE]
      
      tryCatch({
        result <- run_bMIND_for_genes(mat_use, props_use, min_samples = 30)
        if (is.null(result) || is.null(result$mean)) return(NULL)
        list(method = method, result = result)
      }, error = function(e) {
        warning(sprintf("bMIND failed for method %s: %s", method, e$message))
        return(NULL)
      })
    })
    
    # Process bMIND raw results
    for (res in bmind_raw_results) {
      if (!is.null(res)) {
        method <- res$method
        result <- res$result
        spec <- compute_specificity(result$mean)
        
        genes <- rownames(spec)
        cells <- colnames(spec)
        bin_mat <- matrix("no", nrow = length(genes), ncol = length(cells), dimnames = list(genes, cells))
        for (cell in cells) bin_mat[spec[, cell] > 0.5, cell] <- "yes"
        
        colnames(spec) <- paste0(cells, "_specificity")
        colnames(bin_mat) <- paste0(cells, "_marker")
        combined <- cbind(spec[rownames(bin_mat), , drop = FALSE], bin_mat)
        marker_table <- data.frame(Gene = rownames(combined), combined, row.names = NULL)
        
        bMIND_results[[method]] <- result$mean
        bMIND_marker_table[[method]] <- marker_table
      }
    }
    message("bMIND deconvolution finished for all valid transformations.")
    
  } else if (!run_bmind) {
    message("Skipping bMIND analysis because run_bmind = FALSE.")
  } else {
    warning("Package 'MIND' is not installed. Skipping bMIND deconvolution.")
  }
  
  
  # === Spearman correlation between Original data and cell type proportions ===
  # This part of the code correctly identifies gold standard markers
  # based on correlation with cell type proportions.
  if (!is.null(cell_proportion) && nrow(cell_proportion) > 0) {
    common_samples <- intersect(colnames(data), rownames(cell_proportion))
    if (length(common_samples) >= 3) {
      expr_aligned <- data[, common_samples, drop = FALSE]
      cell_props_aligned <- cell_proportion[common_samples, , drop = FALSE]
      
      cor_mat <- matrix(NA, nrow = nrow(expr_aligned), ncol = ncol(cell_props_aligned),
                        dimnames = list(rownames(expr_aligned), colnames(cell_props_aligned)))
      pval_mat <- matrix(NA, nrow = nrow(expr_aligned), ncol = ncol(cell_props_aligned),
                         dimnames = list(rownames(expr_aligned), colnames(cell_props_aligned)))
      
      for (i in seq_len(nrow(expr_aligned))) {
        gene_exp_vec <- as.numeric(expr_aligned[i, ])
        for (j in seq_len(ncol(cell_props_aligned))) {
          cell_prop_vec <- as.numeric(cell_props_aligned[, j])
          valid <- complete.cases(gene_exp_vec, cell_prop_vec)
          if (sum(valid) >= 3) {
            cor_test <- suppressWarnings(cor.test(gene_exp_vec[valid], cell_prop_vec[valid], method = "spearman"))
            cor_mat[i, j] <- cor_test$estimate
            pval_mat[i, j] <- cor_test$p.value
          }
        }
      }
      
      col_names_props <- colnames(cell_props_aligned)
      
      combined_mat <- do.call(cbind, lapply(col_names_props, function(cell) {
        out <- cbind(cor_mat[, cell], pval_mat[, cell])
        colnames(out) <- c(paste0(cell, "_cor"), paste0(cell, "_pval"))
        return(out)
      }))
      rownames(combined_mat) <- rownames(expr_aligned)
      
      for (cell_type in col_names_props) {
        cor_col_name <- paste0(cell_type, "_cor")
        pval_col_name <- paste0(cell_type, "_pval")
        
        if (cor_col_name %in% colnames(combined_mat) && pval_col_name %in% colnames(combined_mat)) {
          current_cor_vals <- combined_mat[, cor_col_name]
          current_pval_vals <- combined_mat[, pval_col_name]
          
          gold_standard_for_cell <- rownames(combined_mat)[
            !is.na(current_cor_vals) & !is.na(current_pval_vals) &
              current_cor_vals > cor_threshold_gold_standard &
              current_pval_vals < pval_threshold_gold_standard
          ]
          if (length(gold_standard_for_cell) > 0) {
            gold_standard_markers[[cell_type]] <- gold_standard_for_cell
          } else {
            gold_standard_markers[[cell_type]] <- character(0)
          }
        } else {
          gold_standard_markers[[cell_type]] <- character(0)
        }
      }
      message("Gold standard marker identification complete.")
      
    } else {
      warning("Not enough overlapping samples (at least 3) to compute Spearman correlation.")
    }
  } else {
    message("No cell proportion data provided for Spearman correlation calculation.")
  }
  
  
  # --- Calculate Recall and Generate Heatmaps ---
  recall_heatmap_files <- list()
  
  if (!is.null(cell_proportion) && nrow(cell_proportion) > 0 &&
      length(gold_standard_markers) > 0 && any(sapply(gold_standard_markers, length) > 0)) {
    
    all_cell_types <- colnames(cell_proportion)
    all_methods <- names(transformed_list)
    
    # EDec Recall
    if (length(EDec_marker_table) > 0 && !all(sapply(EDec_marker_table, is.null))) {
      message("Calculating EDec recall...")
      EDec_recall_list <- calculate_recall_matrix(EDec_marker_table, gold_standard_markers, all_cell_types, all_methods)
      if (!is.null(EDec_recall_list) && !all(is.na(EDec_recall_list$recall))) {
        recall_heatmap_files$EDec <- plot_recall_heatmap(EDec_recall_list$recall, "EDec", EDec_recall_list$tp, EDec_recall_list$total)
      }
    } else {
      message("EDec marker table is empty or contains no valid results. Skipping EDec recall calculation.")
    }
    
    # Rodeo Recall
    if (length(rodeo_marker_table) > 0 && !all(sapply(rodeo_marker_table, is.null))) {
      message("Calculating Rodeo recall...")
      rodeo_recall_list <- calculate_recall_matrix(rodeo_marker_table, gold_standard_markers, all_cell_types, all_methods)
      if (!is.null(rodeo_recall_list) && !all(is.na(rodeo_recall_list$recall))) {
        recall_heatmap_files$rodeo <- plot_recall_heatmap(rodeo_recall_list$recall, "Rodeo", rodeo_recall_list$tp, rodeo_recall_list$total)
      }
    } else {
      message("Rodeo marker table is empty or contains no valid results. Skipping Rodeo recall calculation.")
    }
    
    # csSAM Recall
    if (length(csSAM_marker_table) > 0 && !all(sapply(csSAM_marker_table, is.null))) {
      message("Calculating csSAM recall...")
      csSAM_recall_list <- calculate_recall_matrix(csSAM_marker_table, gold_standard_markers, all_cell_types, all_methods)
      if (!is.null(csSAM_recall_list) && !all(is.na(csSAM_recall_list$recall))) {
        recall_heatmap_files$csSAM <- plot_recall_heatmap(csSAM_recall_list$recall, "csSAM", csSAM_recall_list$tp, csSAM_recall_list$total)
      }
    } else {
      message("csSAM marker table is empty or contains no valid results. Skipping csSAM recall calculation.")
    }
    
    # bMIND Recall
    if (run_bmind && length(bMIND_marker_table) > 0 && !all(sapply(bMIND_marker_table, is.null))) {
      message("Calculating bMIND recall...")
      bMIND_recall_list <- calculate_recall_matrix(bMIND_marker_table, gold_standard_markers, all_cell_types, all_methods)
      if (!is.null(bMIND_recall_list) && !all(is.na(bMIND_recall_list$recall))) {
        recall_heatmap_files$bMIND <- plot_recall_heatmap(bMIND_recall_list$recall, "bMIND", bMIND_recall_list$tp, bMIND_recall_list$total)
      }
    } else if (!run_bmind) {
      message("Skipping bMIND recall because run_bmind = FALSE.")
    } else {
      message("bMIND marker table is empty or contains no valid results. Skipping bMIND recall calculation.")
    }
    
  } else {
    message("Skipping recall heatmap generation (no cell proportion data, or no gold standard markers identified).")
  }
  
  
  # --- Generate CV Ridge Plot ---
  cv_ridge_plot_file <- tempfile(fileext = ".png")
  if (nrow(cv_df) > 0) { # Only plot if data exists
    suppressWarnings({
      p <- ggplot(cv_df, aes(x = CV, y = Method, fill = Method)) +
        ggridges::geom_density_ridges(scale = 1.5, alpha = 0.65, rel_min_height = 0.01,
                                      color = "gray20", linewidth = 0.4) +
        ggplot2::facet_grid(. ~ Group, scales = "fixed") +
        ggplot2::coord_cartesian(xlim = c(0, 4)) +
        ggplot2::scale_fill_brewer(palette = "Set2") +
        ggridges::theme_ridges(font_size = 12, grid = TRUE) +
        ggplot2::labs(title = "CV Distribution by Method",
                      x = "Coefficient of Variation (CV)", y = NULL) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16),
          strip.text = ggplot2::element_text(size = 12, face = "bold", color = "#333333"),
          axis.text.y = ggplot2::element_text(size = 10, face = "bold", color = "black"),
          axis.text.x = ggplot2::element_text(size = 10, face = "bold", color = "black"),
          axis.title.x = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"),
          legend.position = "none",
          panel.spacing = ggplot2::unit(1.2, "lines")
        )
    })
    ggplot2::ggsave(cv_ridge_plot_file, p, width = 6.5, height = 4, dpi = 150)
  } else {
    message("No CV data available for plotting ridge plot.")
    cv_ridge_plot_file <- NULL # Set to NULL if no plot generated
  }
  
  
  # --- Generate Marker Barplots ---
  marker_barplot_files <- list()
  marker_barplot_files$EDec <- save_marker_barplot(EDec_marker_table, "EDec")
  marker_barplot_files$rodeo <- save_marker_barplot(rodeo_marker_table, "Rodeo")
  marker_barplot_files$csSAM <- save_marker_barplot(csSAM_marker_table, "csSAM")
  if (run_bmind) {
    marker_barplot_files$bMIND <- save_marker_barplot(bMIND_marker_table, "bMIND")
  }
  
  
  # --- Prepare and Render R Markdown Report ---
  # Only write HTML report if output_html is not NULL (and not an empty string)
  # Only write HTML report if output_html is not NULL (and not an empty string)
  if (!is.null(output_html) && output_html != "") {
    tmp_rmd <- tempfile(fileext = ".Rmd")
    # message("Writing Rmd content to: ", tmp_rmd)
    # Base RMarkdown content
    rmd_content <- c(
      "---",
      "title: \"Summary Report\"",
      "output: html_document",
      "---",
      "<span style='font-size:22px;'>This report summarizes the coefficient of variation (CV) and visualizes ridge density distributions across different
      transformation methods. Additionally, it includes a heatmap illustrating the recall rates based on cancer-type-specific gold-standard proteins,
      and the expression of the top five proteins identified by the min-score method. The aim is to evaluate the suitability of various mass
      spectrometry (MS) proteomics quantification formats for deconvolution in identifying cell-type-specific protein expression.</span>",
      "```{r setup, include=FALSE}",
      "library(ggplot2)",
      "library(ggridges)",
      "library(dplyr)",
      "library(knitr)",
      "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
      "```",
      "",
      "## CV Summary Table",
      "<span style='font-size:22px;'> This table presents a statistical summary of the coefficient of variation (CV) across various transformation methods.
      It includes the mean, median, range of CV values, and the proportion of CVs greater than 0.25 for each method, enabling a comparative assessment of
      variability and stability introduced by each method.</span>",
      "```{r}",
      "if (nrow(cv_summary) > 0) {",
      "  knitr::kable(cv_summary, digits = 3, row.names = FALSE)",
      "} else {",
      "  cat(\"No CV summary data available.\")",
      "}",
      "```",
      "",
      "## CV Ridge Plot",
      "",
      "<span style='font-size:22px; color:#444;'>",
      "This figure shows the density ridge plot of the CV for each transformation method, allowing for a visual comparison
      of CV distributions across methods.",
      "</span>",
      "",
      "```{r cv_ridge_plot, echo=FALSE, out.width='100%'}",
      if (!is.null(cv_ridge_plot_file) && file.exists(cv_ridge_plot_file)) {
        paste0("knitr::include_graphics('", cv_ridge_plot_file, "')")
      } else {
        "cat(\"No CV ridge plot generated.\")"
      },
      "```"
    )
    # Append Recall Heatmaps
    if (!is.null(recall_heatmap_files$EDec) && file.exists(recall_heatmap_files$EDec)) {
      rmd_content <- c(rmd_content,
                       "\n## Recall Heatmap \n",
                       "<span style='font-size:22px;'>This plot presents the recall rates of gold-standard proteins of each cell type.</span>\n",
                       "```{r EDec_heatmap, echo=FALSE, out.width='100%'}\n",
                       paste0("knitr::include_graphics('", recall_heatmap_files$EDec, "')\n"),
                       "```\n"
      )
    }
    if (!is.null(recall_heatmap_files$rodeo) && file.exists(recall_heatmap_files$rodeo)) {
      rmd_content <- c(rmd_content,
                       "\n##\n",
                       "```{r rodeo_heatmap, echo=FALSE, out.width='100%'}\n",
                       paste0("knitr::include_graphics('", recall_heatmap_files$rodeo, "')\n"),
                       "```\n"
      )
    }
    if (!is.null(recall_heatmap_files$csSAM) && file.exists(recall_heatmap_files$csSAM)) {
      rmd_content <- c(rmd_content,
                       "\n## \n",
                       "```{r csSAM_heatmap, echo=FALSE, out.width='100%'}\n",
                       paste0("knitr::include_graphics('", recall_heatmap_files$csSAM, "')\n"),
                       "```\n"
      )
    }
    if (!is.null(recall_heatmap_files$bMIND) && file.exists(recall_heatmap_files$bMIND)) {
      rmd_content <- c(rmd_content,
                       "\n## \n",
                       "```{r bMIND_heatmap, echo=FALSE, out.width='100%'}\n",
                       paste0("knitr::include_graphics('", recall_heatmap_files$bMIND, "')\n"),
                       "```\n"
      )
    }
    # Append Marker Barplots
    if (!is.null(marker_barplot_files$EDec) && file.exists(marker_barplot_files$EDec)) {
      rmd_content <- c(rmd_content,
                       "\n## Marker Barplot (Min_score Method)\n",
                       "<span style='font-size:22px;'>This plot shows the specificity of the top five proteins in each cell type, identified by the Min_score.</span>\n",
                       "```{r EDec_barplot, echo=FALSE, out.width='100%'}\n",
                       paste0("knitr::include_graphics('", marker_barplot_files$EDec, "')\n"),
                       "```\n"
      )
    }
    if (!is.null(marker_barplot_files$rodeo) && file.exists(marker_barplot_files$rodeo)) {
      rmd_content <- c(rmd_content,
                       "\n## \n",
                       "```{r rodeo_barplot, echo=FALSE, out.width='100%'}\n",
                       paste0("knitr::include_graphics('", marker_barplot_files$rodeo, "')\n"),
                       "```\n"
      )
    }
    if (!is.null(marker_barplot_files$csSAM) && file.exists(marker_barplot_files$csSAM)) {
      rmd_content <- c(rmd_content,
                       "\n## \n",
                       "```{r csSAM_barplot, echo=FALSE, out.width='100%'}\n",
                       paste0("knitr::include_graphics('", marker_barplot_files$csSAM, "')\n"),
                       "```\n"
      )
    }
    if (!is.null(marker_barplot_files$bMIND) && file.exists(marker_barplot_files$bMIND)) {
      rmd_content <- c(rmd_content,
                       "\n## \n",
                       "```{r bMIND_barplot, echo=FALSE, out.width='100%'}\n",
                       paste0("knitr::include_graphics('", marker_barplot_files$bMIND, "')\n"),
                       "```\n"
      )
    }
    writeLines(rmd_content, tmp_rmd)
    
    # Try rendering, catch any errors
    
    tryCatch({
      
      rmarkdown::render(tmp_rmd,
                        
                        output_file = basename(output_html),
                        
                        output_dir = dirname(output_html),
                        
                        quiet = TRUE)
      
      message("HTML Report Finished: ", output_html)
      
    }, error = function(e) {
      
      warning("Failed to render HTML report: ", conditionMessage(e))
      
    }, finally = {
      
      if (file.exists(tmp_rmd)) {
        
        unlink(tmp_rmd) # Clean up temporary Rmd file
        
        # message("Temporary Rmd file removed: ", tmp_rmd)
        
      }
      
      # Also clean up temporary plot files if they exist
      
      all_temp_plots <- c(cv_ridge_plot_file, unlist(recall_heatmap_files), unlist(marker_barplot_files))
      
      for (f in all_temp_plots) {
        
        if (!is.null(f) && file.exists(f)) {
          
          unlink(f)
          
        }
        
      }
      
    })
    
  } else {
    
    message("output_html is NULL or empty. Skipping HTML report generation.")
    
  }
  
  




  
  # Return results
  return(invisible(list(
    transformed_list = transformed_list,
    cv_summary = cv_summary,
    cv_plot_data = cv_df,
    EDec = list(EDec_results = EDec_results, EDec_marker_table = EDec_marker_table),
    rodeo = list(rodeo_results = rodeo_results, rodeo_marker_table = rodeo_marker_table),
    csSAM = list(csSAM_result = csSAM_results, csSAM_marker_table = csSAM_marker_table),
    bMIND = list(bMIND_result = bMIND_results, bMIND_marker_table = bMIND_marker_table),
    gene_cell_correlation = combined_mat,
    gold_standard_markers = gold_standard_markers
  )))
}