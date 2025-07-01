#' Plot Recall Heatmap
#'
#' Generates a heatmap visualizing recall rates for deconvolution markers.
#' Includes options to display true positives and total gold standard markers.
#'
#' @param recall_matrix A matrix of recall values (cell types x methods).
#' @param title_suffix Suffix for the plot title (e.g., "EDec").
#' @param tp_matrix Optional matrix of true positive counts.
#' @param total_matrix Optional matrix of total gold standard markers.
#' @return The file path to the saved PNG image of the heatmap, or NULL if no plot generated.
#' @export
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom grDevices png dev.off
#' @importFrom ggplot2 ggsave
plot_recall_heatmap <- function(recall_matrix, title_suffix, tp_matrix = NULL, total_matrix = NULL) {
  if (is.null(recall_matrix) || all(is.na(recall_matrix)) || nrow(recall_matrix) == 0 || ncol(recall_matrix) == 0) {
    message(paste0("No data to plot recall heatmap for ", title_suffix, "."))
    return(NULL)
  }
  
  df_long <- reshape2::melt(recall_matrix, varnames = c("CellType", "Method"), value.name = "Recall")
  
  if (!is.null(tp_matrix) && !is.null(total_matrix) &&
      all(dim(tp_matrix) == dim(recall_matrix)) && all(dim(total_matrix) == dim(recall_matrix))) {
    df_long$TP     <- reshape2::melt(tp_matrix)[, 3]
    df_long$Total  <- reshape2::melt(total_matrix)[, 3]
    df_long$Label  <- ifelse(
      is.na(df_long$Recall) | is.na(df_long$TP) | is.na(df_long$Total) | df_long$Total == 0,
      "-",
      sprintf("%.2f\n(%d/%d)", df_long$Recall, df_long$TP, df_long$Total)
    )
  } else {
    df_long$Label <- sprintf("%.2f", df_long$Recall)
  }
  
  # Ensure levels are in order for plotting consistency
  df_long$Method   <- factor(df_long$Method, levels = colnames(recall_matrix))
  df_long$CellType <- factor(df_long$CellType, levels = rev(rownames(recall_matrix)))
  
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Method, y = CellType, fill = Recall)) +
    ggplot2::geom_tile(color = "grey80", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = Label), size = 2.8, fontface = "bold", lineheight = 0.95) +
    ggplot2::scale_fill_gradient(low = "white", high = "skyblue", limits = c(0, 1), na.value = "grey90") +
    ggplot2::labs(
      title = paste(title_suffix, "- Recall Heatmap"),
      x = NULL,
      y = NULL,
      fill = "Recall"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold", size = 13),
      axis.text.x     = ggplot2::element_text(angle = 0, hjust = 0.5, face = "bold", color = "black", size = 8.5),
      axis.text.y     = ggplot2::element_text(size = 10, face = "bold", color = "black"),
      axis.ticks      = ggplot2::element_blank(),
      panel.grid      = ggplot2::element_blank(),
      legend.position = "right",
      legend.title    = ggplot2::element_text(face = "bold", size = 8),
      legend.text     = ggplot2::element_text(size = 8),
      legend.key.height = ggplot2::unit(0.6, "cm"),
      legend.key.width  = ggplot2::unit(0.4, "cm")
    )
  
  file_path <- tempfile(fileext = ".png")
  ggplot2::ggsave(filename = file_path, plot = p, width = 6, height = 1.7, dpi = 150)
  return(file_path)
}


#' Save Marker Barplot
#'
#' Generates and saves a bar plot showing the specificity of top marker genes
#' for each cell type, typically identified by the "Min_score" method.
#'
#' @param marker_tables A named list of data frames, where each data frame
#'   contains marker information for a deconvolution method. Expects a "Min_score" entry.
#' @param method_name The name of the deconvolution method (e.g., "EDec", "Rodeo").
#' @return The file path to the saved PNG image of the barplot, or NULL if no plot generated.
#' @export
#' @import ggplot2
#' @import dplyr
#' @importFrom reshape2 melt
#' @importFrom purrr map_dfr
#' @importFrom ggplot2 ggsave
save_marker_barplot <- function(marker_tables, method_name) {
  if (is.null(marker_tables) || length(marker_tables) == 0) return(NULL)
  # Ensure "Min_score" method exists in the provided marker_tables
  if (!"Min_score" %in% names(marker_tables)) {
    warning(paste("Min_score marker table not found for method:", method_name, "Skipping barplot."))
    return(NULL)
  }
  
  tbl <- marker_tables[["Min_score"]]
  if (is.null(tbl) || nrow(tbl) == 0 || !"Gene" %in% colnames(tbl)) {
    warning(paste("Min_score marker table is empty or invalid for method:", method_name, "Skipping barplot."))
    return(NULL)
  }
  
  spec_cols <- grep("_specificity$", colnames(tbl), value = TRUE)
  if (length(spec_cols) == 0) {
    warning(paste("No specificity columns found in Min_score marker table for method:", method_name, "Skipping barplot."))
    return(NULL)
  }
  
  df_long <- reshape2::melt(tbl[, c("Gene", spec_cols)], id.vars = "Gene")
  colnames(df_long) <- c("Gene", "CellType", "Specificity")
  df_long$CellType <- sub("_specificity", "", df_long$CellType)
  df_long$Specificity <- as.numeric(df_long$Specificity)
  df_long <- df_long[!is.na(df_long$Specificity) & df_long$Specificity >= 0 & df_long$Specificity <= 1, ]
  
  unique_celltypes <- unique(df_long$CellType)
  
  all_plot_data <- purrr::map_dfr(unique_celltypes, function(cell) {
    top_genes <- df_long %>%
      dplyr::filter(CellType == cell) %>%
      dplyr::arrange(desc(Specificity)) %>%
      dplyr::pull(Gene) %>%
      unique() %>%
      utils::head(5) # Use utils::head for clarity
    
    df_long %>%
      dplyr::filter(Gene %in% top_genes) %>%
      dplyr::mutate(MarkerSourceCellType = cell)
  })
  
  if (nrow(all_plot_data) == 0) {
    warning(paste("No marker genes to plot for method:", method_name, "after filtering for top 5."))
    return(NULL)
  }
  
  p <- ggplot2::ggplot(all_plot_data, ggplot2::aes(x = CellType, y = Specificity, fill = CellType)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::facet_wrap(~ MarkerSourceCellType + Gene, scales = "free_y", ncol = 5) + # Changed nrow to ncol for better layout with 5 genes per celltype
    
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::scale_fill_brewer(palette = "Accent") +
    
    ggplot2::labs(
      title = paste0(method_name, " - Cell Type Specificity (Min_score)"),
      y = "Specificity", x = "Cell Type"
    ) +
    
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 10, face = "bold", color = "black"),
      axis.text.y = ggplot2::element_text(size = 8, face = "bold", color = "black"),
      axis.title.y = ggplot2::element_text(size = 12),
      panel.spacing = ggplot2::unit(0.8, "lines"),
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "none",
      strip.text = ggplot2::element_text(size = 8, face = "bold") # Make facet labels slightly smaller
    )
  
  
  file_path <- tempfile(fileext = ".png")
  ggplot2::ggsave(file_path, p, width = 8, height = 6, dpi = 150) # Increased dpi for better quality
  return(file_path)
}