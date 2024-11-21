#' Genotype Frequency Analysis Across TCGA Cancer Projects
#'
#' This function analyzes the frequency of genotypes associated with a specific HLA allele across TCGA cancer projects.
#' It can compute genotype frequencies for all projects or a specific project, and optionally focuses on a specific genotype
#' (when a second allele is provided). The output includes a table of genotype frequencies and a customizable plot.
#'
#' @import dplyr
#' @import ggplot2
#' @import scales
#' @import grDevices
#' @import randomcoloR
#'
#' @param allele1 A string specifying the primary HLA allele to analyze (e.g., "A*02:01", "B*51:07", "C*20:04").
#' @param data A data frame containing HLA allele data with the following structure:
#'   - Columns for HLA alleles (e.g., "A_allele1", "A_allele2", "B_allele1", "B_allele2", "C_allele1", "C_allele2").
#'   - A project_id column indicating TCGA project identifiers.
#'   The dataset named `data` is included in the package as an example dataset.
#' @param allele2 (Optional) A string specifying a second allele to analyze a specific genotype
#'   (e.g., "A*11:01", "B*27:05", "C*07:02"). Defaults to NULL.
#' @param project_filter (Optional) A string specifying a TCGA project ID to filter results (e.g., "TCGA-BRCA").
#'   If provided, a vertical bar plot of genotype frequencies for the specified project is generated. Defaults to NULL.
#'
#' @return A list containing:
#'   - results_table: A data frame of genotype frequencies, including:
#'       - project_id: TCGA project identifiers.
#'       - genotype: The specific genotype (e.g., "A*02:01 - A*11:01").
#'       - count: The number of samples with the given genotype in the project.
#'       - total_count: The total number of samples in the project.
#'       - frequency: The frequency of the genotype within the project.
#'   - plot: A ggplot2 object visualizing the genotype frequencies.
#'     - For all projects: A stacked bar plot of frequencies by genotype across projects.
#'     - For a specific project: A vertical bar plot of frequencies for each genotype in the project.
#'
#' @examples
#' # Analyze all genotypes associated with "A*11:01" across projects
#' result <- genotype_freq("A*11:01", data)
#' print(result$results_table)
#'
#' # Analyze the specific genotype "A*11:01 - A*02:01" across projects
#' result_specific <- genotype_freq("A*11:01", data, allele2 = "A*02:01")
#' print(result_specific$results_table)
#'
#' # Analyze all genotypes associated with "A*11:01" in a specific project
#' result_project <- genotype_freq("A*11:01", data, project_filter = "TCGA-LIHC")
#' print(result_project$results_table)
#'
#' @export
#' 
genotype_freq <- function(allele1, data, allele2 = NULL, project_filter = NULL) {
  # Validate allele format
  if (!grepl("^[A-C]\\*\\d{2}:\\d{2}$", allele1)) {
    stop("Invalid allele1 format. Use formats like A*02:01, B*51:07, or C*20:04.")
  }
  if (!is.null(allele2) && !grepl("^[A-C]\\*\\d{2}:\\d{2}$", allele2)) {
    stop("Invalid allele2 format. Use formats like A*02:01, B*51:07, or C*20:04.")
  }
  
  # Identify allele type and related columns
  allele_type <- substr(allele1, 1, 1)
  allele_columns <- paste0(allele_type, c("_allele1", "_allele2"))
  
  # Filter data for allele1
  filtered_data <- data %>%
    dplyr::filter(
      (!is.na(.data[[allele_columns[1]]]) & .data[[allele_columns[1]]] == allele1) |
        (!is.na(.data[[allele_columns[2]]]) & .data[[allele_columns[2]]] == allele1)
    ) %>%
    dplyr::mutate(
      genotype = dplyr::if_else(
        .data[[allele_columns[1]]] < .data[[allele_columns[2]]],
        paste(.data[[allele_columns[1]]], .data[[allele_columns[2]]], sep = " - "),
        paste(.data[[allele_columns[2]]], .data[[allele_columns[1]]], sep = " - ")
      )
    )
  
  # Specific genotype (if allele2 is provided)
  if (!is.null(allele2)) {
    filtered_data <- filtered_data %>%
      dplyr::filter(
        (.data[[allele_columns[1]]] == allele1 & .data[[allele_columns[2]]] == allele2) |
          (.data[[allele_columns[1]]] == allele2 & .data[[allele_columns[2]]] == allele1)
      )
    if (nrow(filtered_data) == 0) {
      message(paste("The genotype", paste(allele1, allele2, sep = " - "), "is not found in the dataset."))
      return(NULL)
    }
  }
  
  # Count genotypes by project_id
  total_samples <- data %>%
    dplyr::group_by(project_id) %>%
    dplyr::summarise(total_count = dplyr::n(), .groups = "drop")
  
  genotype_counts <- filtered_data %>%
    dplyr::group_by(project_id, genotype) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::left_join(total_samples, by = "project_id") %>%
    dplyr::mutate(frequency = count / total_count)
  
  # Generate a color palette with distinct colors
  unique_genotypes <- unique(genotype_counts$genotype)
  num_colors <- length(unique_genotypes)
  set.seed(123) # Ensure reproducibility
  color_palette <- randomcoloR::distinctColorPalette(num_colors) # Generate distinct colors
  
  # Generate plot
  if (!is.null(project_filter)) {
    genotype_counts <- genotype_counts %>%
      dplyr::filter(project_id == project_filter)
    plot <- ggplot2::ggplot(genotype_counts, ggplot2::aes(x = genotype, y = frequency)) +
      ggplot2::geom_bar(stat = "identity", fill = "blue", alpha = 0.9, color = "black") +
      ggplot2::labs(title = paste("Genotype Frequency in", project_filter),
                    x = "Genotype", y = "Frequency") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  } else {
    plot <- ggplot2::ggplot(genotype_counts, ggplot2::aes(x = project_id, y = frequency, fill = genotype)) +
      ggplot2::geom_bar(stat = "identity", position = "stack", color = "black") +
      ggplot2::scale_fill_manual(values = color_palette) +
      ggplot2::labs(title = paste("Frequencies of Genotypes Associated with",allele1, "Across TCGA Projects"),
                    x = "", y = "Frequency") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }
  
  print(plot)
  return(list(results_table = genotype_counts, plot = plot))
}
