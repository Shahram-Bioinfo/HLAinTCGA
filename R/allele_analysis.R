utils::globalVariables(c(
  "project_id", "Sex", "AJCC", "number_of_allele", 
  "number_of_sample", "allele_frequency", "Race Category",
  "Ethnicity Category", "genotype", "frequency", "total_count"
))

#' Allele Frequency and Expression Analysis Across TCGA Cancer Projects
#'
#' This function analyzes the frequency and expression levels of a specific HLA allele across TCGA cancer projects.
#' It generates a summary table and visualizations to explore allele distributions.
#'
#' @param allele A string specifying the HLA allele to analyze (e.g., "A*02:01", "B*51:07", "C*20:04").
#' @param data A data frame containing allele-specific columns, expression data, and TCGA project IDs.
#' @return A list with two elements:
#'   - `results_table`: A data frame with allele frequency and expression statistics across TCGA projects.
#'   - `plot`: A combined visualization showing allele frequency (bar plot) and expression levels (box plot).
#'
#' @details
#' This function performs the following steps:
#' 1. Validates the input allele format.
#' 2. Filters the dataset to identify rows containing the specified allele.
#' 3. Calculates allele frequencies and expression statistics across TCGA projects.
#' 4. Generates two visualizations:
#'    - A bar plot showing allele frequencies across projects.
#'    - A box plot visualizing expression levels for the specified allele.
#' 5. Combines the visualizations into a single layout.
#'
#' @examples
#' # Example usage:
#' result <- allele_freq("A*02:01", data = data)
#' print(result$results_table)
#' print(result$plot)
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom cowplot plot_grid
#' @importFrom stats var sd reorder
#' @export

allele_freq <- function(allele, data) {
  # Validate allele format
  if (!grepl("^[A-C]\\*\\d{2}:\\d{2}$", allele)) {
    stop("Invalid allele format. Use formats like A*02:01, B*51:07, or C*20:04.")
  }
  
  # Identify allele type and related columns
  allele_type <- substr(allele, 1, 1)  # Extract allele type (A, B, or C)
  allele_columns <- paste0(allele_type, c("_allele1", "_allele2"))
  expression_columns <- paste0(allele_type, c("1_exp", "2_exp"))
  
  # Check if the allele exists in the data
  allele_exists <- any(
    (!is.na(data[[allele_columns[1]]]) & data[[allele_columns[1]]] == allele) |
      (!is.na(data[[allele_columns[2]]]) & data[[allele_columns[2]]] == allele)
  )
  if (!allele_exists) {
    message(paste("The allele", allele, "is not found in the dataset."))
    return(NULL)
  }
  
  # Count total samples by project_id
  total_samples <- dplyr::summarise(
    dplyr::group_by(data, project_id),
    number_of_sample = dplyr::n()
  )
  
  # Filter data for the specified allele
  filtered_data <- dplyr::filter(
    data,
    (!is.na(.data[[allele_columns[1]]]) & .data[[allele_columns[1]]] == allele) |
      (!is.na(.data[[allele_columns[2]]]) & .data[[allele_columns[2]]] == allele)
  )
  
  # Check if filtered data is empty
  if (nrow(filtered_data) == 0) {
    message(paste("No data found for allele", allele, "in the dataset."))
    return(NULL)
  }
  
  # Extract expression levels for the allele
  filtered_data <- dplyr::mutate(
    filtered_data,
    expression = dplyr::case_when(
      .data[[allele_columns[1]]] == allele ~ .data[[expression_columns[1]]],
      .data[[allele_columns[2]]] == allele ~ .data[[expression_columns[2]]],
      TRUE ~ NA_real_
    )
  )
  
  # Calculate summary statistics for expression levels
  expression_summary <- dplyr::summarise(
    dplyr::group_by(filtered_data, project_id),
    mean_expression = mean(expression, na.rm = TRUE),
    min_expression = min(expression, na.rm = TRUE),
    max_expression = max(expression, na.rm = TRUE),
    variance_expression = var(expression, na.rm = TRUE),
    std_error_expression = sd(expression, na.rm = TRUE) / sqrt(dplyr::n()),
    number_of_allele = sum(!is.na(expression)),
    .groups = "drop"
  )
  
  # Merge total samples with expression summary and calculate allele frequency
  combined_results <- dplyr::mutate(
    dplyr::left_join(total_samples, expression_summary, by = "project_id"),
    allele_frequency = ifelse(is.na(number_of_allele), 0, number_of_allele) / (2 * number_of_sample)
  )
  
  # Generate frequency plot (horizontal bar chart)
  frequency_plot <- ggplot2::ggplot(combined_results, ggplot2::aes(
    x = allele_frequency, 
    y = reorder(project_id, allele_frequency), 
    fill = project_id
  )) +
    ggplot2::geom_bar(stat = "identity", alpha = 0.8, width = 0.6, color = "black") +
    ggplot2::scale_fill_manual(values = scales::hue_pal()(length(unique(combined_results$project_id)))) +
    ggplot2::labs(title = paste("Allele Frequency of", allele), x = "", y = "") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_blank(),
      legend.position = "none",  # Remove legend
      panel.grid.major = ggplot2::element_line(color = "grey80"),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  # Generate expression plot
  expression_plot <- ggplot2::ggplot(filtered_data, ggplot2::aes(
    x = project_id, 
    y = expression, 
    fill = project_id
  )) +
    ggplot2::geom_boxplot(outlier.colour = "red", outlier.shape = 16, notch = FALSE, alpha = 0.8) +
    ggplot2::scale_fill_manual(values = scales::hue_pal()(length(unique(filtered_data$project_id)))) +
    ggplot2::labs(
      title = paste("Allele-Specific Expression of", allele, "Across TCGA Cancer Projects"),
      x = "",
      y = "Expression Levels (TPM)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
      legend.position = "none",  # Remove legend
      panel.grid.major = ggplot2::element_line(color = "grey80"),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  # Combine the two plots using cowplot
  combined_plot <- cowplot::plot_grid(
    frequency_plot,
    expression_plot,
    ncol = 2,
    rel_widths = c(1, 2)
  )
  
  # Print the combined plot
  print(combined_plot)
  
  # Return results and the combined plot
  return(list(results_table = combined_results, plot = combined_plot))
}

#' Allele Frequency and Expression Analysis by Gender Across TCGA Cancer Projects
#'
#' This function analyzes the frequency and expression levels of a specific HLA allele across TCGA cancer projects,
#' stratified by gender (Male/Female). It generates a summary table and visualizations to explore allele distributions.
#'
#' @param allele A string specifying the HLA allele to analyze (e.g., "A*02:01", "B*51:07", "C*20:04").
#' @param data A data frame containing allele-specific columns, expression data, and TCGA project IDs.
#' @param project_filter (Optional) A string specifying a TCGA project ID to filter the analysis (e.g., "TCGA-BRCA").
#'   If NULL, the analysis is conducted across all projects. Defaults to NULL.
#' @param show_expression Logical; if TRUE, includes box plots of expression levels in the output. Defaults to FALSE.
#'
#' @return A list with two elements:
#'   - `results_table`: A data frame with allele frequency and expression statistics stratified by gender.
#'   - `plot`: A combined visualization showing allele frequency and expression levels.
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom cowplot plot_grid
#' @importFrom stats var sd reorder
#' @export

allele_freq_gender <- function(allele, data, project_filter = NULL, show_expression = FALSE) {
  # Validate allele format
  if (!grepl("^[A-C]\\*\\d{2}:\\d{2}$", allele)) {
    stop("Invalid allele format. Use formats like A*02:01, B*51:07, or C*20:04.")
  }
  
  # Identify allele type and related columns
  allele_type <- substr(allele, 1, 1)  # Extract the first character (A, B, or C)
  allele_columns <- paste0(allele_type, c("_allele1", "_allele2"))  # Columns for alleles
  expression_columns <- paste0(allele_type, c("1_exp", "2_exp"))  # Columns for expressions
  
  # Check if the allele exists in the dataset
  allele_exists <- any(
    (!is.na(data[[allele_columns[1]]]) & data[[allele_columns[1]]] == allele) |
      (!is.na(data[[allele_columns[2]]]) & data[[allele_columns[2]]] == allele)
  )
  if (!allele_exists) {
    message(paste("The allele", allele, "is not found in the dataset."))
    return(NULL)
  }
  
  # Filter data for the specified allele and gender categories
  filtered_data <- data %>%
    filter(
      (!is.na(.data[[allele_columns[1]]]) & .data[[allele_columns[1]]] == allele) |
        (!is.na(.data[[allele_columns[2]]]) & .data[[allele_columns[2]]] == allele),
      Sex %in% c("Male", "Female")
    ) %>%
    mutate(
      expression = case_when(
        .data[[allele_columns[1]]] == allele ~ .data[[expression_columns[1]]],
        .data[[allele_columns[2]]] == allele ~ .data[[expression_columns[2]]],
        TRUE ~ NA_real_
      )
    )
  
  # Count total samples grouped by project_id and gender categories
  total_samples <- data %>%
    filter(Sex %in% c("Male", "Female")) %>%
    group_by(project_id, Sex) %>%
    summarise(number_of_sample = n(), .groups = "drop")
  
  # Compute allele frequency and summary statistics for expression
  expression_summary <- filtered_data %>%
    group_by(project_id, Sex) %>%
    summarise(
      mean_expression = mean(expression, na.rm = TRUE),
      number_of_allele = sum(!is.na(expression)),
      .groups = "drop"
    )
  
  # Combine total samples with expression summary
  combined_results <- total_samples %>%
    left_join(expression_summary, by = c("project_id", "Sex")) %>%
    mutate(
      allele_frequency = ifelse(is.na(number_of_allele), 0, number_of_allele) / (2 * number_of_sample)
    )
  
  # Plot for all projects if no specific project is filtered
  if (is.null(project_filter)) {
    frequency_plot <- ggplot(combined_results, aes(x = project_id, y = allele_frequency, fill = Sex)) +
      geom_bar(stat = "identity", position = "stack", alpha = 0.8, colour = "black") +
      scale_fill_manual(values = c("Male" = "#1f77b4", "Female" = "#ff7f0e")) +
      labs(
        title = paste("Allele Frequency of", allele, "by Gender Across All TCGA Projects"),
        x = "",
        y = "Allele Frequency"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_blank()
      )
    
    print(frequency_plot)
    return(list(results_table = combined_results, plot = frequency_plot))
  }
  
  # Filter results for the specific project
  filtered_results <- combined_results %>%
    filter(project_id == project_filter)
  
  # Generate frequency plot for the filtered project
  frequency_plot <- ggplot(filtered_results, aes(x = Sex, y = allele_frequency, fill = Sex)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8, colour = "black") +
    scale_fill_manual(values = c("Male" = "#1f77b4", "Female" = "#ff7f0e")) +
    labs(
      title = paste("Allele Frequency of", allele, "in", project_filter, "by Gender"),
      x = "Gender",
      y = "Allele Frequency"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid.major = element_line(color = "grey80"),
      panel.grid.minor = element_blank()
    )
  
  # Generate expression plot for the filtered project
  expression_plot <- ggplot(filtered_data %>% filter(project_id == project_filter),
                            aes(x = Sex, y = expression, fill = Sex)) +
    geom_boxplot(outlier.colour = "red", notch = FALSE, alpha = 0.8) +
    scale_fill_manual(values = c("Male" = "#1f77b4", "Female" = "#ff7f0e")) +
    labs(
      title = paste("Expression Levels of", allele, "in", project_filter, "by Gender"),
      x = "Gender",
      y = "Expression Levels (TPM)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid.major = element_line(color = "grey80"),
      panel.grid.minor = element_blank()
    )
  
  # Combine the frequency and expression plots into a single visualization using cowplot
  combined_plot <- cowplot::plot_grid(frequency_plot, expression_plot, ncol = 2, rel_widths = c(1, 1.5))
  
  # Print the combined plot
  print(combined_plot)
  
  # Return results table and combined plot
  return(list(results_table = filtered_results, plot = combined_plot))
}

#' Allele Frequency and Expression Analysis by Race Across TCGA Cancer Projects
#'
#' This function analyzes the frequency and expression levels of a specific HLA allele across TCGA cancer projects,
#' stratified by race categories. It generates a summary table and visualizations to explore allele distributions.
#'
#' @param allele A string specifying the HLA allele to analyze (e.g., "A*02:01", "B*51:07", "C*20:04").
#' @param data A data frame containing allele-specific columns, expression data, and TCGA project IDs.
#' @param project_filter (Optional) A string specifying a TCGA project ID to filter the analysis (e.g., "TCGA-BRCA").
#'   If NULL, the analysis is conducted across all projects. Defaults to NULL.
#' @param show_expression Logical; if TRUE, includes box plots of expression levels in the output. Defaults to FALSE.
#'
#' @return A list with two elements:
#'   - `results_table`: A data frame with allele frequency and expression statistics stratified by race.
#'   - `plot`: A combined visualization showing allele frequency and expression levels.
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom cowplot plot_grid
#' @importFrom stats var sd reorder
#' @export

allele_freq_race <- function(allele, data, project_filter = NULL, show_expression = FALSE) {
  # Validate allele format
  if (!grepl("^[A-C]\\*\\d{2}:\\d{2}$", allele)) {
    stop("Invalid allele format. Use formats like A*02:01, B*51:07, or C*20:04.")
  }
  
  # Identify allele type and columns
  allele_type <- substr(allele, 1, 1)
  allele_columns <- paste0(allele_type, c("_allele1", "_allele2"))
  expression_columns <- paste0(allele_type, c("1_exp", "2_exp"))
  
  # Check if allele exists in the data
  allele_exists <- any(
    (!is.na(data[[allele_columns[1]]]) & data[[allele_columns[1]]] == allele) |
      (!is.na(data[[allele_columns[2]]]) & data[[allele_columns[2]]] == allele)
  )
  if (!allele_exists) {
    message(paste("The allele", allele, "is not found in the dataset."))
    return(NULL)
  }
  
  # Filter data for the allele and race categories
  filtered_data <- data %>%
    filter(
      (!is.na(.data[[allele_columns[1]]]) & .data[[allele_columns[1]]] == allele) |
        (!is.na(.data[[allele_columns[2]]]) & .data[[allele_columns[2]]] == allele),
      `Race Category` %in% c("White", "Asian", "Black or African American")
    ) %>%
    mutate(
      expression = case_when(
        .data[[allele_columns[1]]] == allele ~ .data[[expression_columns[1]]],
        .data[[allele_columns[2]]] == allele ~ .data[[expression_columns[2]]],
        TRUE ~ NA_real_
      )
    )
  
  # Count total samples by project_id and race categories
  total_samples <- data %>%
    filter(`Race Category` %in% c("White", "Asian", "Black or African American")) %>%
    group_by(project_id, `Race Category`) %>%
    summarise(number_of_sample = n(), .groups = "drop")
  
  # Compute allele expression summary
  expression_summary <- filtered_data %>%
    group_by(project_id, `Race Category`) %>%
    summarise(
      mean_expression = mean(expression, na.rm = TRUE),
      number_of_allele = sum(!is.na(expression)),
      .groups = "drop"
    )
  
  # Combine total samples with expression summary
  combined_results <- total_samples %>%
    left_join(expression_summary, by = c("project_id", "Race Category")) %>%
    mutate(
      allele_frequency = ifelse(is.na(number_of_allele), 0, number_of_allele) / (2 * number_of_sample)
    )
  
  # Plot for all projects if no specific project is provided
  if (is.null(project_filter)) {
    frequency_plot <- ggplot(combined_results, aes(x = project_id, y = allele_frequency, fill = `Race Category`)) +
      geom_bar(stat = "identity", position = "stack", alpha = 0.8, colour = "black") +
      scale_fill_manual(values = c("White" = "#ff6347", "Asian" = "#4682b4", "Black or African American" = "#32cd32")) +
      labs(
        title = paste("Allele Frequency of", allele, "by Race Across All TCGA Projects"),
        x = "",
        y = "Allele Frequency"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_blank()
      )
    
    print(frequency_plot)
    return(list(results_table = combined_results, plot = frequency_plot))
  }
  
  # Plot for a specific project
  filtered_results <- combined_results %>%
    filter(project_id == project_filter)
  
  # Generate frequency plot for the specific project
  frequency_plot <- ggplot(filtered_results, aes(x = `Race Category`, y = allele_frequency, fill = `Race Category`)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8, colour = "black") +
    scale_fill_manual(values = c("White" = "#ff6347", "Asian" = "#4682b4", "Black or African American" = "#32cd32")) +
    labs(
      title = paste("Allele Frequency of", allele, "in", project_filter, "by Race"),
      x = "Race Category",
      y = "Allele Frequency"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid.major = element_line(color = "grey80"),
      panel.grid.minor = element_blank()
    )
  
  # Generate expression plot for the specific project
  expression_plot <- ggplot(filtered_data %>% filter(project_id == project_filter),
                            aes(x = `Race Category`, y = expression, fill = `Race Category`)) +
    geom_boxplot(outlier.colour = "red", notch = FALSE, alpha = 0.8) +
    scale_fill_manual(values = c("White" = "#ff6347", "Asian" = "#4682b4", "Black or African American" = "#32cd32")) +
    labs(
      title = paste("Expression Levels of", allele, "in", project_filter, "by Race"),
      x = "Race Category",
      y = "Expression Levels (TPM)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid.major = element_line(color = "grey80"),
      panel.grid.minor = element_blank()
    )
  
  # Combine the frequency and expression plots using cowplot
  combined_plot <- cowplot::plot_grid(frequency_plot, expression_plot, ncol = 2, rel_widths = c(1, 1.5))
  
  # Print the combined plot
  print(combined_plot)
  
  # Return results table and combined plot
  return(list(results_table = filtered_results, plot = combined_plot))
}

#' Allele Frequency and Expression Analysis by Ethnicity Across TCGA Cancer Projects
#'
#' This function analyzes the frequency and expression levels of a specific HLA allele across TCGA cancer projects,
#' stratified by ethnicity categories. It generates a summary table and visualizations to explore allele distributions.
#'
#' @param allele A string specifying the HLA allele to analyze (e.g., "A*02:01", "B*51:07", "C*20:04").
#' @param data A data frame containing allele-specific columns, expression data, and TCGA project IDs.
#' @param project_filter (Optional) A string specifying a TCGA project ID to filter the analysis (e.g., "TCGA-BRCA").
#'   If NULL, the analysis is conducted across all projects. Defaults to NULL.
#' @param show_expression Logical; if TRUE, includes box plots of expression levels in the output. Defaults to FALSE.
#'
#' @return A list with two elements:
#'   - `results_table`: A data frame with allele frequency and expression statistics stratified by ethnicity.
#'   - `plot`: A combined visualization showing allele frequency and expression levels.
#'
#' @importFrom ggplot2 ggplot geom_bar aes labs scale_fill_manual theme_minimal element_text element_blank
#' @importFrom dplyr filter group_by summarise mutate left_join case_when
#' @importFrom cowplot plot_grid
#' @import ggplot2
#' @import dplyr
#' @importFrom cowplot plot_grid
#' @importFrom stats var sd reorder
#' @export

allele_freq_ethnicity <- function(allele, data, project_filter = NULL, show_expression = FALSE) {
  
  # Validate allele format
  if (!grepl("^[A-C]\\*\\d{2}:\\d{2}$", allele)) {
    stop("Invalid allele format. Use formats like A*02:01, B*51:07, or C*20:04.")
  }
  
  # Identify allele type and related columns
  allele_type <- substr(allele, 1, 1)  # Extract the first character (A, B, or C)
  allele_columns <- paste0(allele_type, c("_allele1", "_allele2"))  # Columns for alleles
  expression_columns <- paste0(allele_type, c("1_exp", "2_exp"))  # Columns for expressions
  
  # Check if the allele exists in the dataset
  allele_exists <- any(
    (!is.na(data[[allele_columns[1]]]) & data[[allele_columns[1]]] == allele) |
      (!is.na(data[[allele_columns[2]]]) & data[[allele_columns[2]]] == allele)
  )
  if (!allele_exists) {
    message(paste("The allele", allele, "is not found in the dataset."))
    return(NULL)
  }
  
  # Filter data for the specified allele and ethnicity categories
  filtered_data <- data %>%
    dplyr::filter(
      (!is.na(.data[[allele_columns[1]]]) & .data[[allele_columns[1]]] == allele) |
        (!is.na(.data[[allele_columns[2]]]) & .data[[allele_columns[2]]] == allele),
      `Ethnicity Category` %in% c("Hispanic Or Latino", "Not Hispanic Or Latino")
    ) %>%
    dplyr::mutate(
      expression = dplyr::case_when(
        .data[[allele_columns[1]]] == allele ~ .data[[expression_columns[1]]],
        .data[[allele_columns[2]]] == allele ~ .data[[expression_columns[2]]],
        TRUE ~ NA_real_
      )
    )
  
  # Count total samples grouped by project_id and ethnicity categories
  total_samples <- data %>%
    dplyr::filter(`Ethnicity Category` %in% c("Hispanic Or Latino", "Not Hispanic Or Latino")) %>%
    dplyr::group_by(project_id, `Ethnicity Category`) %>%
    dplyr::summarise(number_of_sample = dplyr::n(), .groups = "drop")
  
  # Compute allele frequency and summary statistics for expression
  expression_summary <- filtered_data %>%
    dplyr::group_by(project_id, `Ethnicity Category`) %>%
    dplyr::summarise(
      mean_expression = mean(expression, na.rm = TRUE),
      number_of_allele = sum(!is.na(expression)),
      .groups = "drop"
    )
  
  # Combine total samples with expression summary
  combined_results <- total_samples %>%
    dplyr::left_join(expression_summary, by = c("project_id", "Ethnicity Category")) %>%
    dplyr::mutate(
      allele_frequency = ifelse(is.na(number_of_allele), 0, number_of_allele) / (2 * number_of_sample)
    )
  
  # Plot for all projects if no specific project is filtered
  if (is.null(project_filter)) {
    frequency_plot <- ggplot2::ggplot(combined_results, ggplot2::aes(x = project_id, y = allele_frequency, fill = `Ethnicity Category`)) +
      ggplot2::geom_bar(stat = "identity", position = "stack", alpha = 0.8, colour = "black") +
      ggplot2::scale_fill_manual(values = c("Hispanic Or Latino" = "#ff6347", "Not Hispanic Or Latino" = "#4682b4")) +
      ggplot2::labs(
        title = paste("Allele Frequency of", allele, "by Ethnicity Across All TCGA Projects"),
        x = "",
        y = "Allele Frequency"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        legend.position = "top",
        panel.grid.major = ggplot2::element_line(color = "grey80"),
        panel.grid.minor = ggplot2::element_blank()
      )
    
    print(frequency_plot)
    return(list(results_table = combined_results, plot = frequency_plot))
  }
  
  # Filter results for the specific project
  filtered_results <- combined_results %>%
    dplyr::filter(project_id == project_filter)
  
  # Generate frequency plot for the filtered project
  frequency_plot <- ggplot2::ggplot(filtered_results, ggplot2::aes(x = `Ethnicity Category`, y = allele_frequency, fill = `Ethnicity Category`)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge", alpha = 0.8, colour = "black") +
    ggplot2::scale_fill_manual(values = c("Hispanic Or Latino" = "#ff6347", "Not Hispanic Or Latino" = "#4682b4")) +
    ggplot2::labs(
      title = paste("Allele Frequency of", allele, "in", project_filter, "by Ethnicity"),
      x = "Ethnicity Category",
      y = "Allele Frequency"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major = ggplot2::element_line(color = "grey80"),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  # Generate expression plot for the filtered project
  expression_plot <- ggplot2::ggplot(filtered_data %>% dplyr::filter(project_id == project_filter),
                                     ggplot2::aes(x = `Ethnicity Category`, y = expression, fill = `Ethnicity Category`)) +
    ggplot2::geom_boxplot(outlier.colour = "red", notch = FALSE, alpha = 0.8) +
    ggplot2::scale_fill_manual(values = c("Hispanic Or Latino" = "#ff6347", "Not Hispanic Or Latino" = "#4682b4")) +
    ggplot2::labs(
      title = paste("Expression Levels of", allele, "in", project_filter, "by Ethnicity"),
      x = "Ethnicity Category",
      y = "Expression Levels (TPM)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major = ggplot2::element_line(color = "grey80"),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  # Combine the frequency and expression plots into a single visualization using cowplot
  combined_plot <- cowplot::plot_grid(frequency_plot, expression_plot, ncol = 2, rel_widths = c(1, 1.5))
  
  # Print the combined plot
  print(combined_plot)
  
  # Return results table and combined plot
  return(list(results_table = filtered_results, plot = combined_plot))
}



#' Allele Frequency and Expression Analysis by Disease Stage Across TCGA Cancer Projects
#'
#' This function analyzes the frequency and expression levels of a specific HLA allele across TCGA cancer projects,
#' stratified by disease stages. It generates a summary table and visualizations to explore allele distributions.
#'
#' @param allele A string specifying the HLA allele to analyze (e.g., "A*02:01", "B*51:07", "C*20:04").
#' @param data A data frame containing allele-specific columns, expression data, TCGA project IDs, and disease stages (AJCC column).
#' @param project_filter (Optional) A string specifying a TCGA project ID to filter the analysis (e.g., "TCGA-BRCA").
#'   If NULL, the analysis is conducted across all projects. Defaults to NULL.
#' @param show_expression Logical; if TRUE, includes box plots of expression levels in the output. Defaults to FALSE.
#' @param bar_colors (Optional) A vector of colors to use for the bar plot fills. Defaults to using the rainbow palette.
#'
#' @return A list with two elements:
#'   - `results_table`: A data frame with allele frequency and expression statistics stratified by disease stage.
#'   - `plot`: A combined visualization showing allele frequency and expression levels.
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom cowplot plot_grid
#' @importFrom stats var sd reorder
#' @export

allele_freq_stage <- function(allele, data, project_filter = NULL, show_expression = FALSE, bar_colors = NULL) {
  # Validate allele format
  if (!grepl("^[A-C]\\*\\d{2}:\\d{2}$", allele)) {
    stop("Invalid allele format. Use formats like A*02:01, B*51:07, or C*20:04.")
  }
  
  # Define the order of disease stages
  stage_order <- c("STAGE I", "STAGE IA", "STAGE IB",
                   "STAGE II", "STAGE IIA", "STAGE IIB",
                   "STAGE III", "STAGE IIIA", "STAGE IIIB", "STAGE IIIC",
                   "STAGE IV", "STAGE X")
  
  # Normalize stage names to uppercase and apply stage ordering
  data <- data %>%
    mutate(AJCC = factor(toupper(AJCC), levels = stage_order)) %>%
    filter(AJCC != "NA")
  
  # Identify allele type and related columns
  allele_type <- substr(allele, 1, 1)
  allele_columns <- paste0(allele_type, c("_allele1", "_allele2"))
  expression_columns <- paste0(allele_type, c("1_exp", "2_exp"))
  
  # Check if allele exists in the data
  allele_exists <- any(
    (!is.na(data[[allele_columns[1]]]) & data[[allele_columns[1]]] == allele) |
      (!is.na(data[[allele_columns[2]]]) & data[[allele_columns[2]]] == allele)
  )
  if (!allele_exists) {
    message(paste("The allele", allele, "is not found in the dataset."))
    return(NULL)
  }
  
  # Filter data for the allele
  filtered_data <- data %>%
    filter(
      (!is.na(.data[[allele_columns[1]]]) & .data[[allele_columns[1]]] == allele) |
        (!is.na(.data[[allele_columns[2]]]) & .data[[allele_columns[2]]] == allele)
    ) %>%
    mutate(
      expression = case_when(
        .data[[allele_columns[1]]] == allele ~ .data[[expression_columns[1]]],
        .data[[allele_columns[2]]] == allele ~ .data[[expression_columns[2]]],
        TRUE ~ NA_real_
      )
    )
  
  # Remove NA values from filtered data for expression plot
  filtered_data <- filtered_data %>% filter(!is.na(expression))
  
  # Count total samples grouped by project and disease stages
  total_samples <- data %>%
    filter(!is.na(AJCC)) %>%
    group_by(project_id, AJCC) %>%
    summarise(number_of_sample = n(), .groups = "drop")
  
  # Compute allele frequency and expression summary
  expression_summary <- filtered_data %>%
    group_by(project_id, AJCC) %>%
    summarise(
      mean_expression = mean(expression, na.rm = TRUE),
      number_of_allele = sum(!is.na(expression)),
      .groups = "drop"
    )
  
  # Combine total samples with expression summary
  combined_results <- total_samples %>%
    left_join(expression_summary, by = c("project_id", "AJCC")) %>%
    mutate(
      allele_frequency = ifelse(is.na(number_of_allele), 0, number_of_allele) / (2 * number_of_sample)
    )
  
  # Set bar colors
  if (is.null(bar_colors)) {
    bar_colors <- rainbow(length(stage_order))
  }
  
  # Plot for all projects if no project filter is specified
  if (is.null(project_filter)) {
    frequency_plot <- ggplot(combined_results, aes(x = project_id, y = allele_frequency, fill = AJCC)) +
      geom_bar(stat = "identity", position = "stack", alpha = 0.8, colour = "black") +
      scale_fill_manual(values = bar_colors) +
      labs(
        title = paste("Allele Frequency of", allele, "Across All TCGA Projects"),
        x = "Project ID",
        y = "Allele Frequency"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_blank()
      )
    
    print(frequency_plot)
    return(list(results_table = combined_results, plot = frequency_plot))
  }
  
  # Filter results for the specific project
  filtered_results <- combined_results %>%
    filter(project_id == project_filter)
  
  # Generate frequency plot for the specific project
  frequency_plot <- ggplot(filtered_results, aes(x = AJCC, y = allele_frequency, fill = AJCC)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8, colour = "black") +
    scale_fill_manual(values = bar_colors) +
    labs(
      title = paste("Allele Frequency of", allele, "in", project_filter, "by Disease Stage"),
      x = "Disease Stage",
      y = "Allele Frequency"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.grid.major = element_line(color = "grey80"),
      panel.grid.minor = element_blank()
    )
  
  # Generate expression plot for the specific project
  expression_plot <- ggplot(filtered_data %>% filter(project_id == project_filter),
                            aes(x = AJCC, y = expression, fill = AJCC)) +
    geom_boxplot(outlier.colour = "red", notch = FALSE, alpha = 0.8) +
    scale_fill_manual(values = bar_colors) +
    labs(
      title = paste("Expression Levels of", allele, "in", project_filter, "by Disease Stage"),
      x = "Disease Stage",
      y = "Expression Levels (TPM)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.grid.major = element_line(color = "grey80"),
      panel.grid.minor = element_blank()
    )
  
  # Combine frequency and expression plots using cowplot
  combined_plot <- cowplot::plot_grid(frequency_plot, expression_plot, ncol = 2, rel_widths = c(1, 1.5))
  
  # Print the combined plot
  print(combined_plot)
  
  # Return results table and combined plot
  return(list(results_table = filtered_results, plot = combined_plot))
}
