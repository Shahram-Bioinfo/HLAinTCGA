#' @title HLA Typing and Expression Data
#' @name data
#' @docType data
#' @description HLA allele frequency and expression data from TCGA projects.
#' @format A data frame with columns:
#' \describe{
#'   \item{A1_exp}{Expression level of the first HLA-A allele (TPM).}
#'   \item{A2_exp}{Expression level of the second HLA-A allele (TPM).}
#'   \item{B1_exp}{Expression level of the first HLA-B allele (TPM).}
#'   \item{B2_exp}{Expression level of the second HLA-B allele (TPM).}
#'   \item{C1_exp}{Expression level of the first HLA-C allele (TPM).}
#'   \item{C2_exp}{Expression level of the second HLA-C allele (TPM).}
#'   \item{A_allele1}{The first HLA-A allele (e.g., "A*02:01").}
#'   \item{A_allele2}{The second HLA-A allele (e.g., "A*32:01").}
#'   \item{B_allele1}{The first HLA-B allele (e.g., "B*03:01").}
#'   \item{B_allele2}{The second HLA-B allele (e.g., "B*15:01").}
#'   \item{C_allele1}{The first HLA-C allele (e.g., "C*07:01").}
#'   \item{C_allele2}{The second HLA-C allele (e.g., "C*07:02").}
#'   \item{AJCC}{Cancer stage classification based on AJCC.}
#'   \item{Ethnicity Category}{Ethnicity of the patient.}
#'   \item{Race Category}{Race of the patient.}
#'   \item{Sex}{Sex of the patient (e.g., "Male", "Female").}
#'   \item{patient_id}{Unique identifier for the patient.}
#'   \item{project_id}{The TCGA project identifier (e.g., "TCGA-GBM").}
#'   \item{sample_id}{Sample identifier associated with the patient.}
#' }
#' @source Extracted from TCGA supplementary files.
"data"
