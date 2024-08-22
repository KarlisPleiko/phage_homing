# Load necessary libraries
required_packages <- c("Rcpp", "readxl", "edgeR", "tidyverse", "ggplot2", "dplyr", "UpSetR", "ggpubr", "purrr")
lapply(required_packages, require, character.only = TRUE)

#Import datasets. File paths have to be changed to where the files are located. 
#Sample data files can be found here under Supplementary Data: https://academic.oup.com/nar/article/49/7/e38/6097687#235655415
#
biopanning_r2 <- read_excel("Downloads/gkaa1279_supplemental_files/S3-Kaarel Kurm - round2_allpeptides.xlsx", skip = 5)
biopanning_r2_2 <- read_excel("Downloads/gkaa1279_supplemental_files/S4-Kaarel Kurm - 20191114_Maarja_allpeptides.xlsx", skip = 5)

# Combine datasets
biopanning_R2 <- full_join(biopanning_r2, biopanning_r2_2, by = "Peptides")
View(biopanning_R2)

# Replace NA values with 0 for further data processing
biopanning_R2[, -1] <- lapply(biopanning_R2[, -1], function(x) ifelse(is.na(x), 0, x))

# Create datasets for differential homing analysis
datasets <- list(
  br2_br_vs_li = biopanning_R2[, c(1,5,11,17,3,9,15)],
  br2_br_vs_ki = biopanning_R2[, c(1,6,12,18,3,9,15)],
  br2_br_vs_lu = biopanning_R2[, c(1,4,10,16,3,9,15)],
  br2_br_vs_mus = biopanning_R2[, c(1,7,13,19,3,9,15)]
)

# Create group
group_biopanning <- factor(c(1,1,1,2,2,2))

# Create model matrix
biopanning_design <- model.matrix(~group_biopanning)

# Define a function to process each dataset and update column names
process_dataset <- function(data, design, output_file_prefix) {
  # Create DGEList and filter low abundance peptides
  dge <- DGEList(counts = data[, 2:7], genes = data[, 1], group = group_biopanning)
  keep <- rowSums(dge$counts > 5) >= 2
  dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]
  
  # Normalize samples and estimate dispersions
  dge_filtered_norm <- calcNormFactors(dge_filtered)
  dge_de <- dge_filtered_norm %>%
    estimateGLMCommonDisp(design) %>%
    estimateGLMTrendedDisp(design) %>%
    estimateGLMTagwiseDisp(design)
  
  # Fit model and perform statistical test
  dge_fit <- glmQLFit(dge_de, design)
  dge_test <- glmQLFTest(dge_fit)
  
  # Extract results
  res <- topTags(dge_test, sort.by = "logFC", n = 200000, p.value = 0.05)$table
  
  # Apply logFC threshold for filtered results
  res_lfc2 <- res %>% filter(logFC > 2)
  
  # Update column names to reflect comparison
  colnames(res)[2:4] <- paste0(colnames(res)[2:4], "_", output_file_prefix)
  colnames(res_lfc2)[2:4] <- paste0(colnames(res_lfc2)[2:4], "_", output_file_prefix)
  
  # Save both filtered and unfiltered results
  write.csv(res, file = paste0(output_file_prefix, "_full_results.csv"))
  write.csv(res_lfc2, file = paste0(output_file_prefix, "_lfc2_results.csv"))
  
  # Return results for viewing or further processing
  return(list(
    full_results = res,
    lfc2_results = res_lfc2
  ))
}

# Apply the processing function to each dataset
results_list <- lapply(names(datasets), function(name) {
  data <- datasets[[name]]
  process_dataset(data, biopanning_design, name)
})

# Full join all full results into one data frame by "Peptides"
combined_full_results <- reduce(lapply(results_list, `[[`, "full_results"), full_join, by = "Peptides")
write.csv(combined_full_results, "combined_full_results.csv")

# Full join all lfc2 results into one data frame by "Peptides"
combined_lfc2_results <- reduce(lapply(results_list, `[[`, "lfc2_results"), full_join, by = "Peptides")
write.csv(combined_lfc2_results, "combined_lfc2_results.csv")

# Optionally, view the combined results
View(combined_full_results)
View(combined_lfc2_results)

# Create a binary presence/absence matrix for UpSetR visualization
presence_absence_matrix <- combined_lfc2_results %>%
  mutate(across(starts_with("logFC"), ~ ifelse(!is.na(.), 1, 0))) %>%
  select(Peptides, starts_with("logFC")) %>%
  column_to_rownames("Peptides")

# Rename the columns for the UpSet plot
colnames(presence_absence_matrix) <- c(
  "brain_vs_liver",
  "brain_vs_kidney",
  "brain_vs_lung",
  "brain_vs_muscle"
)

upset(
  presence_absence_matrix, 
  sets = colnames(presence_absence_matrix), 
  keep.order = TRUE, 
  order.by = "freq", 
  main.bar.color = "black", 
  matrix.color = "red", 
  set_size.show = TRUE  # Hide the single set bars
)

