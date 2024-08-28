# Load necessary libraries
library(readr)
library(janitor)
library(dplyr)
library(readxl)
library(edgeR)
library(ggplot2)
library(ggrepel)

#import the dataset (change the file location to where you have saved it)
play_off_cpm_filt <- read_csv("~/Library/CloudStorage/play_off_all_cpm_filt-book-chapter-input.csv")
play_off_cpm_filt <- play_off_cpm_filt[,-1]

# Calculate amplification speed correction
play_off_cpm_filt$amp_speed_ratio <- ((play_off_cpm_filt$amplified_input_pool_2_s41_r1_001 + play_off_cpm_filt$amplified_input_pool_3_s42_r1_001 + play_off_cpm_filt$amplified_input_pool_1_s56_r1_001) / 3) / 
  ((play_off_cpm_filt$play_off_starting_pool_1_s43_r1_001 + play_off_cpm_filt$play_off_starting_pool_2_s44_r1_001 + play_off_cpm_filt$play_off_starting_pool_3_s45_r1_001) / 3)

# Apply amplification speed correction
play_off_cpm_filt_amp_speed_corr <- play_off_cpm_filt %>%
  mutate(across(3:17, ~ .x / amp_speed_ratio)) %>%
  mutate(across(22:56, ~ .x / amp_speed_ratio))

# Recalculate CPM after amplification correction
columns_to_sum <- c(1:56)
total_counts_2 <- colSums(play_off_cpm_filt_amp_speed_corr[,columns_to_sum])
play_off_cpm_filt_cpm <- play_off_cpm_filt_amp_speed_corr
play_off_cpm_filt_cpm[,columns_to_sum] <- sweep(play_off_cpm_filt_amp_speed_corr[,columns_to_sum], 2, total_counts_2, FUN = "/") * 1e6
View(play_off_cpm_filt_cpm)

# Sample-specific analysis: ctx-both-sx-ctrl-vs-blast
pept_group <- factor(c(1,1,1,1,1,2,2,2,2,2))
play_off_cpm_filt_amp_speed_corr_both_ctx <- play_off_cpm_filt_cpm[,c(58,5,10,15, 39, 44,24,29,34,49,54)]
View(play_off_cpm_filt_amp_speed_corr_both_ctx)
play_off_cpm_filt_amp_speed_corr_ctx_both_dge <- DGEList(counts = play_off_cpm_filt_amp_speed_corr_both_ctx[,-1], genes = play_off_cpm_filt_amp_speed_corr_both_ctx[,1], group = pept_group)

# Normalize factors and estimate dispersion
play_off_cpm_filt_amp_speed_corr_ctx_both_dge <- calcNormFactors(play_off_cpm_filt_amp_speed_corr_ctx_both_dge)
pept_design <- model.matrix(~pept_group)
play_off_cpm_filt_amp_speed_corr_ctx_both_dge_de <- estimateGLMCommonDisp(play_off_cpm_filt_amp_speed_corr_ctx_both_dge, pept_design)
play_off_cpm_filt_amp_speed_corr_ctx_both_dge_de <- estimateGLMTagwiseDisp(play_off_cpm_filt_amp_speed_corr_ctx_both_dge_de, pept_design)

# Fit GLM model and perform differential expression analysis
play_off_cpm_filt_amp_speed_corr_ctx_both_fit <- glmQLFit(play_off_cpm_filt_amp_speed_corr_ctx_both_dge_de, pept_design)
play_off_cpm_filt_amp_speed_corr_ctx_both_comp <-glmQLFTest(play_off_cpm_filt_amp_speed_corr_ctx_both_fit)
play_off_cpm_filt_amp_speed_corr_ctx_both_res <- topTags(play_off_cpm_filt_amp_speed_corr_ctx_both_comp, sort.by = "logFC", n = Inf)$table
View(play_off_cpm_filt_amp_speed_corr_ctx_both_res)
# Label and annotate results
play_off_cpm_filt_amp_speed_corr_ctx_both_res$Differential_homing <- "NO"
play_off_cpm_filt_amp_speed_corr_ctx_both_res$Differential_homing[play_off_cpm_filt_amp_speed_corr_ctx_both_res$logFC > 1] <- "logFC > 1 in blast"
play_off_cpm_filt_amp_speed_corr_ctx_both_res$Differential_homing[play_off_cpm_filt_amp_speed_corr_ctx_both_res$logFC < -1] <- "logFC < -1 in ctrl"
play_off_cpm_filt_amp_speed_corr_ctx_both_res$delabel <- ifelse(play_off_cpm_filt_amp_speed_corr_ctx_both_res$Differential_homing != "NO", play_off_cpm_filt_amp_speed_corr_ctx_both_res$brand_name, NA)

# Plot results
ggplot(data = play_off_cpm_filt_amp_speed_corr_ctx_both_res, aes(x = logFC, y = logCPM, col = Differential_homing, label = delabel)) +
  geom_text_repel(max.overlaps = 10) +
  geom_point() +
  theme_minimal() +
  geom_vline(xintercept = c(1,-1), col = "black") +
  geom_hline(yintercept = 11, col = "black") +
  scale_y_continuous() +
  xlab("log2 Fold Change") +
  ylim(0,20) +
  xlim(-2.5,2.5)

# Save results to CSV
write.csv(play_off_cpm_filt_amp_speed_corr_ctx_both_res, "Downloads/play_off_cpm_filt_amp_speed_corr_ctx_both_res.csv")
