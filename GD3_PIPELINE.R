# ============================================================
# GD3 Plasma and CSF Proteomics Analysis Pipeline
# ============================================================
# Authors: Pavithra Krishnaswami 
# Manuscript: "Plasma and CSF Proteomics: Neuroinflammatory 
#              and Lysosomal Dysregulation in Gaucher Disease Type 3"
# Journal: Movement Disorders 
#
# Note: This script represents the consolidated analysis pipeline
# used to generate the results in this manuscript. Original 
# analyses were performed across multiple scripts and consolidated
# for clarity and reproducibility.
# ============================================================
# ==============================================================================

rm(list = ls())

library(readxl)
library(arrow)
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(ggplot2)
library(ggrepel)
library(ggsignif)
library(ggvenn)
library(patchwork)
library(gridExtra)
library(grid)
library(scales)
library(openxlsx)
library(OlinkAnalyze)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

##Merge TDU-TDR with Metadata to fill in Age, Gender and Race
TDU_plasma <- read_excel("Training/Training/GD3 data analysis v2/TDU-TDR compiled.xlsx")
Metadata_plasma <- read_excel("Training/Training/GD3 data analysis v2/Metadata_plasma_1.xlsx")
Metadata_plasma_merged <- merge(TDU_plasma, Metadata_plasma, by = "Screen ID")

##Read Plasma files
Metadata_plasma <- read_excel("Training/Training/GD3 data analysis v2/Metadata_plasma_1.xlsx")
OLINKdata_Plasma <- read_parquet("Training/Training/GD3 data analysis v2/GD3_Plasma_NPX.parquet")  

##Read CSF files
Metadata_CSF <- read_excel("Training/Training/GD3 data analysis v2/Metadata_CSF.xlsx")
OLINKdata_CSF <- read_parquet("Training/Training/GD3 data analysis_symposium 2025/AMETHIST__Venglustat_CSF_NPX_2024-03-13.parquet") 
output_dir <- "Training/Training/GD3 data analysis v5"

##count number of samples in Plasma and CSF 
##Plasma
sample_counts_plasma <- Metadata_plasma %>%
  group_by(sample_type) %>%
  dplyr::summarize(total_samples = n())
print(sample_counts_plasma)

##CSF
sample_counts_CSF <- Metadata_CSF %>%
  group_by(sample_type) %>%
  dplyr::summarize(total_samples = n())
print(sample_counts_CSF)
##Merge Metadata and parquet
Plasma_merged_data <- merge(Metadata_plasma, OLINKdata_Plasma, by.x = "Sample_ID", by.y = "SampleID")
head(Plasma_merged_data)
CSF_merged_data <- merge(Metadata_CSF, OLINKdata_CSF, by.x = "Sample_ID", by.y = "SampleID")
head(CSF_merged_data)

##Plasma Volcano plot p<0.05
# Rename column if needed
colnames(Plasma_merged_data)[colnames(Plasma_merged_data) == "Sample_ID"] <- "SampleID"
colnames(CSF_merged_data)[colnames(CSF_merged_data) == "Sample_ID"] <- "SampleID"
# Standardize the group label used for plotting/coloring everywhere below
recode_group <- function(df) {
  df %>% mutate(Group = case_when(
    sample_type == "healthy control" ~ "Control",
    sample_type == "disease baseline" ~ "GD3",
    TRUE ~ sample_type
  ), Group = factor(Group, levels = c("Control", "GD3")))
}
Plasma_merged_data <- recode_group(Plasma_merged_data)
CSF_merged_data <- recode_group(CSF_merged_data)

group_colors <- c("GD3" = "#E69F00", "Control" = "#3B7FB6")  # orange / blue, per Figure 1 legend text

# ==============================================================================
# 2. SHARED T-TEST HELPER
# ==============================================================================

run_ttest <- function(merged_data) {
  res <- olink_ttest(df = merged_data, variable = "sample_type", alternative = "two.sided")
  res <- na.omit(res)
  res$log2_fc <- res$estimate
  res$neg_log10_p <- -log10(res$p.value)
  res
}

ttest_plasma <- run_ttest(Plasma_merged_data)
ttest_csf <- run_ttest(CSF_merged_data)

# ==============================================================================
# FIGURE 1A/1B - VOLCANO PLOTS
# ==============================================================================

make_volcano <- function(ttest_results, title_text, highlight_genes = c("CD63", "GPNMB")) {
  ttest_results$significance <- case_when(
    ttest_results$p.value < 0.05 & abs(ttest_results$log2_fc) > 0.5 ~ "Significant & Large Effect",
    ttest_results$p.value < 0.05 ~ "Significant",
    abs(ttest_results$log2_fc) > 0.5 ~ "Large Effect Only",
    TRUE ~ "Not Significant"
  )

  top_by_pvalue <- ttest_results %>% arrange(p.value) %>% slice_head(n = 12)
  highlight_data <- ttest_results %>% filter(Assay %in% highlight_genes)
  genes_to_label <- bind_rows(top_by_pvalue, highlight_data) %>% distinct(Assay, .keep_all = TRUE)

  colors <- c("Significant & Large Effect" = "#d73027",
              "Significant" = "#fc8d59",
              "Large Effect Only" = "#91bfdb",
              "Not Significant" = "#878787")

  ggplot(ttest_results, aes(x = log2_fc, y = neg_log10_p)) +
    geom_point(aes(color = significance), alpha = 0.8, size = 2) +
    geom_text_repel(data = genes_to_label, aes(label = Assay), max.overlaps = Inf,
                     size = 6, fontface = "bold", box.padding = 0.5, point.padding = 0.3,
                     force = 2, min.segment.length = 0.1) +
    scale_color_manual(values = colors, name = "Classification") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.8, linewidth = 0.8) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black", alpha = 0.8, linewidth = 0.8) +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12),
          legend.position = "bottom",
          panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
          plot.margin = margin(20, 20, 20, 20)) +
    labs(title = title_text, x = "Log\u2082 Fold Change", y = "\u2212Log\u2081\u2080 (p-value)")
}

volcano_plasma <- make_volcano(ttest_plasma, "Volcano Plot: GD3 Plasma")
volcano_csf <- make_volcano(ttest_csf, "Volcano Plot: GD3 CSF")

ggsave(file.path(output_dir, "Plasma volcano plot-8x8.tif"), volcano_plasma,
       width = 8, height = 8, dpi = 300, compression = "lzw", bg = "white")
ggsave(file.path(output_dir, "CSF volcano plot-8x8.tif"), volcano_csf,
       width = 8, height = 8, dpi = 300, compression = "lzw", bg = "white")

# ==============================================================================
# FIGURE 1C/1D - TOP 10 SIGNIFICANT PROTEIN BOXPLOTS
# ==============================================================================

make_top10_boxplot <- function(merged_data, ttest_results, title_text) {
  top10 <- ttest_results %>%
    filter(p.value < 0.05) %>%
    arrange(p.value) %>%
    slice_head(n = 10) %>%
    mutate(sig_label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ "ns"
    ))

  plot_data <- merged_data %>%
    filter(Assay %in% top10$Assay) %>%
    mutate(Assay = factor(Assay, levels = top10$Assay))

  y_pos <- plot_data %>%
    group_by(Assay) %>%
    summarise(y = max(NPX, na.rm = TRUE) * 1.08, .groups = "drop")

  annot_data <- top10 %>%
    left_join(y_pos, by = "Assay") %>%
    mutate(x_start = 1, x_end = 2, x_text = 1.5, y_text = y * 1.03)

  ggplot(plot_data, aes(x = Group, y = NPX, fill = Group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
    facet_wrap(~ Assay, scales = "free_y", ncol = 5) +
    scale_fill_manual(values = group_colors) +
    geom_segment(data = annot_data, aes(x = x_start, xend = x_end, y = y, yend = y), inherit.aes = FALSE) +
    geom_text(data = annot_data, aes(x = x_text, y = y_text, label = sig_label), inherit.aes = FALSE, size = 5) +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          strip.text = element_text(face = "bold")) +
    labs(title = title_text, x = NULL, y = "NPX")
}

top10_plasma <- make_top10_boxplot(Plasma_merged_data, ttest_plasma, "Top 10 Significant Proteins - GD3 Plasma")
top10_csf <- make_top10_boxplot(CSF_merged_data, ttest_csf, "Top 10 Significant Proteins - GD3 CSF")

ggsave(file.path(output_dir, "boxplots_top10_gd3_plasma.tif"), top10_plasma,
       width = 14, height = 8, dpi = 300, compression = "lzw", bg = "white")
ggsave(file.path(output_dir, "boxplots_top10_gd3_csf.tif"), top10_csf,
       width = 14, height = 8, dpi = 300, compression = "lzw", bg = "white")

# ==============================================================================
# FIGURE 2A/2B - GSEA DOT PLOTS
# ==============================================================================

run_gsea <- function(ttest_results) {
  gene_list <- ttest_results$statistic
  names(gene_list) <- ttest_results$Assay
  gene_list <- sort(gene_list, decreasing = TRUE)
  gene_list <- gene_list[!duplicated(names(gene_list))]

  gseGO(geneList = gene_list, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
        pvalueCutoff = 0.05, pAdjustMethod = "BH", scoreType = "std", eps = 0)
}

make_gsea_plot <- function(gsea_results, title_text) {
  gsea_data <- as.data.frame(gsea_results@result)
  if (nrow(gsea_data) == 0) return(NULL)

  plot_data <- gsea_data %>%
    slice_head(n = 10) %>%
    mutate(Gene_Count = setSize,
           GeneRatio_calc = setSize / max(setSize),
           Description_wrapped = str_wrap(Description, width = 45),
           Description_factor = fct_reorder(Description_wrapped, GeneRatio_calc))

  ggplot(plot_data, aes(x = GeneRatio_calc, y = Description_factor)) +
    geom_point(aes(size = Gene_Count, color = p.adjust), alpha = 0.85) +
    scale_color_gradient(low = "#d73027", high = "#313695", name = "Adjusted\np-value") +
    scale_size_continuous(name = "Gene\nCount", range = c(4, 10)) +
    theme_classic(base_size = 14) +
    theme(axis.text.y = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          legend.position = "right") +
    labs(title = title_text, x = "Gene Ratio", y = "GO Biological Process")
}

gsea_plasma <- run_gsea(ttest_plasma)
gsea_csf <- run_gsea(ttest_csf)

gsea_plot_plasma <- make_gsea_plot(gsea_plasma, "GSEA Results: GD3 Plasma")
gsea_plot_csf <- make_gsea_plot(gsea_csf, "GSEA Results: GD3 CSF")

ggsave(file.path(output_dir, "gsea_plot_gd3_plasma.tif"), gsea_plot_plasma,
       width = 12, height = 10, dpi = 300, compression = "lzw", bg = "white")
ggsave(file.path(output_dir, "gsea_plot_gd3_csf.tif"), gsea_plot_csf,
       width = 12, height = 10, dpi = 300, compression = "lzw", bg = "white")

# ==============================================================================
# FIGURE 2C/2D - NEUROINFLAMMATORY PATHWAYS + TARGET MOLECULE TABLE
# ==============================================================================

neuro_keywords <- c("neuro", "glial", "microglial", "axon", "neuron", "brain",
                     "central nervous", "neuroinflam", "astrocyte", "oligodendrocyte",
                     "synapt", "dendrit", "myelin", "neural")
inflammatory_keywords <- c("inflamm", "cytokine", "interleukin", "tumor necrosis",
                            "chemokine", "interferon", "immune", "TNF", "IL-",
                            "activation", "response", "leukocyte", "macrophage")

score_pathway <- function(pathway) {
  pathway_lower <- tolower(pathway)
  score <- sum(sapply(neuro_keywords, function(k) grepl(k, pathway_lower))) * 3 +
    sum(sapply(inflammatory_keywords, function(k) grepl(k, pathway_lower))) * 2
  if (grepl("neuro.*inflam|glial.*activ|microglial|axon.*injur|response.*axon", pathway_lower)) score <- score + 5
  score
}

get_top_neuro_pathways <- function(gsea_results) {
  gsea_data <- as.data.frame(gsea_results@result)
  gsea_data$neuro_score <- sapply(gsea_data$Description, score_pathway)
  gsea_data %>% filter(neuro_score > 0) %>% arrange(desc(neuro_score), p.adjust) %>% slice_head(n = 10)
}

make_neuro_bar <- function(neuro_pathways, title_text) {
  plot_data <- neuro_pathways %>%
    mutate(neg_log_p = -log10(p.adjust),
           Pathway_short = ifelse(nchar(Description) > 50, paste0(substr(Description, 1, 47), "..."), Description),
           Pathway_factor = fct_reorder(Pathway_short, neg_log_p))

  ggplot(plot_data, aes(x = neg_log_p, y = Pathway_factor)) +
    geom_bar(stat = "identity", fill = "#4472C4", alpha = 0.8, width = 0.7) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "orange", linewidth = 1.2) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, max(plot_data$neg_log_p) * 1.15)) +
    theme_classic(base_size = 11) +
    theme(axis.text.y = element_text(size = 10, face = "bold"),
          plot.title = element_text(size = 13, face = "bold", hjust = 0)) +
    labs(title = title_text, x = "-log(adjusted p-value)", y = NULL)
}

make_target_table <- function(neuro_pathways, title_text) {
  gene_table <- data.frame(Pathway = character(), Target_Molecules = character())
  for (i in seq_len(min(5, nrow(neuro_pathways)))) {
    pw <- neuro_pathways$Description[i]
    pw_short <- ifelse(nchar(pw) > 45, paste0(substr(pw, 1, 42), "..."), pw)
    genes <- strsplit(as.character(neuro_pathways$core_enrichment[i]), "/")[[1]]
    gene_table <- rbind(gene_table, data.frame(Pathway = pw_short, Target_Molecules = paste(genes, collapse = ", ")))
  }
  tg <- tableGrob(gene_table, rows = NULL,
                   theme = ttheme_minimal(core = list(fg_params = list(hjust = 0, x = 0.02, fontsize = 9)),
                                          colhead = list(fg_params = list(fontface = "bold", fontsize = 10))))
  title_grob <- textGrob(title_text, gp = gpar(fontsize = 11, fontface = "bold"), just = "left", x = 0.02)
  combined <- gtable_add_rows(tg, heights = unit(0.5, "cm"), pos = 0)
  gtable_add_grob(combined, title_grob, t = 1, l = 1, r = ncol(combined))
}

neuro_plasma <- get_top_neuro_pathways(gsea_plasma)
neuro_csf <- get_top_neuro_pathways(gsea_csf)

fig2c <- (make_neuro_bar(neuro_plasma, "Neuroinflammatory Pathways - Plasma")) /
  wrap_elements(make_target_table(neuro_plasma, "Target Molecules - Plasma")) +
  plot_layout(heights = c(2, 1.5))
fig2d <- (make_neuro_bar(neuro_csf, "Neuroinflammatory Pathways - CSF")) /
  wrap_elements(make_target_table(neuro_csf, "Target Molecules - CSF")) +
  plot_layout(heights = c(2, 1.5))

ggsave(file.path(output_dir, "Neuro GSEA and table-Plasma.tif"), fig2c,
       width = 10, height = 10, dpi = 300, compression = "lzw", bg = "white")
ggsave(file.path(output_dir, "Neuro GSEA and table-CSF.tif"), fig2d,
       width = 10, height = 10, dpi = 300, compression = "lzw", bg = "white")

# ==============================================================================
# FIGURE 3A (left) - VENN DIAGRAM OF DE PROTEINS, p < 0.01
# (venn_diagram_common_proteins_fixed.tif was NOT reproducible from any of the
#  three source scripts - the only Venn code found there was for a different
#  comparison, overlap of proteins within specific shared pathways. This is a
#  new, from-scratch build of the 2-set DE-protein Venn the caption describes.)
# ==============================================================================

sig_plasma_01 <- ttest_plasma %>% filter(p.value < 0.01) %>% pull(Assay)
sig_csf_01 <- ttest_csf %>% filter(p.value < 0.01) %>% pull(Assay)
common_proteins <- intersect(sig_plasma_01, sig_csf_01)

venn_plot <- ggvenn(
  list(Plasma = sig_plasma_01, CSF = sig_csf_01),
  fill_color = c("#4472C4", "#ED7D31"),
  stroke_size = 0.5, text_size = 5, set_name_size = 6, show_percentage = FALSE
) + ggtitle("Differentially Expressed Proteins (p < 0.01)")

ggsave(file.path(output_dir, "venn_diagram_common_proteins_fixed.tif"), venn_plot,
       width = 8, height = 8, dpi = 300, compression = "lzw", bg = "white")

# ==============================================================================
# FIGURE 3A (right) - UP/DOWN TABLE FOR CO-DETECTED PROTEINS
# ==============================================================================

reg_table <- data.frame(Protein = common_proteins) %>%
  left_join(ttest_plasma %>% transmute(Protein = Assay, Plasma = ifelse(log2_fc > 0, "Up", "Down")), by = "Protein") %>%
  left_join(ttest_csf %>% transmute(Protein = Assay, CSF = ifelse(log2_fc > 0, "Up", "Down")), by = "Protein") %>%
  pivot_longer(cols = c(Plasma, CSF), names_to = "Compartment", values_to = "Regulation") %>%
  mutate(Compartment = factor(Compartment, levels = c("Plasma", "CSF")))

updown_plot <- ggplot(reg_table, aes(x = Compartment, y = fct_rev(Protein), fill = Regulation)) +
  geom_tile(color = "white", linewidth = 1) +
  scale_fill_manual(values = c("Up" = "#7B2D8E", "Down" = "#4CAF50"), na.value = "grey90") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(face = "bold")) +
  labs(title = "Co-detected Protein Regulation", x = NULL, y = NULL)

ggsave(file.path(output_dir, "common_proteins_updown_table.tif"), updown_plot,
       width = 6, height = max(4, length(common_proteins) * 0.3), dpi = 300, compression = "lzw", bg = "white")

# ==============================================================================
# FIGURE 3B - TREM2 AND CTSL BOXPLOTS, PLASMA & CSF
# ==============================================================================

make_protein_pair_boxplot <- function(merged_data, ttest_results, assays, title_text, ncol = NULL) {
  p_values <- ttest_results %>% filter(Assay %in% assays) %>%
    transmute(Assay, sig_label = case_when(
      p.value < 0.001 ~ "***", p.value < 0.01 ~ "**", p.value < 0.05 ~ "*", TRUE ~ "ns"
    ))

  plot_data <- merged_data %>% filter(Assay %in% assays) %>%
    mutate(Assay = factor(Assay, levels = assays))

  y_pos <- plot_data %>% group_by(Assay) %>% summarise(y = max(NPX, na.rm = TRUE) * 1.08, .groups = "drop")
  annot_data <- p_values %>% left_join(y_pos, by = "Assay") %>%
    mutate(x_start = 1, x_end = 2, x_text = 1.5, y_text = y * 1.03)

  ggplot(plot_data, aes(x = Group, y = NPX, fill = Group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.5) +
    facet_wrap(~ Assay, scales = "free_y", ncol = ncol) +
    scale_fill_manual(values = group_colors) +
    geom_segment(data = annot_data, aes(x = x_start, xend = x_end, y = y, yend = y), inherit.aes = FALSE) +
    geom_text(data = annot_data, aes(x = x_text, y = y_text, label = sig_label), inherit.aes = FALSE, size = 6) +
    theme_classic(base_size = 13) +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          strip.text = element_text(face = "bold", size = 13)) +
    labs(title = title_text, x = NULL, y = "NPX")
}

trem2_ctsl_plasma <- make_protein_pair_boxplot(Plasma_merged_data, ttest_plasma, c("TREM2", "CTSL"),
                                                "TREM2 and CTSL - Plasma, Control vs GD3")
trem2_ctsl_csf <- make_protein_pair_boxplot(CSF_merged_data, ttest_csf, c("TREM2", "CTSL"),
                                             "TREM2 and CTSL - CSF, Control vs GD3")

ggsave(file.path(output_dir, "TREM2 and CTSL-Plasma-Control vs GD3.tif"), trem2_ctsl_plasma,
       width = 8, height = 5, dpi = 300, compression = "lzw", bg = "white")
ggsave(file.path(output_dir, "TREM2 and CTSL-CSF-Control vs GD3.tif"), trem2_ctsl_csf,
       width = 8, height = 5, dpi = 300, compression = "lzw", bg = "white")

# ==============================================================================
# FIGURE 4B - PUBLISHED PD BIOMARKER BOXPLOTS, CSF
# ==============================================================================

pd_biomarkers_csf <- c("BECN1", "CTSD", "DDC", "GPNMB", "LAMP2", "MAP1LC3B2")

pd_biomarker_boxplot_csf <- make_protein_pair_boxplot(
  CSF_merged_data, ttest_csf, pd_biomarkers_csf,
  "Published PD Biomarkers - CSF, Control vs GD3", ncol = 3
)

ggsave(file.path(output_dir, "pd_biomarkers_boxplot_csf.tif"), pd_biomarker_boxplot_csf,
       width = 12, height = 8, dpi = 300, compression = "lzw", bg = "white")



cat("\nDone. Compare the regenerated files in '", output_dir, "' against the originals before submission.\n", sep = "")
