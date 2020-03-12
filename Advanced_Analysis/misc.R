#################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Lab: Havrada 
# Author: Owen Wilkins, Deepanshi Shokeen

#################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


library(Seurat)
library(cowplot)
library(jackstraw)
library(lfa)
library(readr)
library(rtracklayer)
library(ggplot2)
library(tidyverse)
library(hdf5r)

dir_raw <- "/Users/deepanshishokeen/Desktop/scrna/scrna-seq/01_pre_processing/raw_data/"

read_10x_multi <- function(dir_in, group){
  #dir_in <- control[6]
  x <- Read10X(data.dir = paste0(dir_in, "/filtered_feature_bc_matrix/"))
  x <- CreateSeuratObject(counts = x, project = "PD_PBMC", min.cells = 10, min.features = 100)
  x[["sample"]] <- paste0(group)
  x
}



# set names for files to read in 
control <- c("3401HC", "3402HC", "MH_2", "PD125", "PD133", "PD135")
pd <- c("MH_1", "MH_3", "PD123", "PD126", "PD134")
# add directories to these 
control <- paste0(dir_raw, control)
pd <- paste0(dir_raw, pd)

# read in files 
bc_mat_cont <- lapply(control, read_10x_multi, "control")
bc_mat_pd <- lapply(pd, read_10x_multi, "pd")

# Create a merged Seurat object
merged_seurat <- merge(x = bc_mat_cont[[1]], 
                       y = c(bc_mat_cont[[2]], bc_mat_cont[[3]], bc_mat_cont[[4]], bc_mat_cont[[5]], bc_mat_cont[[6]],
                             bc_mat_pd[[1]], bc_mat_pd[[2]], bc_mat_pd[[3]], bc_mat_pd[[4]], bc_mat_pd[[5]]), 
                       add.cell.id = c(control, pd))

# add mito qc column 
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

# add log10 gene/UMI (complexity measure of dataset)
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

# quick visualization of basic QC stats 
p <- VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"), ncol = 3)
save_plot(paste0(dir2, "merged-data_basic-qc-stats-boxplot.png"), p, base_height = 8.5, base_width = 15)

# Visualize the number of cell counts per cell
p <- merged_seurat@meta.data %>%
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
save_plot(paste0(dir2, "cell_counts.png"), p, base_height = 6, base_width = 4)

# Visualize the number UMIs/transcripts per cell
p <- merged_seurat@meta.data %>% 
  ggplot(aes(color=sample, x=nCount_RNA, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500) + ggtitle("UMI counts")
save_plot(paste0(dir2, "nUMI-per-cell.png"), p, base_height = 6, base_width = 4)

# Visualize the distribution of genes detected per cell via histogram
p <- merged_seurat@meta.data %>% 
  ggplot(aes(color=sample, x=nFeature_RNA, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300) + ggtitle("No. of genes detected /cell")
save_plot(paste0(dir2, "genes-detected-per-cell_hist.png"), p, base_height = 6, base_width = 4)

# Visualize the distribution of genes detected per cell via boxplot
p <- merged_seurat@meta.data %>% 
  ggplot(aes(x=sample, y=log10(nFeature_RNA), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")
save_plot(paste0(dir2, "genes-deteected-per-cell_boxplot.png"), p, base_height = 6, base_width = 4)

# Visualize the correlation between genes detected and number of UMIs and determine 
# whether strong presence of cells with low numbers of genes/UMIs
p <- merged_seurat@meta.data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) + 
  geom_point() + 
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)
save_plot(paste0(dir2, "genes-deteected-vs-nUMI_per-cell2.png"), p, base_height = 6, base_width = 9)

# Visualize the distribution of mitochondrial gene expression detected per cell
p <- merged_seurat@meta.data %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)
save_plot(paste0(dir2, "mito-ratio-per-cell.png"), p, base_height = 6, base_width = 4)

# Visualize the overall novelty of the gene expression by visualizing the genes detected per UMI
p <- merged_seurat@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
save_plot(paste0(dir2, "expression-novelty.png"), p, base_height = 6, base_width = 4)

# Filter out low quality reads using selected thresholds 
filtered_seurat <- subset(x = merged_seurat, 
                          subset= (nCount_RNA >= 500) & 
                            (nFeature_RNA >= 250) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))

# check how many cells were filtered
nrow(merged_seurat@meta.data)
nrow(filtered_seurat@meta.data)
saveRDS(filtered_seurat, file = paste0(dir1, "merged_seurat_filtered_counts.rds"))

# save raw counts 
exp_raw <- filtered_seurat@assays$RNA@counts
saveRDS(exp_raw, file = paste0(dir1, "merged_seurat_filtered_counts_raw.rds"))

####
bc_mat_cont <- lapply(control, read_10x_multi, "control")
bc_mat_cont <- lapply(pd, read_10x_multi, "pd")
