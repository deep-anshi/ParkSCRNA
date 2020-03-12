########## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Dot Plot
# 
# Lab: Havrada Lab, Dartmouth College
# Author: Deepanshi Shokeen

########## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(Seurat)
library(cowplot)
library(jackstraw)
library(lfa)
library(readr)
library(rtracklayer)
library(SCINA)
library(preprocessCore)
library(stringr)
library(future)
library(MAST)
library(AnnotationHub)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
rm(list=ls())


# Setting up directories
dir1 <- "/dartfs-hpc/rc/lab/H/HavrdaM/scrna-seq/01_pre_processing/files/"
dir2 <- "/dartfs-hpc/rc/lab/H/HavrdaM/scrna-seq/02_cell-type-classification/files/"
dir_raw <- "/dartfs-hpc/rc/lab/H/HavrdaM/scrna-seq/01_pre_processing/raw_data/"
dir_out_files <- "/dartfs-hpc/rc/lab/H/HavrdaM/scrna-seq/02_cell-type-classification/files/"
dir_out_figures <- "/dartfs-hpc/rc/lab/H/HavrdaM/scrna-seq/03_differential_expression/figures/"

# change the below directory to save the output in different directory.
dir_out_figs <- "/dartfs-hpc/rc/lab/H/HavrdaM/scrna-seq/03_differential_expression/figures/"
 


labels <- read.csv(paste0(dir_out_files, 'SCINA_cell_type_labels.csv'))

# read in Seurat dataset
pd.combined <- readRDS(paste0(dir1, "integrated_seurat_SCTransform_clustered-snn-0.4.rds"))
pd.combined$cellname <- sapply(rownames(pd.combined@meta.data), function(x) strsplit(x, "raw_data/")[[1]][2])
names(Idents(pd.combined)) <- sapply(names(Idents(pd.combined)), function(x) strsplit(x, "raw_data/")[[1]][2])

# check cells are in same order in Seurat dataset and SCINA labels before relabelling cell types 
all(pd.combined@meta.data$cellname==labels[,2])


# rename identities of cells with cell type classifications 
Idents(pd.combined) <- labels[,1]
ind_unk <- which(pd.combined@meta.data$celltype=="unknown")
pd.combined_sub <- pd.combined[, !rownames(pd.combined@meta.data) %in% rownames(pd.combined@meta.data)[ind_unk]]

master_res_list <- readRDS(paste0(dir_out_files, "DEG_results_SCINA_all_cell-types-list.rds"))

DefaultAssay(pd.combined_sub) <- "SCT"

#table(pd.combined_sub@meta.data$celltype.sample)
pd.combined_sub@meta.data$celltype.sample[pd.combined_sub@meta.data$celltype.sample=="b_cell_control"] <- "B-cell (C)"
pd.combined_sub@meta.data$celltype.sample[pd.combined_sub@meta.data$celltype.sample=="b_cell_pd"] <- "B-cell (PD)"
pd.combined_sub@meta.data$celltype.sample[pd.combined_sub@meta.data$celltype.sample=="CD14_monocyte_control"] <- "CD14 monocyte (C)"
pd.combined_sub@meta.data$celltype.sample[pd.combined_sub@meta.data$celltype.sample=="CD14_monocyte_pd"] <- "CD14 monocyte (PD)"
pd.combined_sub@meta.data$celltype.sample[pd.combined_sub@meta.data$celltype.sample=="CD16_monocyte_control"] <- "CD16 monocyte (C)"
pd.combined_sub@meta.data$celltype.sample[pd.combined_sub@meta.data$celltype.sample=="CD16_monocyte_pd"] <- "CD16 monocyte (PD)"
pd.combined_sub@meta.data$celltype.sample[pd.combined_sub@meta.data$celltype.sample=="CD4_t_cell_control"] <- "CD4 T-cell (C)"
pd.combined_sub@meta.data$celltype.sample[pd.combined_sub@meta.data$celltype.sample=="CD4_t_cell_pd"] <- "CD4 T-cell (PD)"
pd.combined_sub@meta.data$celltype.sample[pd.combined_sub@meta.data$celltype.sample=="cytotoxic_t_cell_control"] <- "Cytotoxic T-cell (C)"
pd.combined_sub@meta.data$celltype.sample[pd.combined_sub@meta.data$celltype.sample=="cytotoxic_t_cell_pd"] <- "Cytotoxic T-cell (PD)"
pd.combined_sub@meta.data$celltype.sample[pd.combined_sub@meta.data$celltype.sample=="dendritic_cell_control"] <- "Dendritic cell (C)"
pd.combined_sub@meta.data$celltype.sample[pd.combined_sub@meta.data$celltype.sample=="dendritic_cell_pd"] <- "Dendritic cell (PD)"
pd.combined_sub@meta.data$celltype.sample[pd.combined_sub@meta.data$celltype.sample=="megakaryocyte_control"] <- "Megakaryocyte (C)"
pd.combined_sub@meta.data$celltype.sample[pd.combined_sub@meta.data$celltype.sample=="megakaryocyte_pd"] <- "Megakaryocyte (PD)"
pd.combined_sub@meta.data$celltype.sample[pd.combined_sub@meta.data$celltype.sample=="nk_cell_control"] <- "NK-cell (C)"
pd.combined_sub@meta.data$celltype.sample[pd.combined_sub@meta.data$celltype.sample=="nk_cell_pd"] <- "NK-cell (PD)"
pd.combined_sub@meta.data$celltype.sample[pd.combined_sub@meta.data$celltype.sample=="plasma_cell_control"] <- "Plasma cell (C)"
pd.combined_sub@meta.data$celltype.sample[pd.combined_sub@meta.data$celltype.sample=="plasma_cell_pd"] <- "Plasma cell (PD)"
pd.combined_sub@meta.data$celltype.sample[pd.combined_sub@meta.data$celltype.sample=="plasmacytoid_dendritic_cell_control"] <- "Plasmacytoid dendritic cell (C)"
pd.combined_sub@meta.data$celltype.sample[pd.combined_sub@meta.data$celltype.sample=="plasmacytoid_dendritic_cell_pd"] <- "Plasmacytoid dendritic cell (PD)"

pd.combined_sub@meta.data$celltype[pd.combined_sub@meta.data$celltype=="b_cell"] <- "B cell"
pd.combined_sub@meta.data$celltype[pd.combined_sub@meta.data$celltype=="CD14_monocyte"] <- "CD14 monocyte"
pd.combined_sub@meta.data$celltype[pd.combined_sub@meta.data$celltype=="CD16_monocyte"] <- "CD16 monocyte"
pd.combined_sub@meta.data$celltype[pd.combined_sub@meta.data$celltype=="CD4_t_cell"] <- "CD4 T-cell"
pd.combined_sub@meta.data$celltype[pd.combined_sub@meta.data$celltype=="cytotoxic_t_cell"] <- "Cytotoxic T-cell"
pd.combined_sub@meta.data$celltype[pd.combined_sub@meta.data$celltype=="dendritic_cell"] <- "Dendritic cell"
pd.combined_sub@meta.data$celltype[pd.combined_sub@meta.data$celltype=="megakaryocyte"] <- "Megakaryocyte"
pd.combined_sub@meta.data$celltype[pd.combined_sub@meta.data$celltype=="nk_cell"] <- "NK-cell"
pd.combined_sub@meta.data$celltype[pd.combined_sub@meta.data$celltype=="plasma_cell"] <- "Plasma cell"
pd.combined_sub@meta.data$celltype[pd.combined_sub@meta.data$celltype=="plasmacytoid_dendritic_cell"] <- "Plasmacytoid dendritic cell"

pd.combined_sub@meta.data$sample[pd.combined_sub@meta.data$sample=="pd"] <- "(PD)"
pd.combined_sub@meta.data$sample[pd.combined_sub@meta.data$sample=="control"] <- "(Control)"



# Specifying gene names to be included in Dot Plot
markers.to.plot <- c("ASAH1", "CASP1", "CTSB", "GSDMD", "HLA-DQA1", "LRRK2",
                      "PYCARD", "ATP13A2", "BST1", "CTSD", "IL18", "SCARB2", "SNCA", "VPS35")

 
p <-DotPlot(pd.combined_sub, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8,
        split.by = "sample") + RotatedAxis() 
save_plot(paste0(dir_out_figs, "DotPlot_ParkGenes_1.png"), p, base_height = 7, base_width = 15)


p <-DotPlot(pd.combined, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8,
            split.by = "sample") + RotatedAxis() 
save_plot(paste0(dir_out_figs, "DotPlot_ParkGenes_2.png"), p, base_height = 7, base_width = 15)

