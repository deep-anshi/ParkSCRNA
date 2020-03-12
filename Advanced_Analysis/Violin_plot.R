#################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Plotting violin plots from RNA assay for specific genes to interpret it for single cell rna sequencing data. 


# Lab: Havrada Lab, Dartmouth College
# Author: Deepanshi Shokeen

#################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Loading modules
library(Seurat)
library(cowplot)
library(jackstraw)
library(lfa)
library(readr)
library(rtracklayer)
library(SCINA)
library(preprocessCore)

# setting up directories
dir1 <- "/dartfs-hpc/rc/lab/H/HavrdaM/scrna-seq/01_pre_processing/files/"
dir_out_files <- "/dartfs-hpc/rc/lab/H/HavrdaM/scrna-seq/02_cell-type-classification/files/"
dir_out_figures <- "/dartfs-hpc/rc/lab/H/HavrdaM/scrna-seq/deepanshi/trial/"

# read cell labels from .csv file 
labels <- read.csv(paste0(dir_out_files, 'SCINA_cell_type_labels.csv'))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# visualize cell type classifications
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# read in Seurat dataset
pd.combined <- readRDS(paste0(dir1, "integrated_seurat_SCTransform_clustered-snn-0.4.rds"))
pd.combined$cellname <- sapply(rownames(pd.combined@meta.data), function(x) strsplit(x, "raw_data/")[[1]][2])
names(Idents(pd.combined)) <- sapply(names(Idents(pd.combined)), function(x) strsplit(x, "raw_data/")[[1]][2])

# check cells are in same order in Seurat dataset and SCINA labels before relabelling cell types 
all(pd.combined@meta.data$cellname==labels[,2])


# rename identities of cells with cell type classifications 
Idents(pd.combined) <- labels[,1]

# visualize w/ unknown cell types 
p <- DimPlot(pd.combined, reduction = "umap", label = TRUE)
save_plot(paste0(dir_out_figures, "all-groups_Seurat-v3_split-by-dataset_SCINA-cell-types_w-unknown.png"), 
          p, base_height = 7, base_width = 16)



######
seurat_integrated.2 <- subset(seurat_integrated, cells = names(Idents(seurat_integrated)[-unknown_ind]))


# getr indicies of cells of unknown cell type 
unknown_ind <- which(Idents(pd.combined)=="unknown")

# subset dataset to remove cells of unknwon cell type 
pd.combined.2 <- subset(pd.combined, cells = names(Idents(pd.combined)[-unknown_ind]))

# visualize w/o unknown cell types 
p <- DimPlot(pd.combined.2, reduction = "umap", split.by = "status", label = TRUE)
save_plot(paste0(dir_out_figures, "all-groups_Seurat-v3_split-by-dataset_SCINA-cell-types_wo-unknown.png"), 
          p, base_height = 7, base_width = 16)

# visualize LRRK2
p <- FeaturePlot(pd.combined.2, features = c("LRRK2"), min.cutoff = "q1", split.by = "status")
#p <- FeaturePlot(seurat_integrated, features = c("LRRK2"))
#p <- FeaturePlot(seurat_integrated, features = c("LRRK2"), split.by = "sample")
save_plot(paste0(dir1, "LRRK2.png"), p, base_height = 7, base_width = 15)

#seurat_integrated[["RNA"]]@assays["LRRK2"]



# Vln plot 
DefaultAssay(pd.combined) <- "RNA"
#DefaultAssay(seurat_integrated.2) <- "RNA"
#plots <- VlnPlot(pd.combined, features = c("LRRK2"), 
#                split.by = "sample", group.by = "celltype", 
#               pt.size = 0)
#save_plot(paste0(dir_out_figures, "Seurat-v3_SCINA-cell-types_LRRK2-Violin.png"), plots, base_height = 7, base_width = 15)

#LRRK2
plots <- VlnPlot(pd.combined, features = c("LRRK2"), 
                split.by = "sample", 
                 pt.size = 0)
save_plot(paste0(dir_out_figures, "Seurat-v3_SCINA-cell-types_LRRK2-Violin.png"), plots, base_height = 7, base_width = 15)

#IL18
plots <- VlnPlot(pd.combined, features = c("IL18"), 
                split.by = "sample", 
                 pt.size = 0)
save_plot(paste0(dir_out_figures, "Seurat-v3_SCINA-cell-types_IL18-Violin.png"), plots, base_height = 7, base_width = 15)

#NEK7
plots <- VlnPlot(pd.combined, features = c("NEK7"), 
                split.by = "sample", 
                 pt.size = 0)
save_plot(paste0(dir_out_figures, "Seurat-v3_SCINA-cell-types_NEK7-Violin.png"), plots, base_height = 7, base_width = 15)

#NLRP3
plots <- VlnPlot(pd.combined, features = c("NLRP3"), 
                split.by = "sample", 
                 pt.size = 0)
save_plot(paste0(dir_out_figures, "Seurat-v3_SCINA-cell-types_NLRP3-Violin.png"), plots, base_height = 7, base_width = 15)

#PYCARD
plots <- VlnPlot(pd.combined, features = c("PYCARD"), 
                split.by = "sample", 
                 pt.size = 0)
save_plot(paste0(dir_out_figures, "Seurat-v3_SCINA-cell-types_PYCARD-Violin.png"), plots, base_height = 7, base_width = 15)

#CASP1
plots <- VlnPlot(pd.combined, features = c("CASP1"), 
                split.by = "sample", 
                 pt.size = 0)
save_plot(paste0(dir_out_figures, "Seurat-v3_SCINA-cell-types_CASP1-Violin.png"), plots, base_height = 7, base_width = 15)

#GSDMD 
plots <- VlnPlot(pd.combined, features = c("GSDMD"), 
                split.by = "sample", 
                 pt.size = 0)
save_plot(paste0(dir_out_figures, "Seurat-v3_SCINA-cell-types_GSDMD-Violin.png"), plots, base_height = 7, base_width = 15)

#IL1B
plots <- VlnPlot(pd.combined, features = c("IL1B"), 
                split.by = "sample", 
                 pt.size = 0)
save_plot(paste0(dir_out_figures, "Seurat-v3_SCINA-cell-types_LRRK2-Violin.png"), plots, base_height = 7, base_width = 15)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# with dots
plots <- VlnPlot(pd.combined, features = c("NEK7"), 
                split.by = "sample", 
                 pt.size = 0.25)
save_plot(paste0(dir_out_figures, "Seurat-v3_SCINA-cell-types_NEK7-Violin.png"), plots, base_height = 7, base_width = 15)

# nlrp3
plots <- VlnPlot(pd.combined, features = c("NLRP3"), 
                split.by = "sample", 
                 pt.size = 0.25)
save_plot(paste0(dir_out_figures, "Seurat-v3_SCINA-cell-types_NLRP3-Violin.png"), plots, base_height = 7, base_width = 15)


# IL1B
plots <- VlnPlot(pd.combined, features = c("IL1B"), 
                split.by = "sample", 
                 pt.size = 0.25)
save_plot(paste0(dir_out_figures, "Seurat-v3_SCINA-cell-types_IL1B-Violin.png"), plots, base_height = 7, base_width = 15)


#CASP1
plots <- VlnPlot(pd.combined, features = c("CASP1"), 
                split.by = "sample", 
                 pt.size = 0.25)
save_plot(paste0(dir_out_figures, "Seurat-v3_SCINA-cell-types_CASP1-Violin.png"), plots, base_height = 7, base_width = 15)


#IL18
plots <- VlnPlot(pd.combined, features = c("IL18"), 
                split.by = "sample", 
                 pt.size = 0.25)
save_plot(paste0(dir_out_figures, "Seurat-v3_SCINA-cell-types_IL18-Violin.png"), plots, base_height = 7, base_width = 15)


#GSDMD
plots <- VlnPlot(pd.combined, features = c("GSDMD"), 
                split.by = "sample", 
                 pt.size = 0.25)
save_plot(paste0(dir_out_figures, "Seurat-v3_SCINA-cell-types_GSDMD-Violin.png"), plots, base_height = 7, base_width = 15)

