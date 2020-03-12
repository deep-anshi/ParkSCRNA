
#################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Perform PBMC cell type classification using prior knowledge-based approach, SCINA 
# (http://lce.biohpc.swmed.edu/scina/download.php#usage), using PBMC cell type marker 
# annotations from Ding et al, 2019, BioRxiv
# 
# Author: Owen Wilkins, Deepanshi Shokeen

#################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(Seurat)
library(cowplot)
library(jackstraw)
library(lfa)
library(readr)
library(rtracklayer)
library(SCINA)
library(preprocessCore)

# set up directories 
dir1 <- "/dartfs-hpc/rc/lab/H/HavrdaM/scrna-seq/scrna_sex_qc_4/files/"
dir_out_files <- "/dartfs-hpc/rc/lab/H/HavrdaM/scrna-seq/scrna_sex_qc_4/files/"
dir_out_figures <- "/dartfs-hpc/rc/lab/H/HavrdaM/scrna-seq/scrna_sex_qc_4/figures/"

 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read in each 10X xaset using Seurat functions 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# read in cell type annotations (signatures)
signatures=preprocess.signatures(paste0(dir_out_files, 'Ding-et-al_2019_PMBC_annotations_positive-markers-SCINA.csv'))

# check how many labels used more than once 
tab1 <- table(unlist(signatures))
tab1[tab1>1]

# read in expression matrix of raw counts 
exp_raw <- readRDS(file=paste0(dir1, "merged_seurat_filtered_counts_raw.rds"))
dim(exp_raw)

# Log scale and quantile normalization
exp_raw=log(exp_raw+1)
exp=normalize.quantiles(as.matrix(exp_raw))

# rename rows and columns of normalized data 
colnames(exp) <- colnames(exp_raw)
rownames(exp) <- rownames(exp_raw)

# run SCINA
results = SCINA(exp, signatures, 
                max_iter = 100, 
                convergence_n = 10, 
                convergence_rate = 0.999, 
                sensitivity_cutoff = 0.9, 
                rm_overlap=FALSE, 
                allow_unknown=TRUE,
                log_file='SCINA.log')

# save labels with cell IDs
cell_names <- sapply(colnames(exp), function(x) strsplit(x, "raw_data/")[[1]][2])
labels <- data.frame(results[[1]], cell_names)
labels$results..1.. <- as.character(labels$results..1..)
labels$cell_names <- as.character(labels$cell_names)
dim(labels)
#28087     2

# plot heatmap of cells by marker genes used for classification 
#ppi=300
#png(paste0(dir_out_figures, "plotheat_SCINA.png"), width=6*ppi, height=5*ppi, res=ppi)
#p <- plotheat.SCINA(exp, results, signatures)
#print(p)
#dev.off()

#save_plot(paste0(dir_out_figures, "SCINA-heatmap.png"), 
#          p, base_height = 7, base_width = 16)

# write cell labels to csv file 
write.csv(labels, file=paste0(dir_out_files, 'SCINA_cell_type_labels.csv'), row.names = FALSE)
labels <- read.csv(paste0(dir_out_files, 'SCINA_cell_type_labels.csv'))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# visualize cell type classifications
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# read in Seurat dataset
seurat_integrated <- readRDS(paste0(dir1, "integrated_seurat_SCTransform_clustered-snn-0.4.rds"))
#seurat_integrated <- readRDS(paste0(dir1, "integrated_seurat_SCTransform.rds"))
seurat_integrated$cellname <- sapply(rownames(seurat_integrated@meta.data), function(x) strsplit(x, "raw_data/")[[1]][2])
names(Idents(seurat_integrated)) <- sapply(names(Idents(seurat_integrated)), function(x) strsplit(x, "raw_data/")[[1]][2])
# check cells are in same order in Seurat dataset and SCINA labels before relabelling cell types 
all(seurat_integrated@meta.data$cellname==labels[,2])

# remove labels for cells filtered during QC (if they exist)
#table(labels[,2] %in% seurat_integrated@meta.data$cellname)
#ind1 <- which(!labels[,2] %in% rownames(seurat_integrated@meta.data))
#labels <- labels[-ind1,]
#all(names(Idents(seurat_integrated))==labels[,2])

# rename identities of cells with cell type classifications 
Idents(seurat_integrated) <- labels[,1]

# visualize w/ unknown cell types 
p <- DimPlot(seurat_integrated, reduction = "umap", label = TRUE)
p
save_plot(paste0(dir_out_figures, "all-groups_Seurat-v3_SCINA-cell-types_w-unknown.png"), p, base_height = 7, base_width = 8)

# getr indicies of cells of unknown cell type 
unknown_ind <- which(Idents(seurat_integrated)=="unknown")

# subset dataset to remove cells of unknwon cell type 
seurat_integrated.2 <- subset(seurat_integrated, cells = names(Idents(seurat_integrated)[-unknown_ind]))

# visualize w/o unknown cell types 
p <- DimPlot(seurat_integrated.2, reduction = "umap", label = TRUE)
p
save_plot(paste0(dir_out_figures, "all-groups_Seurat-v3_SCINA-cell-types_wo-unknown.png"), 
          p, base_height = 7, base_width = 8)


Idents(seurat_integrated.2) <- labels[,1]

# visualize w/o unknown cell types w/ sample types split
p <- DimPlot(seurat_integrated.2, reduction = "umap", label = TRUE)+ NoLegend()
p
save_plot(paste0(dir_out_figures, "all-groups_Seurat-v3_split-by-dataset_SCINA-cell-types_wo-unknown2.png"), 
          p, base_height = 7, base_width = 8)


seurat_integrated.2@meta.data$subject <- NA
seurat_integrated.2@meta.data$cellname2 <- sapply(seurat_integrated.2@meta.data$cellname, function(x) gsub("MH_1", "MH1", x))
seurat_integrated.2@meta.data$cellname2 <- sapply(seurat_integrated.2@meta.data$cellname2, function(x) gsub("MH_2", "MH2", x))
seurat_integrated.2@meta.data$cellname2 <- sapply(seurat_integrated.2@meta.data$cellname2, function(x) gsub("MH_3", "MH3", x))

seurat_integrated.2@meta.data$subject <- sapply(seurat_integrated.2@meta.data$cellname2, function(x) strsplit(x, "_")[[1]][1])

head(seurat_integrated@meta.data)
seurat_integrated.2@meta.data$cellname


table(seurat_integrated.2@meta.data$subject)
which(seurat_integrated.2@meta.data$subject=="MH")
seurat_integrated.2@meta.data[c(6198, 6199),]
seurat_integrated.2@meta.data$subject[6888]
seurat_integrated.2@meta.data$cellname2[c(6198, 6199)]

# rename identities of cells with cell type classifications 
Idents(seurat_integrated.2) <- seurat_integrated.2@meta.data$subject
p <- DimPlot(seurat_integrated.2, reduction = "umap", label = F)
p
save_plot(paste0(dir_out_figures, "integrated2.png"), 
          p, base_height = 7, base_width = 8)



p <- DimPlot(seurat_integrated.2, reduction = "umap", label = FALSE)

p
save_plot(paste0(dir_out_figures, "integrated.png"), 
          p, base_height = 7, base_width = 15)


seurat_integrated.2 <- RunUMAP(seurat_integrated.2, reduction = "pca", dims = 1:40)
# Plot the UMAP
p <- DimPlot(seurat_integrated.2,
             reduction = "umap",
             label = TRUE,
             label.size = 6,
             plot.title = "", split.by = "sample")
p
save_plot(paste0(dir_out_figures, "blah.png"), 
          p, base_height = 7, base_width = 8)



seurat_integrated.2.cont <- subset(seurat_integrated.2, cells = rownames(seurat_integrated.2@meta.data)[seurat_integrated.2$sample=="control"])
seurat_integrated.2.pd <- subset(seurat_integrated.2, cells = rownames(seurat_integrated.2@meta.data)[seurat_integrated.2$sample=="pd"])

  
mixmet_cont<-MixingMetric(seurat_integrated.2.cont, "subject")
mixmet_pd<-MixingMetric(seurat_integrated.2.pd, "subject")
mixingmets <- c(mixmet_cont, mixmet_pd)
mixingmets_df <- data.frame(mixingmets)
mixingmets_df$cell <- c(rep("Control", length(mixmet_cont)), rep("PD", length(mixmet_pd)))
mixingmets_df$col <- c(rep("#00AFBB", length(mixmet_cont)), rep("#E7B800", length(mixmet_pd)))

p<-ggviolin(mixingmets_df, "cell", "mixingmets", fill = "col",
         palette = c("#00AFBB", "#E7B800"),
         add = c("mean_sd"), add.params = list(fill = "white"))
save_plot(paste0(dir_out_figures, "mixingmets.png"), 
          p, base_height = 5.5, base_width = 4)


localmet <- LocalStruct(seurat_integrated.2, "subject")
plot()


##### highlight individual cell types ##### 

# write function to plot UMAP and highlight a given cell type
plot1 <- function(x){
  p <- DimPlot(seurat_integrated.2, reduction = "umap",
               split.by = "sample", label = TRUE,
               cells.highlight = rownames(seurat_integrated.2@meta.data)[which(Idents(seurat_integrated.2)==x)])
  save_plot(paste0(dir_out_figures, "Seurat-v3_SCINA_", x, ".png"), p, base_height = 7, base_width = 15)
}

# apply for each cell type 
sapply(names(table(Idents(seurat_integrated.2))), plot1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# evaluate Seurat-defined clusters vs cell types (are SCINA ID'd cell types defined by Seurat ID'Dd clusters?,
# if not, indicates we may need to change resolution/dimensionality for clustering, as some cell types are 
# being clustered w/ other cell types, and won'r be distingued in any down stream analyses using these dimensions)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


DefaultAssay(seurat_integrated) <- "RNA"
DefaultAssay(seurat_integrated.2) <- "RNA"
# Normalize RNA data for visualization purposes
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)
FeaturePlot(seurat_integrated.2, c("LRRK2"))

p <- FeaturePlot(seurat_integrated.2, features = c("LRRK2"))
p <- FeaturePlot(seurat_integrated.2, features = c("LRRK2"), split.by = "sample")
p
save_plot(paste0(dir_out_figures, "Seurat-v3_SCINA-cell-types_LRRK2-UMAP.png"), p, base_height = 7, base_width = 15)

seurat_integrated.2@meta.data$celltype <- Idents(seurat_integrated.2)


plots <- VlnPlot(seurat_integrated.2, features = c("LRRK2"), 
                 split.by = "sample", group.by = "celltype", 
                 pt.size = 0)
plots
save_plot(paste0(dir_out_figures, "Seurat-v3_SCINA-cell-types_LRRK2-Violin.png"), plots, base_height = 7, base_width = 15)



# Creating feature plots
p <- FeaturePlot(seurat_integrated, features = c("LRRK2"))
p <- FeaturePlot(seurat_integrated, features = c("LRRK2"), split.by = "sample")
p
save_plot(paste0(dir1, "LRRK2.png"), p, base_height = 7, base_width = 15)

seurat_integrated[["RNA"]]@assays["LRRK2"]
