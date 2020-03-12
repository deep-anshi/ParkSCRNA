
########## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Heatmaps for Differentially expressed genes 
# Lab: Havrada lab, Dartmouth College
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
library(ComplexHeatmap)
library(RColorBrewer)
library(stats)
library(viridis)
library(colorspace)
library(circlize)

dir1 <- "/dartfs-hpc/rc/lab/H/HavrdaM/scrna-seq/01_pre_processing/files/"
dir2 <- "/dartfs-hpc/rc/lab/H/HavrdaM/scrna-seq/02_cell-type-classification/files/"

data_dir <- "/dartfs-hpc/rc/lab/H/HavrdaM/scrna-seq/03_differential_expression/"

# read in data sets 
mono_14 <- read.csv(paste0(data_dir, "CD14_monocyte.interferon.response.csv"), stringsAsFactors = F)
mono_16 <- read.csv(paste0(data_dir, "CD16_monocyte.interferon.response.csv"), stringsAsFactors = F)
deg<- read.csv(paste0(data_dir, "interferon.response.csv"), stringsAsFactors = F)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
labels <- read.csv(paste0(dir2, 'SCINA_cell_type_labels.csv'))
pd.combined <- readRDS(paste0(dir1, "integrated_seurat_SCTransform_clustered-snn-0.4.rds"))
pd.combined$cellname <- sapply(rownames(pd.combined@meta.data), function(x) strsplit(x, "raw_data/")[[1]][2])
names(Idents(pd.combined)) <- sapply(names(Idents(pd.combined)), function(x) strsplit(x, "raw_data/")[[1]][2])
all(pd.combined@meta.data$cellname==labels[,2])
Idents(pd.combined) <- labels[,1]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# setting the slot of data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
i.c <- pd.combined
pd <- GetAssayData(i.c, slot = "scale.data")

# sample -> PD/Control
sample <- i.c@meta.data$sample

# cellname -> different type of cell names clustered
cellname<- labels$results..1..
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# selecting top n number of genes from the DEG analysis
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## selecting the list based on first n number of gene expression
#list <- head(mono_14$X,60)
list <- head(deg$X,300)
mat <- pd[list,]
coln<- colnames(mat)
#mat_con <- pd[list,1:24611]
#coln_con<- colnames(mat_con)
#mat_pd <- pd[list,24612:42480]
#coln_pd<- colnames(mat_pd)

# ~~~~~~~~~~~~~~~~~~~

## selecting the list for gene expression based on the thresholded p value adjusted
### get variables with p value less than 0.05
siggenes <- deg$X[deg$p_val_adj < 0.05]
mat <- pd[list,]
coln<- colnames(mat)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ComplexHeatmaps
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Simple Heatmap 
mycols <- colorRamp2(breaks = c(-2, 0, 2), colors = c("blue", "white", "red"))
png(paste0(data_dir, "Heatmap_MonoCD14_top60.png"), width=2000,height=800)
Heatmap(mat, name = "Heatmap of Top 60 gene in CD14+ Mono", col = mycols)
dev.off()

## Splitting the heatmap by row - dividing into 2 groups
new_mat <- t(mat)
mycols <- colorRamp2(breaks = c(-2, 0, 2), colors = c("blue", "white", "red"))
png(paste0(data_dir, "HM_mono14_split_1.png"), width=2000,height=800)
set.seed(2)
Heatmap(new_mat, name = "Heatmap of Top 60 gene in CD14+ Mono", col = mycols, k = 2)# control, pd
dev.off()

## Split by a vector specifying rowgroups 
mycols <- colorRamp2(breaks = c(-2, 0, 2),  colors = c("blue", "white", "red"))
pa <- read.csv(paste0(dir_out_files, 'patients.csv'))
png(paste0(data_dir, "HM_mono14_split_2.png"), width=2000,height=800)
Heatmap(new_mat, name = "Heatmap of Top 60 gene in CD14+ Mono", col = mycols, split = pa$x, row_names_gp = gpar(fontsize = 7))
dev.off()

## Not clustering columns automatically 
mycols <- colorRamp2(breaks = c(-2, 0, 2), colors = c("blue", "white", "red"))
png(paste0(data_dir, "Heatmap_MonoCD14_top60_supervised_sample_clustering.png"), width=2000,height=800)
Heatmap(mat, name = "Heatmap of Top 60 gene in CD14+ Mono", column_order = coln, cluster_columns = FALSE, col = mycols)
dev.off()


## Heatmap for control and pd separately
### Control Heatmap for CD+14
mycols <- colorRamp2(breaks = c(-2, 0, 2), colors = c("blue", "white", "red"))
png(paste0(data_dir, "HM_mono14_Control_0.png"), width=2000,height=800)
Heatmap(mat_con, name = "Heatmap of Top 60 gene in CD14+ Mono",column_order = coln_con, cluster_columns = FALSE, col = mycols)
dev.off()
###PD Heatmap for CD+14
mycols <- colorRamp2(breaks = c(-2, 0, 2), colors = c("blue", "white", "red"))
png(paste0(data_dir, "HM_mono14_PD_0.png"), width=2000,height=800)
Heatmap(mat_pd, name = "Heatmap of Top 60 gene in CD14+ Mono", column_order = coln_pd, cluster_columns = FALSE, col = mycols)
dev.off()


## With Annotation Bar - Sample
#### With clustering columns automatically 
ha = HeatmapAnnotation(sample = sample, annotation_name_side = "left")
ht_list = Heatmap(mat, name = "Heatmap of Top 300 differenciated gene expressions", 
                  col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
                  top_annotation = ha, 
                  show_column_names = FALSE, row_title = NULL, show_row_dend = FALSE) 
png(paste0(data_dir, "Heatmap_deg_top300_Ann_sample.png"),width=1000,height=1500)
draw(ht_list, row_title = "Genes")
dev.off()

#### Without clustering columns automatically 
ha = HeatmapAnnotation(sample = sample, annotation_name_side = "left")
ht_list = Heatmap(mat, name = "expression", 
                  col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
                  top_annotation = ha, 
                  show_column_names = FALSE, row_title = NULL, show_row_dend = FALSE, column_order = coln, cluster_columns = FALSE) 
png(paste0(data_dir, "Heatmap_MonoCD14_top60_Ann_supervised_sample_clustering.png"),width=1000,height=1400)
draw(ht_list, row_title = "Genes")
dev.off()

## With Annotation Bar - CellName
#### With clustering columns automatically 
ha = HeatmapAnnotation(cellname= cellname, annotation_name_side = "left")
ht_list = Heatmap(mat, name = "Heatmap of Top 300 differenciated gene expressions", 
                  col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
                  top_annotation = ha, 
                  show_column_names = FALSE, row_title = NULL, show_row_dend = FALSE) 
png(paste0(data_dir, "Heatmap_deg_top300_Ann_cellname.png"),width=1000,height=1500)
draw(ht_list, row_title = "Genes")
dev.off()


## With Annotation Bar - Sample, base_mean
ha = HeatmapAnnotation(sample = sample, annotation_name_side = "left")
base_mean = rowMeans(mat)
ht_list = Heatmap(mat_scaled, name = "expression", 
                  col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
                  top_annotation = ha, 
                  show_column_names = FALSE, row_title = NULL, show_row_dend = FALSE) +
  Heatmap(base_mean, name = "base mean", 
          top_annotation = HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2:6), 
                                                                    height = unit(2, "cm"))),
          width = unit(15, "mm"))
 
png(paste0(data_dir, "Heatmap_of_Top60_gene_CD14Mono_0.png"),width=2000,height=800)
draw(ht_list, row_title = "Genes")
dev.off()

## With Annotation Bar - Sample, cell name

#Set annotation
ann <- data.frame(sample, cellname)
colnames(ann) <- c("Sample", "Cell Name")
#colours <- list("Type"=c(aquamarine","N"="royalblue"), "Type2"=c("AT"="limegreen","PT"="gold"))
#colAnn <- HeatmapAnnotation(df=ann, which="col", col=colours, annotation_width=unit(c(1, 4), "cm"), gap=unit(1, "mm"))
colAnn <- HeatmapAnnotation(df=ann, annotation_width=unit(c(1, 4), "cm"), gap=unit(1, "mm"))
hmap <- Heatmap(
  data,
  name = "expression",
  col = greenred(75), 
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_column_dend = TRUE,
  show_row_dend = TRUE,
  row_dend_reorder = TRUE,
  column_dend_reorder = TRUE,
#  clustering_method_rows = "ward.D2",
 # clustering_method_columns = "ward.D2",
  width = unit(100, "mm"),
  top_annotation_height=unit(1.0,"cm"), top_annotation=colAnn)
  
  #bottom_annotation_height=unit(3, "cm"), bottom_annotation=boxplotCol)
png(paste0(data_dir, "Heatmap_of_Top60_gene_CD14Mono_Ann_sample_cellname.png"))
draw(hmap, heatmap_legend_side="left", annotation_legend_side="right", row_title = "Genes")
dev.off()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


ann <- data.frame(sample, cellname)
colnames(ann) <- c("Sample", "Cell Name")
colAnn <- HeatmapAnnotation(df=ann, annotation_width=unit(c(1, 4), "cm"))
hmap <- Heatmap(mat,
  name = "expression",
  col = greenred(75), 
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_column_dend = TRUE,
  show_row_dend = TRUE,
  row_dend_reorder = TRUE,
  column_dend_reorder = TRUE,
  width = unit(100, "mm"),
  top_annotation_height=unit(1.0,"cm"), top_annotation=colAnn)
png(paste0(data_dir, "Heatmap_of_Top60_gene_CD14Mono_Ann_sample_cellname.png"))
draw(hmap, heatmap_legend_side="left", annotation_legend_side="right", row_title = "Genes")
dev.off()
  
# ~~~~~~~~~~~~~~~~~~
ha = HeatmapAnnotation(ann, annotation_name_side = "left",  gap=unit(1, "mm"))
ht_list = Heatmap(mat, name = "Heatmap of Top 60 gene in CD14+ Mono", 
                  col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
                  top_annotation = ha, 
                  show_column_names = FALSE, row_title = NULL, show_row_dend = FALSE) 
png(paste0(data_dir, "Heatmap_MonoCD14_top60_Ann_sample_cellname.png"),width=1000,height=1500)
draw(ht_list, row_title = "Genes")
dev.off()

ha = HeatmapAnnotation(sample = sample, annotation_name_side = "left")
ht_list = Heatmap(mat, name = "Heatmap of Top 60 gene in CD14+ Mono", 
                  col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
                  top_annotation = ha, 
                  show_column_names = FALSE, row_title = NULL, show_row_dend = FALSE) 
png(paste0(data_dir, "Heatmap_MonoCD14_top60_Ann_sample.png"),width=1000,height=1500)
draw(ht_list, row_title = "Genes")
dev.off()





#pd_list <- subset (mono_16, mono_16$avg_logFC > 0 & mono_16$p_val_adj < 0.05)
#write.csv(pd_list, paste0(dir_out_files,"monocd16_pd_specific.csv"))
