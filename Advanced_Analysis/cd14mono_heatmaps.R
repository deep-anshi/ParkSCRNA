########## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Heatmaps for top n differentially expressed genes for specefic cell type, example - cd14+ monocytes 
# Lab: Havrada 
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
cellname<- as.character(labels$results..1..)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 


### extracting the rows ###
# selecting top n number of genes from the DEG analysis
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## selecting the list based on first n number of gene expression
list <- head(mono_14$X,100)
mat <- pd[list,]
coln<- colnames(mat)
 

## selecting the list for gene expression based on the thresholded p value adjusted
### get variables with p value less than 0.05
siggenes <- mono_14$X$X[mono_14$X$p_val_adj < 0.05]
#siggenes <- deg$X[deg$p_val_adj < 0.05]
mat <- pd[list,]
coln<- colnames(mat)  
 
# Creating new matrix after subsetting the cell barcodes
new_mat <-t(rbind(mat, cellname, sample))
mat_for_specific_cellname<- subset(new_mat, cellname=='CD14_monocyte', select = -cellname)
specific_sample <- mat_for_specific_cellname[,101]
matr<- t(mat_for_specific_cellname[,-101])


### SD ###
#SD <- apply(matr,1,sd) #SD <- apply(matr,1,sd)[200]
#matr<- cbind(matr,SD)
#matr <- matr[order(matr[,12025], decreasing = TRUE),] #######

### Extracting the column names - cellname barcodes ###
# ->rowlist <- head(rownames(matr),60)
collist <- colnames(matr)
# ->matt <- mat[rowlist,] 
matt <- mat[,collist]


### heatmap ###
ha = HeatmapAnnotation(sample = specific_sample, annotation_name_side = "left")
ht_list = Heatmap(matt, name = "Heatmap of Top 100 differenciated gene expressions in CD14+ Monocytes", 
                  col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
                  top_annotation = ha, 
                  show_column_names = FALSE, row_title = NULL, show_row_dend = FALSE) 
png(paste0(data_dir, "Heatmap_monocd14_top100_Ann_sample.png"),width=1000,height=1500)
draw(ht_list, row_title = "Genes")
dev.off()

 
colnn<- colnames(matt)
ha = HeatmapAnnotation(sample = specific_sample, annotation_name_side = "left")
ht_list = Heatmap(matt, name = "Heatmap of Top 100 differenciated gene expressions in CD14+ Monocytes", 
                  col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
                  top_annotation = ha, 
                  show_column_names = FALSE, row_title = NULL, show_row_dend = FALSE, column_order = colnn, cluster_columns = FALSE) 
png(paste0(data_dir, "Heatmap_MonoCD14_top100_Ann_supervised_sample_clustering.png"),width=1000,height=1400)
draw(ht_list, row_title = "Genes")
dev.off()
