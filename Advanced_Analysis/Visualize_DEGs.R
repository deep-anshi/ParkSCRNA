########## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Lab: Havrada 
# Author: Deepanshi Shokeen

########## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


library(ggplot2)
library(ggrepel)
library(Seurat)
library(cowplot)
library(jackstraw)
library(lfa)
library(readr)
library(rtracklayer)
library(SCINA)
library(preprocessCore)
library(stringr)


dir1 <- "/dartfs-hpc/rc/lab/H/HavrdaM/scrna-seq/01_pre_processing/files/"
dir2 <- "/dartfs-hpc/rc/lab/H/HavrdaM/scrna-seq/02_cell-type-classification/files/"
dir3 <- "/dartfs-hpc/rc/lab/H/HavrdaM/scrna-seq/03_differential_expression/files/"
dir_out_figs <- "/dartfs-hpc/rc/lab/H/HavrdaM/scrna-seq/03_differential_expression/figures/Boxplot/"

# load custom volcano plotting function 
source(paste0("/dartfs-hpc/rc/lab/H/HavrdaM/scrna-seq/misc/box_plot.R"))

# load SEurat dataset 
pd.combined <- readRDS(paste0(dir1, "integrated_seurat_SCTransform_clustered-snn-0.4.rds"))
#DefaultAssay(pd.combined) <- "SCT"
head(pd.combined@meta.data)

# add individual cell names as variable to metaddata 
pd.combined$cellname <- sapply(rownames(pd.combined@meta.data), function(x) strsplit(x, "raw_data/")[[1]][2])

# read in cell type labels 
labels <- read.csv(paste0(dir2, 'SCINA_cell_type_labels.csv'))

# check cells are in same order 
all(pd.combined@meta.data$cellname==labels[,2])

# add column to metaddata for cell types 
pd.combined$celltype <- as.character(labels[,1])

# set Idents as cell tyep labels 
Idents(pd.combined) <- as.character(labels[,1])

# extract subjects that each cell came from from cellnames 
### need to fix MH samples as they are named with _ before 1,2,3 for 3 samples 
a = str_replace_all(pd.combined@meta.data$cellname, "MH_1", "MH1")
a = str_replace_all(a, "MH_2", "MH2")
a = str_replace_all(a, "MH_3", "MH3")
pd.combined@meta.data$subject = sapply(a, function(x) strsplit(x, "_")[[1]][1])

# add subject sex to metadata 
pd.combined@meta.data$sex <- NA
pd.combined@meta.data$sex[pd.combined@meta.data$subject == "MH1"] <- "F"
pd.combined@meta.data$sex[pd.combined@meta.data$subject == "MH3"] <- "F"
pd.combined@meta.data$sex[pd.combined@meta.data$subject == "PD125"] <- "F"
pd.combined@meta.data$sex[pd.combined@meta.data$subject == "PD133"] <- "F"
pd.combined@meta.data$sex[pd.combined@meta.data$subject == "PD135"] <- "F"
pd.combined@meta.data$sex[pd.combined@meta.data$subject == "3401HC"] <- "M"
pd.combined@meta.data$sex[pd.combined@meta.data$subject == "3402HC"] <- "M"
pd.combined@meta.data$sex[pd.combined@meta.data$subject == "MH2"] <- "M"
pd.combined@meta.data$sex[pd.combined@meta.data$subject == "PD123"] <- "M"
pd.combined@meta.data$sex[pd.combined@meta.data$subject == "PD126"] <- "M"
pd.combined@meta.data$sex[pd.combined@meta.data$subject == "PD134"] <- "M"

# add subject as a numeric variable (in case needed)
pd.combined@meta.data$subject2 <- NA
pd.combined@meta.data$subject2[pd.combined@meta.data$subject == "MH1"] <- "PD1"
pd.combined@meta.data$subject2[pd.combined@meta.data$subject == "MH3"] <- "PD2"
pd.combined@meta.data$subject2[pd.combined@meta.data$subject == "PD125"] <- "C1"
pd.combined@meta.data$subject2[pd.combined@meta.data$subject == "PD133"] <- "C2"
pd.combined@meta.data$subject2[pd.combined@meta.data$subject == "PD135"] <- "C3"
pd.combined@meta.data$subject2[pd.combined@meta.data$subject == "3401HC"] <- "C4"
pd.combined@meta.data$subject2[pd.combined@meta.data$subject == "3402HC"] <- "C5"
pd.combined@meta.data$subject2[pd.combined@meta.data$subject == "MH2"] <- "C6"
pd.combined@meta.data$subject2[pd.combined@meta.data$subject == "PD123"] <- "PD3"
pd.combined@meta.data$subject2[pd.combined@meta.data$subject == "PD126"] <- "PD4"
pd.combined@meta.data$subject2[pd.combined@meta.data$subject == "PD134"] <- "PD5"

pd.combined$celltype.sample <- paste(Idents(pd.combined), pd.combined$sample, sep = "_")

# set these as cell idendities 
Idents(pd.combined) <- "celltype.sample"


# trial


genes <- read.csv(paste0(dir3, "gene_cd14_cd16.csv"), stringsAsFactors = FALSE)
genes[genes==""] <- NA
Uni_dataframe <- data.frame(d = unlist(genes, use.names = FALSE))
Uni_dataframe[!is.na(Uni_dataframe$d), ]
genes <- Uni_dataframe$d
genes <- unique(genes)
genes <- as.character(genes)
genes <- genes[!is.na(genes) ]


for(i in genes){
  gene_name<- i
  data<- pd.combined
  out_dir<- dir_out_figs
  out_name<- paste0(gene_name, "_boxplot_sex_sample.png")
  title<- paste0(gene_name)
  
  box_plot_sex_sample(paste0(gene_name),
                      data,
                      out_dir,
                      out_name,
                      title)
}


