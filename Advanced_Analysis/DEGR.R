########## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Differentially expressed genes 
# Lab: Havrada, Dartmouth College
# Author: Deepanshi Shokeen, Owen Wilkins

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
library(scater)

rm(list=ls())


dir1 <- "/dartfs-hpc/rc/lab/H/HavrdaM/scrna-seq/01_pre_processing/files/"
dir2 <- "/dartfs-hpc/rc/lab/H/HavrdaM/scrna-seq/02_cell-type-classification/files/"
dir_raw <- "/dartfs-hpc/rc/lab/H/HavrdaM/scrna-seq/01_pre_processing/raw_data/"
# change the below directory to save the output in different directory.
dir_out_files <- "/dartfs-hpc/rc/lab/H/HavrdaM/scrna-seq/03_differential_expression/files/"
dir_out_figs <- "/dartfs-hpc/rc/lab/H/HavrdaM/scrna-seq/03_differential_expression/figures/"

# load SEurat dataset 
pd.combined <- readRDS(paste0(dir1, "integrated_seurat_SCTransform_clustered-snn-0.4.rds"))


DefaultAssay(pd.combined) <- "RNA"

pd.combined.sce <- as.SingleCellExperiment(pbmc)
saveRDS(pd.combined.sce, file = paste0(dir1, "integrated_seurat_SCTransform_clustered-snn-0.4-SCE.rds"))

saveRDS(pd.combined.sce, file = paste0(dir1, "integrated_seurat_SCTransform_clustered-snn-0.4-SCE.rds"))


pd.combined <- NormalizeData(pd.combined)
pd.combined <- FindVariableFeatures(pd.combined, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)
# Scale the counts
pd.combined <- ScaleData(pd.combined)
pd.combined <- RunPCA(object = pd.combined)

p <- PCAPlot(pd.combined, split.by = "sex", group.by="sample")  
save_plot(paste0(dir_out_figs, "PCA_sex_RNA.png"), p, base_height = 7, base_width = 15)
p <- PCAPlot(pd.combined, split.by = "sex")  
save_plot(paste0(dir_out_figs, "PCA_sex_RNA_cell-types.png"), p, base_height = 7, base_width = 15)




DefaultAssay(pd.combined) <- "integrated"
pd.combined <- RunPCA(object = pd.combined)
p <- PCAPlot(pd.combined, split.by = "sex", group.by="sample")  
save_plot(paste0(dir_out_figs, "PCA_sex_integrated.png"), p, base_height = 7, base_width = 15)
p <- PCAPlot(pd.combined, split.by = "sex")  
save_plot(paste0(dir_out_figs, "PCA_sex_integrated2.png"), p, base_height = 7, base_width = 15)






# Run UMAP
pd.combined <- RunUMAP(pd.combined, dims = 1:40)                          
p <- DimPlot(pd.combined, reduction = "umap", group.by = "sample")
save_plot(paste0(dir_out_figs, "UMAP_sample_RNA.png"), p, base_height = 7, base_width = 8)

pd.combined2 <- RunUMAP(pd.combined, dims = 1:40)                            
p <- DimPlot(pd.combined2, reduction = "umap", group.by = "sample")
save_plot(paste0(dir_out_figs, "UMAP_sample_RNA_2.png"), p, base_height = 7, base_width = 8)






# Plot the PCA colored by ell cycle phase
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")


DefaultAssay(pd.combined) <- "SCT"
head(pd.combined@meta.data)

# Run PCA
pd.combined <- RunPCA(object = pd.combined)

# Plot PCA
PCAPlot(seurat_integrated,
        split.by = "sample")  


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set up objects for DE analysis of each cell type 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

pd.combined@meta.data$sex2 <- NA
pd.combined@meta.data$sex2[pd.combined@meta.data$sex == "F"] <- 1
pd.combined@meta.data$sex2[pd.combined@meta.data$sex == "M"] <- 0

# add subject age to metadata 
pd.combined@meta.data$age <- NA
pd.combined@meta.data$age[pd.combined@meta.data$subject == "MH1"] <- 77
pd.combined@meta.data$age[pd.combined@meta.data$subject == "MH3"] <- 58
pd.combined@meta.data$age[pd.combined@meta.data$subject == "PD125"] <- 82
pd.combined@meta.data$age[pd.combined@meta.data$subject == "PD133"] <- 78
pd.combined@meta.data$age[pd.combined@meta.data$subject == "PD135"] <- 79
pd.combined@meta.data$age[pd.combined@meta.data$subject == "3401HC"] <- 61
pd.combined@meta.data$age[pd.combined@meta.data$subject == "3402HC"] <- 60
pd.combined@meta.data$age[pd.combined@meta.data$subject == "MH2"] <- 74
pd.combined@meta.data$age[pd.combined@meta.data$subject == "PD123"] <- 73
pd.combined@meta.data$age[pd.combined@meta.data$subject == "PD126"] <- 83
pd.combined@meta.data$age[pd.combined@meta.data$subject == "PD134"] <- 57

# add subject as a numeric variable (in case needed)
pd.combined@meta.data$subject2 <- NA
pd.combined@meta.data$subject2[pd.combined@meta.data$subject == "MH1"] <- 1
pd.combined@meta.data$subject2[pd.combined@meta.data$subject == "MH3"] <- 2
pd.combined@meta.data$subject2[pd.combined@meta.data$subject == "PD125"] <- 3
pd.combined@meta.data$subject2[pd.combined@meta.data$subject == "PD133"] <- 4
pd.combined@meta.data$subject2[pd.combined@meta.data$subject == "PD135"] <- 5
pd.combined@meta.data$subject2[pd.combined@meta.data$subject == "3401HC"] <- 6
pd.combined@meta.data$subject2[pd.combined@meta.data$subject == "3402HC"] <- 7
pd.combined@meta.data$subject2[pd.combined@meta.data$subject == "MH2"] <- 8
pd.combined@meta.data$subject2[pd.combined@meta.data$subject == "PD123"] <- 9
pd.combined@meta.data$subject2[pd.combined@meta.data$subject == "PD126"] <- 10
pd.combined@meta.data$subject2[pd.combined@meta.data$subject == "PD134"] <- 11

# combine cell type with sample type so that we can specify cells from case/control in DE analysis 
pd.combined$celltype.sample <- paste(Idents(pd.combined), pd.combined$sample, sep = "_")

# set these as cell idendities 
Idents(pd.combined) <- "celltype.sample"



library(ggplot2)
# Change colors
p+scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
# Add box plots
ggplot(ToothGrowth, aes(x=dose, y=len, color=supp)) +
  geom_boxplot(color="black")+
  geom_jitter(position=position_jitter(0.2))
# Change the position
ggplot(ToothGrowth, aes(x=dose, y=len, color=supp)) +
  geom_boxplot(position=position_dodge(0.8))+
  geom_jitter(position=position_dodge(0.8))

p<-ggplot(ToothGrowth, aes(x=dose, y=len, color=supp, shape=supp)) +
  geom_jitter(position=position_dodge(0.8))
save_plot(paste0(dir_out_figs, "test.png"), p, base_height = 7, base_width = 6)



cd14_counts <- pd.combined@assays$SCT@counts["S100A9", rownames(pd.combined@meta.data)[pd.combined@meta.data$celltype==paste0("CD14_monocyte")]]
cd14_counts[cd14_counts==0] <- NA
mean(cd14_counts, na.rm=TRUE)

cd14_S100A9 <- cd14_counts
cd14_cells <- rownames(pd.combined@meta.data)[pd.combined@meta.data$celltype==paste0("CD14_monocyte")]
cd14_sex <- pd.combined@meta.data$sex[pd.combined@meta.data$celltype==paste0("CD14_monocyte")]
cd14_group <- pd.combined@meta.data$sample[pd.combined@meta.data$celltype==paste0("CD14_monocyte")]
cd14_subject <- pd.combined@meta.data$subject[pd.combined@meta.data$celltype==paste0("CD14_monocyte")]

df <- data.frame(cd14_S100A9, cd14_cells, cd14_sex, cd14_group, cd14_subject)

# Add violin plot
p<-ggplot(df, aes(x=cd14_group, y=cd14_S100A9, color=cd14_sex)) + 
  geom_violin(trim = FALSE)+
  #geom_jitter(position=position_jitter(0.01)) +
  geom_boxplot(width=0.3, position=position_dodge(0.9)) + 
  ggtitle("CD14 Monocytes - S100A9")+
  geom_jitter(position=position_dodge(0.9))+
  labs(x="Sample group", y="SCT normalized counts/cell")
  #geom_jitter(position=position_dodge(0.9))
save_plot(paste0(dir_out_figs, "S100A9_all_counts-by_group.png"), p, base_height = 7, base_width = 6)

table(pd.combined@meta.data$subject, pd.combined@meta.data$sample)

subject <- c("MH1", "MH3", "PD125", "PD133", "PD135", "3401HC", "3402HC", "MH2", "PD123", "PD126", "PD134")
sex <- c("F", "F", "F", "F", "F", "M", "M", "M", "M", "M", "M")
sample <- c("pd", "pd", "control", "control", "control", "control", "control", "control", "pd", "pd", "pd")
df_2 <- data.frame(subject, sex, sample)
counts_goi <- c()
for(i in 1:length(sex)){
  counts_goi[i] <- mean(df$cd14_S100A9[df$cd14_subject==subject[i]], na.rm=TRUE)
}
df_2$counts <- counts_goi
df_2$subject <- as.character(df_2$subject)
df_2$sex <- as.character(df_2$sex)
df_2$sample <- as.character(df_2$sample)
df_2$subject <- as.character(df_2$subject)

# Add violin plot
p<-ggplot(df_2, aes(x=sample, y=counts, color=sex)) + 
  #geom_jitter(position=position_jitter(0.01)) +
  #geom_boxplot(width=0.3, position=position_dodge(0.9)) + 
  #geom_violin(trim=F) + 
  geom_jitter(position=position_dodge(0.9)) +
  ggtitle("CD14 Monocytes - S100A9")+ geom_text(aes(label=subject))+
  labs(x="Sample group", y="Mean SCT Counts (across expressing cells per subject)")
  #geom_boxplot(position=position_dodge(0.9), width=0.3)+ 
  #stat_summary(fun.y=median, geom="point", shape=18, size=3, color="red", position=position_dodge(0.9))
  #geom_jitter(position=position_dodge(0.9))
save_plot(paste0(dir_out_figs, "S100A9_mean_counts-by_group.png"), p, base_height = 7, base_width = 6)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run DE analysis on each cell type 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



##### IS MAST APPROPRIATE TO USE ONF SCT DATA 





# write function to do DE analysis and annotation for each cell type recursively 
run_DE <- function(cell_type, seurat_object){
  #seurat_object <- pd.combined
  #cell_type <- uniq_cells[9]
  
  # run DE test 
  results <- FindMarkers(seurat_object, assay = "SCT", 
                        slot = "data", 
                        ident.1 = paste0(cell_type, "_pd"), ident.2 = paste0(cell_type, "_control"),
                        test.use = "MAST", 
                        min.pct = 0.1,
                        logfc.threshold = 0.0,
                        latent.vars = c("sex2","age"), 
                        verbose = FALSE)
  # add gene symbol as column 
  results$symbol <- rownames(results)
  # add ENSG IDs from features output of CellRanger 
  genes <- read.table(paste0(dir_raw, "3401HC/filtered_feature_bc_matrix/features.tsv.gz"), stringsAsFactors = FALSE)
  # get just relevant columns 
  genes <- genes[,c(1,2)]
  #table(rownames(b_cell) %in% genes$V2)
  #table(genes$V2 %in% rownames(b_cell))
  # get gene IDs for those in results 
  genes_sub <- genes[genes$V2 %in% rownames(results),]
  # drop duplicates 
  ind1_dups <- which(duplicated(genes_sub$V2))
  genes_sub <- genes_sub[-ind1_dups,]
  colnames(genes_sub) <- c("ENSG", "symbol")
  # merge with Seurat results 
  merge1 <- merge(results, genes_sub, by = "symbol")
  
  # add gene annotation (CHR, start site)
  edb <- EnsDb.Hsapiens.v86
  anno <- select(edb, keys=as.character(merge1$ENSG), 
                 columns=c("GENEBIOTYPE", "SEQNAME", "GENESEQSTART"),
                 keytype="GENEID")
  colnames(anno)[1] <- "ENSG"
  merge2 <- merge(merge1, anno, by="ENSG", all.x = TRUE)
  
  # add FDR, log2FC, reorder columns 
  merge2$fdr <- p.adjust(merge2$p_val, method = "fdr")
  merge2$log2FC <- log2(exp(abs(merge2$avg_logFC)))
  merge2$log2FC[merge2$avg_logFC<0] <- (merge2$log2FC[merge2$avg_logFC<0])* -1
  merge2 <- merge2[order(merge2$fdr),]
  merge3 <- merge2[,c("ENSG", "symbol", "symbol", "SEQNAME", "GENESEQSTART", 
                      "pct.1", "pct.2", "p_val", "p_val_adj", "fdr", "avg_logFC", "log2FC")]
  # add mean SCT counts per cell to results for each gene 
  genes <- merge3$symbol
  merge3$mean_SCT_count <- NA
  for(i in 1:nrow(merge3)){
    SCT_count <- seurat_object@assays$SCT@data[genes[i], rownames(seurat_object@meta.data)[seurat_object@meta.data$celltype==paste0(cell_type)]]
    SCT_count[SCT_count==0] <- NA
    mean_SCT_count <- mean(SCT_count, na.rm=TRUE)
    merge3$mean_SCT_count[i] <- mean_SCT_count
  }
  merge3$mean_SCT_countlog <- log(merge3$mean_SCT_count)
  
  # write results to file 
  write.csv(merge3, paste0(dir_out_files, cell_type, "_MAST_SCT_DEGs.csv"), row.names = FALSE)
  write.csv(merge3[merge3$fdr<0.05,], paste0(dir_out_files, cell_type, "_MAST_SCT_DEGs_FDR_1%.csv"), row.names = FALSE)
  return(merge3)
}

# apply function to each cell type get get annotated DEG results tables 
uniq_cells <- unique(pd.combined@meta.data$celltype)
#res_all_types <- lapply(as.list(uniq_cells), run_DE, pd.combined)

r1 <- run_DE(uniq_cells[1], pd.combined)
r2 <- run_DE(uniq_cells[2], pd.combined)
r3 <- run_DE(uniq_cells[3], pd.combined)
r4 <- run_DE(uniq_cells[4], pd.combined)
r5 <- run_DE(uniq_cells[5], pd.combined)
r6 <- run_DE(uniq_cells[6], pd.combined)
r7 <- run_DE(uniq_cells[7], pd.combined)
#r8 <- run_DE(uniq_cells[8], pd.combined)
r9 <- run_DE(uniq_cells[9], pd.combined) #r9 <- merge3
r10 <- run_DE(uniq_cells[10], pd.combined)
r11 <- run_DE(uniq_cells[11], pd.combined)

# merge into list and save 
master_res_list <- list(r1, r2, r3, r4, r5, r6, r7, r9, r10, r11)
names(master_res_list) <- uniq_cells[-8] # remove unknown name
saveRDS(master_res_list, paste0(dir_out_files, "DEG_results_SCINA_all_cell-types-list.rds"))

# Single cell heatmap of feature expression
x$symbol[x$fdr<0.01]

markers.to.plot <- unlist(lapply(master_res_list[[3]], function(x) {
  x2 <- x[x$fdr<0.01,]
  x3 <- x2[order(abs(x2$log2FC)),]
  x3$symbol[1:20]
}
  ))
markers.to.plot[["CD14_monocyte"]]
names(master_res_list)


#DefaultAssay(pd.combined_sub) <- "RNA"
#Idents(pd.combined_sub) <- pd.combined_sub@meta.data$celltype
#all.genes <- rownames(pd.combined_sub)
#pd.combined_sub <- ScaleData(pd.combined_sub, features = all.genes)
#p <- DoHeatmap(subset(pd.combined_sub, downsample = 100), rev(markers.to.plot), size = 3, group.by = "ident")
#save_plot(paste0(dir_out_figs, "heatmap_SCINA-celltypes_DEGs_FDR-1.png"), p, base_height = 7, base_width = 15)





table(pd.combined_sub@meta.data$sample)


ind_unk <- which(pd.combined@meta.data$celltype=="unknown")
pd.combined_sub <- pd.combined[, !rownames(pd.combined@meta.data) %in% rownames(pd.combined@meta.data)[ind_unk]]


master_res_list <- readRDS(paste0(dir_out_files, "DEG_results_SCINA_all_cell-types-list.rds"))

DefaultAssay(pd.combined_sub) <- "SCT"

table(pd.combined_sub@meta.data$celltype.sample)
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



uniq_cells

cell_type <- "CD14_monocyte"

plot_dimplot <- function(cell_type, cell_type_name){
res_celltype <- master_res_list[[cell_type]]
res_celltype <- res_celltype[order(abs(res_celltype$log2FC), decreasing=TRUE),]
res_celltype_sub <- res_celltype[1:15,]

p <- DotPlot(pd.combined_sub, features = res_celltype_sub$symbol, 
  cols = rep(c("blue", "red"), 10), 
  dot.scale = 8, 
  group.by = "celltype", 
  split.by = "sample") + RotatedAxis() + ggtitle(paste0(cell_type_name, " DEGs (top DEGs based on Fold change and FDR<0.01)"))
save_plot(paste0(dir_out_figs, "DimPlot_", cell_type_name, "_top-15_DEGs_log2FC.png"), p, base_height = 7, base_width = 15)

}
plot_dimplot(uniq_cells[1], "CD4 T-cell")
plot_dimplot(uniq_cells[2], "NK-cell")
plot_dimplot(uniq_cells[3], "CD14 monocyte")
plot_dimplot(uniq_cells[4], "B-cell")
plot_dimplot(uniq_cells[5], "CD16 monocyte")
plot_dimplot(uniq_cells[6], "Cytotoxic T-cell")
plot_dimplot(uniq_cells[7], "Dendritic cell")
plot_dimplot(uniq_cells[9], "Megakaryocyte")
plot_dimplot(uniq_cells[10], "Plasma_cell")
plot_dimplot(uniq_cells[11], "Plasmacytoid_dendritic_cell")




gene <- "ISG15"
cell_type <- "CD14_monocyte"
p <- FeaturePlot(pd.combined_sub, features = c(gene), split.by = "sample", 
                 min.cutoff = 'q1', max.cutoff = "q10",
                 cols = c("grey", "red"))
save_plot(paste0(dir_out_figs, "FeaturePlot_", cell_type, "_", gene, ".png"), p, base_height = 7, base_width = 15)


p <- VlnPlot(pd.combined_sub, features = gene, split.by = "sample", pt.size=0)
save_plot(paste0(dir_out_figs, "VlnPlot_", cell_type, "_", gene, ".png"), p, base_height = 7, base_width = 15)






