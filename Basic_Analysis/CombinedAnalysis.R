library(Seurat)
library(cowplot)
library(dplyr)
library(hdf5r)
library(xlsx)

immune.combined<- readRDS("/Users/deepanshishokeen/desktop/SingleCellRNASeq/SingleCellData/withoutMT_immunecombined_rna.rds")

i.c <- RenameIdents(immune.combined, `0` = "CD14+ Monocytes", `1` = "CD4+ T", `2` = "CD4+ T", 
                    `3` = "CD8+ T", `4` = "NK", `5` = "CD8+ T", `6` = "CD16+ Monocytes", `7` = "CD14+ Monocytes", `8` = "B", `9` = "Dendritic", 
                    `10` = "NK", `11` = "CD14+ Monocytes", `12` = "Megakaryocyte", `13` = "B", `14`="Dendritic", 
                    `15`="NK",`16`="B", `17`="B")
DimPlot(i.c, label = TRUE)

#Showing violin plots for the following markers to show their impact on Monocytes cell 
#in comparison to other cells
plots <- VlnPlot(i.c, features = c("NLRP3", "PYCARD", "CASP1","GSDMD","IL1B"), split.by = "stim", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)

VlnPlot(i.c, features = "NLRP3", split.by = "stim", group.by = "celltype", 
        pt.size = 0, combine = T)
VlnPlot(i.c, features = "NLRP3", split.by = "stim", group.by = "celltype")

VlnPlot(i.c, features = "PYCARD", split.by = "stim", group.by = "celltype", 
        pt.size = 0, combine = FALSE)
VlnPlot(i.c, features = "PYCARD", split.by = "stim", group.by = "celltype")

VlnPlot(i.c, features = "CASP1", split.by = "stim", group.by = "celltype", 
        pt.size = 0, combine = TRUE)
VlnPlot(i.c, features = "CASP1", split.by = "stim", group.by = "celltype")

VlnPlot(i.c, features = "GSDMD", split.by = "stim", group.by = "celltype", 
        pt.size = 0, combine = FALSE)
VlnPlot(i.c, features = "GSDMD", split.by = "stim", group.by = "celltype")

VlnPlot(i.c, features = "IL1B", split.by = "stim", group.by = "celltype", 
        pt.size = 0, combine = FALSE)
VlnPlot(i.c, features = "IL1B", split.by = "stim", group.by = "celltype")



# DotPlot
#! get genes that are different in first plot and view those conserved cell type markers across conditions, showing both the expression level and the percentage of cells in a cluster expressing
Idents(immune.combined) <- factor(Idents(immune.combined), levels = rev(unique(c("CD14+ Mono", "CD4 Naive T", "CD4 Memory T", "CD8+ T", 
                                                                                 "CD16+ Mono", "Natural killer", "CD8+ T", "CD14+ Mono", 
                                                                                 "CD8+ T", "B", "CD14+ Mono", "Dendritic", "Megakaryocyte", "Megakaryocyte", "CD4+ T/ CD8+ T", 
                                                                                 "Megakaryocyte", "B cell subtype", "MT", "Natural killer", "Megakaryocyte", "B cell subtype")))
)

#markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5", 
#                    "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1", 
#                   "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ")
markers.to.plot <- c("FCGR2A", "VAMP4", "KCNS3", "KCNIP3", "LINC00693", "KPNA1", "MED12L",
                     "SPTSSB", "LCORL", "CLCN3", "PAM", "C5orf24", "TRIM40", "FYN", "RPS12",
                     "GS1-124K5.11", "FAM49B", "UBAP2", "GBF1", "RNF141", "SCAF11", "FBRSL1",
                     "CAB39L", "MBNL2", "MIPOL1", "RPS6KL1", "CD19", "NOD2", "CNOT1", "CHRNB1",
                     "UBTF", "FAM171A2", "BRIP1", "DNAH17", "ASXL3", "MEX3C", "CRLS1", "DYRK1A")

#DotPlot(pbmc, features = features) + RotatedAxis()
DotPlot(i.c, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, 
        split.by = "stim") + RotatedAxis()

# Single cell heatmap of feature expression
#DoHeatmap(subset(pbmc, downsample = 100), features = features, size = 3)
DoHeatmap(subset(immune.combined, downsample = 100), rev(markers.to.plot), size = 3)
# dataset$variable <- factor(dataset$variable, levels = rev(unique(dataset$variable)), ordered=TRUE)





# Finding differenciated cells

i.c$celltype.stim <- paste(Idents(i.c), i.c$stim, sep = "_")
i.c$celltype <- Idents(i.c)
Idents(i.c) <- "celltype.stim"
b.interferon.response <- FindMarkers(i.c, ident.1 = "B_STIM", ident.2 = "B_CTRL", verbose = FALSE)
head(b.interferon.response, n = 200)
dim(b.interferon.response)
#[1] 32  5

cd14mono.interferon.response <- FindMarkers(i.c, ident.1 = "CD14+ Monocytes_STIM", ident.2 = "CD14+ Monocytes_CTRL", verbose = FALSE)
head(cd14mono.interferon.response, n = 400)
dim(cd14mono.interferon.response)
#[1] 41  5

cd4T.interferon.response <- FindMarkers(i.c, ident.1 = "CD4+ T_STIM", ident.2 = "CD4+ T_CTRL", verbose = FALSE)
head(cd4T.interferon.response, n = 200)
dim(cd4T.interferon.response)
#[1] 11  5


cd8T.interferon.response <- FindMarkers(i.c, ident.1 = "CD8+ T_STIM", ident.2 = "CD8+ T_CTRL", verbose = FALSE)
head(cd8T.interferon.response, n = 200)
dim(cd8T.interferon.response)
#[1] 47  5

cd16mono.interferon.response <- FindMarkers(i.c, ident.1 = "CD16+ Monocytes_STIM", ident.2 = "CD16+ Monocytes_CTRL", verbose = FALSE)
head(cd16mono.interferon.response, n = 200)
dim(cd16mono.interferon.response)
#[1] 33  5

nk.interferon.response <- FindMarkers(i.c, ident.1 = "NK_STIM", ident.2 = "NK_CTRL", verbose = FALSE)
head(nk.interferon.response, n = 200)
dim(nk.interferon.response)
#[1] 22  5

dendritic.interferon.response <- FindMarkers(i.c, ident.1 = "Dendritic_STIM", ident.2 = "Dendritic_CTRL", verbose = FALSE)
head(dendritic.interferon.response, n = 200)
dim(dendritic.interferon.response)
#[1] 27  5

mega.interferon.response <- FindMarkers(i.c, ident.1 = "Megakaryocyte_STIM", ident.2 = "Megakaryocyte_CTRL", verbose = FALSE)
head(mega.interferon.response, n = 1000)
dim(mega.interferon.response)
#[1] 783   5

##part of dendritic
#unknown.interferon.response <- FindMarkers(i.c, ident.1 = "Unknown_STIM", ident.2 = "Unknown_CTRL", verbose = FALSE)
#head(unknown.interferon.response, n = 200)
#dim(unknown.interferon.response)
#[1] 296   5


#bsubtype.interferon.response <- FindMarkers(i.c, ident.1 = "B cell subtype_STIM", ident.2 = "B cell subtype_CTRL", verbose = FALSE)
#head(bsubtype.interferon.response, n = 200)
#dim(bsubtype.interferon.response)
#[1] 116   5


FeaturePlot(i.c, features = c("RPS4Y1","RPS4X", "HLA-DRB5","XIST"), split.by = "stim", max.cutoff = 3, 
            cols = c("grey", "red"))

