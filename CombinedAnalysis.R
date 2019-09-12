library(Seurat)
library(cowplot)
library(dplyr)
library(hdf5r)
library(xlsx)

immune.combined<- readRDS("/Users/deepanshishokeen/desktop/SingleCellRNASeq/SingleCellData/immunecombined_rna.rds")

#Find markers/conserved markers for each Clusters
markers.0 <- FindConservedMarkers(immune.combined, ident.1 = 0, grouping.var = "stim", verbose = FALSE)
head(markers.0)
markers.1 <- FindConservedMarkers(immune.combined, ident.1 = 1, grouping.var = "stim", verbose = FALSE)
head(markers.1)
markers.2 <- FindConservedMarkers(immune.combined, ident.1 = 2, grouping.var = "stim", verbose = FALSE)
head(markers.2)
markers.3 <- FindConservedMarkers(immune.combined, ident.1 = 3, grouping.var = "stim", verbose = FALSE)
head(markers.3)
markers.4 <- FindConservedMarkers(immune.combined, ident.1 = 4, grouping.var = "stim", verbose = FALSE)
head(markers.4)
markers.5 <- FindConservedMarkers(immune.combined, ident.1 = 5, grouping.var = "stim", verbose = FALSE)
head(markers.5)
markers.6 <- FindConservedMarkers(immune.combined, ident.1 = 6, grouping.var = "stim", verbose = FALSE)
head(markers.6)
markers.7 <- FindConservedMarkers(immune.combined, ident.1 = 7, grouping.var = "stim", verbose = FALSE)
head(markers.7)
markers.8 <- FindConservedMarkers(immune.combined, ident.1 = 8, grouping.var = "stim", verbose = FALSE)
head(markers.8)
markers.9 <- FindConservedMarkers(immune.combined, ident.1 = 9, grouping.var = "stim", verbose = FALSE)
head(markers.9)
markers.10 <- FindConservedMarkers(immune.combined, ident.1 = 10, grouping.var = "stim", verbose = FALSE)
head(markers.10)
markers.11 <- FindConservedMarkers(immune.combined, ident.1 = 11, grouping.var = "stim", verbose = FALSE)
head(markers.11)
markers.12 <- FindConservedMarkers(immune.combined, ident.1 = 12, grouping.var = "stim", verbose = FALSE)
head(markers.12)
markers.13 <- FindConservedMarkers(immune.combined, ident.1 = 13, grouping.var = "stim", verbose = FALSE)
head(markers.13)
markers.14 <- FindConservedMarkers(immune.combined, ident.1 = 14, grouping.var = "stim", verbose = FALSE)
head(markers.14)
markers.15 <- FindConservedMarkers(immune.combined, ident.1 = 15, grouping.var = "stim", verbose = FALSE)
head(markers.15)
markers.16 <- FindConservedMarkers(immune.combined, ident.1 = 16, grouping.var = "stim", verbose = FALSE)
head(markers.16)
markers.17 <- FindConservedMarkers(immune.combined, ident.1 = 17, grouping.var = "stim", verbose = FALSE)
head(markers.17)
markers.18 <- FindConservedMarkers(immune.combined, ident.1 = 18, grouping.var = "stim", verbose = FALSE)
head(markers.18)
markers.19 <- FindConservedMarkers(immune.combined, ident.1 = 19, grouping.var = "stim", verbose = FALSE)
head(markers.19)
markers.20 <- FindConservedMarkers(immune.combined, ident.1 = 20, grouping.var = "stim", verbose = FALSE)
head(markers.20)

#Cluster Tree
pbmc33k <- BuildClusterTree(object=immune.combined, slot = "scale.data",dims=1:20)
PlotClusterTree(pbmc33k)

#We can explore these marker genes for each cluster and use them to annotate our clusters as specific cell types.
FeaturePlot(immune.combined, features = c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A", 
                                          "CCL2", "PPBP"), min.cutoff = "q9")

#Feature plots
FeaturePlot(immune.combined, features = c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A", 
                                          "CCL2", "PPBP"), min.cutoff = "q9")
FeaturePlot(immune.combined, features = c("IL7R", "CD4", "CD14", "LYZ", "MS4A1","MS4A7", "NKG7", "FCER1A", "CST3"), min.cutoff = "q9")
FeaturePlot(immune.combined, features = c("CLEC7A",  "CLEC4C", "NRP2", "IL3RA","CD1C", "S100A4"), min.cutoff = "q9")


#vln plots
#general
VlnPlot(object = immune.combined, features = c("CD3D",  "IL7R", "CD4", "SELL", "CREM", "CD14", "LYZ", "MS4A1", "FCGR3A", "MS4A7",  
                                               "NKG7", "FCER1A", "CST3", "CLEC7A",  "CX3CR1", "CLEC4C", "NRP2", "IL3RA", 
                                               "CD1C", "CD8A", "CD8B", "GNLY", "CD79A","PPBP"))

#CD 4 T Cells
VlnPlot(object = immune.combined, features= c("IL7R", "CD40LG", "CD4", "LTA", "IFNG", "CD40",
                                              "LCK", "ITK", "LRRN3", "RRM2", "KIAA0101", "CST7",
                                              "IL2RB", "CD96", "CDCA7", "PTPRC"))

#MK
VlnPlot(object = immune.combined, features = c("PPBP", "PF4", "FLI1", "MYH9"))

#NK Cells? - ImmGen report + YN lyons et al 
VlnPlot(object = immune.combined, features= c("GNLY", "NKG7","PRF1", "CD8B", "PRF1", "CD8B", "GZMH", 
                                              "FCGR3A", "KLRB1", "GZMA","GZMB", "KLRC1", "TBX21"))

#CD8 T cells - may need to change
VlnPlot(object = immune.combined, features.plot= c("CD8A", "IL2RB",  "CCL5",  "KLRG1", "GZMH", 
                                                   "TBX21", "PRF1", "GNLY", "ADRB2", "CST7", "CD8B", "FCGR3A"))


#Effector memory vs CTL CD8 T Cells - Hidalgo et al 2008
VlnPlot(object = immune.combined, features.plot= c("GNLY", "GZMB", "GZMA", "CD160", "EOMES", 
                                                   "PRF1", "IL2RB", "CD8A", "RRM2", "LCK",
                                                   "NCR3","TBX21"))

#Tregs - may need to change
VlnPlot(object = immune.combined, features.plot= c("SHMT2", "CALM2",  "ITGAE", "TNFRSF4", "HLA-DMA", "HLA-DPA1", "HLA-DPB1", "HLA-DQB1",
                                                   "ENTPD1", "HLA-DRA", "GBP2", "GBP5", "CITED2", "HLA-DRB5", "GIMAP1",
                                                   "LAIR2", "IL2RA", "FOXP3", "CD5"))

#Th1 cells - Hamalainen et al 2001 + rndsystems t cell subsets image
VlnPlot(object= immune.combined, features.plot = c("CD4", "CCL3", "IL12RB2", "IL7R", "CCR2", "NFKBIA", 
                                                   "TIMP1", "JUN", "IRF1", "CASP1", "CLU", "CCL5", "CCL4", "IFNG", "TNF",
                                                   "CCR1", "CCR5", "CXCR3", "IFNGR2", "IL12RB2", "IL27RA",
                                                   "STAT1", "STAT4"))

#intermediate monocytes - may need to change
VlnPlot(object = immune.combined, features.plot= c("TIMP1", "PLAC8", "CSK",  "MGLL", "FGD2",
                                                   "NKG7", "SNX5", "MARCKSL1", "COTL1"))

#nonclassical monocytes
VlnPlot(object = immune.combined, features.plot= c("INSIG1", "RHOC", "EVL", "ABCC3", "IFITM1", 
                                                   "HSPB1", "LTB", "ABI3", "IFITM3", "FCGR3A"))

#monocyte derived dendritic cells 
VlnPlot(object = immune.combined, features.plot= c("CCND2", "DUSP5", "CD1A", "PRKACB", "CD1C", "TRIB2", 
                                                   "CD36", "TLR4", "TLR7", "CX3CR1", "FCGR3A"))

#CD4 CTL 
VlnPlot(object = immune.combined, features.plot= c("GNLY", "GZMB", "GZMA", "NKG7", "LCK", "ITK", "LRRN3", 
                                                   "RRM2", "PRF1", "IL2RB", "KLRK1", "CDCA7", "KIAA0101"))

#dcs
VlnPlot(object= immune.combined, features.plot = c("FCER1A", "CST3", "CLEC7A",  "NRP2", "IL3RA", "CD209", "CLEC4C",
                                                   "CD1C"))


i.c <- RenameIdents(immune.combined, `0` = "CD14+ Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T", 
                    `3` = "CD8+ T", `4` = "CD16+ Mono", `5` = "Natural killer", `6` = "CD8+ T", `7` = "CD14+ Mono", `8` = "CD8+ T", `9` = "B", 
                    `10` = "CD14+ Mono", `11` = "Dendritic", `12` = "Megakaryocyte", `13` = "Megakaryocyte", `14`="CD4+ T/ CD8+ T", 
                    `15`="Megakaryocyte",`16`="B cell subtype", `17`="MT", `18`="Natural killer", `19`="Megakaryocyte", `20`="B cell subtype")

DimPlot(i.c, label = TRUE)


FeaturePlot(i.c, features = c("NLRP3", "PYCARD", "CASP1","GSDMD","IL1B"), split.by = "stim", max.cutoff = 3, 
            cols = c("grey", "red"))


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



#DotPlot
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



i.c$celltype.stim <- paste(Idents(i.c), i.c$stim, sep = "_")
i.c$celltype <- Idents(i.c)
Idents(i.c) <- "celltype.stim"
b.interferon.response <- FindMarkers(i.c, ident.1 = "B_STIM", ident.2 = "B_CTRL", verbose = FALSE)
head(b.interferon.response, n = 20)

MT.interferon.response <- FindMarkers(i.c, ident.1 = "MT_STIM", ident.2 = "MT_CTRL", verbose = FALSE)
head(MT.interferon.response, n = 20)

cd14mono.interferon.response <- FindMarkers(i.c, ident.1 = "CD14+ Mono_STIM", ident.2 = "CD14+ Mono_CTRL", verbose = FALSE)
head(cd14mono.interferon.response, n = 20)

cd4naiveT.interferon.response <- FindMarkers(i.c, ident.1 = "CD4 Naive T_STIM", ident.2 = "CD4 Naive T_CTRL", verbose = FALSE)
head(cd4naiveT.interferon.response, n = 20)

cd4memoryT.interferon.response <- FindMarkers(i.c, ident.1 = "CD4 Memory T_STIM", ident.2 = "CD4 Memory T_CTRL", verbose = FALSE)
head(cd4memoryT.interferon.response, n = 20)

cd8T.interferon.response <- FindMarkers(i.c, ident.1 = "CD8+ T_STIM", ident.2 = "CD8+ T_CTRL", verbose = FALSE)
head(cd8T.interferon.response, n = 20)

cd16mono.interferon.response <- FindMarkers(i.c, ident.1 = "CD16+ Mono_STIM", ident.2 = "CD16+ Mono_CTRL", verbose = FALSE)
head(cd16mono.interferon.response, n = 20)

nk.interferon.response <- FindMarkers(i.c, ident.1 = "Natural killer_STIM", ident.2 = "Natural killer_CTRL", verbose = FALSE)
head(nk.interferon.response, n = 20)

dendritic.interferon.response <- FindMarkers(i.c, ident.1 = "Dendritic_STIM", ident.2 = "Dendritic_CTRL", verbose = FALSE)
head(dendritic.interferon.response, n = 20)

mega.interferon.response <- FindMarkers(i.c, ident.1 = "Megakaryocyte_STIM", ident.2 = "Megakaryocyte_CTRL", verbose = FALSE)
head(mega.interferon.response, n = 20)

cd4Torcd8T.interferon.response <- FindMarkers(i.c, ident.1 = "CD4+ T/ CD8+ T_STIM", ident.2 = "CD4+ T/ CD8+ T_CTRL", verbose = FALSE)
head(cd4Torcd8T.interferon.response, n = 20)

bsubtype.interferon.response <- FindMarkers(i.c, ident.1 = "B cell subtype_STIM", ident.2 = "B cell subtype_CTRL", verbose = FALSE)
head(bsubtype.interferon.response, n = 20)

 #"CD14+ Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T", 
#"CD8+ T", `4` = "CD16+ Mono", `5` = "Natural killer", `6` = "CD8+ T", `7` = "CD14+ Mono", `8` = "CD8+ T", `9` = "B", 
#"CD14+ Mono", `11` = "Dendritic", `12` = "Megakaryocyte", `13` = "Megakaryocyte", `14`="CD4+ T/ CD8+ T", 
#"Megakaryocyte",`16`="B cell subtype", `17`="MT", `18`="Natural killer", `19`="Megakaryocyte", `20`="B cell subtype")

FeaturePlot(i.c, features = c("RPS4Y1","RPS4X", "HLA-DRB5","XIST"), split.by = "stim", max.cutoff = 3, 
            cols = c("grey", "red"))