library(Seurat)
library(cowplot)
library(dplyr)
library(hdf5r)
library(xlsx)


immune.combined<- readRDS("/Users/deepanshishokeen/desktop/SingleCellRNASeq/SingleCellData/withoutMT_immunecombined_rna.rds")

# Cluster Tree
pbmc33k <- BuildClusterTree(object=immune.combined, slot = "scale.data",dims=1:20)
PlotClusterTree(pbmc33k)


# Find markers/conserved markers for each Clusters
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

#CD 4 T Cells => CLUSTER 1;2
VlnPlot(object = immune.combined, features= c("IL7R", "CD40LG", "CD4", "LTA", "IFNG", "CD40",
                                              "LCK", "ITK", "LRRN3", "RRM2", "KIAA0101", "CST7",
                                              "IL2RB", "CD96", "CDCA7", "PTPRC"))

#MK
VlnPlot(object = immune.combined, features = c("PPBP", "PF4", "FLI1", "MYH9"))

#NK Cells? - ImmGen report + YN lyons et al  => CLUSTER 4
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



i.c <- RenameIdents(immune.combined, `0` = "CD14+ Monocytes", `1` = "CD4+ T", `2` = "CD4+ T", 
                    `3` = "CD8+ T", `4` = "NK", `5` = "CD8+ T", `6` = "CD16+ Monocytes", `7` = "CD14+ Monocytes", `8` = "B", `9` = "Dendritic", 
                    `10` = "NK", `11` = "CD14+ Monocytes", `12` = "Megakaryocyte", `13` = "B", `14`="Dendritic", 
                    `15`="NK",`16`="B", `17`="B")
DimPlot(i.c, label = TRUE)


