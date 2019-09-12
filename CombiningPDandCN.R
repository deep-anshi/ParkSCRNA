library(Seurat)
library(cowplot)
library(dplyr)
library(hdf5r)
library(xlsx)
#### Merging PD ####
pd_1.data <- Read10X_h5("/Users/deepanshishokeen/desktop/QBS195/SingleCellData/pd/115/filtered_feature_bc_matrix.h5")
pd_1 <- CreateSeuratObject(counts = pd_1.data, project = "PD_1", min.cells = 5)
pd_2.data <- Read10X_h5("/Users/deepanshishokeen/desktop/QBS195/SingleCellData/pd/117/filtered_feature_bc_matrix.h5")
pd_2 <- CreateSeuratObject(counts = pd_2.data, project = "PD_2", min.cells = 5)
pd_3.data <- Read10X_h5("/Users/deepanshishokeen/desktop/QBS195/SingleCellData/pd/123/filtered_feature_bc_matrix.h5")
pd_3 <- CreateSeuratObject(counts = pd_3.data, project = "PD_3", min.cells = 5)
pd_4.data <- Read10X_h5("/Users/deepanshishokeen/desktop/QBS195/SingleCellData/pd/126/filtered_feature_bc_matrix.h5")
pd_4 <- CreateSeuratObject(counts = pd_4.data, project = "PD_4", min.cells = 5)
pd_5.data <- Read10X_h5("/Users/deepanshishokeen/desktop/QBS195/SingleCellData/pd/134/filtered_feature_bc_matrix.h5")
pd_5 <- CreateSeuratObject(counts = pd_5.data, project = "PD_5", min.cells = 5)

##Merging 2 seurat objects
pd.combined12 <- merge(pd_1,y=pd_2, add.cell.ids = c("PD1", "PD2"), project = "PD12")
pd.combined123 <- merge(pd.combined12,y=pd_3, add.cell.ids = c("PD12", "PD3"), project = "PD123")
pd.combined1234 <- merge(pd.combined123,y=pd_4, add.cell.ids = c("PD123", "PD4"), project = "PD1234")
pd.combined12345 <- merge(pd.combined1234,y=pd_5, add.cell.ids = c("PD1234", "PD5"), project = "PD12345")

#### Merging Control ####
cn_1.data <- Read10X_h5("/Users/deepanshishokeen/desktop/QBS195/SingleCellData/control/116/filtered_feature_bc_matrix.h5")
cn_1 <- CreateSeuratObject(counts = cn_1.data, project = "CN_1", min.cells = 5)
cn_2.data <- Read10X_h5("/Users/deepanshishokeen/desktop/QBS195/SingleCellData/control/124/filtered_feature_bc_matrix.h5")
cn_2 <- CreateSeuratObject(counts = cn_2.data, project = "CN_2", min.cells = 5)
cn_3.data <- Read10X_h5("/Users/deepanshishokeen/desktop/QBS195/SingleCellData/control/125/filtered_feature_bc_matrix.h5")
cn_3 <- CreateSeuratObject(counts = cn_3.data, project = "CN_3", min.cells = 5)
cn_4.data <- Read10X_h5("/Users/deepanshishokeen/desktop/QBS195/SingleCellData/control/133/filtered_feature_bc_matrix.h5")
cn_4 <- CreateSeuratObject(counts = cn_4.data, project = "CN_4", min.cells = 5)
cn_5.data <- Read10X_h5("/Users/deepanshishokeen/desktop/QBS195/SingleCellData/control/135/filtered_feature_bc_matrix.h5")
cn_5 <- CreateSeuratObject(counts = cn_5.data, project = "CN_5", min.cells = 5)

##Merging 2 seurat objects
cn.combined12 <- merge(cn_1,y=cn_2, add.cell.ids = c("cn1", "cn2"), project = "cn12")
cn.combined123 <- merge(cn.combined12,y=cn_3, add.cell.ids = c("cn12", "cn3"), project = "cn123")
cn.combined1234 <- merge(cn.combined123,y=cn_4, add.cell.ids = c("cn123", "cn4"), project = "cn1234")
cn.combined12345 <- merge(cn.combined1234,y=cn_5, add.cell.ids = c("cn1234", "cn5"), project = "cn12345")


# Set up control object
#ctrl <- CreateSeuratObject(counts = cn.combined12345, project = "IMMUNE_CTRL", min.cells = 5)
ctrl <- cn.combined12345
ctrl$stim <- "CTRL"
ctrl <- subset(ctrl, subset = nFeature_RNA > 500)
ctrl <- NormalizeData(ctrl, verbose = FALSE)
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000)

# Set up stimulated object
#stim <- CreateSeuratObject(counts = stim.data, project = "IMMUNE_STIM", min.cells = 5)
stim <- pd.combined12345
stim$stim <- "STIM"
stim <- subset(stim, subset = nFeature_RNA > 500)
stim <- NormalizeData(stim, verbose = FALSE)
stim <- FindVariableFeatures(stim, selection.method = "vst", nfeatures = 2000)

#performing integration
immune.anchors <- FindIntegrationAnchors(object.list = list(ctrl, stim), dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

##perform an integrated analysis
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)

#Finding the optimal number of dims / PCA's
ElbowPlot(immune.combined)

# t-SNE and Clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

DimPlot(immune.combined, reduction = "umap", split.by = "stim")

#saveRDS(immune.combined, file = "/Users/deepanshishokeen/desktop/QBS195/SingleCellData/immunecombined_integrated.rds")

##indentify conserved cell types
DefaultAssay(immune.combined) <- "RNA"

#Saving in RDS format to load object in effective way
saveRDS(immune.combined, file = "/Users/deepanshishokeen/desktop/QBS195/SingleCellData/immunecombined_rna.rds")
