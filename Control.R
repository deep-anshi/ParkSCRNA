library(dplyr)
library(Seurat)
library(hdf5r)


cn_1.data <- Read10X_h5("/Users/deepanshishokeen/desktop/QBS195/SingleCellData/control/116/filtered_feature_bc_matrix.h5")
cn_1 <- CreateSeuratObject(counts = cn_1.data, project = "CN_1")
cn_1
head(colnames(cn_1))
table(cn_1$orig.ident)

cn_2.data <- Read10X_h5("/Users/deepanshishokeen/desktop/QBS195/SingleCellData/control/124/filtered_feature_bc_matrix.h5")
cn_2 <- CreateSeuratObject(counts = cn_2.data, project = "CN_2")
cn_2
head(colnames(cn_2))
table(cn_2$orig.ident)

cn_3.data <- Read10X_h5("/Users/deepanshishokeen/desktop/QBS195/SingleCellData/control/125/filtered_feature_bc_matrix.h5")
cn_3 <- CreateSeuratObject(counts = cn_3.data, project = "CN_3")
cn_3
head(colnames(cn_3))
table(cn_3$orig.ident)

cn_4.data <- Read10X_h5("/Users/deepanshishokeen/desktop/QBS195/SingleCellData/control/133/filtered_feature_bc_matrix.h5")
cn_4 <- CreateSeuratObject(counts = cn_4.data, project = "CN_4")
cn_4
head(colnames(cn_4))
table(cn_4$orig.ident)

cn_5.data <- Read10X_h5("/Users/deepanshishokeen/desktop/QBS195/SingleCellData/control/135/filtered_feature_bc_matrix.h5")
cn_5 <- CreateSeuratObject(counts = cn_5.data, project = "CN_5")
cn_5
head(colnames(cn_5))
table(cn_5$orig.ident)


##Merging 2 seurat objects
cn.combined12 <- merge(cn_1,y=cn_2, add.cell.ids = c("cn1", "cn2"), project = "cn12")
cn.combined12
# notice the cell names now have an added identifier
head(colnames(cn.combined12))
table(cn.combined12$orig.ident)

cn.combined123 <- merge(cn.combined12,y=cn_3, add.cell.ids = c("cn12", "cn3"), project = "cn123")
cn.combined123
head(colnames(cn.combined123))
table(cn.combined123$orig.ident)

cn.combined1234 <- merge(cn.combined123,y=cn_4, add.cell.ids = c("cn123", "cn4"), project = "cn1234")
cn.combined1234
head(colnames(cn.combined1234))
table(cn.combined1234$orig.ident)

cn.combined12345 <- merge(cn.combined1234,y=cn_5, add.cell.ids = c("cn1234", "cn5"), project = "cn12345")
cn.combined12345
head(colnames(cn.combined12345))
table(cn.combined12345$orig.ident)

#####    INDENTIFYING THE cn GENES--> CELL TYPES  #####
### Initialize the Seurat object with the raw (non-normalized data).
#pbmc <- CreateSeuratObject(counts = data, project = "pbmc3k", min.cells = 3, min.features = 200)
#pbmc

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
cn.combined12345[["percent.mt"]] <- PercentageFeatureSet(cn.combined12345, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(cn.combined12345, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(cn.combined12345, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(cn.combined12345, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

cn.combined12345 <- subset(cn.combined12345, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


##normalizing the data
cn.combined12345 <- NormalizeData(cn.combined12345, normalization.method = "LogNormalize", scale.factor = 10000)
# OR pbmc <- NormalizeData(pbmc)

##Identification of highly variable features
cn.combined12345 <- FindVariableFeatures(cn.combined12345, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(cn.combined12345), 10)

# plot variable features with and without labels   
plot1 <- VariableFeaturePlot(cn.combined12345)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#CombinePlots(plots = list(plot1, plot2))

##Scaling the data
all.genes <- rownames(cn.combined12345)
cn.combined12345 <- ScaleData(cn.combined12345, features = all.genes)

##perform linear dimensional reduction 
cn.combined12345 <- RunPCA(cn.combined12345, features = VariableFeatures(object = cn.combined12345))


# Examine and visualize PCA results a few different ways
print(cn.combined12345[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(cn.combined12345, dims = 1:2, reduction = "pca")

DimPlot(cn.combined12345, reduction = "pca")

DimHeatmap(cn.combined12345, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(cn.combined12345, dims = 1:15, cells = 500, balanced = TRUE)##---> error margins too large

##determing dimenionality of dataset

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
cn.combined12345 <- JackStraw(cn.combined12345, num.replicate = 100)
cn.combined12345 <- ScoreJackStraw(cn.combined12345, dims = 1:20)

JackStrawPlot(cn.combined12345, dims = 1:15)

ElbowPlot(cn.combined12345)



##cluster the cells
cn.combined12345 <- FindNeighbors(cn.combined12345, dims = 1:10)
cn.combined12345 <- FindClusters(cn.combined12345, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(cn.combined12345), 5)

##run non linear dimensional reduction(UMAP/tSNE)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')

cn.combined12345 <- RunUMAP(cn.combined12345, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters

DimPlot(cn.combined12345, reduction = "umap")

##Finding differentially expressed features
# find all markers of cluster 1

cluster1.markers <- FindMarkers(cn.combined12345, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3

cluster5.markers <- FindMarkers(cn.combined12345, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones

pbmc.markers <- FindAllMarkers(cn.combined12345, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)


cluster1.markers <- FindMarkers(cn.combined12345, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(cn.combined12345, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
VlnPlot(cn.combined12345, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(cn.combined12345, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                                           "CD8A"))

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(cn.combined12345, features = top10$gene) + NoLegend()

##CHANGE THIS--> GET THE CELL TYPE
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(cn.combined12345)
pbmc <- RenameIdents(cn.combined12345, new.cluster.ids)
DimPlot(cn.combined12345, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


