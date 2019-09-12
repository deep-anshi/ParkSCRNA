library(dplyr)
library(Seurat)
library(hdf5r)

pd_1.data <- Read10X_h5("/Users/deepanshishokeen/desktop/SingleCellData/pd/115/filtered_feature_bc_matrix.h5")
pd_1 <- CreateSeuratObject(counts = pd_1.data, project = "PD_1")
pd_1
head(colnames(pd_1))
table(pd_1$orig.ident)

pd_2.data <- Read10X_h5("/Users/deepanshishokeen/desktop/SingleCellData/pd/117/filtered_feature_bc_matrix.h5")
pd_2 <- CreateSeuratObject(counts = pd_2.data, project = "PD_2")
pd_2
head(colnames(pd_2))
table(pd_2$orig.ident)

pd_3.data <- Read10X_h5("/Users/deepanshishokeen/desktop/SingleCellData/pd/123/filtered_feature_bc_matrix.h5")
pd_3 <- CreateSeuratObject(counts = pd_3.data, project = "PD_3")
pd_3
head(colnames(pd_3))
table(pd_3$orig.ident)

pd_4.data <- Read10X_h5("/Users/deepanshishokeen/desktop/SingleCellData/pd/126/filtered_feature_bc_matrix.h5")
pd_4 <- CreateSeuratObject(counts = pd_4.data, project = "PD_4")
pd_4
head(colnames(pd_4))
table(pd_4$orig.ident)

pd_5.data <- Read10X_h5("/Users/deepanshishokeen/desktop/SingleCellData/pd/134/filtered_feature_bc_matrix.h5")
pd_5 <- CreateSeuratObject(counts = pd_5.data, project = "PD_5")
pd_5
head(colnames(pd_5))
table(pd_5$orig.ident)


##Merging 2 seurat objects
pd.combined12 <- merge(pd_1,y=pd_2, add.cell.ids = c("PD1", "PD2"), project = "PD12")
pd.combined12
# notice the cell names now have an added identifier
head(colnames(pd.combined12))
table(pd.combined12$orig.ident)


pd.combined123 <- merge(pd.combined12,y=pd_3, add.cell.ids = c("PD12", "PD3"), project = "PD123")
pd.combined123
head(colnames(pd.combined123))
table(pd.combined123$orig.ident)

pd.combined1234 <- merge(pd.combined123,y=pd_4, add.cell.ids = c("PD123", "PD4"), project = "PD1234")
pd.combined1234
head(colnames(pd.combined1234))
table(pd.combined1234$orig.ident)

pd.combined12345 <- merge(pd.combined1234,y=pd_5, add.cell.ids = c("PD1234", "PD5"), project = "PD12345")
pd.combined12345
head(colnames(pd.combined12345))
table(pd.combined12345$orig.ident)

#####    INDENTIFYING THE PD GENES--> CELL TYPES  #####

### Initialize the Seurat object with the raw (non-normalized data).
#pbmc <- CreateSeuratObject(counts = data, project = "pbmc3k", min.cells = 3, min.features = 200)
#pbmc


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pd.combined12345[["percent.mt"]] <- PercentageFeatureSet(pd.combined12345, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pd.combined12345, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pd.combined12345, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pd.combined12345, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

pd.combined12345 <- subset(pd.combined12345, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


##normalizing the data
pd.combined12345 <- NormalizeData(pd.combined12345, normalization.method = "LogNormalize", scale.factor = 10000)
# OR pbmc <- NormalizeData(pbmc)

##Identification of highly variable features
pd.combined12345 <- FindVariableFeatures(pd.combined12345, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pd.combined12345), 10)

# plot variable features with and without labels   ------->> DIDNT RUN
plot1 <- VariableFeaturePlot(pd.combined12345)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

##Scaling the data
all.genes <- rownames(pd.combined12345)
pd.combined12345 <- ScaleData(pd.combined12345, features = all.genes)

##perform linear dimensional reduction 
pd.combined12345 <- RunPCA(pd.combined12345, features = VariableFeatures(object = pd.combined12345))

# Examine and visualize PCA results a few different ways
print(pd.combined12345[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pd.combined12345, dims = 1:2, reduction = "pca")

DimPlot(pd.combined12345, reduction = "pca")

DimHeatmap(pd.combined12345, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(pd.combined12345, dims = 1:15, cells = 500, balanced = TRUE)##---> error margins too large


##determing dimenionality of dataset

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pd.combined12345 <- JackStraw(pd.combined12345, num.replicate = 100)
pd.combined12345 <- ScoreJackStraw(pd.combined12345, dims = 1:20)

JackStrawPlot(pd.combined12345, dims = 1:15)

ElbowPlot(pd.combined12345)

##cluster the cells
pd.combined12345 <- FindNeighbors(pd.combined12345, dims = 1:10)
pd.combined12345 <- FindClusters(pd.combined12345, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pd.combined12345), 5)

##run non linear dimensional reduction(UMAP/tSNE)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pd.combined12345 <- RunUMAP(pd.combined12345, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pd.combined12345, reduction = "umap")

##Finding differentially expressed features
# find all markers of cluster 1
cluster1.markers <- FindMarkers(pd.combined12345, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pd.combined12345, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pd.combined12345, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

cluster1.markers <- FindMarkers(pd.combined12345, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(pd.combined12345, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
VlnPlot(pd.combined12345, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(pd.combined12345, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                               "CD8A"))

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pd.combined12345, features = top10$gene) + NoLegend()

##CHANGE THIS--> GET THE CELL TYPE
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pd.combined12345)
pbmc <- RenameIdents(pd.combined12345, new.cluster.ids)
DimPlot(pd.combined12345, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()





