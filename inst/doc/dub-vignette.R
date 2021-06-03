## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE, warning=FALSE---------------------------------------------
if (!require(DUBStepR))
  install.packages("DUBStepR")

## ----warning=FALSE------------------------------------------------------------
library(DUBStepR)

## ----message=FALSE, warning=FALSE---------------------------------------------
# install.packages(c("Seurat", "hdf5r"), repos = "https://cloud.r-project.org")
library(Seurat)
library(dplyr)

## ----warning=FALSE------------------------------------------------------------
seuratObj <- CreateSeuratObject(counts = Read10X_h5("pbmc_1k_v2_filtered_feature_bc_matrix.h5"), assay = "RNA", project = "10k_PBMC")
seuratObj

## -----------------------------------------------------------------------------
seuratObj <- NormalizeData(object = seuratObj, normalization.method = "LogNormalize")

## ----echo=TRUE, message=TRUE, results=FALSE, warning=FALSE--------------------
dubstepR.out <- DUBStepR(input.data = seuratObj@assays$RNA@data, min.cells = 0.05*ncol(seuratObj), optimise.features = TRUE, k = 10, num.pcs = 20, error = 0)
seuratObj@assays$RNA@var.features <- dubstepR.out$optimal.feature.genes
seuratObj

## ----echo=TRUE, message=TRUE, warning=FALSE-----------------------------------
seuratObj <- ScaleData(seuratObj, features = rownames(seuratObj))
seuratObj <- RunPCA(seuratObj, features = VariableFeatures(object = seuratObj), npcs = 30)
ElbowPlot(seuratObj, ndims = 30)

## ----echo=TRUE, fig.height=15, fig.width=15, message=TRUE---------------------
seuratObj <- RunUMAP(seuratObj, dims = 1:10, n.components = 2, seed.use = 2019)
FeaturePlot(seuratObj, features = VariableFeatures(object = seuratObj)[1:9], cols = c("lightgrey", "magenta"))

## ----echo=TRUE, fig.height=15, fig.width=15, message=TRUE---------------------
FeaturePlot(seuratObj, features = c("MS4A1", "NKG7", "CD3E", "IL7R", "CD8A", "CD14", "CST3", "FCGR3A", "PPBP"))

## ----message=TRUE-------------------------------------------------------------
seuratObj <- FindNeighbors(seuratObj, reduction = "pca", dims = 1:10)
seuratObj <- FindClusters(seuratObj)
DimPlot(seuratObj, reduction = "umap", label = TRUE, pt.size = 0.5, repel = T, label.size = 5)

## ----fig.height=10, fig.width=15, message=TRUE--------------------------------
top.10.markers <- FindAllMarkers(object = seuratObj, assay = "RNA", logfc.threshold = 0.5, min.pct = 0.5, only.pos = TRUE) %>% filter(p_val_adj < 0.1) %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(object = seuratObj, features = unique(top.10.markers$gene), size = 5)

## ----echo=TRUE, fig.height=7, fig.width=10, message=TRUE----------------------
cell.types <- c("0" = "CD14+ Monocytes", "5" = "Inflammatory CD14+ Monocytes", "1" = "Naive CD4+ T cells", "3" = "Memory CD4+ T cells", "4" = "Naive CD8+ T cells", "2" = "B cells", "6" = "NK cells", "7" = "CD16+ Monocytes", "8" = "Platelets")
seuratObj <- RenameIdents(seuratObj, cell.types)
DimPlot(seuratObj, reduction = "umap", label = TRUE, pt.size = 1, repel = TRUE, label.size = 5) + NoLegend()

