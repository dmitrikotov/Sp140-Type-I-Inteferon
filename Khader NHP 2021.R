library(Seurat)
library(tidyverse)

load(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Khader NHP data/GSE149758_CD45_expData.Rda")
khader.meta <- read.delim("/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Khader NHP data/GSE149758_CD45_meta_data.txt")
row.names(khader.meta) <- khader.meta$barcode

khader.cd45 <- CreateSeuratObject(counts = expData, metadata = khader.meta)
khader.cd45 <- AddMetaData(object = khader.cd45, metadata = khader.meta)
k.umap <- cbind(matrix(khader.cd45$UMAP_1),matrix(khader.cd45$UMAP_2))
colnames(k.umap) <- c("UMAP_1","UMAP_2")
row.names(k.umap) <- row.names(khader.meta)

khader.cd45[["umap"]] <- CreateDimReducObject(embeddings = k.umap, key = "UMAP_", assay = DefaultAssay(khader.cd45))
Idents(khader.cd45) <- "Cluster_new"
DimPlot(khader.cd45, label = T)

#Myeloid 6 (mainly) and pDC express IFNb1 specifically in Active Mtb lungs.
FeaturePlot(khader.cd45, features = c("IFNB1"), order = T, split.by = "Condition")
ifnb.pos.khader <- subset(khader.cd45, IFNB1 > 0)
VlnPlot(ifnb.pos.khader, features = "IFNB1", group.by = "Condition", pt.size = 3)
VlnPlot(ifnb.pos.khader, features = "IFNB1", idents = c("Myeloid 6","pDC"), pt.size = 3)

#For further defining myeloid clusters compare flow markers from this paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3959879/ to scRNA-seq data:
#This is how we defined "Myeloid 6" as the IMs and Myeloid cells 1, 3, 4, 5, 8 as AMs. Myeloid 2 is likely the monocytes.
FeaturePlot(khader.cd45, features = c("CD163","MRC1","ITGAM","ITGAX","CD14","CD68"), order = T)

pDCs <- subset(khader.cd45, idents = "pDC")
Idents(pDCs) <- "Condition"
pDC.active.control <- FindMarkers(pDCs, ident.1 = "active", ident.2 = "control")


#Fig. 3D:
FeaturePlot(khader.cd45, features = c("IFNB1"), order = T, split.by = "Condition")
