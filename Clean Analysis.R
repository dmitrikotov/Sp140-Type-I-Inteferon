library(Seurat)
library(tidyverse)
library(patchwork)
library(nichenetr)
library(UCell)

# Load the infected cell dataset - use the 30000 expected cell HTO and ADT mapping -remap HTO for B6 and Sp140 using 5 tags.
infected.rna <- Read10X(data.dir = "/Users/DmitriKotov/Documents/Sequencing Data/scRNA-seq 10-9-20/RVDK002A/outs/filtered_feature_bc_matrix", strip.suffix = TRUE)
infected.hto <- Read10X(data.dir = "/Users/DmitriKotov/Documents/Sequencing Data/scRNA-seq 10-9-20/ADT and HTO for 30000 cells/HTO for RVDK002G/umi_count", gene.column = 1)
infected.adt <- Read10X(data.dir = "/Users/DmitriKotov/Documents/Sequencing Data/scRNA-seq 10-9-20/ADT and HTO for 30000 cells/ADT for RVDK002D/umi_count", gene.column = 1)
B6.rna <- Read10X(data.dir = "/Users/DmitriKotov/Documents/Sequencing Data/scRNA-seq 10-9-20/RVDK002B/outs/filtered_feature_bc_matrix", strip.suffix = TRUE)
B6.hto <- Read10X(data.dir = "/Users/DmitriKotov/Documents/Sequencing Data/scRNA-seq 10-9-20/ADT and HTO for 30000 cells/Fixed_HTO for RVDK002H/umi_count", gene.column = 1)
B6.adt <- Read10X(data.dir = "/Users/DmitriKotov/Documents/Sequencing Data/scRNA-seq 10-9-20/ADT and HTO for 30000 cells/ADT for RVDK002E/umi_count", gene.column = 1)
Sp140.rna <- Read10X(data.dir = "/Users/DmitriKotov/Documents/Sequencing Data/scRNA-seq 10-9-20/RVDK002C/outs/filtered_feature_bc_matrix", strip.suffix = TRUE)
Sp140.hto <- Read10X(data.dir = "/Users/DmitriKotov/Documents/Sequencing Data/scRNA-seq 10-9-20/ADT and HTO for 30000 cells/Fixed_HTO for RVDK002I/umi_count", gene.column = 1)
Sp140.adt <- Read10X(data.dir = "/Users/DmitriKotov/Documents/Sequencing Data/scRNA-seq 10-9-20/ADT and HTO for 30000 cells/ADT for RVDK002F/umi_count", gene.column = 1)

#Keep cells with data for RNA, ADT, and HTO - infected
joint.infected = Reduce("intersect", list(colnames(infected.rna), colnames(infected.hto), colnames(infected.adt)))
infected.rna = infected.rna[, joint.infected]
infected.hto = as.matrix(infected.hto[-nrow(infected.hto), joint.infected])
infected.adt = as.matrix(infected.adt[-nrow(infected.adt), joint.infected])

#Keep cells with data for RNA, ADT, and HTO - B6
joint.B6 = Reduce("intersect", list(colnames(B6.rna), colnames(B6.hto), colnames(B6.adt)))
B6.rna = B6.rna[, joint.B6]
B6.hto = as.matrix(B6.hto[-nrow(B6.hto), joint.B6])
B6.adt = as.matrix(B6.adt[-nrow(B6.adt), joint.B6])

#Keep cells with data for RNA, ADT, and HTO - Sp140
joint.Sp140 = Reduce("intersect", list(colnames(Sp140.rna), colnames(Sp140.hto), colnames(Sp140.adt)))
Sp140.rna = Sp140.rna[, joint.Sp140]
Sp140.hto = as.matrix(Sp140.hto[-nrow(Sp140.hto), joint.Sp140])
Sp140.adt = as.matrix(Sp140.adt[-nrow(Sp140.adt), joint.Sp140])

# Confirm that the HTO have the correct names
rownames(infected.hto)
rownames(infected.hto) <- c("HTO_A", "HTO_B", "HTO_C", "HTO_D","HTO_E", "HTO_F")
rownames(B6.hto) <- c("HTO_A", "HTO_B", "HTO_D","HTO_E", "HTO_F")
rownames(Sp140.hto) <- c("HTO_A", "HTO_B", "HTO_C", "HTO_D","HTO_E")

# Confirm that the ADT have the correct names
rownames(infected.adt)
rownames(infected.adt) <- c("Ly6C", "CD44", "Siglec-F", "PD-L1","CD169", "CSF-1R","Ly6G","CD11b","CD86","MHCII","CX3CR1")
rownames(B6.adt) <- c("Ly6C", "CD44", "Siglec-F", "PD-L1","CD169", "CSF-1R","Ly6G","CD11b","CD86","MHCII","CX3CR1")
rownames(Sp140.adt) <- c("Ly6C", "CD44", "Siglec-F", "PD-L1","CD169", "CSF-1R","Ly6G","CD11b","CD86","MHCII","CX3CR1")

# Initialize the Seurat object with the raw (non-normalized data).
infected <- CreateSeuratObject(counts = infected.rna, project = "infected")
infected #check seurat object

# Add HTO data as a new assay independent from RNA
infected[["HTO"]] <- CreateAssayObject(counts = infected.hto)

# Add ADT data as a new assay independent from RNA
infected[["ADT"]] <- CreateAssayObject(counts = infected.adt)

#Create Seurat object for B6 and Sp140
B6 <- CreateSeuratObject(counts = B6.rna, project = "B6")
B6[["HTO"]] <- CreateAssayObject(counts = B6.hto)
B6[["ADT"]] <- CreateAssayObject(counts = B6.adt)
Sp140 <- CreateSeuratObject(counts = Sp140.rna, project = "Sp140")
Sp140[["HTO"]] <- CreateAssayObject(counts = Sp140.hto)
Sp140[["ADT"]] <- CreateAssayObject(counts = Sp140.adt)

#check for dying cells based on mitochondrial percentage
infected[["percent.mt"]] <- PercentageFeatureSet(infected, pattern = "^mt-")
B6[["percent.mt"]] <- PercentageFeatureSet(B6, pattern = "^mt-")
Sp140[["percent.mt"]] <- PercentageFeatureSet(Sp140, pattern = "^mt-")

#standard filters are 
VlnPlot(infected, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# standard log-normalization
infected <- NormalizeData(infected)
B6 <- NormalizeData(B6)
Sp140 <- NormalizeData(Sp140)

# choose ~1k variable features
infected <- FindVariableFeatures(infected)
B6 <- FindVariableFeatures(B6)
Sp140 <- FindVariableFeatures(Sp140)

# standard scaling (no regression)
infected <- ScaleData(infected)
B6 <- ScaleData(B6)
Sp140 <- ScaleData(Sp140)

#HTO demultiplexing

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
infected <- NormalizeData(infected, assay = "HTO", normalization.method = "CLR")
B6 <- NormalizeData(B6, assay = "HTO", normalization.method = "CLR")
Sp140 <- NormalizeData(Sp140, assay = "HTO", normalization.method = "CLR")

#Run HTO demux
infected <- HTODemux(infected, assay = "HTO", positive.quantile = 0.99)
B6 <- HTODemux(B6, assay = "HTO", positive.quantile = 0.99)
Sp140 <- HTODemux(Sp140, assay = "HTO", positive.quantile = 0.99)

# Global classification results
table(infected$HTO_classification.global)

# Group cells based on the max HTO signal
Idents(infected) <- "HTO_maxID"

#Visualize multiplets
RidgePlot(infected, assay = "HTO", features = rownames(infected[["HTO"]])[1:2], ncol = 2)
FeatureScatter(infected, feature1 = "hto_HTO-A", feature2 = "hto_HTO-B")
Idents(infected) <- "HTO_classification.global"
VlnPlot(infected, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

# First, we will remove negative cells from the object
infected.subset <- subset(infected, idents = "Negative", invert = TRUE)
B6.subset <- subset(B6, idents = "Negative", invert = TRUE)
Sp140.subset <- subset(Sp140, idents = "Negative", invert = TRUE)

# Calculate a tSNE embedding of the HTO data
DefaultAssay(infected.subset) <- "HTO"
infected.subset <- ScaleData(infected.subset, features = rownames(infected.subset), 
                             verbose = FALSE)
infected.subset <- RunPCA(infected.subset, features = rownames(infected.subset), approx = FALSE)
infected.subset <- RunTSNE(infected.subset, dims = 1:6, perplexity = 100, check_duplicates = FALSE) #check why there are duplicates
DimPlot(infected.subset)

# Extract the singlets
infected.singlet <- subset(infected, idents = "Singlet")
B6.singlet <- subset(B6, idents = "Singlet")
Sp140.singlet <- subset(Sp140, idents = "Singlet")

#Filter to exclude cells that have unique feature counts more than 4500 (multiplets) and less than 200 (empty drops) and removing cells that have >5% mitochondrial counts (dead/dying cells)
#Macrophages were heavily removed when removing cells with feature counts over 2500, hence why I raised it to 4500.
infected.singlet <- subset(infected.singlet, subset = nFeature_RNA < 4500 & nFeature_RNA > 200 & percent.mt < 5)
B6.singlet <- subset(B6.singlet, subset = nFeature_RNA < 4500 & nFeature_RNA > 200 & percent.mt < 5)
Sp140.singlet <- subset(Sp140.singlet, subset = nFeature_RNA < 4500 & nFeature_RNA > 200 & percent.mt < 5)

#for dividing the infected into B6 or Sp140
infected.Sp140 <- subset(infected.singlet, subset = HTO_maxID == "HTO-A" | HTO_maxID == "HTO-B" | HTO_maxID == "HTO-C")
infected.B6 <- subset(infected.singlet, subset = HTO_maxID == "HTO-D" | HTO_maxID == "HTO-E" | HTO_maxID == "HTO-F")

#Add metadata about genotype and infection status (Mtb+ = cell was infected; Mtb- = cell was not infected but came from an infected mouse)
infected.Sp140[["state"]] <- "Mtb+"
infected.B6[["state"]] <- "Mtb+"
infected.Sp140[["cell_type"]] <- "Sp140"
infected.B6[["cell_type"]] <- "B6"

B6.naive <- subset(B6.singlet, subset = HTO_maxID == "HTO-A" | HTO_maxID == "HTO-B")
B6.immune <- subset(B6.singlet, subset = HTO_maxID == "HTO-D" | HTO_maxID == "HTO-E" | HTO_maxID == "HTO-F")
B6.naive[["state"]] <- "naive"
B6.immune[["state"]] <- "Mtb-"
B6.naive[["cell_type"]] <- "B6"
B6.immune[["cell_type"]] <- "B6"

Sp140.naive <- subset(Sp140.singlet, subset = HTO_maxID == "HTO-D" | HTO_maxID == "HTO-E")
Sp140.immune <- subset(Sp140.singlet, subset = HTO_maxID == "HTO-A" | HTO_maxID == "HTO-B" | HTO_maxID == "HTO-C")
Sp140.naive[["state"]] <- "naive"
Sp140.immune[["state"]] <- "Mtb-"
Sp140.naive[["cell_type"]] <- "Sp140"
Sp140.immune[["cell_type"]] <- "Sp140"

#continue with integrated analysis here: #what is appropriate dims for integration?
object.list <-c(infected.Sp140, infected.B6, B6.naive, B6.immune, Sp140.naive, Sp140.immune)
immune.anchors <- FindIntegrationAnchors(object.list = object.list, dims = 1:30) #resulted in cells with duplicat naming. Look into subsetting by hash.ID instead of HTO_maxID
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:30)
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)

# t-SNE and Clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.8)

#Visualize integrated dataset
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "state")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
p3 <- DimPlot(immune.combined, reduction = "umap", group.by = "cell_type")
plot_grid(p1, p2, p3)

DimPlot(immune.combined, reduction = "umap", split.by = "cell_type")
DimPlot(immune.combined, reduction = "umap", split.by = "state")

#Find markers for each cluster
DefaultAssay(immune.combined) <- "RNA"
markers.1 <- FindConservedMarkers(immune.combined, ident.1 = 1, grouping.var = "state", verbose = FALSE)
markers.1b <- FindConservedMarkers(immune.combined, ident.1 = 1, grouping.var = "cell_type", verbose = FALSE)
head(markers.1)

FeaturePlot(immune.combined, features = c("adt_Ly6G", "adt_Ly6C", "adt_CX3CR1", "adt_CSF-1R", "adt_MHCII", "adt_Siglec-F",
                                          "Ly6g", "Ly6c1", "Cx3cr1","Csf1r","H2-Ab1","Siglecf"), min.cutoff = "q05", max.cutoff = "q95", ncol = 6)

#Improve clustering through Weighted nearest neighbor analysis

# Process ADT data and  set a dimensional reduction name to avoid overwriting the RNA PCA
DefaultAssay(immune.combined) <- "ADT"
VariableFeatures(immune.combined) <- rownames(immune.combined[["ADT"]])
immune.combined <- NormalizeData(immune.combined, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')

immune.combined <- FindMultiModalNeighbors(
  immune.combined, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:10), modality.weight.name = "RNA.weight"
)

immune.combined <- RunUMAP(immune.combined, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
immune.combined <- FindClusters(immune.combined, graph.name = "wsnn", algorithm = 3, resolution = 1.5, verbose = FALSE)

plot <- DimPlot(immune.combined, reduction = "wnn.umap", label = TRUE)

#Name the various clusters
immune.combined <- RenameIdents(immune.combined, `0` = "Nos2 Neut", `1` = "Act Neut", `2` = "Neut", 
                                `3` = "Neut", `4` = "AM", `5` = "Aged Neut", `6` = "Aged Neut", `7` = "Mature Neut", `8` = "Mature Neut", `9` = "IM", 
                                `10` = "Mono", `11` = "ISG IM", `12` = "Inflam Mono", `13` = "ISG Neut",`14` = "Siglec F Neut",`15` = "CD16-2 Mono",
                                `16` = "CD63 Neut",`17` = "CD63 Neut",`18` = "Mmp8 Neut",`19` = "Int Mono",`20` = "Lyve1 IM",`21` = "Int Mono",`22` = "Stfa Neut",
                                `23` = "DC",`24` = "B cell",`25` = "Basophil",`26` = "Cycling AM",`27`="Aged Neut")

#Re-order clusters to group related cell types
levels(x = immune.combined) <- c("Basophil","B cell","DC","AM","Cycling AM","Mono","Int Mono","Inflam Mono","CD16-2 Mono","Lyve1 IM","IM","IFN IM",
                                 "Neut","Siglec F Neut","Mmp8 Neut","Stfa Neut","Mature Neut","Act Neut","Aged Neut","ISG Neut","Nos2 Neut","CD63 Neut")

plot <- DimPlot(immune.combined, reduction = "wnn.umap", label = TRUE, repel = TRUE)

##Optional - Select cDC1
#select.cells <- CellSelector(plot = plot)
#Idents(immune.combined, cells = select.cells) <- "cDC1"

#add metadata to immune.combined to specify infection status and cluster in one value.
immune.combined$celltype.state <- paste(Idents(immune.combined), immune.combined$state, sep = "_")
immune.combined$celltype <- Idents(immune.combined)

#Visualize RNA weights for each cluster
VlnPlot(immune.combined, features = "RNA.weight", sort = TRUE, pt.size = 0.1) +
  NoLegend()

#cell cycle score
#function to convert human gene lists to mouse

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

m.s.genes <-convert_human_to_mouse_symbols(s.genes)
m.g2m.genes <- convert_human_to_mouse_symbols(g2m.genes)

immune.combined <- CellCycleScoring(immune.combined, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)
# view cell cycle scores and phase assignments
head(immune.combined[[]])
#visualize different cell cycle phases
DimPlot(immune.combined, reduction = "wnn.umap")
#switch ident back to clusters
Idents(immune.combined) <- "old.ident"

#Visualize contribution of RNA and ADT to clustering
immune.combined <- RunUMAP(immune.combined, reduction = 'pca', dims = 1:30, assay = 'RNA', 
                           reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
immune.combined <- RunUMAP(immune.combined, reduction = 'apca', dims = 1:10, assay = 'ADT', 
                           reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

p2 <- DimPlot(immune.combined, reduction = 'rna.umap', label = TRUE, 
              repel = TRUE, label.size = 2.5) + NoLegend()
p3 <- DimPlot(immune.combined, reduction = 'adt.umap', label = TRUE, 
              repel = TRUE, label.size = 2.5) + NoLegend()
p2 + p3 + p1

#Add scores for IFNg,Type 1 Inteferon sig.
cytokine <- read.csv(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/B6 vs Sp140 scRNA-seq 10-9-20/exvivo_markers.csv")
cytokine.ifn <- na.omit(cytokine) %>% filter(avgExpr > 2) %>% filter(gsi > 1.5) #gsi > 1.3 for TGFb and TNF?
cytokine.ifnb <- filter(cytokine.ifn, group == "IFNB")  
cytokine.ifng <- filter(cytokine.ifn, group == "IFNG")  
cytokine.tnf <- filter(cytokine.ifn, group == "TNFA")  
cytokine.tgfb <- filter(cytokine.ifn, group == "TGFB") 
cytokine.null <- filter(cytokine.ifn, group == "null")

cyto.sig <- list()
cyto.sig$ifnb <- setdiff(cytokine.ifnb$mouse_symbol, cytokine.ifng$mouse_symbol) %>% setdiff(cytokine.tnf$mouse_symbol) %>% setdiff(cytokine.tgfb$mouse_symbol) %>% setdiff(cytokine.null$mouse_symbol)
cyto.sig$ifng <- setdiff(cytokine.ifng$mouse_symbol, cytokine.ifnb$mouse_symbol) %>% setdiff(cytokine.tnf$mouse_symbol) %>% setdiff(cytokine.tgfb$mouse_symbol) %>% setdiff(cytokine.null$mouse_symbol)

immune.combined <- AddModuleScore_UCell(immune.combined, features = cyto.sig)

Idents(immune.combined) <- "celltype"
VlnPlot(immune.combined, features = c("ifnb_UCell","ifng_UCell"), split.by = "cell_type", ncol = 2)
FeaturePlot(immune.combined, reduction = "wnn.umap", features = c("ifnb_UCell","ifng_UCell"), split.by = "cell_type")

#Label cells based on cytokine signature expression
ifnb.pos <- subset(immune.combined, ifnb_UCell > 0.35)
DimPlot(ifnb.pos, reduction = "wnn.umap")
ifng.pos <- subset(immune.combined, ifng_UCell > 0.4)
DimPlot(ifng.pos, reduction = "wnn.umap")

immune.combined$ifn_sig <- ifelse(immune.combined$ifnb_UCell >= 0.35 & immune.combined$ifng_UCell >= 0.4, "both", 
                                  ifelse(immune.combined$ifnb_UCell >= 0.35 & immune.combined$ifng_UCell < 0.4, "ifnb", 
                                         ifelse(immune.combined$ifnb_UCell < 0.35 & immune.combined$ifng_UCell >= 0.4, "ifng","neither")))

#Add Mouse Identity Data for infected mice
immune.combined$cell_type.state.hash <- paste(immune.combined$cell_type, immune.combined$state, immune.combined$hash.ID, sep= "_")
combos <- unique(immune.combined$cell_type.state.hash)
super.cells <- lapply(X = combos, FUN = function(x){
  Idents(immune.combined) <- "cell_type.state.hash"
  WhichCells(immune.combined, idents = x)
})
names(super.cells) <- combos

immune.combined$mouse_id <- ifelse(colnames(immune.combined) %in% c(super.cells$`B6_Mtb+_HTO-E`,super.cells$`B6_Mtb-_HTO-E`,
                                                                    super.cells$`Sp140_Mtb+_HTO-A`,super.cells$`Sp140_Mtb-_HTO-A`), 1,
                                   ifelse(colnames(immune.combined) %in% c(super.cells$`B6_Mtb+_HTO-F`,super.cells$`B6_Mtb-_HTO-F`,
                                                                           super.cells$`Sp140_Mtb+_HTO-B`,super.cells$`Sp140_Mtb-_HTO-B`), 2,
                                          ifelse(colnames(immune.combined) %in% c(super.cells$`B6_Mtb+_HTO-D`,super.cells$`B6_Mtb-_HTO-D`,
                                                                                  super.cells$`Sp140_Mtb+_HTO-C`,super.cells$`Sp140_Mtb-_HTO-C`), 3,
                                                 ifelse(colnames(immune.combined) %in% c(super.cells$`B6_naive_HTO-A`, super.cells$`Sp140_naive_HTO-D`), 4,
                                                        ifelse(colnames(immune.combined) %in% c(super.cells$`B6_naive_HTO-B`, super.cells$`Sp140_naive_HTO-E`), 5, NA)))))

#Save latest version of immune combined
saveRDS(immune.combined, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/B6 vs Sp140 scRNA-seq 10-9-20/immunecombined4")

