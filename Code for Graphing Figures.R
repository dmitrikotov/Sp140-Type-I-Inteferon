library(Seurat)
library(tidyverse)
library(patchwork)
library(EnhancedVolcano)

immune.combined <- readRDS(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/B6 vs Sp140 scRNA-seq 10-9-20/immunecombined4")

#Graphing the data for Figures.

#Fig. 2A
DimPlot(immune.combined, reduction = "wnn.umap", label = T, repel = T)

#Fig. 2B
levels(immune.combined) <- c("ISG IM","IM","Lyve1 IM","Inflam Mono","Int Mono","Mono","CD16-2 Mono","Cycling AM","AM","CD63 Neut","Nos2 Neut", "ISG Neut","Aged Neut",
                             "Act Neut","mNeut","Stfa Neut", "SigF Neut","Neut","Mmp8 Neut","Basophil","DC","B cell")
DotPlot(immune.combined, features = c("adt_Ly6G","adt_Siglec-F","Mertk","Fcgr1","Ccr2","Fcgr4","Cd19","Fcer1a","Isg15","Lyve1"))

#Fig. 2C
DimPlot(immune.combined, reduction = "wnn.umap", split.by = "cell_type")

#Fig. 2D
DimPlot(immune.combined, reduction = "wnn.umap", split.by = "state")

#Sup. Fig. 3A
#clustering based on RNA
DimPlot(immune.combined, reduction = "rna.umap")
#Clustering based on protein
DimPlot(immune.combined, reduction = "adt.umap")
#Clustering based on RNA and protein
DimPlot(immune.combined, reduction = "wnn.umap")

#Sup. Fig. 3B
VlnPlot(immune.combined, features = "RNA.weight", sort = T)

#Sup. Fig. 3C
Idents(immune.combined) <- "state"
infected <- subset(immune.combined, idents = c("Mtb+","Mtb-"))
naive <- subset(immune.combined, idents = "naive")
Idents(infected) <- "mouse_id"
#plot for infected mice
DimPlot(infected, reduction = "wnn.umap", split.by = "cell_type")
Idents(naive) <- "mouse_id"
#plot for naive mice
DimPlot(naive, reduction = "wnn.umap", split.by = "cell_type")

#Sup. Fig. 4A
Neut.plot <- subset(immune.combined, idents = c("Mmp8 Neut","Neut","CD63 Neut","ISG Neut","SigF Neut","Stfa Neut","mNeut","Act Neut","Aged Neut","Nos2 Neut"))
DotPlot(Neut.plot, features = c("adt_CD11b","adt_Ly6G","Ly6g","Cxcr2","Cxcr4","Mmp8","Sell","Cd63","Isg15","Nos2","Stfa1","adt_Siglec-F"))

#Sup. Fig. 4B
Mac.plot <- subset(immune.combined, idents = c("AM","Cycling AM","CD16-2 Mono","Mono","Int Mono","Inflam Mono","IM","Lyve1 IM","ISG IM"))
levels(Mac.plot) <- c("Lyve1 IM","ISG IM","IM","Inflam Mono","Int Mono","Mono","CD16-2 Mono","Cycling AM","AM")
DotPlot(Mac.plot, features = c("Csf1r","adt_CD11b","adt_Siglec-F","Fcgr1","Mertk","Itgax","Ccr2","adt_CX3CR1",
                               "Fcgr4","Ly6c2","Sell","adt_MHCII","Nos2","Cd63","Isg15","adt_PD-L1","adt_CD86","Mki67"))

#Sup. Fig. 5A
simplified.immune <- immune.combined
DimPlot(simplified.immune, reduction = "wnn.umap", label = TRUE)
simplified.immune <- RenameIdents(simplified.immune, 'SigF Neut' = "Neut", 'Mmp8 Neut' = "Neut", 'Stfa Neut' = "Neut",
                                  'mNeut' = "Neut",'Act Neut' = "Neut",'Aged Neut' = "Neut",'ISG Neut' = "Neut",'Nos2 Neut' = "Neut",'CD63 Neut' = "Neut", 'Int Mono' = "Mono",
                                  'Inflam Mono' = "Mono",'CD16-2 Mono' = "Mono",'Lyve1 IM' = "IM", 'ISG IM' = "IM")
simplified.immune <- RenameIdents(simplified.immune, 'Neut' = "Neutrophil")
DimPlot(simplified.immune, reduction = "wnn.umap", label = TRUE, split.by = "cell_type")
simplified.immune$old.ident.2 <- Idents(simplified.immune)
simplified.immune$celltype.state.condition <- paste(Idents(simplified.immune), simplified.immune$state, simplified.immune$cell_type, sep = "_")
simplified.immune$state.condition <- paste(simplified.immune$state, simplified.immune$cell_type, sep = "_")
Idents(simplified.immune) <- "state"
super.simp <- subset(simplified.immune, idents = "naive")
Idents(super.simp) <- "old.ident.2"
DimPlot(super.simp, reduction = "wnn.umap", split.by = "cell_type")

#Sup. Fig. 5B
Idents(super.simp) <- "celltype.state.condition"
simp.neut.markers <- FindMarkers(super.simp, ident.1 = "Neut_naive_B6", ident.2 = "Neut_naive_Sp140", logfc.threshold = 0, min.pct = 0)
simp.neut.volcano <- simp.neut.markers[c(2,5)]
EnhancedVolcano(simp.neut.volcano,
                lab = rownames(simp.neut.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Naive B6 vs Sp140 KO Neutrophils",
                subtitle = "",
                ylim = c(0,-log10(10e-250)),
                xlim = c(-4,4))

simp.mono.markers <- FindMarkers(super.simp, ident.1 = "Mono_naive_B6", ident.2 = "Mono_naive_Sp140", logfc.threshold = 0, min.pct = 0)
simp.mono.volcano <- simp.mono.markers[c(2,5)]
EnhancedVolcano(simp.mono.volcano,
                lab = rownames(simp.mono.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Naive B6 vs Sp140 KO Monocytes",
                subtitle = "",
                ylim = c(0,-log10(10e-250)),
                xlim = c(-4,4))

simp.IM.markers <- FindMarkers(super.simp, ident.1 = "IM_naive_B6", ident.2 = "IM_naive_Sp140", logfc.threshold = 0, min.pct = 0)
simp.IM.volcano <- simp.IM.markers[c(2,5)]
EnhancedVolcano(simp.IM.volcano,
                lab = rownames(simp.IM.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Naive B6 vs Sp140 KO IM",
                subtitle = "",
                ylim = c(0,-log10(10e-250)),
                xlim = c(-4,4))

simp.AM.markers <- FindMarkers(super.simp, ident.1 = "AM_naive_B6", ident.2 = "AM_naive_Sp140", logfc.threshold = 0, min.pct = 0)
simp.AM.volcano <- simp.AM.markers[c(2,5)]
EnhancedVolcano(simp.AM.volcano,
                lab = rownames(simp.AM.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Naive B6 vs Sp140 KO AM",
                subtitle = "",
                ylim = c(0,-log10(10e-250)),
                xlim = c(-4,4))

#Fig. 3A
FeaturePlot(immune.combined, reduction = "wnn.umap", features = "Ifnb1", order = T, split.by = "state")

#Fig. 3B
FeaturePlot(immune.combined, reduction = "wnn.umap", features = "Ifnb1", order = T, split.by = "cell_type")

#Fig. 3C
Ifn.pos <- subset(immune.combined, Ifnb1 > 0)
Idents(Ifn.pos) <- "cell_type"
levels(Ifn.pos) <- c("B6","Sp140")
VlnPlot(Ifn.pos, features = c("Ifnb1"), pt.size = 2) 

#Fig. 6A
infected <- subset(immune.combined, idents = c("Mtb+","Mtb-"))
FeaturePlot(infected, reduction = "wnn.umap", features = c("Ifnar1","Ifnar2"), order = T)

#Fig. 6B
Idents(simplified.immune) <- "celltype.state.condition"
super.simp2 <- subset(simplified.immune, idents = c("Mtb-","Mtb+"))

#generate volcano plot for neutrophils comparing B6 and Sp140
Idents(super.simp2) <- "cell_type"
b6.sp140.neutrophil.markers <- FindMarkers(simplified.immune, ident.1 = "Neutrophil_Mtb+_B6", ident.2 = "Neutrophil_Mtb+_Sp140")
b6.sp140.neutrophil.volcano <- b6.sp140.neutrophil.markers[c(2,5)]
EnhancedVolcano(b6.sp140.neutrophil.volcano,
                lab = rownames(b6.sp140.neutrophil.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "B6 vs Sp140 KO Neutrophils",
                subtitle = "",)

#generate volcano plot for IMs comparing B6 and Sp140
b6.sp140.IM.markers <- FindMarkers(simplified.immune, ident.1 = "IM_Mtb+_B6", ident.2 = "IM_Mtb+_Sp140")
b6.sp140.IM.volcano <- b6.sp140.IM.markers[c(2,5)]
EnhancedVolcano(b6.sp140.IM.volcano,
                lab = rownames(b6.sp140.IM.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "B6 vs Sp140 KO IM",
                subtitle = "",)

#Fig. 6C
Idents(immune.combined) <- "ifn_sig"
DimPlot(immune.combined, reduction = "wnn.umap", split.by = "cell_type",cols = c("grey","purple","red","blue"), order = c("ifnb","ifng","both","neither"))

#Fig. 6D
Idents(immune.combined) <- "state"
infected <- subset(immune.combined, idents = c("Mtb+","Mtb-"))
Idents(infected) <- "ifn_sig"
DimPlot(infected, reduction = "wnn.umap", split.by = "cell_type",cols = c("grey","purple","red","blue"), order = c("ifnb","ifng","both","neither"))

#Fig. 6E - generating data that was then plotted and statistically analyzed in Prism.

super.simp3 <- subset(super.simp2, idents = c("Neutrophil","Mono","IM","B cell","DC","AM"))
super.simp3$ifn_sig <- ifelse(super.simp3$ifnb_UCell >= 0.35 & super.simp3$ifng_UCell >= 0.4, "both", 
                              ifelse(super.simp3$ifnb_UCell >= 0.35 & super.simp3$ifng_UCell < 0.4, "ifnb", 
                                     ifelse(super.simp3$ifnb_UCell < 0.35 & super.simp3$ifng_UCell >= 0.4, "ifng","neither")))
super.simp3$ifn_sig.condition <- paste(super.simp3$ifn_sig, super.simp3$cell_type, sep = "_")
super.simp3$ifn_sig_combined <- ifelse(super.simp3$ifn_sig =="ifnb" | super.simp3$ifn_sig =="both", "ifnb",
                                       ifelse(super.simp3$ifn_sig == "ifng", "ifng", "neither"))
super.simp.neut <- subset(super.simp3, idents = "Neutrophil")
super.simp.IM <- subset(super.simp3, idents = "IM")
super.simp.mono <- subset(super.simp3, idents ="Mono")
super.simp.b <- subset(super.simp3, idents ="B cell")
super.simp.DC <- subset(super.simp3, idents = "DC")
super.simp.AM <- subset(super.simp3, idents = "AM")

Idents(super.simp.AM) <- "cell_type"
simp.AM.B6 <- subset(super.simp.AM, idents = "B6")
simp.AM.sp140 <- subset(super.simp.AM, idents = "Sp140")
Idents(super.simp.IM) <- "cell_type"
simp.IM.B6 <- subset(super.simp.IM, idents = "B6")
simp.IM.sp140 <- subset(super.simp.IM, idents = "Sp140")
Idents(super.simp.mono) <- "cell_type"
simp.mono.B6  <- subset(super.simp.mono, idents = "B6")
simp.mono.sp140 <- subset(super.simp.mono, idents = "Sp140")
Idents(super.simp.DC) <- "cell_type"
simp.DC.B6 <- subset(super.simp.DC, idents = "B6")
simp.DC.sp140 <- subset(super.simp.DC, idents = "Sp140")
Idents(super.simp.b) <- "cell_type"
simp.b.B6 <- subset(super.simp.b, idents = "B6")
simp.b.sp140 <- subset(super.simp.b, idents = "Sp140")
Idents(super.simp.neut) <- "cell_type"
simp.neut.B6 <- subset(super.simp.neut, idents = "B6")
simp.neut.sp140 <- subset(super.simp.neut, idents = "Sp140")

cyt.objects <- c(simp.AM.B6,simp.AM.sp140,simp.IM.B6,simp.IM.sp140,simp.mono.B6,simp.mono.sp140,simp.DC.B6,simp.DC.sp140,simp.b.B6,simp.b.sp140,simp.neut.B6,simp.neut.sp140)
cyt.sig.total <- lapply(cyt.objects, function(x) {(prop.table(table(Idents(x), x$ifn_sig)))})
cyt.sig.total <- lapply(cyt.sig.total, function(x) {as.data.frame.matrix(x)})
cyt.sig.total <- data.table:::rbindlist(cyt.sig.total, fill=TRUE)
cyt.sig.total$genotype <- rep(c("B6","Sp140"),6)
cyt.sig.total$celltype <- c("AM","AM","IM","IM","Monocyte","Monocyte","DC","DC","B cell","B cell","Neutrophil","Neutrophil")
write.csv(cyt.sig.total, file = "cyt.sig.total.csv")

#Sup. Fig. 9B
#Ifn-b gene signature expression
FeaturePlot(immune.combined, reduction = "wnn.umap", features = "ifnb_UCell", split.by = "state")
#Ifn-g gene signature expression
FeaturePlot(immune.combined, reduction = "wnn.umap", features = "ifng_UCell", split.by = "state")

#Sup. Fig. 9C
#Ifn-b gene signature expression
FeaturePlot(immune.combined, reduction = "wnn.umap", features = "ifnb_UCell", split.by = "cell_type")
#Ifn-g gene signature expression
FeaturePlot(immune.combined, reduction = "wnn.umap", features = "ifng_UCell", split.by = "cell_type")
