#Load libraries
library(tidyverse)
library(ggrepel)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(DESeq2)

#Load in raw counts for mouse BMDMs stimulated with various cytokines.
mouse.stim <- read.csv(file = "/Users/dmitrikotov/Documents/Sequencing Data/Things to Upload to GEO/geo_submission_May17_BMDM/raw_counts.csv")
mouse.unstim <- mouse.stim %>% dplyr::select("X","But1","But2","But3")
mouse.IFNb1 <- mouse.stim %>% dplyr::select("X","Bib1","Bib2","Bib3")
mouse.TNF <- mouse.stim %>% dplyr::select("X","Btn1","Btn2","Btn3")
mouse.IFNg <- mouse.stim %>% dplyr::select("X","Big1","Big2","Big3")
mouse.TGFb <- mouse.stim %>% dplyr::select("X","Btg1","Btg2","Btg3")

mouse.counts <- left_join(mouse.unstim,mouse.IFNb1, by = "X") %>% left_join(mouse.IFNg, by = "X") %>% left_join(mouse.TNF, by = "X") %>% left_join(mouse.TGFb, by = "X")
  
#convert ENSG gene names to gene symbol
mouse.counts$gene <- mapIds(org.Mm.eg.db, mouse.counts$X, keytype="ENSEMBL", column="SYMBOL", multiVals = "first")
mouse.counts <- mouse.counts[!is.na(mouse.counts$gene),]
mouse.counts$gene <- make.unique(mouse.counts$gene, sep="-")
row.names(mouse.counts) <- mouse.counts$gene
mouse.counts <- mouse.counts[,2:(length(mouse.counts)-1)]

#Generate metadata
samples <- colnames(mouse.counts)
coldata <- data.frame(matrix(ncol = 1, nrow = 15))
rownames(coldata) <- samples
colnames(coldata) <- "stim"
coldata$stim <- c(rep("unstim",3),rep("IFNb1",3),rep("IFNg",3),rep("TNF",3),rep("TGFb",3))

# create dds object
dds <- DESeqDataSetFromMatrix(countData = mouse.counts,
                              colData = coldata,
                              design = ~ stim)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$stim <- relevel(dds$stim, ref = "unstim")

# run DESeq
dds <- DESeq(dds)

#get average normalized counts for the various groups
mouse.normalized.counts <- counts(dds, normalized = T)
mouse.normalized.counts <- as.data.frame(mouse.normalized.counts)
mouse.normalized.counts$mean <- rowMeans(mouse.normalized.counts)
mouse.normalized.counts$Unstim_mean <- rowMeans(subset(mouse.normalized.counts, select = c(grep("ut", names(mouse.normalized.counts), value=TRUE))), na.rm = TRUE)
mouse.normalized.counts$IFNB1_mean <- rowMeans(subset(mouse.normalized.counts, select = c(grep("ib", names(mouse.normalized.counts), value=TRUE))), na.rm = TRUE)
mouse.normalized.counts$IFNG_mean <- rowMeans(subset(mouse.normalized.counts, select = c(grep("ig", names(mouse.normalized.counts), value=TRUE))), na.rm = TRUE)
mouse.normalized.counts$TGFB1_mean <- rowMeans(subset(mouse.normalized.counts, select = c(grep("tg", names(mouse.normalized.counts), value=TRUE))), na.rm = TRUE)
mouse.normalized.counts$TNF_mean <- rowMeans(subset(mouse.normalized.counts, select = c(grep("tn", names(mouse.normalized.counts), value=TRUE))), na.rm = TRUE)
mouse.normalized.counts$max <- apply(subset(mouse.normalized.counts, select = c("IFNB1_mean","IFNG_mean","TGFB1_mean","TNF_mean")), 1, max)
maxn <- function(n) function(x) order(x, decreasing = TRUE)[n]
mouse.normalized.counts$second_max <- apply(subset(mouse.normalized.counts, select = c("IFNB1_mean","IFNG_mean","TGFB1_mean","TNF_mean")), 1, function(x) x[maxn(2)(x)])
mouse.normalized.counts$mstr <- mouse.normalized.counts$max/mouse.normalized.counts$second_max
mouse.normalized.counts$gene <- row.names(mouse.normalized.counts)
mouse.simp.counts <- subset(mouse.normalized.counts, select = c("IFNB1_mean","IFNG_mean","TGFB1_mean","TNF_mean","mstr","gene"))

# aggregate results
resultsNames(dds) 
m.b <- lfcShrink(dds, coef = 2, type = "apeglm")
m.b_df <- as.data.frame(m.b@listData)
m.b_df$condition <- "IFNB1"
m.b_df$gene <- rownames(m.b_df)
m.b_df <- inner_join(m.b_df, mouse.simp.counts, by = "gene")

m.g <- lfcShrink(dds, coef = 3, type = "apeglm")
m.g_df <- as.data.frame(m.g@listData)
m.g_df$condition <- "IFNG"
m.g_df$gene <- rownames(m.g_df)

m.tgf <- lfcShrink(dds, coef = 4, type = "apeglm")
m.tgf_df <- as.data.frame(m.tgf@listData)
m.tgf_df$condition <- "TGFB1"
m.tgf_df$gene <- rownames(m.tgf_df)

m.tnf <- lfcShrink(dds, coef = 5, type = "apeglm")
m.tnf_df <- as.data.frame(m.tnf@listData)
m.tnf_df$condition <- "TNF"
m.tnf_df$gene <- rownames(m.tnf_df)

m.tgf.filt <- m.tgf_df %>% dplyr::filter(log2FoldChange > 1.5)
m.tnf.filt <- m.tnf_df %>% dplyr::filter(log2FoldChange > 1.5)

#Combine datasets and exclude genes induced by TNF or TGF-beta.
mouse.IFNs <- inner_join(m.b_df, m.g_df, by = "gene") %>% filter(log2FoldChange.x >1 | log2FoldChange.y > 1)
mouse.IFNs.specific <- mouse.IFNs[!mouse.IFNs$gene %in% m.tgf.filt$gene,]
mouse.IFNs.specific <- mouse.IFNs.specific[!mouse.IFNs.specific$gene %in% m.tnf.filt$gene,]
mouse.IFNs.specific$FCratio <- mouse.IFNs.specific$log2FoldChange.y / mouse.IFNs.specific$log2FoldChange.x

#Filter data to find genes primarily induced by type I (ifnb) or II (ifng) interferons.
specific.ifnb.mouse <- mouse.IFNs.specific %>% filter(FCratio < 0.66  & log2FoldChange.x > 4 & IFNB1_mean > 4000 & mstr > 3)
specific.ifng.mouse <- mouse.IFNs.specific %>% filter(FCratio > 1.5 & log2FoldChange.y > 2 & IFNG_mean > 500 & mstr > 3)