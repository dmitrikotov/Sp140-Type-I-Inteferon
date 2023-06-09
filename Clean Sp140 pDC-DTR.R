library(tidyverse)
library(DESeq2)
library(mygene)
library(biomaRt)
library(stringr)
library(fuzzyjoin)

#Read in counts
pDC.counts <- read_csv("/Users/dmitrikotov/Documents/Sequencing Data/Sp140 KO x pDC-DTR bulk RNA-seq 11-30-22/30-749975653/hit-counts/raw_counts.csv")
ID.gene <- pDC.counts[c(1,26)]
rownames <- pDC.counts$ID
pDC.counts <- pDC.counts[2:25]
rownames(pDC.counts) <- rownames

#Generate metadata
samples <- colnames(pDC.counts)
coldata <- data.frame(matrix(ncol = 3, nrow = 24))
rownames(coldata) <- samples
colnames(coldata) <- c("batch","genotype","pDC.DTR")
coldata$genotype <- c(rep("Sp140",6), rep("B6",6),rep("Sp140",6), rep("B6",6))
coldata$batch <- c(rep("Exp1",12), rep("Exp2",12))
coldata$pDC.DTR <- c(rep("minus",3), rep("plus",3),rep("plus",3),rep("minus",3),rep("minus",3),rep("plus",3),rep("minus",3),rep("plus",3))
#DK22 was mis-genotyped and is Sp140 KO not B6
coldata[22,] <- c("Exp2","Sp140","plus")
coldata$geno.pDC <- paste(coldata$genotype, coldata$pDC.DTR, sep = "_")

#Analyze counts with DESeq2
dds <- DESeqDataSetFromMatrix(countData = pDC.counts, colData = coldata, design = ~batch + geno.pDC)
dds <- DESeq(dds)
res <- results(dds)

#Comapre Sp140 KO
summary(res<-results(dds, contrast=c("geno.pDC","Sp140_minus","Sp140_plus"),alpha=0.05))
res <- as.data.frame(res)
res$gene_id <- row.names(res)
res <- left_join(res, ID.gene, by = c("gene_id" = "ID"))

#compare pDC-DTR
summary(pdc.dtr<-results(dds, contrast=c("geno.pDC","B6_plus","Sp140_plus"),alpha=0.05))
pdc.dtr <- as.data.frame(pdc.dtr)
pdc.dtr$gene_id <- row.names(pdc.dtr)
pdc.dtr <- left_join(pdc.dtr, ID.gene, by = c("gene_id" = "ID"))

#compare B6 to Sp140 KO
summary(b6.sp140<-results(dds, contrast=c("geno.pDC","B6_minus","Sp140_minus"),alpha=0.05))
b6.sp140 <- as.data.frame(b6.sp140)
b6.sp140$gene_id <- row.names(b6.sp140)
b6.sp140 <- left_join(b6.sp140, ID.gene, by = c("gene_id" = "ID"))

#compare B6 pDC-DTR
summary(b6.dtr<-results(dds, contrast=c("geno.pDC","B6_minus","B6_plus"),alpha=0.05))
b6.dtr <- as.data.frame(b6.dtr)
b6.dtr$gene_id <- row.names(b6.dtr)
b6.dtr <- left_join(b6.dtr, ID.gene, by = c("gene_id" = "ID"))

#Graph plots prettier just by entering ENSMUSG id
graph <- function(x){
  d <- plotCounts(dds, gene= x, intgroup="geno.pDC", returnData=TRUE)
  ggplot(d, aes(x = geno.pDC, y = count, color = geno.pDC)) + 
    geom_boxplot(coef=0, outlier.shape = NA) +
    geom_jitter(size = 3, alpha = 0.9) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(size=11)
    ) +
    ggtitle(x) +
    xlab("")
}

#Example for plotting graph
graph("ENSMUSG00000025330")

