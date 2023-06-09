# Sp140 KO as a model for type I inteferon-driven TB susceptibility.
Code in this repository was used to analyze single cell RNA-seq data in the manuscript titled "Cellular sources and targets of type I interferons that drive susceptibility to tuberculosis"

File descriptions:

Clean Analysis.R - describes how the data was analyzed using Seurat. This code also includes defining and classifying cells based on the type I interferon and type II interferon gene signature.

Code for Graphing Figures.R - contains the code used to generate the figures in the manuscript.

Clean Sp140 pDC-DTR.R - script for using DESeq2 to analyze bulk RNA-seq data from Sp140 KO pDC-DTR mouse lungs, including graphing gene expression.  

Mouse Interferon Signature.R - script for using DESeq2 to analyze bulk RNA-seq data of mouse bone marrow-derived macrophages stimulated with cytokines to generate signatures of genes predominantly upregulated by type I or II inteferon.  

Khader NHP 2021.R - describes how I reanalyzed data from GSE149758 from Esaulova E., et al. 2021 and contains the code for creating the plot in Figure 3D.
