# IGAN
A method for Inferring cell-cell communication pathways represented by spatial gene associations based on spatial transcriptomic data
![image](overview2.png)  

## Installation
IGAN can be installed in R by devtools:

``devtools::install_github('ZJC/IGAN')``

## Tutorials
``library(IGAN)``  
``library(clusterProfiler)``  
``library(org.Hs.eg.db)``  
``library(org.Mm.eg.db)``  
``library(ggplot2)``  
``library(ggsankey)``  
``library(networkD3)``  
``library(KEGGREST)``  
``library(tidyverse)``  
``library(plotly)``  
``library(magick)``  
``library(CellChat)``  
``library(patchwork)``  
``library(NMF)``  
``library(ggalluvial)``  
``library(ComplexHeatmap)``  

#loading input datas  
##spatial is coordinate information, a dataframe which rows represent spots and 2 cols represent coordinate.  

``spatial = read.csv('spatial.csv')``  

##ident is celltype information, a 1 column dataframe which rows represent spots.

``ident = read.csv('ident.csv')``  

##M is gene expression matrix, a dataframe which cols represent genes and rows represent spots.  

``M = read.csv('M')``

##gene is genes' name of the dataset, a 1 column dataframe which rows represent genes.  

``gene = read.csv('gene.csv')``  

##genelist1 is gene list selected in sending spots, a 1 column dataframe which rows represent genes.  

``genelist1 = read.csv('genelist1')``  

##genelist2 is gene list selected in receiving spots, a 1 column dataframe which rows represent genes.  

``genelist2 = read.csv('genelist2')``  
