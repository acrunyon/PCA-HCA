# PCA on WICA scatterplots
#setwd(WD_plots)
library(ggplot2)
library(plyr)
library(lubridate)
library(dplyr)
library(SPEI)
library(ggrepel)
library(gridExtra)
library(grid)
library(Hmisc)
library(cowplot)

rm(list=ls())
reportwd<-("C:/Users/achildress/OneDrive - DOI/Documents/RSS/Working/WICA/Report/FigUpdates/")

setwd(reportwd)
load("ScatterData.RData")

GCMs<-c("HadGEM2-CC365.rcp45","MRI-CGCM3.rcp45","IPSL-CM5A-MR.rcp85","CSIRO-Mk3-6-0.rcp45")
colors5 <-   c("#12045C","#9A9EE5","#F28995","#E10720","grey")

################### PCA test, based on Annual_Means df
head(Annual_Means)
AM_rownames <- Annual_Means[,-1]
rownames(AM_rownames) <- Annual_Means[,1]

WICA.pca <- prcomp(AM_rownames[,c(1:5)], center = TRUE,scale. = TRUE)

head(WICA.pca$rotation)
head(WICA.pca$x)
summary(WICA.pca)

str(WICA.pca)

library(ggbiplot)
ggbiplot(WICA.pca, labels=rownames(AM_rownames))
