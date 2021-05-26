# PCA on WICA scatterplots
#setwd(WD_plots)
library(ggplot2)
library(plyr)
library(lubrGCMate)
library(dplyr)
library(SPEI)
library(ggrepel)
library(grGCMExtra)
library(grGCM)
library(Hmisc)
library(cowplot)
library(ggbiplot)

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


ggbiplot(WICA.pca, labels=rownames(AM_rownames))

WICA.pca.x<-as.data.frame(WICA.pca$x)

# This method doesn't return the exact models we selected,
# but we may not have selected the most divergent models
rownames(WICA.pca.x)[which.min(WICA.pca.Future_Means$PC1)]
rownames(WICA.pca.x)[which.max(WICA.pca.Future_Means$PC1)]
rownames(WICA.pca.x)[which.min(WICA.pca.Future_Means$PC2)]
rownames(WICA.pca.x)[which.max(WICA.pca.Future_Means$PC2)]


########## TEST W/ DINO DATA ###########
rm(list=ls())
DINOwd<-("C:/Users/achildress/OneDrive - DOI/Documents/RSS/Working/DINO/MACA/S/Figs MACA/")

load(paste0(DINOwd,"DINO-S_T1.RData"))
GCMs<-c("GFDL-ESM2G.rcp85","CSIRO-Mk3-6-0.rcp45","MIROC-ESM.rcp85","HadGEM2-ES365.rcp85")

# Format df
head(Future_Means)
FM<-subset(Future_Means, select = c("GCM","DeltaTavg","DeltaTmax99","DeltaGrowLen","DeltaSp.SM","DeltaCS.SM","DeltaD"))
FM<- FM[,-1]
rownames(FM) <- Future_Means[,1]
head(FM)

DINO.pca <- prcomp(FM[,c(1:6)], center = TRUE,scale. = TRUE)
summary(DINO.pca)

ggbiplot(DINO.pca, labels=rownames(FM))

DINO.pca.x<-as.data.frame(DINO.pca$x)

# This method doesn't return the exact models we selected,
# but we may not have selected the most divergent models
rownames(DINO.pca.x)[which.min(DINO.pca.x$PC1)]
rownames(DINO.pca.x)[which.max(DINO.pca.x$PC1)]
rownames(DINO.pca.x)[which.min(DINO.pca.x$PC2)]
rownames(DINO.pca.x)[which.max(DINO.pca.x$PC2)]

############## TEST VS. LASSO APPROACH #################
 

dualscatter = ggplot(Future_Means, aes(DeltaTavg, DeltaPr*365, 
                                       xmin=quantile(Future_Means$DeltaTavg,0.25), 
                                       xmax=quantile(Future_Means$DeltaTavg,0.75), 
                                       ymin=quantile(Future_Means$DeltaPr,0.25)*365, 
                                       ymax=quantile(Future_Means$DeltaPr,0.75)*365))

dualscatter  + geom_text_repel(aes(label=GCM)) +
  geom_point(colour="black",size=4) +
  theme(axis.text=element_text(size=18),
        axis.title.x=element_text(size=18,vjust=-0.2),
        axis.title.y=element_text(size=18,vjust=0.2),
        plot.title=element_text(size=18,face="bold",vjust=2,hjust=0.5),
        legend.text=element_text(size=18), legend.title=element_text(size=16)) + 
  ###
  labs(title =paste(SiteGCM," Changes in climate means in 2040 by GCM run",sep=""), 
       x = "Changes in annual average temperature (F)", # Change
       y = "Changes in annual average precipitation (in)") + #change
  scale_color_manual(name="Scenarios", values=c("black")) +
  # scale_fill_manual(name="Scenarios",values = c("black")) + 
  theme(legend.position="none") +
  geom_rect(color = "black", alpha=0) + 
  geom_hline(aes(yintercept=mean(DeltaPr*365)),linetype=2) + #change
  geom_vline(aes(xintercept=mean(DeltaTavg)),linetype=2) #change

# 4-corners
lx = min(Future_Means$DeltaTavg)
ux = max(Future_Means$DeltaTavg)
ly = min(Future_Means$DeltaPr)*365
uy = max(Future_Means$DeltaPr)*365

#convert to points
ww = c(lx,uy)
wd = c(lx,ly)
hw = c(ux,uy)
hd = c(ux,ly)

#calc EuclGCMian dist of each point from corners
Future_Means$WW.distance <- sqrt((Future_Means$DeltaTavg - ww[1])^2 + (Future_Means$DeltaPr*365 - ww[2])^2)
Future_Means$WD.distance <- sqrt((Future_Means$DeltaTavg - wd[1])^2 + (Future_Means$DeltaPr*365 - wd[2])^2)
Future_Means$HW.distance <- sqrt((Future_Means$DeltaTavg - hw[1])^2 + (Future_Means$DeltaPr*365 - hw[2])^2)
Future_Means$HD.distance <- sqrt((Future_Means$DeltaTavg - hd[1])^2 + (Future_Means$DeltaPr*365 - hd[2])^2)

Future_Means$GCM[which.min(Future_Means$WW.distance)]
Future_Means$GCM[which.min(Future_Means$WD.distance)]
Future_Means$GCM[which.min(Future_Means$HW.distance)]
Future_Means$GCM[which.min(Future_Means$HD.distance)]

# PCA
FM2<-subset(Future_Means, select = c("GCM","DeltaTavg","DeltaPr"))
FM2<- FM2[,-1]
rownames(FM2) <- Future_Means[,1]
head(FM2)

DINO.pca2 <- prcomp(FM2[,c(1:2)], center = TRUE,scale. = TRUE)
summary(DINO.pca2)

ggbiplot(DINO.pca2, labels=rownames(FM2))

DINO.pca.x2<-as.data.frame(DINO.pca2$x)

rownames(DINO.pca.x2)[which.min(DINO.pca.x2$PC1)]
rownames(DINO.pca.x2)[which.max(DINO.pca.x2$PC1)]
rownames(DINO.pca.x2)[which.min(DINO.pca.x2$PC2)]
rownames(DINO.pca.x2)[which.max(DINO.pca.x2$PC2)]



########## TEST W/ YELL DATA -- selecting 3 models ###########
rm(list=ls())
YELLwd<-("C:/Users/achildress/OneDrive - DOI/Documents/RSS/Completed/YELL/MACA/Blacktail/Figs MACA/")

load(paste0(YELLwd,"Blacktail_T1.RData"))
GCMs<-c("MRI-CGCM3.rcp45","MIROC-ESM-CHEM.rcp85","HadGEM2-CC365.rcp85")

# Format df
head(Future_Means)
FM<-subset(Future_Means, select = c("GCM","DeltaSeverity","DeltaAbove32","DeltaMaxD","DeltaApr1SWE","DeltaOverPrecip99","DeltaAbove0SWE",
                                    "DeltaDuration","DeltaBegGrow","Deltasumsnow","DeltaOverHighQ","DeltaPr"))
FM<- FM[,-1]
rownames(FM) <- Future_Means[,1]
head(FM)

YELL.pca <- prcomp(FM[,c(1:11)], center = TRUE,scale. = TRUE)
summary(YELL.pca)

ggbiplot(YELL.pca, labels=rownames(FM))

YELL.pca.x<-as.data.frame(YELL.pca$x)

rownames(YELL.pca.x)[which.min(YELL.pca.x$PC1)]
rownames(YELL.pca.x)[which.max(YELL.pca.x$PC1)]
rownames(YELL.pca.x)[which.min(YELL.pca.x$PC2)]
rownames(YELL.pca.x)[which.max(YELL.pca.x$PC2)]

#Median of Pc2?
rownames(YELL.pca.x)[with(YELL.pca.x, which.min(abs(PC2 - mean(PC2))))]
