# generate data
set.seed(3)
X <- as.data.frame(matrix(rnorm(40), ncol = 2))
names(X) <- c("x","y")
X$id <- LETTERS[seq( from = 1, to = 20 )]

plot(y~x,data=X, cex = 0.5,pch=19,col="grey")
text(y~x, labels=id,data=X, cex=0.9, font=2)

### 4 CORNERS
  #min bounding bo around oints

lx = min(X$x)
ux = max(X$x)
ly = min(X$y)
uy = max(X$y)

  #convert to points
ww = c(lx,uy)
wd = c(lx,ly)
hw = c(ux,uy)
hd = c(ux,ly)

#calc Euclidian dist of each point from corners
X$WW.distance <- sqrt((X$x - ww[1])^2 + (X$y - ww[2])^2)
X$WD.distance <- sqrt((X$x - wd[1])^2 + (X$y - wd[2])^2)
X$HW.distance <- sqrt((X$x - hw[1])^2 + (X$y - hw[2])^2)
X$HD.distance <- sqrt((X$x - hd[1])^2 + (X$y - hd[2])^2)

X$id[which.min(X$WW.distance)]
X$id[which.min(X$WD.distance)]
X$id[which.min(X$HW.distance)]
X$id[which.min(X$HD.distance)]

library(ggplot2)
ggplot(X, aes(x = x, y = y)) +
  geom_point() + geom_text(aes(label=id),hjust=0, vjust=0) +
  geom_rect(color = "black", alpha=0, mapping=aes(xmin=lx, xmax=ux, ymin=ly, ymax=uy),linetype=2) +
  geom_point(aes(x=x[which.min(WW.distance)], y=y[which.min(WW.distance)]), shape=21, size=9, stroke=2, colour="slateblue1") +
  geom_point(aes(x=x[which.min(WD.distance)], y=y[which.min(WD.distance)]), shape=21, size=9, stroke=2, colour="pink") +
  geom_point(aes(x=x[which.min(HW.distance)], y=y[which.min(HW.distance)]), shape=21, size=9, stroke=2, colour="navy") +
  geom_point(aes(x=x[which.min(HD.distance)], y=y[which.min(HD.distance)]), shape=21, size=9, stroke=2, colour="red")


### JHR METHOD
# Identify extreme points in each quadrant (either convex hull or Euclid dist from mean X & mean Y)
  #Convex hull
hpts <- chull(X[1:2]) #convex hull of points
circle.hpts <- c(hpts, hpts[1]) #create 'circle'
lines(X[circle.hpts, ]) #plot points

X.hpts <- X[hpts,] #subset df to points
row.names(X.hpts) <- X.hpts$id

# Select combinations that are maximally divergent from one another - dissimilarity matrix

n.pts = 2 #identify n points to use

# dissimilarity matrix
library(cluster)
diss <- daisy(X.hpts[,1:2], metric = "euclidean", stand = FALSE) # allows weighting of variables as well
diss<-as.matrix(diss)
diss[upper.tri(diss, diag = FALSE)] <- NA 

top<-tail(sort(diss), 3)


index <- which(diss == top[1], arr.ind=TRUE)
paste(rownames(diss)[index[1]], colnames(diss)[index[2]], sep=", ")

