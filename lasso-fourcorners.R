# Convex hull
X <- as.data.frame(matrix(rnorm(40), ncol = 2))
names(X) <- c("x","y")
plot(X, cex = 0.5)
hpts <- chull(X)
hpts <- c(hpts, hpts[1])
lines(X[hpts, ])


# 4 corners
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

X$id <- LETTERS[seq( from = 1, to = 20 )]
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





