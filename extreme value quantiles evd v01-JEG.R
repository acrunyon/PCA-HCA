#  extreme value quantiles evd vxx.R
#  putzing with ways to ID residuals

# v01 - 17 Dec 2020

#  Want to see vignette for this example: 
#   https://cran.r-project.org/web/packages/evd/vignettes/Multivariate_Extremes.pdf

library(evd)   # extreme value package. 

############ Vignette Example  #######
rm(list=ls())

options(show.signif.stars=FALSE)
nn <- nrow(lossalae)
loss <- lossalae/1e+05
lts <- c(1e-04,100)
plot(loss, log="xy", xlim=lts, ylim=lts)   # plot rescaled data

ula <- apply(loss,2,rank)/(nn +1)    # rescale data to uniform
plot(ula)     

         # want points that lie above a threshold, say, > 75th percentile
         # bvtcplot = Bivariate threshold choice plot, which "assists with threshold
         # choice for bivariate data" (from R help)
k0 <- bvtcplot(loss)$k0     # ID threshold for "large" points. In this case the (ranked) position
bvtcplot(loss, spectral = TRUE)

thresh <- apply(loss, 2, sort, decreasing = TRUE)[(k0+5)/2,]   # est threshold for what is essentially
               # fitting quantile models (k0+5/2) ~ 170, so fit to top 170 (of 1500) points

mar1 <- fitted(fpot(loss[,1], thresh[1]))    # fit pareto to margins, save fitted parameters (only)
mar2 <- fitted(fpot(loss[,2], thresh[2]))    # see text to fit parametric distr
rbind(mar1,mar2)

          # plot dependence function. pot = peaks over threshold; for large values only. Joel will 
          # explain this to us. 
abvnonpar(data = loss, method = "pot", k = k0, epmar = TRUE,plot = TRUE, lty = 3)
   
    # fit parametric models
m1 <- fbvpot(loss, thresh, model = "alog", asy1 = 1)
m2 <- fbvpot(loss, thresh, model = "bilog")
m3 <- fbvpot(loss, thresh, model = "bilog", likelihood = "poisson")
round(rbind(fitted(m2), std.errors(m2)), 3)

abvnonpar(data = loss, method = "pot", k = k0, epmar = TRUE,plot = TRUE, lty = 3)
plot(m1, which = 2, add = TRUE)
plot(m2, which = 2, add = TRUE, lty = 4)
plot(m3, which = 2, add = TRUE, lty = 2)

lts <- c(1e-04, 100)    # xy limits
plot(loss, log = "xy", col = "grey", xlim = lts, ylim = lts)

        # plot quantile curves
plot(m1, which = 3, p = c(0.95,0.975,0.99), tlty = 0, add = TRUE)
abline(v=thresh[1], h=thresh[2])

