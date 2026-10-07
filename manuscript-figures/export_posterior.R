script <- sub("^--file=", "", commandArgs()[grepl("^--file=",commandArgs())][1])
script <- gsub("~+~", " ", script, fixed=TRUE) # Decode spaces encoded by Rscript.
here <- dirname(normalizePath(script))
root <- dirname(here)
args <- commandArgs(trailingOnly=TRUE)
out <- if (length(args)) args[1] else file.path(here,'data')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
load(file.path(root,'data','normal-simplified-reduced.RData'))
xlim <- c(.63,1.13); ylim <- c(.77,1.27)
kde <- function(z){k<-MASS::kde2d(z[,1],z[,2],n=180,lims=c(xlim,ylim)); d<-expand.grid(x=k$x,y=k$y);d$z<-as.vector(k$z);d}
g <- expand.grid(x=seq(xlim[1],xlim[2],length.out=180),y=seq(ylim[1],ylim[2],length.out=180))
g$z <- mvtnorm::dmvnorm(g,common_posterior_mean,common_posterior_cov)
write.csv(g,file.path(out,'posterior_analytic.csv'),row.names=FALSE)
write.csv(kde(common_full_mcmc),file.path(out,'posterior_full_mcmc.csv'),row.names=FALSE)
ys<-individual_subfits[[(5-2)*50+1]]
for(i in seq_along(ys))write.csv(kde(ys[[i]]),file.path(out,paste0('posterior_subset_',i,'.csv')),row.names=FALSE)
f<-readRDS(file.path(root,'revision-work','refresh','wasp_1_5_100_1.rds'))
write.csv(kde(f$fit$support),file.path(out,'posterior_barycenter.csv'),row.names=FALSE)
cat('Exported density grids from saved posterior draws.\n')
