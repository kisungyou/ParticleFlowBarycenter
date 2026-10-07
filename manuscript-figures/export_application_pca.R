args <- commandArgs(trailingOnly = TRUE)
script <- sub("^--file=", "", commandArgs()[grepl("^--file=",commandArgs())][1])
script <- gsub("~+~", " ", script, fixed=TRUE) # Decode spaces encoded by Rscript.
here <- dirname(normalizePath(script))
project <- dirname(here)
out <- if (length(args)) args[1] else file.path(here,"data")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
for (name in c("pbmc", "fashion")) {
  datafile <- file.path(project, "data",
                       if (name == "pbmc") "processed_pbmc.RData" else "processed_fmnist.RData")
  env <- new.env()
  load(datafile, envir=env)
  X <- as.matrix(env$X_reduced)
  Y <- as.integer(env$Y)
  # Reproduce load_processed_embedding in the submitted vis-real-clustering.R.
  set.seed(100 + match(name, c("pbmc", "news", "fashion")))
  idx <- seq_len(nrow(X))
  if (length(idx) > 5000L) idx <- sort(sample(idx, 5000L, replace=FALSE))
  pc <- stats::prcomp(X[idx,,drop=FALSE], center=TRUE, scale.=FALSE, rank.=2)
  df <- data.frame(x=pc$x[,1], y=pc$x[,2], label=Y[idx], original_row=idx)
  write.csv(df, file.path(out,paste0("pca_",name,".csv")),row.names=FALSE)
  labs <- sort(unique(df$label))
  palette <- grDevices::hcl(h=seq(15,375,length.out=length(labs)+1L)[seq_along(labs)], c=100, l=65)
  write.csv(data.frame(label=labs,color=palette),file.path(out,paste0("pca_",name,"_palette.csv")),row.names=FALSE)
  print(list(name=name,n=nrow(df),range=apply(df[,1:2],2,range)))
}
