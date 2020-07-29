#!/usr/bin/env Rscript

# sample cmd
# inside downloaded dir...
# ./prune.R --files $(pwd)/data/ed8a4dbbb1554a5e3227d6dfb2368828.h5ad,$(pwd)/data/7fd04d1aba61c35843dd2eb6a19d2545.h5ad --i ed8a4,7fd04 --annotations $(pwd)/misc/annotations.txt --o $(pwd)/output  

library(docopt)
library(Seurat)
library(SeuratDisk)
library(Matrix)
library(cli)

doc <- 
  "Annotate and prune integrated HuBMAP object.

Usage:
  prune.R [--files <h5ad-files> --i <ids> --annotations <path-to-annotation-file> --o <output-directory> --verbose <verbosity> --threshold <threshold> --nsample <nsample>]

Options:
  --files <h5ad-files>                        Absolute paths to h5ad files to be integrated (unspaced, comma-separated).
  --i <ids>                                   ID names of objects in h5ad files (unspaced, comma-separated).
  --annotations <path-to-annotation-file>     Absolute path to file with cell type markers.
  --o <output-directory>                      Output directory.
  --threshold <threshold>                     Threshold for pruning [default: 0.75].
  --nsample <nsample>                         Number of subsamples for pre- and post-pruned objects [default: 3].
  --verbose <verbosity>                       Verbosity [default: TRUE]."

opts <- docopt(doc)
files <- unlist(strsplit(opts$files,","))
ids <- unlist(strsplit(opts$i,","))
annotations <- opts$annotations
out <- opts$o
v <- as.logical(opts$verbose)
threshold <- as.numeric(opts$threshold)
num.subsamples <- as.numeric(opts$nsample)

cli_h1("Loading objects")
integrated.all <- LoadH5Seurat(file = paste0(out,"/integrated.h5seurat"),verbose=v)
objs <- list()
for (i in 1:length(ids)) {
  objs[[i]] <- LoadH5Seurat(file=paste0(out,"/",ids[[i]],".h5seurat"),verbose=v)
}

############  ANNOTATE ############  
cli_h1("Adding annotations")
first = T
resolution = NULL
cluster = NULL
for (r in readLines(annotations)) { 
  if (length(r)==0) {next}
  if (any(grep(pattern = "RESOLUTION",x=r))) {
    cluster = NULL
    if (first) {
      store <- as.vector(Idents(integrated.all))
      first = F
    } else {
      store <- cbind(store,as.vector(Idents(integrated.all)))
    }
    
    resolution = unlist(strsplit(x=r, split=" "))[2]
    Idents(integrated.all) <- paste0('integrated_snn_res.',resolution)
    next
  }
  if (is.null(cluster)) {
    cluster = r 
  } else {
    names <- unlist(strsplit(x=r,split=","))
    lst <- as.list(rep(cluster,length(names)))
    names(lst) <- names
    integrated.all <- RenameIdents(integrated.all, lst)
    cluster = NULL
  }
}
store <- cbind(store,as.vector(Idents(integrated.all)))

integrated.all <- 
  AddMetaData(
    integrated.all,
    col.name = "my.labels",
    metadata = apply(X = store,
                     MARGIN = 1,
                     FUN = function(row) {
                       row <- row[is.na(as.numeric(row))] # all non-numbers
                       if (length(row) == 0) {return("Unknown")}
                       return(row[length(row)]) # prioritizes highest resolution ID (rightmost column)
                     }))
Idents(integrated.all) <- 'my.labels'
SaveH5Seurat(integrated.all, filename = paste0(out,"/integrated.h5seurat"),overwrite = T,verbose=v)

############ SUBSAMPLE AND PRUNE ############  
cli_h1("Subsampling to obtain transfer scores (pre-pruning)")
anchors.subsampling <- list()
predictions <- list()
integrated.all <- AddMetaData(integrated.all,
                              col.name='rand.sample',
                              metadata=sample(1:num.subsamples,size=length(Cells(integrated.all)),replace=T))
for (i in 1:num.subsamples) {
  sub <- subset(integrated.all,subset=(rand.sample!=i))
  query <- subset(integrated.all,subset=(rand.sample==i))
  DefaultAssay(sub) <- "integrated"
  DefaultAssay(query) <- "integrated"
  sub[['SCT']] <- NULL
  query[['SCT']] <- NULL
  anchors.subsampling[[i]] <- FindTransferAnchors(reference = sub, query = query, npcs=NULL, k.filter=NA)
  predictions[[i]] <- TransferData(anchorset = anchors.subsampling[[i]], refdata = Idents(sub))
}

# collect scores/labels from the sub-transfers
all.scores <- rep(0,length(Cells(integrated.all)))
names(all.scores) <- Cells(integrated.all)
all.labels <- rep(0,length(Cells(integrated.all)))
names(all.labels) <- Cells(integrated.all)
for (i in 1:length(predictions)) {
  sub <- subset(integrated.all,subset=(rand.sample==i))@meta.data
  cells <- rownames(sub)
  all.scores[cells] <- predictions[[i]]$prediction.score.max
  all.labels[cells] <- predictions[[i]]$predicted.id
}
integrated.all <- AddMetaData(integrated.all,col.name='subsampling.score',
                              metadata=all.scores)
integrated.all <- AddMetaData(integrated.all,col.name='subsampling.label',
                              metadata=all.labels)

cli_h1("Pruning based on transfer scores")
# prune out cells with subsampling.score < threshold (default 0.75)
bad.cells <- Cells(integrated.all)[which(integrated.all[['subsampling.score',drop=T]]< threshold)]
integrated.pruned <- subset(integrated.all,
                            cells=Cells(integrated.all)[!Cells(integrated.all) %in% bad.cells])

cli_h1("Re-processing pruned reference")
# re-process pruned
DefaultAssay(integrated.pruned) <- 'integrated'
integrated.pruned <- RunPCA(integrated.pruned)
integrated.pruned <- RunUMAP(integrated.pruned,dims=1:50)

############  TRANSFER ############  
cli_h1("Subsampling to obtain transfer scores + labels (post-pruning)")
anchors.subsampling.postprune <- list()
predictions <- list()
integrated.pruned <- AddMetaData(integrated.pruned,
                                 col.name='rand.sample',
                                 metadata=sample(1:num.subsamples,size=length(Cells(integrated.pruned)),replace=T))
for (i in 1:num.subsamples) {
  sub <- subset(integrated.pruned,subset=(rand.sample!=i))
  query <- subset(integrated.pruned,subset=(rand.sample==i))
  DefaultAssay(sub) <- "integrated"
  DefaultAssay(query) <- "integrated"
  sub[['SCT']] <- NULL
  query[['SCT']] <- NULL
  anchors.subsampling.postprune[[i]] <- FindTransferAnchors(reference = sub, query = query, npcs=NULL, k.filter=NA)
  predictions[[i]] <- TransferData(anchorset = anchors.subsampling.postprune[[i]], refdata = Idents(sub))
}

# get transferred labels for the pruned out cells as well
i <- i + 1
{
  sub <- integrated.pruned
  query <- subset(integrated.all,cells=Cells(integrated.all)[which(!(Cells(integrated.all) %in% Cells(integrated.pruned)))])
  DefaultAssay(sub) <- "integrated"
  DefaultAssay(query) <- "integrated"
  sub[['SCT']] <- NULL
  query[['SCT']] <- NULL
  anchors.subsampling.postprune[[i]] <- FindTransferAnchors(reference = sub, query = query, npcs=NULL, k.filter=NA)
  predictions[[i]] <- TransferData(anchorset = anchors.subsampling.postprune[[i]], refdata = Idents(sub))
}

# make sure you have all cells
should.have <- NULL
for (i in 1:length(objs)) {
  should.have <- c(should.have,paste0(Cells(objs[[i]]),"_",i))
}
have <- NULL
for (i in 1:length(predictions)) {
  have <- c(have,rownames(predictions[[i]]))
}
stopifnot(setdiff(should.have,have)==0)
stopifnot(setdiff(have,should.have)==0)

# rbind the predictions metadata, and add to integrated.pruned object 
predictions.all <- Reduce(rbind,predictions)
names(predictions.all) <- paste0(names(predictions.all),".subsampling")
names(predictions.all)[which(names(predictions.all)=='predicted.id.subsampling')] <-
  'subsampling.label.postprune'
names(predictions.all)[which(names(predictions.all)=='prediction.score.max.subsampling')] <-
  'subsampling.score.postprune'
labs <- predictions.all[Cells(integrated.pruned),]$'subsampling.label.postprune'
scores <- predictions.all[Cells(integrated.pruned),]$'subsampling.score.postprune'
integrated.pruned <- AddMetaData(integrated.pruned,col.name='subsampling.label.postprune',metadata=labs)
integrated.pruned <- AddMetaData(integrated.pruned,col.name='subsampling.score.postprune',metadata=scores)

# also add these to the original objs
for (i in 1:length(objs)) {
  df <- predictions.all[paste0(Cells(objs[[i]]),"_",i), ]
  rownames(df) <- sapply(rownames(df),function(x) {return(Seurat:::ExtractField(x,field=1,delim="_"))})
  objs[[i]] <- AddMetaData(objs[[i]],metadata=df)
  Idents(objs[[i]]) <- 'subsampling.label.postprune'
}

cli_h1("Saving objects")
SaveH5Seurat(integrated.all, filename = paste0(out,"/integrated.h5seurat"), overwrite=T, verbose=v)
SaveH5Seurat(integrated.pruned, filename = paste0(out,"/integrated_pruned.h5seurat"), overwrite=T, verbose=v)
for (i in 1:length(ids)) {
  SaveH5Seurat(objs[[i]], filename=paste0(out,"/",ids[[i]],".h5seurat"),overwrite=T, verbose=v)
}

