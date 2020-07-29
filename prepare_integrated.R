#!/usr/bin/env Rscript

# sample cmd
# inside downloaded dir...
# ./prepare_integrated.R --files $(pwd)/data/ed8a4dbbb1554a5e3227d6dfb2368828.h5ad,$(pwd)/data/7fd04d1aba61c35843dd2eb6a19d2545.h5ad --i ed8a4,7fd04 --ensembl $(pwd)/misc/mart_export.txt --markers $(pwd)/misc/markers.txt --o $(pwd)/output  

library(docopt)
library(Seurat)
library(SeuratDisk)
library(Matrix)
library(cli)

doc <- 
  "Prepare integrated HuBMAP object.

Usage:
  prepare_integrated.R [--files <h5ad-files> --i <ids> --ensembl <path-to-ensembl-file> --markers <path-to-markers-file> --o <output-directory> --verbose <verbosity>]

Options:
  --files <h5ad-files>                Absolute paths to h5ad files to be integrated (unspaced, comma-separated).
  --i <ids>                           ID names of objects in h5ad files (unspaced, comma-separated).
  --ensembl <path-to-ensembl-file>    Absolute path to file with ensemble ID dataframe.
  --markers <path-to-markers-file>    Absolute path to file with cell type markers.
  --o <output-directory>              Output directory.
  --verbose <verbosity>               Verbosity [default: TRUE]."

opts <- docopt(doc)
files <- unlist(strsplit(opts$files,","))
ids <- unlist(strsplit(opts$i,","))
out <- opts$o
markers <- unlist(lapply(FUN=function(x){strsplit(x,",")}, 
                         X=readLines(opts$markers)))
v <- as.logical(opts$verbose)

ens_table <- read.table(opts$ensembl,sep = "\t",header = T,stringsAsFactors = F)
ens_table <- ens_table[!duplicated(ens_table$Gene.stable.ID),]
rownames(ens_table) <- ens_table$Gene.stable.ID

objs <- list()
mito.genes <- vector()
ribo.genes <- vector()
for(i in 1:length(ids)) {
  # read file and make seurat
  id <- ids[i]
  cli_h1(paste0("Preparing ",id))
  file.name <- Seurat:::ExtractField(files[[i]],field=1,delim=".h5ad")
  suppressMessages(suppressWarnings(
    Convert(paste0(file.name,".h5ad"), dest = "h5seurat", overwrite = TRUE,verbose=v)))
  objs[[i]] <- suppressMessages(suppressWarnings(
    LoadH5Seurat(paste0(file.name,".h5seurat"),graphs=F,verbose=v)))
  # give gene names
  ens_names <- sapply(rownames(objs[[i]]),Seurat:::ExtractField,1,"\\.")
  gene_names <- ens_table[ens_names,]$Gene.name
  gene_names[is.na(gene_names)] <- ens_names[is.na(gene_names)] # if NA, keep ensemble name
  gene_names <- gsub(pattern = '_', replacement = '-', x = gene_names)
  take <- !duplicated(gene_names)
  objs[[i]][["RNA"]] <- CreateAssayObject(data = objs[[i]][["RNA"]][take,])
  rownames(objs[[i]][['RNA']]@data) <- gene_names[take]
  # fix meta.features
  objs[[i]][['RNA']]@meta.features <- data.frame(row.names=rownames(objs[[i]]))
  # make counts
  objs[[i]][['RNA']]@counts <- Matrix(exp(objs[[i]][['RNA']]@data)-1, sparse=T)
  cnts <- objs[[i]][['RNA']]@counts
  # percent mito
  mito.genes.1 <- grep(pattern="MT-",x=rownames(objs[[i]]),value=T)
  percent.mito <- colSums(cnts[mito.genes.1,]) / colSums(cnts)
  mito.genes <- union(mito.genes, mito.genes.1)
  ribo.genes <- union(ribo.genes,
                      union(grep(pattern="RPL",x=rownames(objs[[i]]),value=T),
                            grep(pattern="RPS",x=rownames(objs[[i]]),value=T)))
  # subset cells: 200 < nFeature_RNA < 5000, percent_mito < 30%
  message(paste0("pre-filter: ",length(Cells(objs[[i]])), " cells"))
  cells <- which(objs[[i]][['nFeature_RNA']] > 200 &
                   objs[[i]][['nFeature_RNA']] < 5000 &
                   percent.mito < 0.30)
  objs[[i]] <- subset(objs[[i]], cells=cells)
  message(paste0("post-filter: ",length(Cells(objs[[i]])), " cells"))
}
for(i in 1:length(ids)) {
  # run sctransform
  cli_h1(paste0("SCTransforming ",ids[[i]]))
  objs[[i]] <- SCTransform(objs[[i]], ncells=5000, verbose=v)
  # save seurat
  suppressMessages(suppressWarnings(SaveH5Seurat(object = objs[[i]],filename = paste0(out,"/",ids[[i]],".h5seurat"),overwrite = T,verbose=v)))
}

# (save all rownames for later,) and REMOVE mito genes and ribo genes from scale data
orig.rownames <- list()
for (i in 1:length(ids)) {
  orig.rownames[[i]] <- rownames(objs[[i]][['SCT']])
  VariableFeatures(objs[[i]]) <- setdiff(setdiff(VariableFeatures(objs[[i]]), mito.genes),
                                         ribo.genes)
  objs[[i]][['SCT']]@scale.data <- objs[[i]][['SCT']]@scale.data[VariableFeatures(objs[[i]]),]
}
save(orig.rownames, file = paste0(out,"/orig.rownames.rda"))

# PROCESS objs individually
for (i in 1:length(ids)) {
  cli_h1(paste0("Processing ",ids[[i]]," (PCA, UMAP, neighbor-finding, clustering)"))
  objs[[i]] <- RunPCA(objs[[i]],verbose=v) # based on scale.data
  objs[[i]] <- RunUMAP(objs[[i]],dims=1:50,verbose=v) # based on PCA
  objs[[i]] <- FindNeighbors(objs[[i]],dims=1:50,verbose=v) # SNN graph made in PCA space
  objs[[i]] <- FindClusters(objs[[i]],resolution=1,verbose=v) # based on SNN
  suppressMessages(suppressWarnings(SaveH5Seurat(object = objs[[i]],filename = paste0(out,"/",ids[[i]],".h5seurat"),overwrite = T,verbose=v)))
}

# RELOAD objs and ADDBACK misc
objs <- list()
load(paste0(out,"/orig.rownames.rda"))
for (i in 1:length(ids)) {
  objs[[i]] <- LoadH5Seurat(file = paste0(out,"/",ids[[i]],".h5seurat"),verbose=v)
  rownames(objs[[i]][['SCT']]@misc$vst.out$model_pars_fit) <- orig.rownames[[i]]
  colnames(objs[[i]][['SCT']]@misc$vst.out$model_pars_fit) <- c("theta","(Intercept)","log_umi")
  rownames(objs[[i]][['SCT']]@misc$vst.out$cell_attr) <- Cells(objs[[i]])
}

# INTEGRATE objs and process
cli_h1("Finding anchors and integrating")
feats.all <- SelectIntegrationFeatures(objs, nfeatures = 3000,verbose=v)
objs <- PrepSCTIntegration(objs, anchor.features = feats.all,verbose=v)
anchors.all <- FindIntegrationAnchors(objs, normalization.method = "SCT", anchor.features = feats.all,verbose=v)
integrated.all <- IntegrateData(anchorset = anchors.all, normalization.method = "SCT",verbose=v)
cli_h1("Processing integrated object (PCA, UMAP, neighbor-finding, clustering)")
integrated.all <- RunPCA(integrated.all,verbose=v)
integrated.all <- RunUMAP(integrated.all,dims=1:50,verbose=v)
integrated.all <- FindNeighbors(integrated.all,dims=1:50,verbose=v)
integrated.all <- FindClusters(integrated.all,resolution=c(10,5,1.2),verbose=v)
integrated.all <- AddMetaData(integrated.all,
                              metadata = unname(sapply(FUN = Seurat:::ExtractField,
                                                       X = Cells(integrated.all),
                                                       delim = "_", field = 2)),
                              col.name = "sample")
SaveH5Seurat(integrated.all, filename = paste0(out,"/integrated.h5seurat"),overwrite = T,verbose=v)

# Feature plots and resolution plots for annotation
integrated.all <- LoadH5Seurat(file = paste0(out,"/integrated.h5seurat"))
cli_h1("Making plots")
out.plots <- paste0(out,"/plots")
system(paste0("mkdir ",out.plots))
DefaultAssay(integrated.all) <- "SCT"
for (marker in markers) {
  png(paste0(out.plots,"/feature_",marker,".png"))
  print(FeaturePlot(integrated.all,marker))
  dev.off()
}
DefaultAssay(integrated.all) <- "integrated"
for (name in grep(x=names(integrated.all[[]]), pattern='^integrated_snn_res', value=T)) {
  resolution <- Seurat:::ExtractField(name,field=2,delim="res.")
  png(paste0(out.plots,"/clusters_res.",resolution,".png"),width=1000)
  print(DimPlot(integrated.all,group.by=name,label=T))
  dev.off()
}
png(paste0(out.plots,"/BATCH.png"))
print(DimPlot(integrated.all,group.by='sample',label=T))
dev.off()

