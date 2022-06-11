############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
############    --   Data to post in DICE's website    --    ############
############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############


# ---> About the script.
# Version: 0
# Subversion: 2
# Second subversion of the script to prepare all data necessary to post in DICE's website.


### -------------------------- Description -------------------------- ###


### --------------------------- Libraries --------------------------- ###
library(Seurat)
library(COTAN)
library(data.table)
library(tidyr)
library(stringr)
library(gtools)
library(ggplot2)
# source('/home/vfajardo/scripts/functions/R_handy_functions.0.4.R')


### --------------------------- Functions --------------------------- ###


### ----------------------- General Arguments ----------------------- ###
# ---> General definitions.
total.workers <- as.integer(system(command='echo $PBS_NUM_PPN', intern=TRUE))
# ---> Path definitions.
gen.path <- '/mnt/hpcscratch/vfajardo/BISB/BNFO-286'
data.path <- paste0(gen.path, '/data_to_start')
reports.path <- paste0(gen.path, '/cotan_analysis')
if(!dir.exists(reports.path)) dir.create(reports.path)
# ---> File definitions.
seurat.obj.file <- paste0(data.path, '/SeuratObj_PT-6-ALL-SF_Processed.RDS')
feature.info.file <- paste0(data.path, '/FeaturesInfo.tsv')
# To define: Gene signatures and gene markers of interest.
# ---> Check directories and files.
# if(!all(dir.exists(gen.data.path), dir.exists(gen.reports.path))) stop(paste0('Following paths must be already defined:\n', gen.data.path, '\n', gen.reports.path))
# if(!dir.exists(reports.path)) dir.create(reports.path)
essential.files <- c(seurat.obj.file)
essential.files <- essential.files[!unlist(lapply(X=essential.files, FUN=file.exists))]
if(length(essential.files) > 0) stop(paste0('Next essential files are not appropriately defined -make sure you\'ve got their names right-:\n', paste0(essential.files, collapse='\n'), '\n'))


### ------------------------- Data Loading -------------------------- ###

# Seurat obejcts.
seurat.obj <- readRDS(file=seurat.obj.file)
# Features info.
# feature.info <- read.delim(file=feature.info.file, header=FALSE, col.names=c('ensembl', 'name', 'type'), stringsAsFactors=FALSE)


### ---------------------- Data preprocessing ----------------------- ###

# ---> Process metadata first.
seurat.obj@meta.data[, 'cell.barcode'] <- Cells(seurat.obj)
seurat.obj@meta.data[, 'UMAP_1'] <- seurat.obj@reductions$umap@cell.embeddings[, 'UMAP_1']
seurat.obj@meta.data[, 'UMAP_2'] <- seurat.obj@reductions$umap@cell.embeddings[, 'UMAP_2']

# ---> Set appropriate column names.
# tmp.cols <- meta.cols[names(meta.cols) %in% colnames(seurat.obj@meta.data)]
# seurat.obj@meta.data <- seurat.obj@meta.data[, names(tmp.cols)]
# colnames(seurat.obj@meta.data) <- tmp.cols

# ---> Get subset of cells.
tmp.file.name <- str_replace(string=seurat.obj.file, pattern='.RDS', replacement='_Whole_Sample_Subset.RDS')
# tmp.file.name <- str_replace(string=seurat.obj.file, pattern='.RDS', replacement='_Sample_Subset.RDS')
if(!file.exists(tmp.file.name)){
  cells.subset <- sample(x=Cells(seurat.obj), size=65000)
  seurat.obj <- subset(x=seurat.obj, cells=cells.subset)
  saveRDS(file=tmp.file.name, object=seurat.obj)
}else{
  seurat.obj <- readRDS(file=tmp.file.name)
}

### ------------------------- Main program -------------------------- ###

# ---> Split dataset into cross-reactive and non-cross-reactive T cells.
# meta.data <- as.data.table(seurat.obj@meta.data)
Idents(seurat.obj) <- seurat.obj@meta.data$pr.tag
cr.obj <- subset(seurat.obj, idents='pR')
ncr.obj <- subset(seurat.obj, idents='non-pR')
# Get further sample from the whole dataset.
all.cells.subset <- sample(x=Cells(seurat.obj), size=30000)
sample.seurat.obj <- subset(x=seurat.obj, cells=all.cells.subset)
seurat.objs <- list(
  # 'CR'=cr.obj,
  # 'NCR'=ncr.obj,
  'ALL'=sample.seurat.obj
)
objs.of.int <- 'ALL'

# ---> Load COTAN files in case they have been previously loaded. Else, go through the process.
coexpression.path <- paste0(reports.path, '/coexpression_analysis')
cotan.file.names <- paste0(coexpression.path, '/', objs.of.int, '_COTANObject.RDS')
names(cotan.file.names) <- objs.of.int
tmp.check <- all(file.exists(cotan.file.names))
if(tmp.check){
  cotan.objs <- lapply(X=c('CR', 'NCR'), FUN=function(obj){
    readRDS(file=cotan.file.names[[obj]])
  })
}else{
  # ---> Create COTAN object.
  cotan.objs <- list()
  for(obj in objs.of.int){
    obj.title <-
    if(obj=='ALL'){
      'ALL'
    }else{
      if(obj=='CR'){
        'Cross-reactive'
      }else{
        'Non-cross-reactive'
      }
    }
    count.data <- as.data.frame(GetAssayData(object=seurat.objs[[obj]], assay='RNA', slot='counts'))
    cotan.obj <- new('scCOTAN', raw=count.data)
    cotan.obj <- initRaw(cotan.obj, GEO=obj.title, sc.method="10x", cond=obj.title)
    tmp.obj <- list(cotan.obj); names(tmp.obj) <- obj
    cotan.objs <- c(cotan.objs, tmp.obj)
  }

  # ---> Object preprocessing.

  clean.path <- paste0(reports.path, '/cleaning')
  if(!dir.exists(clean.path)) dir.create(clean.path)

  # @ Remove mitochondrial genes.
  for(obj in objs.of.int){
    genes.to.rm <- rownames(cotan.objs[[obj]]@raw[grep("^MT-", rownames(cotan.objs[[obj]]@raw)), ])
    cotan.objs[[obj]]@raw <- cotan.objs[[obj]]@raw[!rownames(cotan.objs[[obj]]@raw) %in% genes.to.rm, ]
  }

  # @ First round of cleaning.
  # Stimation with mu linear method.
  ttm.list <- list()
  for(obj in objs.of.int){
  # for(obj in c('NCR')){
    for(i in 1:5){
      cotan.obj <- cotan.objs[[obj]]
      cat(paste0('Iteration ', i, '\n'))
      ttm <- clean(cotan.obj)
      # Output PCA.
      tmp.file.name <- paste0(clean.path, '/', obj, '_PCA_Iteration-', i, '.pdf')
      pdf(file=tmp.file.name)
      print(ttm$pca.cell.2)
      dev.off()
      # Remove B cells.
      if(length(ttm$cl1) < length(ttm$cl2)){
        to.rem <- ttm$cl1
      }else{
        to.rem <- ttm$cl2
      }
      cotan.objs[[obj]]@raw <- cotan.obj@raw[, !colnames(cotan.obj@raw) %in% to.rem]
    }
    cotan.objs[[obj]] <- ttm$object
    tmp.obj <- list(ttm); names(tmp.obj) <- obj
    ttm.list <- c(ttm.list, tmp.obj)
    rm(tmp.obj)
  }

  # @ Assess cell efficiency.
  for(obj in objs.of.int){
  # for(obj in c('NCR')){
    # Output cell efficiency.
    nu.ests <- round(cotan.objs[[obj]]@nu, digits = 7)
    plot.nu <- ggplot(ttm.list[[obj]]$pca_cells, aes(x=PC1, y=PC2, color=log(nu.ests))) +
      geom_point(size=1, alpha=0.8) +
      scale_color_gradient2(low="#E64B35B2", mid="#4DBBD5B2", high="#3C5488B2", midpoint=log(mean(nu.ests)), name="ln (nu)") +
      ggtitle("Cells PCA coloured by cell efficiency") +
      theme(plot.title=element_text(color="#3C5488FF", size=20), legend.title=element_text(color="#3C5488FF", size=14, face="italic"), legend.text=element_text(color="#3C5488FF", size=11), legend.key.width=unit(2, "mm"), legend.position="right")
    tmp.file.name <- paste0(clean.path, '/', obj, '_PCAByCellEfficiency.pdf')
    pdf(file=tmp.file.name)
    print(plot.nu)
    dev.off()
    # Output cell efficiency for low-efficiency cells.
    low.ests <- nu.ests[log(nu.ests) < (-0.5)]
    plot.nu <- ggplot(ttm.list[[obj]]$pca_cells[log(nu.ests) < (-0.5), ], aes(x=PC1, y=PC2, color=log(low.ests))) +
      geom_point(size=1, alpha=0.8) +
      scale_color_gradient2(low="#E64B35B2", mid="#4DBBD5B2", high="#3C5488B2", midpoint=log(mean(nu.ests)), name="ln (nu)") +
      ggtitle("Cells PCA coloured by cells efficiency") +
      theme(plot.title=element_text(color="#3C5488FF", size=20), legend.title=element_text(color="#3C5488FF", size=14, face="italic"), legend.text=element_text(color="#3C5488FF", size=11), legend.key.width=unit(2, "mm"), legend.position="right")
    tmp.file.name <- paste0(clean.path, '/', obj, '_PCAByCellEfficiency_LowCells.pdf')
    pdf(file=tmp.file.name)
    print(plot.nu)
    dev.off()
    # Output cell efficiency for high-efficiency cells.
    high.ests <- nu.ests[log(nu.ests) >= (0.5)]
    plot.nu <- ggplot(ttm.list[[obj]]$pca_cells[log(nu.ests) >= (0.5), ], aes(x=PC1, y=PC2, color=log(high.ests))) +
      geom_point(size=1, alpha=0.8) +
      scale_color_gradient2(low="#E64B35B2", mid="#4DBBD5B2", high="#3C5488B2", midpoint=log(mean(nu.ests)), name="ln (nu)") +
      ggtitle("Cells PCA coloured by cells efficiency") +
      theme(plot.title=element_text(color="#3C5488FF", size=20), legend.title=element_text(color="#3C5488FF", size=14, face="italic"), legend.text=element_text(color="#3C5488FF", size=11), legend.key.width=unit(2, "mm"), legend.position="right")
    tmp.file.name <- paste0(clean.path, '/', obj, '_PCAByCellEfficiency_HighCells.pdf')
    pdf(file=tmp.file.name)
    print(plot.nu)
    dev.off()
  }

  # @ Save results.
  for(obj in objs.of.int){
    tmp.file.name <- paste0(clean.path, '/', obj, '_ttm_Object.RDS')
    # saveRDS(file=tmp.file.name, object=ttm)
    saveRDS(file=tmp.file.name, object=ttm.list[[obj]])
  }


  # ---> COTAN analysis.

  coexpression.path <- paste0(reports.path, '/coexpression_analysis')
  if(!dir.exists(coexpression.path)) dir.create(coexpression.path)

  # ---> Apply COTAN analysis.
  # for(obj in c('CR', 'NCR')){
  for(obj in objs.of.int){
    # @ COTAN analysis.
    cotan.objs[[obj]] <- cotan_analysis(cotan.objs[[obj]], cores=total.workers)
    # @ Co-expression analysis.
    cotan.objs[[obj]] <- get.coex(object=cotan.objs[[obj]])
    # @ Save results.
    tmp.file.name <- paste0(coexpression.path, '/', obj, '_COTANObject.RDS')
    # saveRDS(file=tmp.file.name, object=cr.cotan.obj)
    saveRDS(file=tmp.file.name, object=cotan.objs[[obj]])
    # p.vals <- get.pval(object=cr.cotan.obj)
    # @ Save P-values.
    p.vals <- get.pval(object=cotan.objs[[obj]])
    tmp.file.name <- paste0(coexpression.path, '/', obj, '_PValMatrix.RDS')
    saveRDS(file=tmp.file.name, object=p.vals)
    # @ Save Co-expression matrix.
    coex.mat <- extract.coex(cotan.objs[[obj]])
    tmp.file.name <- paste0(coexpression.path, '/', obj, '_CoexMatrix.RDS')
    saveRDS(file=tmp.file.name, object=coex.mat)
    # Obtain adjacency matrix (undirected -symmetrical- and weighted)
    coex.mat[p.vals>(0.05/ncol(p.vals))] <- 0
    coex.mat <- coex.mat[colSums(coex.mat)!=0, colSums(coex.mat)!=0]
    tmp.file.name <- paste0(coexpression.path, '/', obj, '_AdjMatrix.csv')
    # saveRDS(file=tmp.file.name, object=cr.cotan.obj)
    fwrite(file=tmp.file.name, x=coex.mat, quote=FALSE)
    # Assess edge weight distirbution.
    edge.weights <- abs(as.numeric(coex.mat))
    edge.weights <- edge.weights[edge.weights>0]
    tmp.data <- data.table(edge.weight=edge.weights)
    tmp.ggplot <- ggplot(data=tmp.data, aes(x=edge.weight)) +
      geom_density(fill='lightblue') +
      labs(x='Weight', y='Density') +
      theme_bw()
    tmp.file.name <- paste0(coexpression.path, '/', obj, '_WeightDistribution.pdf')
    pdf(file=tmp.file.name)
    print(tmp.ggplot)
    dev.off()
    sum(edge.weights>0.1)/length(edge.weights)
    # Gather all associations as an adjacency list.
    tmp.data <- as.data.table(gather(data=as.data.table(coex.mat), key='source', value='coeff'))
    tmp.data[, target:=rep(row.names(coex.mat), times=nrow(coex.mat))]
    tmp.data <- tmp.data[(coeff>0 | coeff<0) & source!=target]
    tmp.file.name <- paste0(coexpression.path, '/', obj, '_AdjList-Full.csv')
    # saveRDS(file=tmp.file.name, object=cr.cotan.obj)
    fwrite(file=tmp.file.name, x=tmp.data, quote=FALSE)
    # Filtered adjacency list.
    tmp.data <- tmp.data[coeff > 0.2]
    tmp.file.name <- paste0(coexpression.path, '/', obj, '_AdjList-Filtered.csv')
    # saveRDS(file=tmp.file.name, object=cr.cotan.obj)
    fwrite(file=tmp.file.name, x=tmp.data, quote=FALSE)
  }
}
