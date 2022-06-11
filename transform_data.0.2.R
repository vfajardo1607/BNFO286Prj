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
library(data.table)
library(stringr)
library(gtools)
source('/home/vfajardo/scripts/functions/R_handy_functions.0.4.R')


### --------------------------- Functions --------------------------- ###


### ----------------------- General Arguments ----------------------- ###
# ---> General definitions.
meta.cols <- c(
  # Cell and QC information.
  cell.barcode='cell.barcode',
  nCount_RNA='umi.count', nFeature_RNA='gene.count', percent.mt='mitochondrial.umi.percent', percent.rb='ribosomal.umi.percent',
  # Run information.
  `RNA_snn_res.0.3`='cluster', UMAP_1='umap.1', UMAP_2='umap.2',
  # Experiment information.
  chrom.batch.tag='facs.sorting.batch', seq.batch.tag='sequencing.batch', virus.tag='virus.peptide.pool',
  # Donor information.
  donor.id.tag='donor.id', age.tag='donor.age', gender.tag='donor.sex', blood.draw.phase.tag='blood.draw.phase', hospitalization.tag='hospitalization'
)
subset.size <- 8000
# ---> Path definitions.
gen.data.path <- '/mnt/hpcscratch/vfajardo/BISB/BNFO-286/data_to_start'
# gen.reports.path <- '/mnt/beegfs/vfajardo/R24/paper_developments/sc_eqtl_analysis/data_for_website'
# reports.path <- paste0(gen.reports.path, '/data_', Sys.Date())
reports.path <- gen.data.path
# ---> File definitions.
seurat.obj.file <- paste0(gen.data.path, '/SeuratObj_PT-6-ALL-SF_Processed.RDS')
feature.info.file <- paste0(gen.data.path, '/FeaturesInfo.tsv')
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
feature.info <- read.delim(file=feature.info.file, header=FALSE, col.names=c('ensembl', 'name', 'type'), stringsAsFactors=FALSE)


### ---------------------- Data preprocessing ----------------------- ###

# ---> Process metadata first.
seurat.obj@meta.data[, 'cell.barcode'] <- Cells(seurat.obj)
seurat.obj@meta.data[, 'UMAP_1'] <- seurat.obj@reductions$umap@cell.embeddings[, 'UMAP_1']
seurat.obj@meta.data[, 'UMAP_2'] <- seurat.obj@reductions$umap@cell.embeddings[, 'UMAP_2']

# ---> Set appropriate column names.
# tmp.cols <- meta.cols[names(meta.cols) %in% colnames(seurat.obj@meta.data)]
# seurat.obj@meta.data <- seurat.obj@meta.data[, names(tmp.cols)]
# colnames(seurat.obj@meta.data) <- tmp.cols


### ------------------------- Main program -------------------------- ###

# ---------- Whole dataset.

# I tried multiple options and this one seemed to have worked best.
# ---> Get subset of cells.
tmp.file.name <- str_replace(string=seurat.obj.file, pattern='.RDS', replacement='_Whole_Sample_Subset.RDS')
if(!file.exists(tmp.file.name)){
  cells.subset <- sample(x=Cells(seurat.obj), size=65000)
  tmp.seurat.obj <- subset(x=seurat.obj, cells=cells.subset)
  saveRDS(file=tmp.file.name, object=tmp.seurat.obj)
}else{
  tmp.seurat.obj <- readRDS(file=tmp.file.name)
}

# Translate gene names into Ensembl IDs accordingly.
tmp.genes <- row.names(tmp.seurat.obj)
tmp.genes <- translate.ids(ids=tmp.genes, ensembl=FALSE)
# I wrote the code chunk below to try to get the proper Ensembl ID for all gene names under the assumption that the suffix ".\\d" was preventing a proper translation. However, it did not work well. More to research to get the translation of all genes. However, the loss of less than 0.2% of the whole set of genes must not be a problem. Further, I checked the gene names and none seems to be a receptor or a ligand gene (except, probably, for a TLR gene).
# all(str_starts(string=names(tmp.genes), pattern='ENSG'), na.rm=T)
# sum(is.na(names(tmp.genes)))
# troubling.genes <- tmp.genes[is.na(names(tmp.genes))]
# troubling.genes <- str_replace(string=troubling.genes, pattern='\\.\\d+$', replacement='')
# translate.ids(ids=troubling.genes, ensembl=FALSE)

# # ---> Gene expression UMI counts.
# Get data as it is and output info.
gex.data <- GetAssayData(object=tmp.seurat.obj, assay='RNA', slot='data')
tmp.file.name <- str_replace(string=seurat.obj.file, pattern='.RDS', replacement='_Whole_LogNormalizedGexData.csv')
write.csv(file=tmp.file.name, x=gex.data, quote=FALSE)
# Translate gene names into Ensembl IDs.
gex.data <- gex.data[!is.na(names(tmp.genes)), ]
rownames(gex.data) <- names(tmp.genes)[!is.na(names(tmp.genes))]
# Output info.
tmp.file.name <- str_replace(string=seurat.obj.file, pattern='.RDS', replacement='_Whole_LogNormalizedGexData.txt')
write.table(file=tmp.file.name, x=gex.data, quote=FALSE, sep='\t')

# # ---> Gene expression, scaled data.
# gex.data <- GetAssayData(object=tmp.seurat.obj, assay='RNA', slot='scale.data')
# gex.data <- as.matrix(gex.data)
# # Output info.
# tmp.file.name <- str_replace(string=seurat.obj.file, pattern='.RDS', replacement='_Whole_ScaledData.csv')
# write.csv(file=tmp.file.name, x=gex.data, quote=FALSE)

# ---> Metadata.
meta.data <- as.data.table(tmp.seurat.obj@meta.data)
# Output info.
tmp.file.name <- str_replace(string=seurat.obj.file, pattern='.RDS', replacement='_Whole_MetaData.csv')
fwrite(file=tmp.file.name, x=meta.data, quote=FALSE, na=NA)
# Input for CellPhoneDB.
tmp.vals <- c('0'='TFH', '1'='CD4-CTL', '2'='TH17', '3'='THIFNR', '4'='TH1', '5'='cytTFH', '6'='Proliferating', '7'='CD4-CTL')
cell.pop <- tmp.vals[meta.data[, RNA_snn_res.0.3]]
tmp.data <- data.table(Cell=Cells(tmp.seurat.obj), cell_type=cell.pop)
tmp.file.name <- str_replace(string=seurat.obj.file, pattern='.RDS', replacement='_Whole_CellPopAnns.csv')
fwrite(file=tmp.file.name, x=tmp.data, quote=FALSE, na=NA)
tmp.data <- data.table(Cell=Cells(tmp.seurat.obj), cell_type=cell.pop)
tmp.file.name <- str_replace(string=seurat.obj.file, pattern='.RDS', replacement='_Whole_CellPopAnns.txt')
fwrite(file=tmp.file.name, x=tmp.data, quote=FALSE, na=NA, sep='\t')


# ---------- Downsampled dataset.

# ---> Get subset of cells.
tmp.file.name <- str_replace(string=seurat.obj.file, pattern='.RDS', replacement='_Sample_Subset.RDS')
if(!file.exists(tmp.file.name)){
  cells.subset <- sample(x=Cells(seurat.obj), size=subset.size)
  tmp.seurat.obj <- subset(x=seurat.obj, cells=cells.subset)
  saveRDS(file=tmp.file.name, object=tmp.seurat.obj)
}else{
  tmp.seurat.obj <- readRDS(file=tmp.file.name)
}

# ---> Gene expression raw UMI counts.
gex.data <- GetAssayData(object=tmp.seurat.obj, assay='RNA', slot='counts')
gex.data <- as.matrix(gex.data)
# Output info.
tmp.file.name <- str_replace(string=seurat.obj.file, pattern='.RDS', replacement='_Sample_RawUMIGexData.csv')
write.csv(file=tmp.file.name, x=gex.data, quote=FALSE)
# I used this chunk to prove that the column names are exactly the same.
# gex.data <- read.csv(file=tmp.file.name)
# tmp.cells <- str_replace(string=colnames(gex.data), pattern='\\.', replacement='-')
# tmp.cells[!tmp.cells %in% Cells(tmp.seurat.obj)]
# Cells(tmp.seurat.obj)[!Cells(tmp.seurat.obj) %in% tmp.cells]

# ---> Gene expression UMI counts.
gex.data <- GetAssayData(object=tmp.seurat.obj, assay='RNA', slot='data')
gex.data <- as.matrix(gex.data)
# Output info.
tmp.file.name <- str_replace(string=seurat.obj.file, pattern='.RDS', replacement='_Sample_LogNormalizedGexData.csv')
write.csv(file=tmp.file.name, x=gex.data, quote=FALSE)
# Translate gene names into Ensembl IDs.
gex.data <- gex.data[!is.na(names(tmp.genes)), ]
rownames(gex.data) <- names(tmp.genes)[!is.na(names(tmp.genes))]
# Output info.
tmp.file.name <- str_replace(string=seurat.obj.file, pattern='.RDS', replacement='_Sample_LogNormalizedGexData.txt')
write.table(file=tmp.file.name, x=gex.data, quote=FALSE, sep='\t')

# ---> Gene expression, scaled data.
gex.data <- GetAssayData(object=tmp.seurat.obj, assay='RNA', slot='scale.data')
gex.data <- as.matrix(gex.data)
# Output info.
tmp.file.name <- str_replace(string=seurat.obj.file, pattern='.RDS', replacement='_Sample_ScaledData.csv')
write.csv(file=tmp.file.name, x=gex.data, quote=FALSE)

# ---> Metadata.
meta.data <- as.data.table(cbind(tmp.seurat.obj@meta.data, tmp.seurat.obj@reductions$umap@cell.embeddings))
# Output info.
tmp.file.name <- str_replace(string=seurat.obj.file, pattern='.RDS', replacement='_Sample_MetaData.csv')
fwrite(file=tmp.file.name, x=meta.data, quote=FALSE, na=NA)
# Input for CellPhoneDB.
tmp.vals <- c('0'='TFH', '1'='CD4-CTL', '2'='TH17', '3'='THIFNR', '4'='TH1', '5'='cytTFH', '6'='Proliferating', '7'='CD4-CTL')
cell.pop <- tmp.vals[meta.data[, RNA_snn_res.0.3]]
tmp.data <- data.table(Cell=Cells(tmp.seurat.obj), cell_type=cell.pop)
# tmp.file.name <- str_replace(string=seurat.obj.file, pattern='.RDS', replacement='_Whole_CellPopAnns.csv')
# fwrite(file=tmp.file.name, x=tmp.data, quote=FALSE, na=NA)
tmp.data <- data.table(Cell=Cells(tmp.seurat.obj), cell_type=cell.pop)
tmp.file.name <- str_replace(string=seurat.obj.file, pattern='.RDS', replacement='_Sample_CellPopAnns.txt')
fwrite(file=tmp.file.name, x=tmp.data, quote=FALSE, na=NA, sep='\t')
