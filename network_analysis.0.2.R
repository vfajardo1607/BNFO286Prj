############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
############    ---------   Network analysis    ---------    ############
############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############


# ---> About the script.
# Version: 0
# Subversion: 1
# Second subversion of the script to prepare all data necessary to post in DICE's website.


### -------------------------- Description -------------------------- ###


### --------------------------- Libraries --------------------------- ###
library(Seurat)
library(data.table)
library(stringr)
# library(gtools)
library(ggplot2)
library(ggpubr)
library(VennDiagram)
coll.version <- 0.1
coll.file <- paste0('/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/BISB/BNFO-286/course_project/jobs_scripts/functions_collection.', coll.version, '.R')
source(coll.file)
# source('/home/vfajardo/scripts/functions/R_handy_functions.0.4.R')


### --------------------------- Functions --------------------------- ###


### ----------------------- General Arguments ----------------------- ###
# ---> General definitions.
blank.complement.1 <- theme(line=element_blank(), text=element_blank(), panel.background=element_blank(), panel.border=element_blank(), legend.position='none', axis.ticks=element_line(size=1.2), axis.line=element_line(size=1.2), axis.ticks.length=unit(0.4, "cm")) # Default blank.
# ---> Path definitions.
gen.reports.path <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/BISB/BNFO-286/course_project'
gen.data.path <- paste0(gen.reports.path, '/data')
# gen.reports.path <- '/mnt/beegfs/vfajardo/R24/paper_developments/sc_eqtl_analysis/data_for_website'
# reports.path <- paste0(gen.reports.path, '/data_', Sys.Date())
reports.path <- paste0(gen.reports.path, '/network_analysis')
# ---> File definitions.
htf.target.data.file <- paste0(gen.data.path, '/tf-target-infomation.txt')
net.results.file <- paste0(gen.data.path, '/WholeNetwork_DiffusionResults_CommunityResults.csv')
seurat.obj.file <- paste0(gen.data.path, '/SeuratObj_PT-6-ALL-SF_Processed_Whole_Sample_Subset.RDS')
# Results from gene set annotation analysis.
ann.res.3.file <- paste0(gen.data.path, '/DAVIDResults_C-3.tsv')
ann.res.3.6.file <- paste0(gen.data.path, '/DAVIDResults_C-3-6.tsv')
ann.res.5.file <- paste0(gen.data.path, '/IPAResults_C-5.tsv')
# PageRank centrality.
prank.results.file <- paste0(gen.data.path, '/PageRankResults.csv')
# Cytotoxicity signature/
cyt.sig.file <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/BISB/BNFO-286/course_project/module_signatures/module_signatures_06-02-2022/cell_cytotoxicity_patil_signature.csv'
# feature.info.file <- paste0(gen.data.path, '/FeaturesInfo.tsv')
# To define: Gene signatures and gene markers of interest.
# ---> Check directories and files.
# if(!all(dir.exists(gen.data.path), dir.exists(gen.reports.path))) stop(paste0('Following paths must be already defined:\n', gen.data.path, '\n', gen.reports.path))
# if(!dir.exists(reports.path)) dir.create(reports.path)
essential.files <- c(seurat.obj.file, htf.target.data.file)
essential.files <- essential.files[!unlist(lapply(X=essential.files, FUN=file.exists))]
if(length(essential.files) > 0) stop(paste0('Next essential files are not appropriately defined -make sure you\'ve got their names right-:\n', paste0(essential.files, collapse='\n'), '\n'))


### ------------------------- Data Loading -------------------------- ###

# TF target information
htf.target.data <- fread(file=htf.target.data.file)
# Network analysis results.
net.results <- fread(file=net.results.file)
# Seurat obejcts.
seurat.obj <- readRDS(file=seurat.obj.file)
# Features info.
# feature.info <- read.delim(file=feature.info.file, header=FALSE, col.names=c('ensembl', 'name', 'type'), stringsAsFactors=FALSE)
# Results from gene set annotation analysis.
ann.res.3 <- fread(file=ann.res.3.file)
ann.res.3.6 <- fread(file=ann.res.3.6.file)
ann.res.5 <- fread(file=ann.res.5.file)
# PageRank results.
prank.results <- fread(file=prank.results.file)
# Cytotoxicity signature.
cyt.sig <- fread(file=cyt.sig.file, header=TRUE)[, feature]


### ---------------------- Data preprocessing ----------------------- ###

# List of tissues in the TF data.
orig.tissue <- unique(as.character(htf.target.data[, str_split(string=tissue, pattern=',', simplify=TRUE)]))
# Network analysis results.
tmp.lvls <- net.results[, as.character(unique(sort(community)))]
net.results$community <- factor(x=as.character(net.results$community), levels=tmp.lvls)
net.results[, jund.community:=FALSE]; net.results[community=='3', jund.community:=TRUE]


### ------------------------- Main program -------------------------- ###


# ---> Cross-reactive DEG enrichment per community.
tmp.data <- net.results[
  community %in% tmp.pops,
  .(
    percent=.SD[dea_status=='Up', .N]/.N
  ),
  by=community
]
tmp.data.2 <- as.matrix(tmp.data[, .(percent)])
row.names(tmp.data.2) <- tmp.data[, as.character(community)]
tmp.file.name <- paste0(reports.path, '/CrossReactiveDEGEnrichmentPerCommunity.pdf')
pdf(file=tmp.file.name, width=3)
pheatmap(
  mat=tmp.data.2,
  color=colorRampPalette(RColorBrewer::brewer.pal(n=7, name='YlOrBr'))(100),
  cluster_rows=FALSE, cluster_cols=FALSE
)
dev.off()


# ---> Definition of communities.

# @ Communities 3 and 6.
# Get data.
tmp.data <- ann.res.3.6[Category=='KEGG_PATHWAY']
setorderv(x=tmp.data, cols='Bonferroni')
tmp.data.2 <- head(tmp.data, 5)
# tmp.data.2 <- rbind(tmp.data.2, tmp.data[Term=='hsa04064:NF-kappa B signaling pathway'])
tmp.data.2 <- tmp.data.2[,
  .(
    term=str_replace_all(string=str_extract(string=Term, pattern=':.+$'), pattern=':', replacement=''),
    significance=-log10(Bonferroni),
    enrichment=`Fold Enrichment`
  )
]
tmp.data.2$term <- factor(x=tmp.data.2$term, levels=rev(tmp.data.2$term))
# Plot.
tmp.ggplot <- ggplot(data=tmp.data.2, aes(x=term, y=significance, fill=enrichment)) +
  geom_bar(stat='identity', width=0.7) +
  scale_fill_gradientn(colors=colorRampPalette(RColorBrewer::brewer.pal(n=7, name='YlOrBr'))(100)) +
  scale_y_continuous(expand=expansion(add=0), breaks=scales::pretty_breaks(n=3)) +
  coord_flip() +
  labs(x='Term', y='Significance (-log10(adj. P-value))', fill='Enrichment')
publish.plot(tmp.ggplot=tmp.ggplot, output.path=reports.path, file.name='Community-3-6_DAVID', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE, height=5)

# @ MAPK pathway markers.
tmp.path <- paste0(reports.path, '/mapk_pathway_markers')
if(!dir.exists(tmp.path)) dir.create(tmp.path)
# tmp.markers <- c('JUND', 'DUSP2', 'DUSP4', 'DUSP5', 'MAP2K2', 'MAP2K3', 'NFKB1', 'NR4A1', 'TGFB1', 'VEGFB', 'CD81', 'ALYREF')
tmp.markers <- c('CD81', 'ALYREF')
for(tmp.marker in tmp.markers){
  tmp.ggplot <- vln.plot(seurat.obj=seurat.obj, feature=tmp.marker, slot='cpm', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=NULL, groups.to.filter=NULL, keep=TRUE, feature.thold=NULL, vln.type='violin', color='burst.frequency', size.thold=0, file.name=NULL, adjust.val=0.8, trim.val=TRUE, plot.limits=NULL, add.points=FALSE, this.color.scale=colorRampPalette(RColorBrewer::brewer.pal(n=7, name='YlOrBr'))(100))
  publish.plot(tmp.ggplot=tmp.ggplot, output.path=tmp.path, file.name=paste0('MAPKPathway_Marker-', tmp.marker, '_Cr-vs-Others'), type='pdf', blank.comp=blank.complement.1, do.legend=TRUE)
}


# @ Community 5
# Get data.
tmp.terms <- c(
  'Th1 Pathway',
  'Granzyme B Signaling',
  'Th1 and Th2 Activation Pathway',
  'Granzyme A Signaling',
  'Natural Killer Cell Signaling'
)
all(tmp.terms %in% ann.res.5[, `Ingenuity Canonical Pathways`])
tmp.data <- ann.res.5[`Ingenuity Canonical Pathways` %in% tmp.terms]
setorderv(x=tmp.data, cols='-log(p-value)', order=-1)
tmp.data.2 <- tmp.data[,
  .(
    term=`Ingenuity Canonical Pathways`,
    significance=`-log(p-value)`,
    enrichment=`Ratio`
  )
]
tmp.data.2$term <- factor(x=tmp.data.2$term, levels=rev(tmp.data.2$term))
# Plot.
tmp.ggplot <- ggplot(data=tmp.data.2, aes(x=term, y=significance, fill=enrichment)) +
  geom_bar(stat='identity', width=0.7) +
  scale_fill_gradientn(colors=colorRampPalette(RColorBrewer::brewer.pal(n=7, name='YlOrBr'))(100)) +
  scale_y_continuous(expand=expansion(add=0), breaks=scales::pretty_breaks(n=3)) +
  coord_flip() +
  labs(x='Term', y='Significance (-log10(adj. P-value))', fill='Enrichment')
publish.plot(tmp.ggplot=tmp.ggplot, output.path=reports.path, file.name='Community-5_IPA', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE, height=5)


# @ Overlap of community 5 genes with cytotoxicity signature.
# Get data.
tmp.data <- list(
  `Established`=cyt.sig,
  `Community 5`=net.results[community==5, `shared name`]
)
fill.cols <- c('#cc0000', '#cc6600'); names(fill.cols) <- names(tmp.data)
tmp.val <- length(cyt.sig)
# Output.
tmp.file.name <- paste0(reports.path, '/Community-5_GeneSharingWSignature.C.tiff')
venn.diagram(x=tmp.data, filename=tmp.file.name, imagetype='tiff', na='remove', hyper.test=TRUE, total.population=tmp.val, fill=fill.cols, col=fill.cols)
tmp.file.name <- paste0(reports.path, '/Community-5_GeneSharingWSignature.B.tiff')
venn.diagram(x=tmp.data, filename=tmp.file.name, imagetype='tiff', fill=fill.cols, col=fill.cols, cex=0, cat.cex=0, ext.text=FALSE)
tmp.files <- list.files(path=reports.path, pattern='.log$', full.names=TRUE)
file.remove(tmp.files)

# ---> PageRank results.
prank.results[, in.cyt.sig:='Absent']
prank.results[gene.id %in% cyt.sig, in.cyt.sig:='Present']

# @ Results for up-regulated genes.
# Retrieve data.
tmp.data <- prank.results[dea_status=='Up', .(gene.id, pagerank, in.cyt.sig)][1:10, ]
tmp.data$gene.id <- factor(x=tmp.data$gene.id, levels=tmp.data$gene.id)
tmp.cols <- c(Present='#cc0000', Absent='#4d4dff')
# Output.
tmp.ggplot <- ggplot(data=tmp.data, aes(x=gene.id, y=pagerank, col=in.cyt.sig)) +
  geom_point(size=7) +
  scale_color_manual(values=tmp.cols) +
  scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
  labs(x='Gene ID', y='PageRank centrality value', col='')
publish.plot(tmp.ggplot=tmp.ggplot, output.path=reports.path, file.name='PageRank_TopUpRegulatedGenes', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE, width=14)

# @ Results for down-regulated genes.
# Retrieve data.
tmp.data <- prank.results[dea_status=='Down', .(gene.id, pagerank, in.cyt.sig)][1:5, ]
tmp.data$gene.id <- factor(x=tmp.data$gene.id, levels=tmp.data$gene.id)
tmp.cols <- c(Present='#cc0000', Absent='#4d4dff')
# Output.
tmp.ggplot <- ggplot(data=tmp.data, aes(x=gene.id, y=pagerank, col=in.cyt.sig)) +
  geom_point(size=7) +
  scale_color_manual(values=tmp.cols) +
  scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
  labs(x='Gene ID', y='PageRank centrality value', col='')
publish.plot(tmp.ggplot=tmp.ggplot, output.path=reports.path, file.name='PageRank_TopDownRegulatedGenes', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE, width=7)


# ---> Network propagation signal across clusters (JUND)
# @ Comparison among large communities.
tmp.ggplot <- ggplot(data=net.results[community %in% tmp.pops & community != 3], aes(x=community, y=diffusion_output_heat_JUND)) +
  geom_boxplot(alpha=0.6, fill='blue', width=0.4) +
  labs(x='Community', y='Diffusion heat output') +
  theme_bw()
publish.plot(tmp.ggplot=tmp.ggplot, output.path=reports.path, file.name='DiffusionHeat-JUND', type='pdf', blank.comp=blank.complement.1, do.legend=FALSE, width=10)
# @ Only for communoity 3.
tmp.ggplot <- ggplot(data=net.results[community == '3'], aes(x=community, y=diffusion_output_heat_JUND)) +
  geom_boxplot(alpha=0.6, fill='blue', width=0.4) +
  labs(x='Community', y='Diffusion heat output') +
  theme_bw()
publish.plot(tmp.ggplot=tmp.ggplot, output.path=reports.path, file.name='DiffusionHeat-JUND_Only3', type='pdf', blank.comp=blank.complement.1, do.legend=FALSE, width=2)


# ---> Network propagation signal across clusters (CD81)
# @ Comparison among large communities.
tmp.ggplot <- ggplot(data=net.results[community %in% tmp.pops & community != 3], aes(x=community, y=diffusion_output_heat_CD81)) +
  geom_boxplot(alpha=0.6, fill='blue', width=0.4) +
  labs(x='Community', y='Diffusion heat output') +
  theme_bw()
publish.plot(tmp.ggplot=tmp.ggplot, output.path=reports.path, file.name='DiffusionHeat-CD81', type='pdf', blank.comp=blank.complement.1, do.legend=FALSE, width=10)
# @ Only for communoity 3.
tmp.ggplot <- ggplot(data=net.results[community == '3'], aes(x=community, y=diffusion_output_heat_CD81)) +
  geom_boxplot(alpha=0.6, fill='blue', width=0.4) +
  labs(x='Community', y='Diffusion heat output') +
  theme_bw()
publish.plot(tmp.ggplot=tmp.ggplot, output.path=reports.path, file.name='DiffusionHeat-CD81_Only3', type='pdf', blank.comp=blank.complement.1, do.legend=FALSE, width=2)

# ---> TF targets.
jund.targets <- htf.target.data[TF=='JUND', target]
tmp.data <- net.results[
  community %in% tmp.pops,
  .(target.enrichment=sum(name %in% jund.targets)/.N),
  by=community
]
tmp.data.2 <- as.matrix(tmp.data[, .(target.enrichment)])
row.names(tmp.data.2) <- tmp.data[, as.character(community)]
tmp.file.name <- paste0(reports.path, '/JUNDTargetAnalysis.pdf')
pdf(file=tmp.file.name, width=3)
print(pheatmap(
  mat=tmp.data.2,
  color=colorRampPalette(RColorBrewer::brewer.pal(n=7, name='YlOrBr'))(100),
  cluster_rows=FALSE, cluster_cols=FALSE
))
dev.off()
