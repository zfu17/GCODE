library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(data.table)
library(scCustomize)
library(future)
library(readxl)
library(Signac)
library(dplyr)
library(harmony)

## loading meta data
meta_series = fread('/dcs04/hongkai/data/zfu17/scMulti/ENCODE/snyder_human_scrna_scatac1.csv') %>% as.data.frame()
meta_data = read_xlsx('/dcs04/hongkai/data/zfu17/scMulti/ENCODE/snyder_encode4_human_multi_meta.xlsx') %>% as.data.frame()
good_samples = read_xlsx('/dcs04/hongkai/data/zfu17/scMulti/ENCODE/cellranger_human/snyder_multi_summary.xlsx',sheet='good_samples') %>% as.data.frame()
rownames(meta_series) = meta_series$mid
rownames(good_samples) = good_samples$`Sample ID`
filenames = good_samples$`Sample ID`
fpath='/dcs04/hongkai/data/zfu17/scMulti/ENCODE/cellranger_human/'
dhspath = '/dcs04/hongkai/data/zfu17/scMulti/ENCODE/atac_DHS_mat/'
objpath = '/dcs04/hongkai/data/zfu17/scMulti/ENCODE/seurat_objects/'
filtered_objpath = '/dcs04/hongkai/data/zfu17/scMulti/ENCODE/seurat_objects_filtered_new/'

## loading DHS
load('/dcs04/hongkai/data/zfu17/tools/ENCODE_human_mouse_DHS_noY.rda')

## get gene annotation
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

####################################
### DO NOT NEED TO RUN EVERY TIME ###
### JUST SUMMARY METRICS          ###
####################################
## get a data summary from all cellranger summaries
# summary_table = NULL
# for (name in meta$mid){
#   path = paste0(name,'/summary.csv')
#   file = fread(path) %>% as.data.frame()
#   summary_table = rbind(summary_table,file)
# }
# # write_csv(summary_table,'snyder_multi_summary.csv')
# {par(mfrow=c(4,1))
# par(mar = c(2,2,2,2))
# hist(summary_multi$`Estimated number of cells`,main='Estimated number of cells',xlab='',breaks=13)
# hist(summary_multi[,10],main='ATAC Fraction of high-quality fragments in cells',xlab='',breaks=15)
# hist(summary_multi[,15],main='ATAC Median high-quality fragments per cell',xlab='',breaks=12)
# hist(summary_multi[,30],main='GEX Median genes per cell',xlab='',breaks=10)
# }

####################################
###    CREATING SEURAT OBJECTS   ###
###      CELL-LEVEL FILTERING    ###
####################################
## create a list of seurat objects of RNA-seq
seurat_rna_list = list()
for (name in filenames){
  count = Read10X_h5(paste0(fpath,name,'/filtered_feature_bc_matrix.h5'))
  temp = CreateSeuratObject(counts=count$`Gene Expression`,assay='RNA',project = name)
  rna_name = meta_series[name,]$scRNA
  tissue = dplyr::filter(meta_data,`Experiment accession`==rna_name)$'Biosample term name'%>%unique()
  temp = AddMetaData(temp,metadata = tissue, col.name = 'tissue')
  gen_tissue = tissue
  tissue = str_split(tissue,pattern=' ') %>% unlist()
  if ('ventricle' %in% tissue | 'atrium' %in% tissue) {gen_tissue='heart'}
  if ('colon' %in% tissue) {gen.tissue='colon'}
  if ('liver' %in% tissue) {gen.tissue='liver'}
  if ('lung' %in% tissue) {gen.tissue='lung'}
  temp = AddMetaData(temp,metadata=gen_tissue,col.name='gen.tissue')
  seurat_rna_list[[name]] = temp
  rm(count,temp,tissue,gen_tissue,rna_name)
}
saveRDS(seurat_rna_list,'/dcs04/hongkai/data/zfu17/scMulti/ENCODE/seurat_gex_list.rds')

# loading the object list with GEX data only
seurat_rna_list = readRDS('/dcs04/hongkai/data/zfu17/scMulti/ENCODE/seurat_gex_list.rds')

# add chromatin assay for each object, remove blacklist regions
dhs_over = as.matrix(findOverlaps(human_regions,blacklist_hg38_unified))[,1]
human_regions = human_regions[setdiff(1:length(human_regions),dhs_over),] 

for (name in filenames){
  atac_count = readRDS(paste0(dhspath,name,'_atac_DHS_features_all.rds'))
  atac_count = atac_count[setdiff(1:dim(atac_count)[1],dhs_over),]
  seurat_rna_list[[name]][['ATAC']] = CreateChromatinAssay(
    counts = atac_count,
    fragments = paste0(fpath,name,'/atac_fragments.tsv.gz'),
    annotation = annotation
  )
  print('Chromatin assay added.')
}

# saveRDS(seurat_rna_list,'/dcs04/hongkai/data/zfu17/scMulti/ENCODE/snyder_multi_seurat_list.rds')
# print('Saved the multiome seurat object list, GEX and ATAC included.')

# seurat_rna_list = readRDS('/dcs04/hongkai/data/zfu17/scMulti/ENCODE/snyder_multi_seurat_list.rds')

## QC PLOTS
qc_plots = list()
seurat_rna_list = lapply(seurat_rna_list,function(x){
  DefaultAssay(x) = 'RNA'
  x[['percent.mt']] = PercentageFeatureSet(x, pattern = "^MT-")
  DefaultAssay(x) = 'ATAC'
  x = NucleosomeSignal(object=x,verbose=F)
  x = TSSEnrichment(object=x,verbose=F,fast=F)
  x
  print('finished NS and TSS.')
  saveRDS(x,paste0('/dcs04/hongkai/data/zfu17/scMulti/ENCODE/seurat_objects/',Project(x),'_seurat_obj.rds'))
  # plotting the violin plot
  p = VlnPlot(object=x,
      features = c("nCount_RNA","nFeature_RNA","percent.mt","nCount_ATAC","TSS.enrichment", 
                    "nucleosome_signal"),
      pt.size = 0.05,
      ncol = 3)
  qc_plots[[Project(x)]] = p
  rm(p)
})

## 
# print('Saving complete for all objects.')
# saveRDS(seurat_rna_list,'/dcs04/hongkai/data/zfu17/scMulti/ENCODE/snyder_multi_complete_list.rds')

saving the qc plots
for (i in 1:length(qc_plots)){
  p = qc_plots[[i]]
  ggsave(paste0('/dcs04/hongkai/data/zfu17/scMulti/ENCODE/QC_plots/',names(qc_plots)[i],'.pdf'),p,width=15,height=15,dpi=300)
  rm(p)
}

## filtering by QC
seurat_rna_list = lapply(seurat_rna_list,function(x){
  # DefaultAssay(x) = 'RNA'
  # x[['percent.mt']] = PercentageFeatureSet(x, pattern = "^MT-")
  x = subset(x,
  subset = 
    nCount_ATAC < 60000 &
    nCount_RNA < 40000 &
    nCount_ATAC > 500 &
    nCount_RNA > 500 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1 &
    percent.mt < 10)
})

## AND WE SAVE THIS FILTERED LIST IN THE END