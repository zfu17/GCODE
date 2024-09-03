library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(data.table)
library(scCustomize)
library(future)
# library(EnsDb.Hsapiens.v86)
# library(BSgenome.Hsapiens.UCSC.hg38)
library(readxl)
library(Signac)
library(dplyr)
library(harmony)

Sys.time() %>% print()

## paralle 
# plan('multisession',workers=10)
# options(future.globals.maxSize = 20*1024^3) # 30Gb of RAM

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
# load('/dcs04/hongkai/data/zfu17/tools/ENCODE_human_mouse_DHS_noY.rda')

## get gene annotation
# annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# seqlevelsStyle(annotation) <- "UCSC"

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

## create a list of seurat objects of RNA-seq
# seurat_rna_list = list()
# for (name in filenames){
#   count = Read10X_h5(paste0(fpath,name,'/filtered_feature_bc_matrix.h5'))
#   temp = CreateSeuratObject(counts=count$`Gene Expression`,assay='RNA',project = name)
#   rna_name = meta_series[name,]$scRNA
#   tissue = dplyr::filter(meta_data,`Experiment accession`==rna_name)$'Biosample term name'%>%unique()
#   temp = AddMetaData(temp,metadata = tissue, col.name = 'tissue')
#   gen_tissue = tissue
#   tissue = str_split(tissue,pattern=' ') %>% unlist()
#   if ('ventricle' %in% tissue | 'atrium' %in% tissue) {gen_tissue='heart'}
#   if ('colon' %in% tissue) {gen.tissue='colon'}
#   if ('liver' %in% tissue) {gen.tissue='liver'}
#   if ('lung' %in% tissue) {gen.tissue='lung'}
#   temp = AddMetaData(temp,metadata=gen_tissue,col.name='gen.tissue')
#   seurat_rna_list[[name]] = temp
#   rm(count,temp,tissue,gen_tissue,rna_name)
# }
# saveRDS(seurat_rna_list,'/dcs04/hongkai/data/zfu17/scMulti/ENCODE/seurat_gex_list.rds')

## loading the object list with GEX data only
# seurat_rna_list = readRDS('/dcs04/hongkai/data/zfu17/scMulti/ENCODE/seurat_gex_list.rds')

## add chromatin assay for each object, remove blacklist regions
# dhs_over = as.matrix(findOverlaps(human_regions,blacklist_hg38_unified))[,1]
# human_regions = human_regions[setdiff(1:length(human_regions),dhs_over),] 

# for (name in filenames){
#   atac_count = readRDS(paste0(dhspath,name,'_atac_DHS_features_all.rds'))
#   atac_count = atac_count[setdiff(1:dim(atac_count)[1],dhs_over),]
#   seurat_rna_list[[name]][['ATAC']] = CreateChromatinAssay(
#     counts = atac_count,
#     fragments = paste0(fpath,name,'/atac_fragments.tsv.gz'),
#     annotation = annotation
#   )
#   print('Chromatin assay added.')
# }

# saveRDS(seurat_rna_list,'/dcs04/hongkai/data/zfu17/scMulti/ENCODE/snyder_multi_seurat_list.rds')
# print('Saved the multiome seurat object list, GEX and ATAC included.')

# seurat_rna_list = readRDS('/dcs04/hongkai/data/zfu17/scMulti/ENCODE/snyder_multi_seurat_list.rds')

## quality control with GEX and ATAC
# qc_plots = list()
# seurat_rna_list = lapply(seurat_rna_list,function(x){
#   DefaultAssay(x) = 'RNA'
#   x[['percent.mt']] = PercentageFeatureSet(x, pattern = "^MT-")
#   DefaultAssay(x) = 'ATAC'
#   x = NucleosomeSignal(object=x,verbose=F)
#   x = TSSEnrichment(object=x,verbose=F,fast=F)
#   x
#   print('finished NS and TSS.')
#   saveRDS(x,paste0('/dcs04/hongkai/data/zfu17/scMulti/ENCODE/seurat_objects/',Project(x),'_seurat_obj.rds'))
#   # plotting the violin plot
#   p = VlnPlot(object=x,
#       features = c("nCount_RNA","nFeature_RNA","percent.mt","nCount_ATAC","TSS.enrichment", 
#                     "nucleosome_signal"),
#       pt.size = 0.05,
#       ncol = 3)
#   qc_plots[[Project(x)]] = p
#   rm(p)
# })


## 
# print('Saving complete for all objects.')
# saveRDS(seurat_rna_list,'/dcs04/hongkai/data/zfu17/scMulti/ENCODE/snyder_multi_complete_list.rds')
# print('big list saved.')

# saving the qc plots
# for (i in 1:length(qc_plots)){
#   p = qc_plots[[i]]
#   ggsave(paste0('/dcs04/hongkai/data/zfu17/scMulti/ENCODE/QC_plots/',names(qc_plots)[i],'.pdf'),p,width=15,height=15,dpi=300)
#   rm(p)
# }


## filtering by subsetting
# seurat_rna_list = lapply(seurat_rna_list,function(x){
#   # DefaultAssay(x) = 'RNA'
#   # x[['percent.mt']] = PercentageFeatureSet(x, pattern = "^MT-")
#   x = subset(x,
#   subset = 
#     nCount_ATAC < 60000 &
#     nCount_RNA < 40000 &
#     nCount_ATAC > 500 &
#     nCount_RNA > 500 &
#     nucleosome_signal < 2 &
#     TSS.enrichment > 1 &
#     percent.mt < 10)
# })

# loading individual filtered seurat object, keep only the GEX portion
# seurat_rna_list = list()
# for (name in filenames){
#   obj = readRDS(paste0(filtered_objpath,name,'_filtered_seurat_obj.rds'))
#   DefaultAssay(obj) = 'RNA'
#   obj[['ATAC']] = NULL ## remove the chromatin assay
#   seurat_rna_list[[name]] = obj
#   seurat_rna_list[[name]]
#   print('object loaded, only GEX was kept')
#   rm(obj)
# }

# Sys.time() %>% print()
# print('Finished loading all filtered objects, GEX kept.')


## merging all objects
# merged_filtered=merge(seurat_rna_list[[1]],as.vector(seurat_rna_list[2:length(seurat_rna_list)]),
#                            add.cell.ids = filenames,
#                            project = 'Snyder Human Multiome')
# Sys.time() %>% print()
# print('Finished merging.')    
# merged_filtered     
# saveRDS(merged_filtered,'/dcs04/hongkai/data/zfu17/scMulti/ENCODE/merged_GEX_filtered.rds')
# rm(seurat_rna_list)         

## loading merged object
# merged_filtered = readRDS('/dcs04/hongkai/data/zfu17/scMulti/ENCODE/merged_GEX_filtered.rds')
# print('Merged object loaded')
# merged_filtered

# # ## standard workflow to the merged object
# DefaultAssay(merged_filtered) = 'RNA'
# merged_filterd = NormalizeData(merged_filtered)
# merged_filtered = FindVariableFeatures(merged_filtered)
# merged_filtered = ScaleData(merged_filtered)
# merged_filtered = RunPCA(merged_filtered)
# elbow = ElbowPlot(merged_filtered,ndims=50)
# merged_filtered = FindNeighbors(merged_filtered,dims=1:20)
# merged_filtered = FindClusters(merged_filtered)
# merged_filtered = RunUMAP(merged_filtered,dim=1:20)

# print('Pipeline finished, begin saving.')

# saveRDS(merged_filtered,'/dcs04/hongkai/data/zfu17/scMulti/ENCODE/merged_GEX_filterd_processed.rds')

# # ## plotting
# p1 = DimPlot(merged_filtered,reduction='umap',group.by = 'orig.ident')+ ggtitle('Sample')
# p2 = DimPlot(merged_filtered,reduction='umap',group.by = 'tissue')+ ggtitle('Tissue type')
# p3 = DimPlot(merged_filtered,reduction='umap',group.by = 'gen.tissue')+ ggtitle('General tissue type')
# # # grid.arrange(p1, p2, ncol = 1, nrow = 2)

# save(elbow,p1,p2,p3,'/dcs04/hongkai/data/zfu17/scMulti/ENCODE/merged_plots.Rdata')

## integrate these objects
# seurat_rna_list = lapply(seurat_rna_list,function(x){
#   x = NormalizeData(x)
#   x = FindVariableFeatures(x,selection.method = "vst", nfeatures = 2000)
# })
# print('Objects loaded, normalized, variable features found in each.')
# features <- SelectIntegrationFeatures(object.list = seurat_rna_list)
# print('Integration features selected. Run harmony.')

# gex.anchors = FindIntegrationAnchors(object.list = seurat_rna_list, anchor.features = features)
# rm(seurat_rna_list)
# print('Begin integration.')
# gex_integrated = IntegrateData(anchorset = gex.anchors)
# saveRDS(gex_integrated,'/dcs04/hongkai/data/zfu17/scMulti/ENCODE/integrated_GEX.rds')
# print('Integration finished. Begin pipeline.')


# gex_integrated = ScaleData(gex_integrated)
# gex_integrated = RunPCA(gex_integrated)
# gex_integrated = FindNeighbors(gex_integrated,dims=1:20)
# gex_integrated = FindClusters(gex_integrated)
# gex_integrated = RunUMAP(gex_integrated,dim=1:20)
# print('Begin saving the processed integrated object')
# saveRDS(gex_integrated,'/dcs04/hongkai/data/zfu17/scMulti/ENCODE/integrated_GEX_processed.rds')

## run harmony, since seurat integration (finding anchors doesn't work)
# merged_filtered = readRDS('/dcs04/hongkai/data/zfu17/scMulti/ENCODE/merged_GEX_filtered_processed.rds') 
# DefaultAssay(merged_filtered) = 'RNA'
# print('data loaded. begin harmony')
# merged_filtered = RunHarmony(merged_filtered,group.by.vars = 'orig.ident') 
# print('harmony finished, saving data')
# saveRDS(merged_filtered,'/dcs04/hongkai/data/zfu17/scMulti/ENCODE/merged_GEX_harmony.rds')

# merged_harmony = readRDS('/dcs04/hongkai/data/zfu17/scMulti/ENCODE/merged_GEX_harmony.rds')
# merged_harmony = RunUMAP(merged_harmony,reduction='harmony',dim=1:20)
# merged_harmony = FindNeighbors(merged_harmony,reduction = "harmony", dims = 1:20)
# merged_harmony = FindClusters(merged_harmony,resolution = 0.5)
# saveRDS(merged_harmony,'/dcs04/hongkai/data/zfu17/scMulti/ENCODE/merged_GEX_harmony_processed.rds')

# merged_filtered = readRDS('/dcs04/hongkai/data/zfu17/scMulti/ENCODE/merged_GEX_filtered_processed.rds')
# merged_filtered = FindClusters(merged_filtered,resolution=1,n.iter=1000)
# saveRDS(merged_filtered,'/dcs04/hongkai/data/zfu17/scMulti/ENCODE/merged_GEX_filtered_cluster.rds')

### get variable features 
## loading individual filtered seurat object, keep only the GEX portion
varlist = list()
for (name in filenames){
  obj = readRDS(paste0(filtered_objpath,name,'_filtered_seurat_obj.rds'))
  DefaultAssay(obj) = 'RNA'
  varlist[[name]] = obj@assays[["RNA"]]@var.features
  print('get 1 var')
  rm(obj)
}
varlist = do.call(cbind,varlist) %>% as.data.frame()
colnames(varlist) = filenames
fwrite(varlist,'/dcs04/hongkai/data/zfu17/scMulti/ENCODE/variable_features.csv')