library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(gridExtra)
library(data.table)
library(scCustomize)
library(future)
library(Signac)
library(GenomicAlignments)
library(GenomicRanges)
library(rtracklayer)

## multi-core processing
plan('multicore',workers=8)
options(future.globals.maxSize = 20 * 1024 ^ 3) ## reserve max RAM = 20 Gb

## loading meta data
meta = fread('/dcs04/hongkai/data/zfu17/scMulti/ENCODE/snyder_human_scrna_scatac.csv') %>% as.data.frame()
rownames(meta) = meta$mid

## loading DHS
load('/dcs04/hongkai/data/zfu17/tools/ENCODE_human_mouse_DHS_noY.rda')

##### remove genome blacklist regions #####
blacklist = import('/dcs04/hongkai/data/zfu17/tools/hg38-blacklist.v2.bed',format='bed') 
dhs_over = as.matrix(findOverlaps(human_regions,blacklist_hg38_unified))[,1] ##7163 regions removed
human_regions = human_regions[setdiff(1:length(human_regions),dhs_over),] ##2723210 features

## get gene annotation
# annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# seqlevelsStyle(annotation) <- "UCSC"

## creating a list of chromatin assays 
## feature matrix obtained by counting the human DHS sites
cpath = '/dcs04/hongkai/data/zfu17/scMulti/ENCODE/cellranger_human/'
setwd(cpath)

for (name in meta$mid[2:2]){
    barcodes = fread(paste0(cpath,name,'/per_barcode_metrics.csv')) %>% as.data.frame()
    atac_cells = dplyr::filter(barcodes,is_cell==1)$barcode
    frag = CreateFragmentObject(path=paste0(cpath,name,'/atac_fragments.tsv.gz'),cells=atac_cells,validate.fragments=T)
    featmat = FeatureMatrix(fragments=frag,features=human_regions,cells=atac_cells,verbose=T)
    print(paste0(name,': feature matrix constructed'))
        saveRDS(featmat,paste0('/dcs04/hongkai/data/zfu17/scMulti/ENCODE/atac_DHS_mat/',name,'_atac_DHS_features_all.rds'))
    }