library(CellChat)
library(celltalker)
library(ggplot2)                  
library(patchwork)
library(ggalluvial)
library(igraph)
library(dplyr)
library(tidyr)
options(stringsAsFactors = FALSE)
library(Seurat)
library(abind)
library(reshape2)
library(CombinePValue)
library(stringr)
library(purrr)
library(remotes)
library(jsonlite)
library(RColorBrewer)
library(DT)
library(MuDataSeurat)
library(SingleCellExperiment)
library(SeuratDisk)
## R path: /home/panyan/anaconda3/envs/seurat4/lib/R/bin
## cellphonedb path: /home/panyan/anaconda3/envs/seurat4/bin/cellphonedb

### Read in datasets
numgene <- 0
numcell <- 0
pvalthre <- 0.05
iter_nums <- 30
num_lrp <- 10
weight1 <- 1
weight2 <- 1
weight3 <- 1
task_id <- 0

weight <- c(weight1,weight2,weight3)
args <- commandArgs(trailingOnly = TRUE)
filetype <- args[1]
genecount_path <- args[2]
genemeta_path <- args[3]

if (filetype=="txt"){
genecount <- as.matrix(read.table(genecount_path, header = TRUE, sep = ','))
genemeta <- read.csv(genemeta_path, row.names=NULL)
rownames(genemeta) <- genemeta$Cell_id
genemeta <- genemeta[, !colnames(genemeta) %in% "Cell_id"] ## the rownames of genemeta must be cell_id (same to the colnames of genecount)
ser.obj <- CreateSeuratObject(genecount, meta.data=genemeta,min.cells = numcell, min.features = numgene)
} else if (filetype=="h5ad"){
ser.obj <- ReadH5AD(genecount_path)
} else if (filetype=="h5Seurat"){
ser.obj <- LoadH5Seurat(genecount_path)
} else {
ser.obj <- readRDS(genecount_path)

}

ser.obj <- subset(ser.obj, subset = nFeature_RNA >= numgene)
ser.obj <- subset(ser.obj, subset = nCount_RNA >= numcell)
ser.obj <- ser.obj[, ser.obj@meta.data$Celltype != "NA"]
print(ser.obj)
genemeta <- ser.obj@meta.data
genecount <- GetAssayData(object = ser.obj, assay = DefaultAssay(ser.obj))


if (!class(ser.obj) == "Seurat") {
  print("The rds file is not Seurat.")
}

### Read in L-R pair database
pairs <- readRDS('../online/collected_pairs_new.rds')

### Read in project name
project = 'test'

start_time <- Sys.time()

### Proprocessing
ser.obj <- NormalizeData(ser.obj)
ser.obj <- FindVariableFeatures(ser.obj)
ser.obj <- ScaleData(ser.obj)
ser.obj <- RunPCA(ser.obj)
ser.obj <- RunUMAP(ser.obj,reduction="pca",dims=1:15)
ser.obj <- FindNeighbors(ser.obj,reduction="pca",dims=1:15)
ser.obj <- FindClusters(ser.obj,resolution=0.5)

### CellChat

cellchat<- createCellChat(object =ser.obj, group.by = "Celltype")
future::plan("multicore", workers = 50) # do parallel
memory_size = 10
options(future.globals.maxSize = memory_size*1000 * 1024^2) # 1024^2 = MB

cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat@idents = droplevels(cellchat@idents, exclude = setdiff(levels(cellchat@idents),unique(cellchat@idents)))
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat) ### rds

### CellphoneDB
anno <- data.frame(Cell = rownames(genemeta),
                   cell_type = genemeta$Celltype)
outputpath <- "./"
write.table(anno, paste(outputpath,"temp_anno.txt",sep=""), col.names = TRUE, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(genecount, paste(outputpath,"temp_matrix.txt",sep=""), col.names = TRUE, sep = "\t", row.names = TRUE, quote = FALSE)
system(paste('cellphonedb method statistical_analysis --threads=12 ',outputpath,'temp_anno.txt ',outputpath,'temp_matrix.txt --counts-data=gene_name --project-name=test --iterations=', iter_nums, sep=''))

### Celltalker

celltalker <- celltalk(input_object=ser.obj, metadata_grouping='Celltype', 
                             ligand_receptor_pairs=pairs,
                             number_cells_required=100,
                              min_expression=numgene,
                              max_expression=20000,
			     scramble_times=iter_nums)

### Combined

cpdb_means <- read.table(paste(outputpath,'../online/out/test/means.txt',sep=""), sep='\t',header=TRUE)
cpdb_pval <- read.table(paste(outputpath,'../online/out/test/pvalues.txt',sep=""), sep='\t',header=TRUE)
cpdb_decon <- read.table(paste(outputpath,'../online/out/test/deconvoluted.txt',sep=""), sep='\t')
colnames(cpdb_decon) <- cpdb_decon[1,]
cpdb_decon <- cpdb_decon[-1,]
celltype_cpdb <- colnames(cpdb_decon)[-c(1:6)]

### convert each to a 3D matrix
num <- dim(cellchat@net$prob)[1]
celltype <- colnames(cellchat@net$weight)
### cellphonedb
num_pairs <- nrow(cpdb_means)
cpdb_means_3d <- array(rep(0, num*num*num_pairs), dim=c(num, num, num_pairs))
cpdb_pval_3d <- array(rep(1, num*num*num_pairs), dim=c(num, num, num_pairs))
idx <- rep(0,num)
for (i in 1:num)
{
    if (celltype_cpdb[i] %in% celltype) idx[i] <- grep(celltype_cpdb[i], celltype, fixed=TRUE)
}
for (i in 1:num_pairs)
{
    for (j in 1:num*num)
    {
        if (idx[j%/%num] == 0) next
        cpdb_means_3d[idx[j%/%num],idx[((j-1)%%num)+1],i] <- cpdb_means[i,j+11]
        cpdb_pval_3d[idx[j%/%num],idx[((j-1)%%num)+1],i] <- cpdb_pval[i,j+11]
    }
}
# cpdb_means_3d <- as.double(cpdb_means_3d)
# cpdb_pval_3d <- as.double(cpdb_pval_3d)
cpdb_pairs <- cpdb_means$interacting_pair


###celltalker
num_pairs <- length(table(celltalker$interaction))
ct_pairs <- names(table(celltalker$interaction))
ct_ligand <- as.character(map(strsplit(ct_pairs, split = '_'), 1))
ct_receptor <- as.character(map(strsplit(ct_pairs, split = '_'), 2))
ct_means_3d <- array(rep(0, num*num*num_pairs), dim=c(num, num, num_pairs))
ct_pval_3d <- array(rep(1, num*num*num_pairs), dim=c(num, num, num_pairs))
num_pairs_temp <- celltalker$interaction[1]
j = 1
for (i in 1:nrow(celltalker))
{
    if (celltalker$interaction[i] != num_pairs_temp) {
        j = j+1
        num_pairs_temp = celltalker$interaction[i]
    }
    if (celltalker$cell_type1[i] %in% celltype && celltalker$cell_type2[i] %in% celltype){ 
        idx1 <- grep(celltalker$cell_type1[i], celltype)
        idx2 <- grep(celltalker$cell_type2[i], celltype)
        ct_means_3d[idx1, idx2, j] <- celltalker$value[i]
        ct_pval_3d[idx1, idx2, j] <- celltalker$p_val[i]
    }
}

cc_prob_3d <- cellchat@net$prob
cc_pval_3d <- cellchat@net$pval

## normalization
ct_means_3d_nml <- (ct_means_3d - min(ct_means_3d)) / (max(ct_means_3d) - min(ct_means_3d))
cpdb_means_3d_nml <- (cpdb_means_3d - min(cpdb_means_3d)) / (max(cpdb_means_3d) - min(cpdb_means_3d))
cc_prob_3d_nml <- (cc_prob_3d - min(cc_prob_3d)) / (max(cc_prob_3d) - min(cc_prob_3d))
cc_pairs <- cellchat@LR$LRsig$interaction_name
cc_ligand <- cellchat@LR$LRsig$ligand
cc_receptor <- cellchat@LR$LRsig$receptor
intersect_1 <- intersect(cc_pairs, ct_pairs)
intersect_2 <- intersect(cc_pairs, cpdb_means$interacting_pair)
cpdb_pairs <- cpdb_means$interacting_pair
cpdb_ligand <- cpdb_means$gene_a
cpdb_receptor <- cpdb_means$gene_b

## integration
for (i in 1:dim(ct_means_3d)[3])
{
    if(ct_pairs[i] %in% intersect_1)
    {
        idx <- grep(ct_pairs[i], intersect_1)
        cc_prob_3d_nml[,,idx] <- (cc_prob_3d_nml[,,idx]*weight[1]+ct_means_3d_nml[,,i]*weight[3])/(weight[1]+weight[3])
        for(m in 1:num)
        {
            for(n in 1:num)
            {
                cc_pval_3d[m,n,idx] <- selfcontained.test(pvalue=c(cc_pval_3d[m,n,idx],ct_pval_3d[m,n,i]))$'significance level for combining pvalues'
            }
        }
    }
}

for (i in 1:dim(cpdb_means_3d)[3])
{
    if(cpdb_pairs[i] %in% intersect_2)
    {
        idx <- grep(cpdb_pairs[i], intersect_2)
        cc_prob_3d_nml[,,idx] <- (cc_prob_3d_nml[,,idx]*weight[1]+cpdb_means_3d_nml[,,i]*weight[2])/(weight[1]+weight[2])
        for(m in 1:num)
        {
            for(n in 1:num)
            {
                cc_pval_3d[m,n,idx] <- selfcontained.test(pvalue=c(cc_pval_3d[m,n,idx],cpdb_pval_3d[m,n,i]))$'significance level for combining pvalues'
            }
        }
    }
}

pval_3d <- abind(cc_pval_3d, cpdb_pval_3d[,,!cpdb_pairs %in% intersect_2], ct_pval_3d[,,!ct_pairs %in% intersect_1], along=3)
means_3d <- abind(cc_prob_3d_nml, cpdb_means_3d_nml[,,!cpdb_pairs %in% intersect_2], ct_means_3d_nml[,,!ct_pairs %in% intersect_1], along=3)
means_3d <- (means_3d - min(means_3d)) / (max(means_3d) - min(means_3d))
print(dim(pval_3d))
print(dim(means_3d))

## save to mysql csv

combined_pairs <- c(cc_pairs, cpdb_pairs[!cpdb_pairs %in% intersect_2], ct_pairs[!ct_pairs %in% intersect_1])
combined_ligand <- c(cc_ligand, cpdb_ligand[!cpdb_pairs %in% intersect_2], ct_ligand[!ct_pairs %in% intersect_1])
combined_receptor <- c(cc_receptor, cpdb_receptor[!cpdb_pairs %in% intersect_2], ct_receptor[!ct_pairs %in% intersect_1])
dim_a = dim(pval_3d)[1]
dim_b = dim(pval_3d)[3]

id = 1
df <- matrix(0, ncol = 10, nrow = dim_a * dim_a * dim_b)
for (i in 1:dim_a){
    for (j in 1:dim_a){
        for (k in 1:dim_b){
            if (pval_3d[i,j,k] <= pvalthre)
            {
                df[id, 1] = id
                df[id, 2] = 'Integrated'
                df[id, 3] = task_id
                df[id, 4] = combined_ligand[k]
                df[id, 5] = combined_receptor[k]
                df[id, 6] = combined_pairs[k]
                df[id, 7] = celltype[i]
                df[id, 8] = celltype[j]
                df[id, 9] = means_3d[i,j,k]
		
		if (pval_3d[i,j,k] < 0.000001) {
  			formatted_num <- '<1e-6'
		}
		else
		{
			formatted_num <- formatC(pval_3d[i,j,k], format = "f", digits = 6)
		}
                df[id, 10] = as.character(formatted_num)
                id = id + 1
            }
        
        }
    }
}
print(id)
df <- as.data.frame(df[1:(id-1),])
colnames(df) <- c('id', 'method', 'project', 'ligand', 'receptor', 'interaction', 'celltype1', 'celltype2', 'score', 'pval')
file <- file(paste(outputpath,task_id,'_combined_3d.csv',sep = ""), encoding="UTF-8")    
write.csv(df, file, quote=F, row.names=FALSE)

                           
end_time <- Sys.time() - start_time
print(paste('Time: ', end_time, ' s', sep=''))
  

