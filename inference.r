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
### this line changed
library(MuDataSeurat)
library(SingleCellExperiment)
library(SeuratDisk)
## R path: /home/panyan/anaconda3/envs/seurat4/lib/R/bin
## cellphonedb path: /home/panyan/anaconda3/envs/seurat4/bin/cellphonedb

### Read in datasets







args <- commandArgs(trailingOnly = TRUE)

task_id = args[1]




inputpath <- paste("/var/www/html/tbcdb/tasks/", substr(task_id, 1, 4), "/", substr(task_id, 5, 8), "/", substr(task_id, 9, 16), "/", sep = "")
outputpath <- paste("/var/www/html/tbcdb/results/", substr(task_id, 1, 4), "/", substr(task_id, 5, 8), "/", substr(task_id, 9, 16), "/", sep = "")


parameterpath = paste(inputpath,'parameter.txt',sep="")
para <- read.table(parameterpath,  sep = " ")
filetype <- strsplit(para[[1]], ":")[[1]][2]
numgene <- as.numeric(strsplit(para[[3]], ":")[[1]][2])
numcell <- as.numeric(strsplit(para[[4]], ":")[[1]][2])
pvalthre <- as.numeric(strsplit(para[[5]], ":")[[1]][2])
iter_nums <- as.numeric(strsplit(para[[6]], ":")[[1]][2])
num_lrp <- as.numeric(strsplit(para[[7]], ":")[[1]][2])
weight1 <- as.numeric(strsplit(para[[8]], ":")[[1]][2])
weight2 <- as.numeric(strsplit(para[[9]], ":")[[1]][2])
weight3 <- as.numeric(strsplit(para[[10]], ":")[[1]][2])

weight <- c(weight1,weight2,weight3)
path_anno <- "/var/www/html/tbcdb/pipelines/data/"

para_file_path <- paste(inputpath, "parameter.txt", sep = "")
infos <- read.table(para_file_path, header = FALSE, sep = " ")

if (filetype=="txt"){
genecount_path = paste(inputpath,task_id,'.txt',sep = "")
genecount <- as.matrix(read.table(genecount_path, header = TRUE, sep = ','))
genemeta_path = paste(inputpath,task_id,'.csv',sep = "")
genemeta <- read.csv(genemeta_path, row.names=NULL)
rownames(genemeta) <- genemeta$Cell_id
genemeta <- genemeta[, !colnames(genemeta) %in% "Cell_id"] ## the rownames of genemeta must be cell_id (same to the colnames of genecount)
ser.obj <- CreateSeuratObject(genecount, meta.data=genemeta,min.cells = numcell, min.features = numgene)
ser.obj <- ser.obj[, ser.obj@meta.data$Celltype != "NA"] ## only recognized celltypes considered
} else if (filetype=="h5ad"){
genecount_path = paste(inputpath,task_id,'.h5ad',sep = "")
ser.obj <- ReadH5AD(genecount_path)
genemeta <- ser.obj@meta.data
ser.obj <- ser.obj[, ser.obj@meta.data$Celltype != "NA"] 
genecount <- GetAssayData(object = ser.obj, assay = DefaultAssay(ser.obj))
} else if (filetype=="h5Seurat"){
genecount_path = paste(inputpath,task_id,'.h5Seurat',sep = "")
ser.obj <- LoadH5Seurat(genecount_path)
genemeta <- ser.obj@meta.data
ser.obj <- ser.obj[, ser.obj@meta.data$Celltype != "NA"]
genecount <- GetAssayData(object = ser.obj, assay = DefaultAssay(ser.obj))
} else {
genecount_path = paste(inputpath,task_id,'.rds',sep = "")
ser.obj <- readRDS(genecount_path)
genecount <- GetAssayData(object = ser.obj, assay = DefaultAssay(ser.obj))
genemeta_path = paste(inputpath,task_id,'.csv',sep = "")
genemeta <- read.csv(genemeta_path, row.names=NULL)
rownames(genemeta) <- genemeta$Cell_id
genemeta <- genemeta[, !colnames(genemeta) %in% "Cell_id"] ## the rownames of genemeta must be cell_id (same to the colnames of genecount)
}





if (!class(ser.obj) == "Seurat") {
  print("The rds file is not Seurat.")
}

### Read in L-R pair database
pairs <- readRDS('/home/panyan/tbcDB/online/collected_pairs_new.rds')

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
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat) ### rds

### CellphoneDB
anno <- data.frame(Cell = rownames(genemeta),
                   cell_type = genemeta$Celltype)
write.table(anno, paste(outputpath,"temp_anno.txt",sep=""), col.names = TRUE, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(genecount, paste(outputpath,"temp_matrix.txt",sep=""), col.names = TRUE, sep = "\t", row.names = TRUE, quote = FALSE)
system(paste('/home/panyan/anaconda3/envs/seurat4/bin/cellphonedb method statistical_analysis --threads=12 ',outputpath,'temp_anno.txt ',outputpath,'temp_matrix.txt --counts-data=gene_name --database=/var/www/html/tbcdb/data/cellphoneDB/cellphonedb_user_2022-05-10-14_51.db --project-name=test --iterations=', iter_nums, sep=''))

### Celltalker

celltalker <- celltalk(input_object=ser.obj, metadata_grouping='Celltype', 
                             ligand_receptor_pairs=pairs,
                             number_cells_required=100,
                              min_expression=numgene,
                              max_expression=20000,
			     scramble_times=iter_nums)

### Combined

cpdb_means <- read.table(paste(outputpath,'results/out/test/means.txt',sep=""), sep='\t',header=TRUE)
cpdb_pval <- read.table(paste(outputpath,'results/out/test/pvalues.txt',sep=""), sep='\t',header=TRUE)
cpdb_decon <- read.table(paste(outputpath,'results/out/test/deconvoluted.txt',sep=""), sep='\t')
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
  


                             
## Draw
invisible({
    num_colors <- 30
    palette <- c(brewer.pal(8, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))

    name_list <- names(table(cellchat@idents))
    value_list <- as.numeric(table(cellchat@idents))
    symbolSize_list <- 10 + value_list*0.05
    symbolSize_list <- lapply(symbolSize_list, function(x) ifelse(x>50, 50, x))
    color_list <- palette[1:length(value_list)]
    nodes_df <- data.frame(
      name = name_list,
      value = value_list,
      symbolSize = symbolSize_list,
      color = color_list
    )
    convertToJSON <- function(csvPath, nodes_df) {
  
        cellchat_df <- read.csv(csvPath, header = TRUE)
        cellchat_df <- cellchat_df[(cellchat_df$pval <= 0.05) & (cellchat_df$score != 0), ]
        print(nrow(cellchat_df))
        cell1_pre <- cellchat_df$celltype1[1]
        cell2_pre <- cellchat_df$celltype2[1]
        count <- 1
        sum <- 0
        links_df <- data.frame()
        attended_nodes <- c()
        for (i in c(1:nrow(cellchat_df))){
            cell1 <- cellchat_df$celltype1[i]
            cell2 <- cellchat_df$celltype2[i]
            if (cell1 != cell1_pre | cell2 != cell2_pre){
                links_df <- rbind(links_df, c(cell1, cell2, count, round(sum/count,3)*80, round(sum/count,3)))
                count <- 1
                sum <- 0
                if (!(cell1 %in% attended_nodes)){
                    attended_nodes <- c(attended_nodes, cell1)
                }
                if (!(cell2 %in% attended_nodes)){
                    attended_nodes <- c(attended_nodes, cell2)
                }
                cell1_pre <- cell1
                cell2_pre <- cell2
            }
            else{
                count <- count + 1
                sum <- sum + cellchat_df$score[i]
            }

        }
        colnames(links_df) <- c("source", "target", "count", "value", "means")
        links_df$value <- lapply(links_df$value, function(x) ifelse(x<0.01, 0, x))
        links_df$value <- lapply(links_df$value, function(x) ifelse(x>7, 7, x))

        links_df$count <- as.integer(links_df$count)
        links_df$means <- as.numeric(links_df$means)

        nodes_df[!(nodes_df$name %in% attended_nodes), 'color'] <- '#808080'

        json_data <- list(
          'nodes' = lapply(unname(split(nodes_df, seq(nrow(nodes_df)))), as.list),
          'links' = lapply(unname(split(links_df, seq(nrow(links_df)))),as.list)
        )

      return(json_data)
    }
                            
                                 
    json_data_combined <- convertToJSON(paste(outputpath,task_id,'_combined_3d.csv',sep = ""), nodes_df)
    json_data = list(
    'combined' = json_data_combined
    )
    

    write(toJSON(json_data, pretty = TRUE, row.names = FALSE, auto_unbox=TRUE), file = paste(outputpath,task_id,'_combined.json',sep = ""))
})

## Draw Statistics
celltype_df <- table(df$celltype1)
out_celltype = data.frame()
for (i in c(1:length(celltype_df))){
    out_celltype = rbind(out_celltype, c(names(celltype_df)[i], celltype_df[i][[1]]))
}
colnames(out_celltype) <- c('label', 'value')
merged_df2 = df[df$ligand!='',]
a = head(sort(table(merged_df2$ligand), decreasing = TRUE),num_lrp)
out_ligand = data.frame()
for (i in c(1:10)){
    out_ligand = rbind(out_ligand, c(names(a)[i], a[i][[1]]))
}
merged_df2 = df[df$receptor!='',]
a = head(sort(table(merged_df2$receptor), decreasing = TRUE),num_lrp)
out_receptor = data.frame()
for (i in c(1:10)){
    out_receptor = rbind(out_receptor, c(names(a)[i], a[i][[1]]))
}
a = head(sort(table(df$interaction), decreasing = TRUE),num_lrp)
out_interaction = data.frame()
for (i in c(1:10)){
    out_interaction = rbind(out_interaction, c(names(a)[i], a[i][[1]]))
}

json_data = list(
    'celltype' = out_celltype,
    'ligand' = unname(out_ligand),
    'receptor' = unname(out_receptor),
    'interaction' = unname(out_interaction)
)
write(toJSON(json_data, pretty = TRUE, row.names = FALSE, auto_unbox=TRUE), file = paste(outputpath,task_id,'_combined_statistics.json',sep = ""))


## Draw heatplot for celltype-celltype
means_2d <- means_3d[,,1]
for (i in 2:dim_b)
{
    means_2d <- means_2d+means_3d[,,i]
}

rownames(means_2d) <- celltype
colnames(means_2d) <- celltype
melt_weight =  melt(as.matrix(means_2d),na.rm=T)
colnames(melt_weight) <- c('source_celltype', 'target_celltype', 'value')
melt_weight$value <- (melt_weight$value - min(melt_weight$value))/(max(melt_weight$value)- min(melt_weight$value))
ggp <- ggplot(melt_weight, aes(source_celltype, target_celltype)) +                           
  geom_tile(aes(fill = value))+ theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1), axis.text.y = element_text(size = 12))
ggsave(paste(outputpath,task_id,"celltype_heatplot.png",sep = ""), plot = ggp, width = 9, height = 6, dpi = 500) 
## Draw heatplot for source cell-top L-R pairs
means_2d <- means_3d[1,,]
for (i in 2:num)
{
    means_2d <- means_2d+means_3d[i,,]
}

colnames(means_2d) <- combined_pairs
melt_weight =  melt(as.matrix(means_2d[,names(a)]),na.rm=T)
colnames(melt_weight) <- c('source_celltype', 'ligand_receptor_pairs', 'value')
ggp <- melt_weight %>% mutate(`normalized-value` = (melt_weight$value - min(melt_weight$value))/(max(melt_weight$value)- min(melt_weight$value))) %>% 
  ggplot(aes(x=source_celltype, y = ligand_receptor_pairs, color = value, size = `normalized-value`)) + 
  geom_point() + theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1), axis.text.y = element_text(size = 12))
ggsave(paste(outputpath,task_id,"source_heatplot.png",sep = ""), plot = ggp, width = 9, height = 6, dpi = 500)   


## Draw heatplot for target cell-top L-R pairs
means_2d <- means_3d[,1,]
for (i in 2:num)
{
    means_2d <- means_2d+means_3d[,i,]
}

colnames(means_2d) <- combined_pairs
melt_weight =  melt(as.matrix(means_2d[,names(a)]),na.rm=T)
colnames(melt_weight) <- c('target_celltype', 'ligand_receptor_pairs', 'value')
ggp <- melt_weight %>% mutate(`normalized-value` = (melt_weight$value - min(melt_weight$value))/(max(melt_weight$value)- min(melt_weight$value))) %>% 
  ggplot(aes(x=target_celltype, y = ligand_receptor_pairs, color = value, size = `normalized-value`)) + 
  geom_point() + theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1), axis.text.y = element_text(size = 12))
ggsave("target_heatplot.png", plot = ggp, width = 9, height = 6, dpi = 500)   