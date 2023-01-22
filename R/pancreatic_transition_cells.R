library(Seurat)
library(dplyr)

# compute transition index
pearson<-abs(pearson)
res_1<-apply(pearson,1,function(x) ks.test(x,pearson[1,],alternative='greater')$statistic)
res_2<-apply(pearson,1,function(x) ks.test(x,pearson[1,],alternative='less')$statistic)
res<-find_ks_d(pearson[which.max(res_1),],pearson[which.max(res_2),])
tmp<-apply(pearson,1,function(x) sum(abs(x)>min(unlist(res)) & abs(x)<max(unlist(res)),na.rm=T)/sum(abs(x)>=0,na.rm=T))
data<-AddMetaData(data,data.frame('pearson'=tmp))

# compare transition index for metaplastic cells
meta<-subset(data,celltype=='Metaplastic')
meta@meta.data %>% filter(time_point!='3M') %>% group_by(time_point) %>% summarize(p_value = wilcox.test(pearson,meta_ctrl$pearson)$p.value) %>% 
mutate(p_adj=p.adjust(p_value,method='BH',n=3))
#  A tibble: 3 x 3
#  time_point  p_value    p_adj
#  <chr>         <dbl>    <dbl>
#1 15M        5.47e- 7 8.21e- 7
#2 5M         5.63e- 1 5.63e- 1
#3 9M         1.04e-46 3.12e-46

# separate transitioning and stable metaplastic cells at 3M PTI
pearson_sub<-subset(meta,time_point %in% c('3M'))$pearson
library(mclust)
fit <- Mclust(pearson_sub,G=2)
pdf('meta_3m_pearson_dist.pdf')
plot(density(pearson_sub),main="Pearson's correlation coefficients distribution for 3M PTI")
rug(jitter(pearson_sub[fit$classification==2 & fit$z[,2]>0.9]),col='red')
rug(jitter(pearson_sub[fit$classification==1 & fit$z[,1]>0.9]),col='blue')
dev.off()
transition_cells<-names(pearson_sub[fit$classification==2 & fit$z[,2]>0.9])
ctrl<-pearson_sub[fit$classification==1 & fit$z[,1]>0.9] %>% names()
trans<-meta[,transition_cells]
trans<-FindVariableFeatures(trans, selection.method = "vst", nfeatures = 2000)
trans<- ScaleData(trans)
trans<- RunPCA(trans, features = VariableFeatures(object = trans))
pdf('elbow.pdf')
print(ElbowPlot(trans,ndim=50))
dev.off()
trans<- FindNeighbors(trans, dims = 1: 13) 
trans<- FindClusters(trans, resolution = 1)
saveRDS(trans,'meta_trans.rds')
for(i in trans$seurat_clusters %>% unique()){
    cells<-subset(trans,seurat_clusters==i) %>% colnames()
    res<-FindMarkers(meta,ident.1=cells,ident.2=ctrl)
    write.csv(res,paste0('meta_deg_cluster',i,'.csv'))
}

# gene expressions violin plot
pdf('trans_cells.pdf',height=3)
VlnPlot(trans, features = c('Muc6','Pgc','Tff1'),pt.size=0.00) # gastric
VlnPlot(trans, features = c('Cd44','Krt7','Fn1'),pt.size=0.00) # tumor
VlnPlot(trans, features = c('Mki67','Cdk1','Cdc20'),pt.size=0.00) #dividing
dev.off()

# transitioning acinar cells at 6W PTI 
acinar<-subset(data,celltype=='Acinar')
pearson_sub<-subset(acinar,time_point %in% c('6W'))$pearson
library(mclust)
fit <- Mclust(pearson_sub,G=2)
pdf('acinar_6W_pearson_dist.pdf')
plot(density(pearson_sub),main="Transition index distribution for 6W PTI")
rug(jitter(pearson_sub[fit$classification==2 & fit$z[,2]==1]),col='red')
rug(jitter(pearson_sub[fit$classification==1 & fit$z[,1]>0.8]),col='blue')
dev.off()
transition_cells<-names(pearson_sub[fit$classification==2 & fit$z[,2]==1])
ctrl<-subset(acinar,time_point %in% c('CTRL')) %>% colnames()
all_cell<-subset(acinar,time_point %in% c('6W')) %>% colnames()
res<-FindMarkers(acinar,ident.1=transition_cells,ident.2=ctrl)
gene_trans<-(res %>% filter(p_val_adj<0.01,avg_log2FC>1))
res<-FindMarkers(acinar,ident.1=all_cell,ident.2=ctrl)
gene_all<-(res %>% filter(p_val_adj<0.01,avg_log2FC>1))
write.csv(gene_trans,'acinar_6w_transition_cell_deg.csv')
write.csv(gene_all,'acinar_6w_all_cell_deg.csv')
