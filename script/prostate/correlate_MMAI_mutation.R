#######umap based on MMAI features##############################################################################



feature=read.csv("/Users/ysong/Desktop/P_Tran_ref/score/to_use/score_111/image_predict_score/MMAI_score_feature_128_96.mod.csv",row.names=2)
mutation=read.table("/Users/ysong/Desktop/P_Tran_ref/score/to_use/score_111/image_predict_score/mutation.169.txt",sep="\t",header=T,check.names=F)
#data=read.csv("/Users/ysong/Desktop/P_Tran_ref/score/to_use/score_111/image_predict_score/raw_score.csv")
colnames(mutation)=gsub("-T","",colnames(mutation))

meta.umap=read.delim2("/Users/ysong/Desktop/P_Tran_ref/score/to_use/score_111/image_predict_score/meta.mutation.numeric.txt",sep="\t",header=T)


feature.meta=merge(feature,meta.umap,by.y="predictive_score",by.x="stadt_predictive_v1_3_0")

rownames(feature.meta)=feature.meta$Tumor_Sample_Barcode



feature_96=feature.meta[c(1,132:227)]
feature_128=feature.meta[,c(3:131)]

library(ggplot2)
library(umap)



my_umap <- as.data.frame(umap(feature_128[-1])$layout)

graph <- bluster::makeSNNGraph(feature_128[-1],k=10)
#graph <- bluster::makeSNNGraph(feature,k=10)
#graph <- bluster::makeSNNGraph(feature,k=5)
#> Warning in (function (to_check, X, clust_centers, clust_info, dtype, nn, :
#> detected tied distances to neighbors, see ?'BiocNeighbors-ties'
clust <- igraph::cluster_louvain(graph)

my_umap$clust <- factor(clust$membership)

my_umap.ann=unique(merge(my_umap,feature.meta[,c(1,228:ncol(feature.meta))],by.x=0,by.y=0))

#write.table(my_umap.ann[,1:69],"meta_w_umap_cluster_predi.txt",sep="\t",quote=F,row.names=F)

my_umap.ann[,-c(2,3,5,10,46,56,58,60)] <- lapply(my_umap.ann[,-c(2,3,5,10,46,56,58,60)], function(x) as.character(x))
my_umap.ann$ARTERA.Score=as.numeric(my_umap.ann$ARTERA.Score)
my_umap.ann$met_survival=as.numeric(my_umap.ann$met_survival)
my_umap.ann$rPFS_time=as.numeric(my_umap.ann$rPFS_time)


write.table(my_umap.ann,"~/Desktop/my_umap.ann.txt",sep="\t",quote=F,row.names=F)


file=merge(my_umap.ann[,1:60],t(mutation),by.y=0,by.x="Row.names")

write.table(file,"~/Desktop/candidate.txt",sep="\t",quote=F)

setwd("~/Desktop/prog_umap")
my_umap.ann=read.delim("my_umap.ann.txt",sep="\t",header=T)


my_umap.ann[62:ncol(my_umap.ann)] <- lapply(my_umap.ann[62:ncol(my_umap.ann)], function(x) as.factor(as.character(x)))
my_umap.ann$clust=as.factor(my_umap.ann$clust)


for (i in 62:ncol(my_umap.ann)){
  tmp=my_umap.ann[,c(1:4,i)]
  colnames(tmp)[5]="gene"
  my_colors <- c("magenta","lightgrey", "cyan")
  ggplot(tmp, aes(V1, V2,shape=clust,color=gene)) + geom_point()+ guides(color=guide_legend(colnames(my_umap.ann[i]))) +scale_color_manual(values = my_colors) +theme_bw()+xlab("umap_1")+ylab("umap_2")
  
    ggsave(paste(colnames(my_umap.ann)[i],".pdf",sep=""),width=5,height=5)
  }








my_umap.ann$Time..2.synchronous..1.metachronous.=as.factor(my_umap.ann$Time..2.synchronous..1.metachronous.)
my_umap.ann$clust=as.factor(my_umap.ann$clust)

my_umap.ann$bin <- rep(NA, nrow(my_umap.ann))

  my_umap.ann[my_umap.ann$ARTERA.Score<median(my_umap.ann$ARTERA.Score), ][, "bin"] <- "less than median"

  my_umap.ann[my_umap.ann$ARTERA.Score>=median(my_umap.ann$ARTERA.Score), ][, "bin"] <- "more than median"

  
  
 # The ggrepel package works great for repelling overlapping text labels away from each other. You can use either geom_label_repel() (draws rectangles around the text) or geom_text_repel() functions.
  
  library(ggplot2)
  library(ggrepel)
  
pdf("~/Desktop/umap.label.pdf",20,20)
sub=my_umap.ann[!grepl("OM",my_umap.ann$Tumor_Sample_Barcode),]
ggplot(my_umap.ann, aes(V1, V2, colour = clust)) +
  geom_point()+  geom_label_repel(aes(label=Tumor_Sample_Barcode,size=1),
                                  box.padding   = 0.05, 
                                  point.padding = 0.05,
                                  segment.color = 'grey50',max.overlaps = Inf)

sub=my_umap.ann[!grepl("OM",my_umap.ann$Tumor_Sample_Barcode),]
ggplot(sub, aes(V1, V2, colour = clust)) +
  geom_point()+  geom_label_repel(aes(label=Tumor_Sample_Barcode,size=1),
                                  box.padding   = 0.05, 
                                  point.padding = 0.05,
                                  segment.color = 'grey50',max.overlaps = Inf)


dev.off()

pdf("~/Desktop/umap.cluster.prog.pdf",10,5)

ggplot(my_umap.ann, aes(V1, V2, colour = Time..2.synchronous..1.metachronous.,shape=clust)) +
  geom_point()


ggplot(my_umap.ann, aes(V1, V2, colour =bin,shape=clust)) +
  geom_point()


ggplot(my_umap.ann, aes(V1, V2)) +
  geom_point(aes(color = ARTERA.Score), size = 2)+  scale_color_gradient2(low = "grey", mid = "lightyellow", high = "red", midpoint = 0.1602506)

ggplot(my_umap.ann, aes(V1, V2)) +
  geom_point(aes(color = met_survival), size = 2)+  scale_color_gradient2(low = "grey", mid = "yellow", high = "red", midpoint = mean(my_umap.ann$met_survival,na.rm=T))

ggplot(my_umap.ann, aes(V1, V2)) +
  geom_point(aes(color = rPFS_time), size = 2)+  scale_color_gradient2(low = "grey", mid = "yellow", high = "red", midpoint = mean(my_umap.ann$rPFS_time,na.rm=T))

dev.off()

write.table(my_umap,"umap_feature_128.txt",sep='\t',quote=F)




library(ggplot2)
library(umap)



my_umap <- as.data.frame(umap(t(feature_128[-1]))$layout)

graph <- bluster::makeSNNGraph(feature_128)
#graph <- bluster::makeSNNGraph(feature,k=10)
#graph <- bluster::makeSNNGraph(feature,k=5)
#> Warning in (function (to_check, X, clust_centers, clust_info, dtype, nn, :
#> detected tied distances to neighbors, see ?'BiocNeighbors-ties'
clust <- igraph::cluster_louvain(graph)

my_umap$clust <- factor(clust$membership)

my_umap.ann=unique(merge(my_umap,feature.meta[,c(1,228:ncol(feature.meta))],by.x=0,by.y=0))

my_umap.ann[,-c(2,3,5,10,46,56,58,60)] <- lapply(my_umap.ann[,-c(2,3,5,10,46,56,58,60)], function(x) as.character(x))
my_umap.ann$ARTERA.Score=as.numeric(my_umap.ann$ARTERA.Score)
my_umap.ann$met_survival=as.numeric(my_umap.ann$met_survival)
my_umap.ann$rPFS_time=as.numeric(my_umap.ann$rPFS_time)

pdf("umap.cluster.pdf",5,5)
ggplot(my_umap.ann, aes(V1, V2, colour = clust)) +
  geom_point()
ggplot(my_umap.ann, aes(V1, V2, colour = bin,shape=clust)) +
  geom_point()
ggplot(my_umap.ann, aes(V1, V2, colour = TP53,shape=clust)) +
  geom_point()
ggplot(my_umap.ann, aes(V1, V2, colour = TP53,shape=ADT..1.yes..0.no.)) +
  geom_point()
ggplot(my_umap.ann, aes(V1, V2)) +
  geom_point(aes(color = ARTERA.Score), size = 2)+  scale_color_gradient2(low = "grey", mid = "yellow", high = "red", midpoint = 0.1602506)

ggplot(my_umap.ann, aes(V1, V2)) +
  geom_point(aes(color = met_survival), size = 2)+  scale_color_gradient2(low = "grey", mid = "yellow", high = "red", midpoint = mean(my_umap.ann$met_survival,na.rm=T))

ggplot(my_umap.ann, aes(V1, V2)) +
  geom_point(aes(color = rPFS_time), size = 2)+  scale_color_gradient2(low = "grey", mid = "yellow", high = "red", midpoint = mean(my_umap.ann$rPFS_time,na.rm=T))

dev.off()



###########################wilcoxin with subfeature##############################
                                                 
feature=read.csv("/Users/ysong/Desktop/P_Tran_ref/score/to_use/score_111/image_predict_score/MMAI_score_feature_128_96.mod.csv",row.names=2)
mutation=read.table("/Users/ysong/Desktop/P_Tran_ref/score/to_use/score_111/image_predict_score/mutation.169.txt",sep="\t",header=T,check.names=F)
mutation[is.na(mutation)]<-0
#data=read.csv("/Users/ysong/Desktop/P_Tran_ref/score/to_use/score_111/image_predict_score/raw_score.csv")
colnames(mutation)=gsub("-T","",colnames(mutation))



meta.umap=read.delim2("/Users/ysong/Desktop/P_Tran_ref/score/to_use/score_111/image_predict_score/meta.mutation.numeric.txt",sep="\t",header=T)


feature.meta.sub=merge(feature,meta.umap,by.y="predictive_score",by.x="stadt_predictive_v1_3_0")

rownames(feature.meta.sub)=feature.meta.sub$Tumor_Sample_Barcode

feature.meta.sub.sub=feature.meta.sub[which(feature.meta.sub$Time..2.synchronous..1.metachronous.==1),]

feature_128=feature.meta.sub.sub[c(3:131)]


library(plotly);library(umap)

library(bluster)

iris.umap = umap::umap((feature_128.f))
my_umap0 <- as.data.frame(umap::umap((feature_128.f[-1]))$layout)





input=feature.meta.sub.sub[,c(4:131,284:399)]

input1=input[!grepl("OM",rownames(input)),]

input=input1

library(dplyr)

t1=data.frame(t(input[129:ncol(input)] %>% summarise_all(n_distinct)))
names(t1)="count"
t1$sample=rownames(t1)

t1f=data.frame(t1[which(t1$count>1),])

t1a=input[which(colnames(input)%in%t1f$sample)]

t1b=input[1:128]

t1ab=cbind(t1a,t1b)

#t1ab.f=t1ab[,-c(36,46:48)]

#t1ab=t1ab.f

t1ab.sub=t1ab[,c(1:95)]

t1ab.sub[is.na(t1ab.sub)]<-0
number_mutation=colSums(t1ab.sub !=0)
write.table(data.frame(number_mutation),"number_mutation_k3.txt",sep="\t")
t1ab=cbind(t1ab.sub,t1b)

y <- data.frame()

for (i in 1:95){
  
  for (j in 96:ncol(t1ab)){
    tmp=t1ab[,c(i,j)]
    
    histdata = tmp[which(tmp[,1]==1),2]
    hordata = tmp[which(tmp[,1]==0),2]
    if(length(histdata)>=3){
      
      Wilcoxson=wilcox.test(histdata,hordata,na.rm=TRUE,p.adjust.method = "fdr")
      Wilcox_P <- data.frame(paste(colnames(t1ab)[i],colnames(t1ab)[j],sep="_"),Wilcoxson$p.value)
      
      y <- rbind.data.frame(y,Wilcox_P)
      
      #write.table(pvalues1,paste(colnames(tmp)[1],colnames(tmp)[2],sep="_"),sep="\t")
    }
  }}
y$adj.p=p.adjust(y$Wilcoxson.p.value,method="fdr")



colnames(y)=c("comp","pval","adj.p")



foo <- data.frame(do.call('rbind', strsplit(as.character(y$comp),'_',fixed=TRUE)))

#y=read.delim("~/Desktop/pval/all_p.f.txt",sep="\t",header=T)

y1=cbind(y,foo)
setwd("~/Desktop/umap_plot_wilcoxin_metachronous.k3.noGhent//")

write.table(y1,"adjp_prog_feature.cluster.k3.v2.txt",sep="\t",quote=F,row.names=F)

y1.f=y1[which(y1$pval<0.05),]


library(reshape2)

y.f2=dcast(X2~X1, data = y1.f, value.var = "pval",fill="NA")


my_umap.ann=merge(my_umap0,y.f2,by.x=0,by.y="X2",all.x=T)

write.table(my_umap.ann,'my_umap.mutation.k2.dup.txt',sep="\t",quote=F)

my_umap.ann=read.table('my_umap.mutation.k2.dup.txt',sep="\t",header=T)

my_umap.ann[-c(1:6)] <- lapply(my_umap.ann[-c(1:6)], function(x) as.numeric(x))

my_umap.ann$count=rowSums(!is.na(my_umap.ann[-c(1:6)]))

my_umap.ann=my_umap.ann[order(my_umap.ann$count,decreasing=T),]

write.table(my_umap.ann,"~/Desktop/umap_plot_wilcoxin_metachronous.k3.noGhent//my_umap.mutation.k3.1.txt",sep="\t",quote=F,row.names=F)


#/Users/ysong/Desktop/umap_plot_wilcoxin_metachronous.k3/RAD51D.pdf

setwd("~/Desktop/umap_plot_wilcoxin_metachronous.k3.noGhent//")
my_umap.ann=read.table("~/Desktop/umap_plot_wilcoxin_metachronous.k3.noGhent/my_umap.mutation.k3.1.txt",sep="\t",header=T)

rownames(my_umap.ann)=my_umap.ann$Row.names
my_umap.ann[-c(1:6)] <- lapply(my_umap.ann[-c(1:6)], function(x) as.numeric(x))


my_umap.ann$clust=as.factor(my_umap.ann$clust)

for( i in 7:34){
  #my_umap.ann[i]=as.character(my_umap.ann[i])
  
  
  ggplot(my_umap.ann, aes(V1, V2,shape=clust)) +
    geom_point(aes(color=my_umap.ann[,i]))+
    scale_colour_gradient2(limits=c(min(my_umap.ann[,i]),
                                            max(my_umap.ann[,i])), high="black",na.value ="ivory3",low="red",midpoint=0.05,mid="blue" )+guides(color=guide_legend(as.character(colnames(my_umap.ann[i]))))

  ggsave(paste(colnames(my_umap.ann)[i],".pdf",sep=""),width=5,height=5)}



adjp=read.table("/Users/ysong/Desktop/umap_plot_wilcoxin_metachronous.k3.noGhent//adjp_prog_feature.cluster.k3.v2.filter.txt",sep="\t",header=T)

my_umap.ann=read.table("~/Desktop/umap_plot_wilcoxin_metachronous.k3.noGhent//my_umap.mutation.k3.1.txt",sep="\t",header=T,row.names=1)

library(circlize)

set.seed(200)
#rownames(my_umap.ann)=my_umap.ann$Row.names
colnames(my_umap.ann)[3]="clust"
mutation_sub_feature=my_umap.ann[order(my_umap.ann$clust),]
mutation_sub_feature1=mutation_sub_feature[-c(1:2)]
mutation_sub_feature1=mutation_sub_feature1[order(match(colnames(mutation_sub_feature1),adjp$comp))]


mutation_sub_feature1$clust=as.factor(mutation_sub_feature1$clust)
ha1=HeatmapAnnotation(cluster=mutation_sub_feature1$clust,col = list(cluster=c("1"="blue","2"="red","3"="yellow","4"="green","5"="pink")))

adjp1=unique(adjp[,c(4,6)])
adjp1=adjp1[order(adjp1$function.),]
#adjp1=adjp1[-22,]
adjp1$function.=as.factor(adjp1$function.)
library(dichromat );library(RColorBrewer)
ttxx <- colorRampPalette(brewer.pal(12,"Paired"))(length(levels(adjp1$function.)))
names(ttxx) <- levels(adjp1$function.)
#ha_row <- HeatmapAnnotation(pathway=data.frame(adjp1$function.), which="row", col = list(pathway = ttxx))




sub1=t(mutation_sub_feature1[,-c(1:3,32)])
sub1=sub1[order(match(rownames(sub1),adjp1$mutation)),]
sub1t=t(sub1)
colnames(my_umap.ann)[5]="cluster3"
sub1t.ann=merge(my_umap.ann[,c(1:2,5)],sub1t,by=0)

data<-replace(data.frame(lapply(sub1t.ann[-c(1:5)], as.character), stringsAsFactors = FALSE),
              !is.na(sub1t.ann[-c(1:5)]), "1")
data[is.na(data)]<-0
data1=cbind(sub1t.ann[1:4],data)


col0 = c("1" = "red","2" = "red","3" = "red","4"="red","5"="red","6"="red")

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "white", col = NA))
  },
  # big red
  
  # big red
  "1" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col0["1"], col = NA))
    
  },
  # big red
  "0" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "grey", col = NA))
    
  },
  # big red
  "3" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col0["3"], col = NA))
    
  },
  "4" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col0["4"], col = NA))
    
  },
  "5" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col0["5"], col = NA))
    
  },
  
  # big red
  "6" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col0["6"], col = NA))
    
  })
rownames(data1)=data1$Row.names
data1[data1==0]<-""

a1=data1[which(data1$clust==1),]
a2=data1[which(data1$clust==2),]
#a3=data1[which(data1$clust==3),]
#a4=data1[which(data1$clust==4),]
data2=cbind(t(a1[-c(1:4)]),t(a2[-c(1:4)]))

#data2=cbind(t(a1[-c(1:4)]),t(a2[-c(1:4)]),t(a3[-c(1:4)]))
col = c("1" = "red","2" = "red","3" = "red","4"="red","5"="red","6"="red")

a1t=t(a1[-c(1:4)])
a2t=t(a2[-c(1:4)])
#a3t=t(a3[-c(1:4)])
#a4t=t(a4[-c(1:4)])
colnames(adjp1)[1]="mutation"

a1t.ann=merge(adjp1,a1t,by.x="mutation",by.y=0)
a2t.ann=merge(adjp1,a2t,by.x="mutation",by.y=0)
#a3t.ann=merge(adjp1,a3t,by.x="mutation",by.y=0)
#a4t.ann=merge(adjp1,a4t,by.x="mutation",by.y=0)

rownames(a1t.ann)=paste(a1t.ann$function.,a1t.ann$mutation,sep="--")
rownames(a2t.ann)=paste(a2t.ann$function.,a2t.ann$mutation,sep="--")
#rownames(a3t.ann)=paste(a3t.ann$function.,a3t.ann$mutation,sep="--")
#rownames(a4t.ann)=paste(a4t.ann$function.,a4t.ann$mutation,sep="--")

a1t.ann=a1t.ann[order(rownames(a1t.ann)),]
a2t.ann=a2t.ann[order(rownames(a2t.ann)),]
#a3t.ann=a3t.ann[order(rownames(a3t.ann)),]
#a4t.ann=a4t.ann[order(rownames(a4t.ann)),]


pct_num1 = data.frame(rowSums(a1t.ann[,-c(1:2)]==1) / ncol(a1t.ann[-c(1:2)]))*100
pct_num2 = data.frame(rowSums(a2t.ann[,-c(1:2)]==1) / ncol(a2t.ann[-c(1:2)]))*100

#pct_num3 = data.frame(rowSums(a3t.ann[,-c(1:2)]==1) / ncol(a3t.ann[-c(1:2)]))*100

#pct_num4 = data.frame(rowSums(a4t.ann[,-c(1:2)]==1) / ncol(a4t.ann[-c(1:2)]))*100
names(pct_num1)="perc_cluster1"
names(pct_num2)="perc_cluster2"
#names(pct_num3)="perc_cluster3"
#names(pct_num4)="perc_cluster4"

pct_all=cbind(pct_num1,pct_num2)

#pct_all=cbind(pct_num1,pct_num2,pct_num3)
p0=""
for (i in 1:nrow(pct_all)){
  wnt=pct_all[i,]
  wnt.no=100-t(wnt)
  wnt.no1=cbind(wnt.no,t(wnt))
  p=chisq.test(wnt.no1)$p.value
  names(p)=rownames(pct_all)[i]
  print(p)
 
}
p=read.table("~/Desktop/umap_plot_wilcoxin_metachronous.k3.noGhent//chiq_test.k2.txt",sep="\t",header=F)

p$adjp=p.adjust(p$V2)

#pct = paste0(round(pct_num * 100, digits = 2), "%")
write.table(p,"chiq_test.adjp.txt",sep="\t",quote=F)


pdf("~/Desktop/mutation.cor.k2..new.pdf",35,10)

t01=oncoPrint(a1t.ann[-c(1:2)], remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col0,row_order=rownames(a1t.ann),
              row_title_gp = gpar(fontsize = 6) ,show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 10),column_title="cluster1")
#t01

t02=oncoPrint(a2t.ann[-c(1:2)], remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col0,row_order=rownames(a2t.ann),
              row_title_gp = gpar(fontsize = 6) ,show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 10),column_title="cluster2")
#t01
#t03=oncoPrint(a3t.ann[-c(1:2)], remove_empty_columns = FALSE, remove_empty_rows = F,
 #             alter_fun = alter_fun, col = col0,row_order=rownames(a3t.ann),
  #            row_title_gp = gpar(fontsize = 6) ,show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 10),column_title="cluster3")
#t01

t01+t02

dev.off()


#################################

sub1[is.na(sub1)]<-1


sub1.ann=merge(adjp1,sub1,by.x="X1",by.y=0)

rownames(sub1.ann)=paste(sub1.ann$function.,sub1.ann$X1,sep="__")

sub1.ann.s=sub1.ann[order(match(sub1.ann$X1,adjp1$X1)),]


ha2=Heatmap(adjp1$function., name = "pathway",  width = unit(1, "cm"),col=ttxx)


p1=Heatmap(sub1.ann.s[-c(1:2)],na_col = "grey",col = colorRamp2(c( 0, 0.05,0.055), c("red","yellow", "black")),name=" ",column_order = colnames(sub1.ann.s[-c(1:2)]),row_order=rownames(sub1.ann.s[-c(1:2)]),top_annotation = ha1,column_names_gp=gpar(fontsize=8))

draw(ha2+p1,heatmap_legend_side="bottom", annotation_legend_side="bottom")
#p2=Heatmap(sub1,na_col = "grey",col = colorRamp2(c( 0, 0.05,0.055), c("red","yellow", "black")),name=" ",column_order = colnames(sub1),top_annotation = ha1,column_names_gp=gpar(fontsize=8))

#draw(p2+ha2,heatmap_legend_side="bottom", annotation_legend_side="bottom")



p2=Heatmap(sub1.ann.s[-c(1:2)],na_col = "grey",col = colorRamp2(c( 0, 0.05,0.055), c("red","yellow", "black")),name=" ",column_order = colnames(sub1.ann.s[-c(1:2)]),top_annotation = ha1,column_names_gp=gpar(fontsize=8))
draw(p2,heatmap_legend_side="right", annotation_legend_side="right")

p3=Heatmap(sub1.ann.s[-c(1:2)],na_col = "grey",name=" ",top_annotation = ha1,column_names_gp=gpar(fontsize=8),col = colorRamp2(c( 0, 0.05,0.055), c("red","yellow", "black")),)

draw(p3,heatmap_legend_side="right", annotation_legend_side="right")

dev.off()
