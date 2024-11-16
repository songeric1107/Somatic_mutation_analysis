clinical=read.delim("/Users/ysong/Desktop/P_Tran_ref/pdac_background/pdac_data_tcga/paad_cptac_2021/data_clinical_patient.txt",sep="\t",header=T)
clinical1=read.delim("/Users/ysong/Desktop/P_Tran_ref/pdac_background/pdac_data_tcga/paad_cptac_2021/paad_cptac_2021_clinical_data_purity.txt",sep="\t",header=T,check.names=F)
cliniclal2=merge(clinical,clinical1,by.x="PATIENT_ID",by.y="Patient ID")
clinical3=cliniclal2[which(cliniclal2$purity=="high purity"),]



marker1=read.delim("/Users/ysong/Desktop/P_Tran_ref/pdac_background/Biomarkers Papers/motiffit.marker.txt",sep="\t",header=F)

marker2=read.delim("/Users/ysong/Desktop/P_Tran_ref/pdac_background/Biomarkers Papers/collisson_marker.txt",sep="\t",header=F)

marker3=read.delim("/Users/ysong/Desktop/P_Tran_ref/pdac_background/Biomarkers Papers/puleo/Puleo_readable_list.tsv",sep="\t",header=T)



marker4=read.delim("/Users/ysong/Desktop/P_Tran_ref/pdac_background/Biomarkers Papers/michelle.txt",sep="\t",header=F)
marker5=read.delim("/Users/ysong/Desktop/P_Tran_ref/pdac_background/Biomarkers Papers/bailey/Bailey_readable_list.tsv",sep="\t",header=F)
marker5.uniq=marker5[which(marker5$V2!="not unique"),]

all_cluster=read.table("/Users/ysong/Desktop/P_Tran_ref/analysis_2024/ajmal_rnaseq/human_pdac/result_321/all_cluster.3.txt",sep="\t",header=T)
rownames(all_cluster)=gsub("[.]","-",rownames(all_cluster))
all_markers=unique(c(marker1$V1,marker2$V1,marker5.uniq$V1))

all_cluster$motiff=gsub("2","classical",   all_cluster$motiff)
all_cluster$motiff=gsub("1","Basal-like",   all_cluster$motiff)
# all_cluster$motiff=gsub("3","na",   all_cluster$motiff)


all_cluster$collisson=gsub("1","QM-PDA",   all_cluster$collisson)
all_cluster$collisson=gsub("2","classical",   all_cluster$collisson)
all_cluster$collisson=gsub("3","exocrine-like",   all_cluster$collisson)
# all_cluster$collisson=gsub("4","NA",   all_cluster$collisson)


all_cluster$Bailey=gsub("1","Squamous",    all_cluster$Bailey)
all_cluster$Bailey=gsub("2","Progenitor",    all_cluster$Bailey)
all_cluster$Bailey=gsub("3","ADEX",    all_cluster$Bailey)

all_cluster$Bailey=gsub("4","Immunogenic",    all_cluster$Bailey)

#all_cluster$Bailey=gsub("5","NA",    all_cluster$Bailey)


all_markers=unique(c(marker1$V1,marker2$V1,marker5.uniq$V1))

ref.sub=s1[which(s1$Group.1%in%all_markers),]


pdf_sub1 <- as.data.frame(t(scale(t(ref.sub[-1]), center = T, scale = T)))
smallx=pdf_sub1
colnames(smallx)=gsub("[.]","-",colnames(smallx))



data3=read.table("/Users/ysong/Desktop/P_Tran_ref/pdac_background/pdac_data_tcga/paad_cptac_2021/data_mrna_seq_v2_rsem.txt",header=T,sep="\t",check.names=F)

s1=aggregate( data3[-c(1:2)], by=list( data3$Hugo_Symbol), FUN=mean)
rownames(s1)=s1$Group.1


s1.log=na.omit(s1[-1])

#s1f=normalize.quantiles(as.matrix(na.omit(s1[-1])))
#rownames(s1f)=rownames(na.omit(s1[-1]))
#colnames(s1f)=colnames(na.omit(s1[-1]))


all_markers=unique(c(marker1$V1,marker2$V1,marker5.uniq$V1))

ref.sub=s1.log[which(rownames(s1.log)%in%all_markers),]


ref.sub1=ref.sub[,which(colnames(ref.sub)%in%clinical3$PATIENT_ID)]



#pdf_sub1 <- as.data.frame(t(scale(t(ref.sub[-1]), center = T, scale = T)))
pdf_sub2 <- as.data.frame(t(scale(t(ref.sub1), center = T, scale = T)))

smallx=pdf_sub2
colnames(smallx)=gsub("[.]","-",colnames(smallx))
twist1=data.frame(t(smallx[which(rownames(smallx)=="TWIST1"),]))


#s1.norm=(as.matrix(s1f-min(s1f)))
#s1.log$zero_count=rowSums( s1.log== 0)

all_cluster.ann=merge(all_cluster,cliniclal2,by.x="X",by.y="PATIENT_ID")
all_cluster.ann$status=as.numeric(all_cluster.ann$status)


all.cluster.annf=all_cluster.ann[order(match(all_cluster.ann$X,colnames(smallx))),]
all.cluster.annf$TUMOR_NECROSIS=as.factor(all.cluster.annf$TUMOR_NECROSIS)
all.cluster.annf$all_cluster=paste("NMF",all.cluster.annf$all_cluster,sep="")
all.cluster.annf$motiff=gsub("classical","Classical",all.cluster.annf$motiff)
all.cluster.annf$Bailey=gsub("classical","Classical",all.cluster.annf$Bailey)
all.cluster.annf$collisson=gsub("classical","Classical",all.cluster.annf$collisson)
all.cluster.annf$collisson=gsub("exocrine-like","Exocrine-like",all.cluster.annf$collisson)

set.seed(10)
myColors<-c("deeppink3","chartreuse4","mediumorchid1","darkblue","darkorange","lightcoral","mediumorchid4","darkseagreen","orangered","cyan","pink","blue","grey","red","orange")

labels=c("166 gene","56 gene","111 gene")


all.cluster.annf$motiff=gsub("classical","Classical",all.cluster.annf$motiff)
all.cluster.annf$collisson=gsub("classical","Classical",all.cluster.annf$collisson)
all.cluster.annf$collisson=gsub("exocrine-like","Exocrine-like",all.cluster.annf$collisson)

set.seed(10)

twist1=data.frame(t(ref.sub1[grep("TWIST1",rownames(ref.sub1)),]))


pdf("fig1.human.pdac_1018.pdf",10,10)
all.cluster.annf=all.cluster.annf[order(all.cluster.annf$all_cluster),]
twist1$sample=rownames(twist1)
twist1=data.frame(twist1[order(match(twist1$sample,all.cluster.annf$X)),])

ha01 = HeatmapAnnotation(  # empty = anno_empty(border = FALSE),
                            raw_cluster1 = anno_block(gp = gpar(fill = c("darkorange","yellow2","antiquewhite4")), labels = c("NFM1","NMF2","NMF3")),
                       #  raw_cluster=all.cluster.annf$all_cluster, 
                       moffit= all.cluster.annf$motiff,
                         collisson= all.cluster.annf$collisson,
                         bailey=all.cluster.annf$Bailey,TWIST1 =twist1$TWIST1,
                         annotation_label = c( "NMF_cluster", "Moffit et al. (2015)", "Collison et al. (2011)","Bailey et al. (2016)","TWIST1"),
                         col=list(#raw_cluster=c("NMF3" = "antiquewhite4","NMF2" = "yellow2","NMF1"="darkorange"),
                                  
                                  moffit=c("Classical" = "purple","Basal-like" = "darkorange","na" = "grey"),
                                  collisson=c("Classical" = "magenta","Exocrine-like" = "cyan","QM-PDA" = "darkorange","NA"="grey"),
                                  bailey=c("ADEX" = "green","Progenitor" = "brown","Squamous" = "darkorange","Immunogenic"="lightgrey","NA"="grey")))



smallx=smallx[,order(match(colnames(smallx),all.cluster.annf$X))]


Heatmap(smallx,name=" ", column_order=all.cluster.annf$X,column_split=all.cluster.annf$all_cluster,top_annotation=ha01,show_column_dend = F,
        show_row_dend = F,row_names_gp = gpar(fontsize = 2),row_gap =unit(3, "mm"), column_title = " ",
        show_column_names = F,show_row_names = F,column_dend_height = unit(20, "mm"),
        row_km=3,row_title = labels,row_dend_reorder = T)


dev.off()

nmf_marker=readRDS("/Users/ysong/Desktop/P_Tran_ref/analysis_2024/ajmal_rnaseq/human_pdac/result_321/nmf.all.marker.ref1.rds")



# Set up the NMF result
raw.nmf3 <- nmf_marker$fit$`3`

cophenetic_scores <- nmf_marker$measures$cophenetic
silhouette_scores <- nmf_marker$measures$silhouette.consensus


# Extract the Cophenetic Index and Silhouette Score from nmf_result

ranks <- nmf_marker$r

# Create a data frame for easy plotting
df <- data.frame(Rank = ranks, Cophenetic = cophenetic_scores, Silhouette = silhouette_scores)

# Load ggplot2 for custom plotting
library(ggplot2)

# Create the plot showing only Cophenetic Index and Silhouette Score
pdf("nmf_eval.pdf")
ggplot(df, aes(x = Rank)) +
  geom_line(aes(y = Cophenetic, color = "Cophenetic Index"), size = 1.2) +
  geom_line(aes(y = Silhouette, color = "Silhouette Score"), size = 1.2) +
  labs(y = "Scores", title = "Cophenetic Index and Silhouette Score for NMF Results") +
  theme_minimal() +
  scale_color_manual("", values = c("Cophenetic Index" = "blue", "Silhouette Score" = "red")) +
  theme(legend.position = "top")
dev.off()






# Increase the size of the PDF output to fit the plot and labels
pdf('metagenes.cluster3.updated.pdf', width = 10, height = 10)  # Adjust width and height

# Optionally adjust margins to give more space around the plot
par(mar = c(5, 5, 4, 1) + 0.1)  # Adjust as needed (bottom, left, top, right margins)

# Generate the consensus map with label size adjustments
xx1 = consensusmap(raw.nmf3,tracks = c('consensus:', 'silhouette:'),
                   annCol = NA,  # Optional: Disable column annotations
                   annRow = NA,  # Optional: Disable row annotations
                   fontsize = 16)  # Increase the font size for both row and column labels
#si1 <- silhouette(raw.nmf3, what = 'consensus')
#plot(si1, cex = 1.5)  # Increase size of plot text and symbols

dev.off()

# Close the PDF device
dev.off()


pdf('metagenes.cluster2.pdf',5,5)
raw.nmf3 <- nmf_marker$fit$`3`


#http://nmf.r-forge.r-project.org/_DEMOS.html
# Increase row and column label size in consensusmap
 xx1 = consensusmap(raw.nmf3, annCol = NA,annRow = NA,  # Optional: Disable column annotations
                    fontsize = 20) 
                                     # Optional: Disable row annotations
                                       # Optional: Disable column annotations
                                      # Adjust the font size for both row and column labels # Adjust the font size for both row and column labels
# Increase row and column label size in coefmap
yy1 = coefmap(fit(raw.nmf3), 
              annCol = NA,annRow = NA,  # Optional: Disable column annotations
              fontsize = 20)  # Adjust the font size for both row and column labels
# Increase row and column label size in coefmap
# Increase the size of silhouette plot
si1 <- silhouette(raw.nmf3, what = 'consensus')
plot(si1, cex = 1.5)  # Increase size of plot text and symbols

dev.off()
  
