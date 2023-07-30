

library(GenVisR);library(maftools)
setwd("/local/scratch/ysong/soma")
path="/local/scratch/ysong/soma"

##maftools read maf do not include synonymous mutation, only use non-synonymous
### add samples name to maf files.
file.names <- dir(path, pattern =".maf$")

out.file<-""
for(i in 1:length(file.names)){
  maf <- data.table::fread(file.names[i])
  maf$Tumor_Sample_Barcode=rep(file.names[i],nrow(maf))
  
  filename <- paste(file.names[i], ".mod.txt", sep="")
  
  write.table(maf, file=filename,sep="\t",quote=FALSE,row.names=F)}


file.names <- dir(path, pattern =".maf.mod.txt")

clinical_tmb=read.table("/local/scratch/ysong/soma/meta_no97.txt",sep="\t",header=T,check.names = F)
library(maftools)


###merge_maf files of SOMA samples



#setwd("/Users/ysong/Desktop/P_Tran_ref/maf/soma/")
#$clinical_tmb=read.table("/Users/ysong/Desktop/P_Tran_ref/onco_508/soma/meta_no97.txt",sep="\t",header=T,check.names = F)

path=getwd()

file.names <- dir(path, pattern =".maf.mod.txt")

## use the default variant classification

maf_2=maftools:::merge_mafs(file.names,clinicalData=clinical_tmb,
                            removeDuplicatedVariants=T)


clinical_maf=getClinicalData(maf_2)
anno_df=clinical_maf

saveRDS(maf_1,"maf_all_soma.rds")


###create onco_matrix as input of  oncoprint

pdf("soma1.pdf")

oncoplot(maf = maf_2,top=nrow(maf_2@data),writeMatrix = TRUE) 
dev.off()   


p2=plotmafSummary(maf = maf_2, rmOutlier = FALSE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE,top=nrow(maf_2@data))



################################


## use the default variant classification, but do not remove duplicate variant

maf_1=maftools:::merge_mafs(file.names,clinicalData=clinical_tmb,
                            removeDuplicatedVariants=FALSE)

clinical_maf=getClinicalData(maf_1)
anno_df=clinical_maf

saveRDS(maf_1,"maf_all_soma_keep.dup.var.rds")



## use defined  classification
maf_3=maftools:::merge_mafs(file.names,clinicalData=clinical_tmb,
                            vc_nonSyn=c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation","3'Flank",
"3'UTR", "5'Flank","5'UTR","Intron","RNA","Silent","IGR","Splice_Region"))
saveRDS(maf_3,"/Users/ysong/Desktop/P_Tran_ref/onco_508/maf_all_.var.515.rds")

maf_3=readRDS("/Users/ysong/Desktop/P_Tran_ref/onco_508/maf_all_.var.515.rds")
all_stat=maftools:::createOncoMatrix(maf_2,g=maf_2@gene.summary$Hugo_Symbol)

all_stat.matrix=all_stat$numericMatrix



pdf("soma3.pdf")

oncoplot(maf = maf_3,top=nrow(maf_3@data),writeMatrix = TRUE) 
dev.off() 


########somatic
maf_2=readRDS("~/Desktop/P_Tran_ref/onco_508/soma/maf_all_soma.rds")

pdf("somatic_interactionpldf")
somaticInteractions(maf = maf_2, top = 25, pvalue = c(0.05, 0.1))
dev.off()


#soma.sig = oncodrive(maf = maf_2, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')

#prog_geneset = survGroup(maf = maf_3, top = 20, geneSetSize = 2, time = "CRPC_from_MDT_(1=yes,_0=no)", Status = "OS_(1=dead,_0=alive)", verbose = FALSE)




fab.ce = clinicalEnrichment(maf = maf_2, clinicalFeature = 'OS_(1=dead,_0=alive)')

pair=fab.ce$groupwise_comparision[p_value < 0.05]


pdf("enrich_515.pdf")

plotEnrichmentResults(enrich_res = fab.ce, pVal = 0.05)
dev.off()

pdf("summary.515.pdf")
plotmafSummary(maf = maf_3, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

plotmafSummary(maf = maf_3, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE, top=30)
dev.off()




onco_mat=read.delim("/Users/ysong/Desktop/P_Tran_ref/maf/soma/onco_matrix_515.txt",sep="\t",header=T,row.names=1,check.names=F)


gene=read.delim("/Users/ysong/Desktop/P_Tran_ref/maf/soma/gene_variant.all_variant.515.txt",sep="\t",header=T,check.names=F)

gene.m=melt(gene)

top53=gene.m[which(gene.m$Hugo_Symbol=="TP53"),]
myColors<-c("deeppink3","chartreuse4","mediumorchid1","darkblue","darkorange","lightcoral","mediumorchid4","darkseagreen","orangered","cyan","black","pink","blue","grey","red","orange","yellow")

pdf("stat_gene_515.pdf")
ggplot(top53,aes(x=as.factor(Hugo_Symbol),y=value,fill=variable))+    geom_bar(position="dodge", stat="identity")+facet_wrap(~variable,ncol=6) +
  theme(legend.position="bottom")+xlab("")

ggplot(top53,aes(x=as.factor(Hugo_Symbol),y=value,fill=variable))+    geom_bar(position="dodge", stat="identity")+
  theme(legend.position="bottom")+xlab("")


dev.off()




anno_df=data.frame(getClinicalData(maf_2))

col = c("Missense_Mutation" = "blue", "In_Frame_Ins" = "red", "Splice_Site" = "green","In_Frame_Del"="purple","Frame_Shift_Ins"="yellow",
        "Frame_Shift_Del"="grey", "Nonsense_Mutation"="brown", "Multi_Hit"="orange","Translation_Start_Site"="pink","Nonstop_Mutation"="cyan",
        "3'Flank"="goldenrod4", "3'UTR"= "aquamarine4","5'Flank" ="lightslateblue",               "5'UTR"="lightcoral","IGR"="black" )


#col = c("Missense_Mutation" = "blue", "In_Frame_Ins" = "red", "Splice_Site" = "#008000","In_Frame_Del"="purple","Frame_Shift_Ins"="yellow",
#       "Frame_Shift_Del"="grey", "Nonsense_Mutation"="brown", "Multi_Hit"="orange","Translation_Start_Site"="pink","Nonstop_Mutation"="cyan")

rate=tmb(maf = maf_2)

foo1 <- data.frame(do.call('rbind', strsplit(as.character( rate$Tumor_Sample_Barcode),'.',fixed=TRUE)))
foo2=data.frame(do.call('rbind',strsplit(as.character(foo1$X1),'_',fix=TRUE)))

rate$index=foo2$X1

anno_df1=merge(anno_df,rate,by="Tumor_Sample_Barcode")


#prepare the topbar
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big blue
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["Missense_Mutation"], col = NA))
  },
  # big red
  In_Frame_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["In_Frame_Ins"], col = NA))
  },
  # small green
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),  
              gp = gpar(fill = col["Splice_Site"], col = NA))},
  
  # small green
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["In_Frame_Del"], col = NA))
  },
  
  Frame_Shift_Ins= function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["Frame_Shift_Ins"], col = NA))},
  
  Frame_Shift_Del= function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["Frame_Shift_Del"], col = NA))},
  Nonsense_Mutation= function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["Nonsense_Mutation"], col = NA))},
  Multi_Hit= function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["Multi_Hit"], col = NA))},
  Translation_Start_Site= function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["Translation_Start_Site"], col = NA))},
  Nonstop_Mutation= function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["Nonstop_Mutation"], col = NA))}
  
)


library(ComplexHeatmap)
library(ComplexHeatmap)
column_title = "OncoPrint"

anno_df1$index.y <- substr(anno_df1$index.y, 0, 2)

#anno_df[c(5,13,17,22,24,26,28,30,32,34,36)] <- data.frame(lapply( anno_df[2:11], factor))

i=c(8,16,20,25,27,29,31,33,35,37,39,41:44)

#any(round(anno_df[,1:ncol(anno_df)]) != anno_df[,1:ncol(anno_df)])
anno_df1[ , i] <- apply(anno_df1[ , i], 2,            # Specify own function within apply
                        function(x) as.numeric(as.character(x)))

anno_df1=anno_df1[order(anno_df1$total_perMB,decreasing=T),]
rownames(anno_df1)=anno_df1$Tumor_Sample_Barcode

ha.top = HeatmapAnnotation(
  dist1 = anno_barplot(
    anno_df1$total_perMB, 
    bar_width = 1, 
    gp = gpar(col = "white", fill = "#FFE200"), 
    border = FALSE,
    
    height = unit(2, "cm")
  ), show_annotation_name = FALSE)


tmb=anno_df1$total_perMB
OS.from.MDT = anno_df1$OS.from.MDT
age = anno_df1$Age
alive=anno_df1$OS..1.dead..0.alive.
Distant.Failure.from.MDT=anno_df1$Distant.Failure.from.MDT

Time.to.Distant.Failure.from.oligomet=anno_df1$Time.to.Distant.Failure.from.oligomet
Post.MDT.Local.Failure..time.from.MDT.to.local.failure=anno_df1$Post.MDT.Local.Failure..time.from.MDT.to.local.failure.

PSA.Failure.time.from.MDT=anno_df1$PSA.Failure.time.from.MDT
PSA.failure.time.from.oligomet  =anno_df1$PSA.failure.time.from.oligomet
location=anno_df1$Location..1.node..2.bone..3.node.bone..visceral.4..visceral.other.5

Imaging.Type=anno_df1$Imaging.Type..1.convention..2.enhanced.
Number.of.mets =anno_df1$Number.of.mets
PSA.at.oligomet =anno_df1$PSA.at.oligomet
last.FU.from.oligo=anno_df1$Last.FU.from.oligomet
Risk.Group=anno_df1$Risk.Group
index=anno_df1$index
Primary.Gleason =anno_df1$Primary.Gleason
CRPC.from.MDT=anno_df1$CRPC.from.MDT

# ha2=HeatmapAnnotation(anno_df,OS.from.MDT = anno_df$OS.from.MDT,cbar = anno_oncoprint_barplot(), age = anno_df$Age,alive=anno_df$OS..1.dead..0.alive.)

ha1 = HeatmapAnnotation(tmb=anno_df1$total_perMB,OS.from.MDT = anno_df1$OS.from.MDT, age = anno_df1$Age,alive=anno_df1$OS..1.dead..0.alive.,Distant.Failure.from.MDT=anno_df1$Distant.Failure.from.MDT,
                        Time.to.Distant.Failure.from.oligomet=anno_df1$Time.to.Distant.Failure.from.oligomet,Post.MDT.Local.Failure..time.from.MDT.to.local.failure=anno_df1$Post.MDT.Local.Failure..time.from.MDT.to.local.failure.,
                        PSA.Failure.time.from.MDT=anno_df1$PSA.Failure.time.from.MDT,PSA.failure.time.from.oligomet  =anno_df1$PSA.failure.time.from.oligomet,location=anno_df1$Location..1.node..2.bone..3.node.bone..visceral.4..visceral.other.5,
                        Imaging.Type=anno_df1$Imaging.Type..1.convention..2.enhanced.,Number.of.mets =anno_df1$Number.of.mets,PSA.at.oligomet =anno_df1$PSA.at.oligomet,last.FU.from.oligo=anno_df1$Last.FU.from.oligomet,
                        Risk.Group=anno_df1$Risk.Group,index=anno_df1$index,Primary.Gleason =anno_df1$Primary.Gleason,CRPC.from.MDT=anno_df1$CRPC.from.MDT,
                        tmb=anno_points(tmb,ylim = c(0, max(tmb, na.rm = TRUE))),
                        OS.from.MDT=anno_oncoprint_barplot(OS.from.MDT, ylim = c(0, max(OS.from.MDT, na.rm = TRUE)), axis = TRUE),
                        last.FU.from.oligo=anno_oncoprint_barplot(last.FU.from.oligo, ylim = c(0, max(last.FU.from.oligo, na.rm = TRUE)), axis = TRUE),
                        CRPC.from.MDT=anno_oncoprint_barplot(CRPC.from.MDT, ylim = c(0, max(CRPC.from.MDT, na.rm = TRUE)), axis = TRUE),
                        Distant.Failure.from.MDT=anno_oncoprint_barplot(anno_df1$Distant.Failure.from.MDT, ylim = c(0, max(Distant.Failure.from.MDT, na.rm = TRUE)), axis = TRUE),
                        Time.to.Distant.Failure.from.oligomet=anno_oncoprint_barplot(anno_df1$Time.to.Distant.Failure.from.oligomet,ylim = c(0, max(anno_df1$Time.to.Distant.Failure.from.oligomet, na.rm = TRUE)), axis = F),
                        Post.MDT.Local.Failure..time.from.MDT.to.local.failure = anno_oncoprint_barplot(anno_df1$Post.MDT.Local.Failure..time.from.MDT.to.local.failure., ylim = c(0, max(anno_df1$Post.MDT.Local.Failure..time.from.MDT.to.local.failure., na.rm = TRUE)), axis = TRUE),
                        alive=anno_oncoprint_barplot(alive,ylim = c(0, max(alive, na.rm = TRUE)), axis = F),
                        
                        age = anno_points(age, ylim = c(0, max(age))),
                        col = list(alive = c("1" = "red", "0" = "blue","NA"="grey")),
                        Imaging.Type = c("1"="red","2"="blue"),Number.of.mets=c("1"="red","2"="blue","3"="green","4"="pink"),
                        Risk.Group=c("FI"="red","H"="blue","Metastatic"="green","Node+"="pink","UI"="yellow","VH"="cyan"),
                        Primary.Gleason=c("3"="red","4"="blue","5"="green"),
                        annotation_height = unit(c(1, 1, 1,1), "cm"))


onco_mats=onco_mat[order(match(colnames(onco_mat),anno_df1$Tumor_Sample_Barcode))]





#gender = anno_df[, "gender"]




#anno_df.sort=anno_df
#onco_mat.sub=onco_mat[which(names(onco_mat)%in%anno_df$Tumor_Sample_Barcode)]
#onco_mat.sub2=onco_mat.sub[order(match(colnames(onco_mat.sub),anno_df.sort$Tumor_Sample_Barcode))]


#onco_mat.sub1=onco_mat.sub[order(match(colnames(onco_mat.sub),anno_df.sort$Tumor_Sample_Barcode))]




setwd("~/Desktop")

amp = ifelse(apply(onco_mats, 1, function(x) sum(grepl("Missense_Mutation", x))/length(x) > 0.1), "high Missense_Mutation events", "low Missense_Mutation events")
amp = factor(amp, levels = c("low Missense_Mutation events", "high Missense_Mutation events"))

gene1=data.frame(maf_2@gene.summary)

gene_plot1=gene1[which(gene1$MutatedSamples>32*.5),]

mat_gs1=onco_mats[which(rownames(onco_mats)%in%gene_plot1$Hugo_Symbol),]

onco_mat.sub1= mat_gs1[order(match(colnames( mat_gs1),anno_df1$Tumor_Sample_Barcode))]

sample_order1=as.character(colnames( onco_mat.sub1))

gene_to_plot=read.table("~/Desktop/P_Tran_ref/matation_to_plot.txt",sep="\t",header=T)

mat_gs=onco_mats[which(rownames(onco_mats)%in%gene_to_plot$name),]


onco_mat.sub2=mat_gs[order(match(colnames(mat_gs),anno_df1$Tumor_Sample_Barcode))]

#onco_mat.sub1= mat_gs1[order(match(colnames( mat_gs1),anno_df1$Tumor_Sample_Barcode))]
sample_order2=as.character(colnames( onco_mat.sub2))


pdf("tmp.pdf")
oncoplot(maf_2,top=10,writeMatrix = TRUE)
dev.off()

library("dplyr"); library(corrplot)

pdf("tmb.cor.v2.pdf",20,20)
sub=select_if(anno_df1, is.numeric)



M = cor(sub[,-c(13)],method="spearman",use="complete.obs")

colour_set <- colorRampPalette(colors = c("#f4ff4d", "#c7d123", "#acb515", "#81890b", "#656c06"))

corrplot(M, tl.col = "blue", bg = "White", tl.srt = 35, 
         title = "\n\n Correlation Plot\n",
         addCoef.col = "black", type = "lower",
         col = colour_set(100))


chart.Correlation(M, histogram=TRUE, pch=19,font.labels=16)

dev.off()



gene=read.table("~/Downloads/gene_mutation_summary.soma.txt",sep="\t",header=T)

gene.m=melt(gene[-c(12:13)])

top53=gene.m[which(gene.m$Hugo_Symbol=="TP53"),]

pdf("stat_gene_513.pdf")
ggplot(top53,aes(x=as.factor(Hugo_Symbol),y=value,fill=variable))+    geom_bar(position="dodge", stat="identity")+facet_wrap(~variable,ncol=4)
dev.off()


###tmb default size captureSize = 50



pdf("tmb.cor.pdf",20,20)
library(GGally)
lowerFn <- function(data, mapping, method = "lm", ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(colour = "black") +
    geom_smooth(method = method, color = "red", ...)
  p
}


ggpairs(sub[,-c(13)], lower=list(continuous = wrap(lowerFn, method = "lm")),
        upper = list(continuous = wrap("cor", size = 8)),
        diag=list(discrete= "barDiag"), 
        axisLabels='none') +
  theme_bw() +
  theme(strip.text.x = element_text(angle=90, hjust=1,size=20),
        strip.background = element_rect(fill = NA),
        strip.text.y = element_text(angle=90, hjust=1,size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()
stat=anno_df1[c(36,43)]
statm=melt(stat)
statm=na.omit(statm)

pdf("stat.pdf")
ggplot(statm,aes(x=as.factor(statm$CRPC.from.MDT..1.yes..0.no.),y=value))+geom_violin()+ geom_jitter(height = 0, width = 0.1)

dev.off()
#total_perMB

gene_to_plot1=read.table("~/Desktop/P_Tran_ref/gene_plot1.txt",sep="\t",header=T)

mat_gs=onco_mats[which(rownames(onco_mats)%in%gene_to_plot1$gene),]


onco_mat.sub3=mat_gs[order(match(colnames(mat_gs),anno_df1$Tumor_Sample_Barcode))]

#onco_mat.sub2= onco_mat.sub1[order(match(colnames( onco_mat.sub1),anno_df1$Tumor_Sample_Barcode))]
sample_order3=as.character(colnames( onco_mat.sub3))

#file:///Users/ysong/Downloads/ComplexHeatmap-supplementary1-4/supplS1_TCGA_OncoPrint/supplS1_OncoPrint.html
pdf("somatic_metastatic_gene_512.v2.pdf",paper="special",  pointsize=10,  width=39.7/2.54,height=30.9/2.54 * 1,pagecentre=T)
column_title = "Mutate in more than 70% of samples "
p1=oncoPrint(  onco_mat.sub1,col=col,
               alter_fun = alter_fun, get_type = function(x) strsplit(x, "\t") [[1]],column_order = sample_order1,
               remove_empty_columns = TRUE, remove_empty_rows = TRUE,
               column_title = column_title,bottom_annotation = ha1,pct_gp = gpar(fontsize = 8))
#,  top_annotation=ha.top)

draw(p1,heatmap_legend_side="bottom", annotation_legend_side="bottom")
column_title = "cancer mutant marker genes "
p=oncoPrint( onco_mat.sub2,col=col,
             alter_fun = alter_fun, get_type = function(x) strsplit(x, "\t") [[1]],column_order = sample_order2,
             remove_empty_columns = TRUE, remove_empty_rows = TRUE,
             column_title = column_title,bottom_annotation = ha1,pct_gp = gpar(fontsize = 8))
#,  top_annotation=ha.top)

draw(p,heatmap_legend_side="bottom", annotation_legend_side="bottom")

pdf("somatic_metastatic_gene_512.v2.pdf",paper="special",  pointsize=10,  width=39.7/2.54,height=50.9/2.54 * 1,pagecentre=T)
column_title = "cancer mutant top10 marker genes and dead/alive associated genes "
p=oncoPrint( onco_mat.sub3,col=col,
             alter_fun = alter_fun, get_type = function(x) strsplit(x, "\t") [[1]],column_order = sample_order3,
             remove_empty_columns = F, remove_empty_rows = TRUE,show_column_names = FALSE, 
             column_names_gp = gpar(cex=0.5, font=1, col= "blue"), 
             column_title = column_title,bottom_annotation = ha1,pct_gp = gpar(fontsize = 8))
#,  top_annotation=ha.top)

draw(p,heatmap_legend_side="bottom", annotation_legend_side="bottom")

p0=oncoPrint( onco_mat.sub3,col=col,
              alter_fun = alter_fun, get_type = function(x) strsplit(x, "\t") [[1]],column_order = sample_order3,
              remove_empty_columns = F, remove_empty_rows = TRUE,show_column_names = T, 
              #column_names_gp = gpar(cex=0.5, font=1, col= "blue"), 
              column_title = column_title,bottom_annotation = ha1,pct_gp = gpar(fontsize = 8))
draw(p0,heatmap_legend_side="bottom", annotation_legend_side="bottom")

dev.off()



#####germ##########
maf_2=readRDS("~/Desktop/P_Tran_ref/onco_508/germ/maf_all_germ.rds")



onco_mat=read.table("~/Desktop/onco_508/germ/onco_matrix.germ.txt",sep="\t",header=T,row.names=1,check.names=F)
anno_df=data.frame(getClinicalData(maf_2))
col = c("Missense_Mutation" = "blue", "In_Frame_Ins" = "red", "Splice_Site" = "green","In_Frame_Del"="purple","Frame_Shift_Ins"="yellow",
        "Frame_Shift_Del"="grey", "Nonsense_Mutation"="brown", "Multi_Hit"="orange","Translation_Start_Site"="pink","Nonstop_Mutation"="cyan",
        "3'Flank"="goldenrod4", "3'UTR"= "aquamarine4","5'Flank" ="lightslateblue",               "5'UTR"="lightcoral","IGR"="black" )

rate=data.frame(tmb(maf = maf_2))

foo1 <- data.frame(do.call('rbind', strsplit(as.character( rate$Tumor_Sample_Barcode),'.',fixed=TRUE)))
foo2=data.frame(do.call('rbind',strsplit(as.character(foo1$X1),'_',fix=TRUE)))

rate$index=foo2$X1

anno_df1=merge(anno_df,rate,by="Tumor_Sample_Barcode")


#prepare the topbar
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big blue
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["Missense_Mutation"], col = NA))
  },
  # big red
  In_Frame_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["In_Frame_Ins"], col = NA))
  },
  # small green
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),  
              gp = gpar(fill = col["Splice_Site"], col = NA))},
  
  # small green
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["In_Frame_Del"], col = NA))
  },
  
  Frame_Shift_Ins= function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["Frame_Shift_Ins"], col = NA))},
  
  Frame_Shift_Del= function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["Frame_Shift_Del"], col = NA))},
  Nonsense_Mutation= function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["Nonsense_Mutation"], col = NA))},
  Multi_Hit= function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["Multi_Hit"], col = NA))},
  Translation_Start_Site= function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["Translation_Start_Site"], col = NA))},
  Nonstop_Mutation= function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["Nonstop_Mutation"], col = NA))}
  
)


library(ComplexHeatmap)
library(ComplexHeatmap)
 column_title = "OncoPrint For Germline"
 
 anno_df1$index.y <- substr(anno_df1$index.y, 0, 2)
 
 #anno_df[c(5,13,17,22,24,26,28,30,32,34,36)] <- data.frame(lapply( anno_df[2:11], factor))
 
 i=c(8,16,20,25,27,29,31,33,35,37,39,41:44)
 
 #any(round(anno_df[,1:ncol(anno_df)]) != anno_df[,1:ncol(anno_df)])
 anno_df1[ , i] <- apply(anno_df1[ , i], 2,            # Specify own function within apply
                     function(x) as.numeric(as.character(x)))
 
 anno_df1=anno_df1[order(anno_df1$total_perMB,decreasing=T),]
 
 rownames(anno_df1)=anno_df1$Tumor_Sample_Barcode
 
 ha.top = HeatmapAnnotation(
   dist1 = anno_barplot(
     anno_df1$total_perMB, 
     bar_width = 1, 
     gp = gpar(col = "white", fill = "#FFE200"), 
     border = FALSE,
     
     height = unit(2, "cm")
   ), show_annotation_name = FALSE)
 
# ha2=HeatmapAnnotation(anno_df,OS.from.MDT = anno_df$OS.from.MDT,cbar = anno_oncoprint_barplot(), age = anno_df$Age,alive=anno_df$OS..1.dead..0.alive.)
 
 ha1 = HeatmapAnnotation(tmb=anno_df1$total_perMB,OS.from.MDT = anno_df1$OS.from.MDT, age = anno_df1$Age,alive=anno_df1$OS..1.dead..0.alive.,Distant.Failure.from.MDT=anno_df1$Distant.Failure.from.MDT,
                         Time.to.Distant.Failure.from.oligomet=anno_df1$Time.to.Distant.Failure.from.oligomet,Post.MDT.Local.Failure..time.from.MDT.to.local.failure=anno_df1$Post.MDT.Local.Failure..time.from.MDT.to.local.failure.,
                         PSA.Failure.time.from.MDT=anno_df1$PSA.Failure.time.from.MDT,PSA.failure.time.from.oligomet  =anno_df1$PSA.failure.time.from.oligomet,location=anno_df1$Location..1.node..2.bone..3.node.bone..visceral.4..visceral.other.5,
                         Imaging.Type=anno_df1$Imaging.Type..1.convention..2.enhanced.,Number.of.mets =anno_df1$Number.of.mets,PSA.at.oligomet =anno_df1$PSA.at.oligomet,
                         Risk.Group=anno_df1$Risk.Group,index=anno_df1$index,Primary.Gleason =anno_df1$Primary.Gleason,
                         tmb=anno_oncoprint_barplot(tmb,ylim = c(0, max(tmb, na.rm = TRUE))),
                         OS.from.MDT=anno_oncoprint_barplot(OS.from.MDT, ylim = c(0, max(OS.from.MDT, na.rm = TRUE)), axis = TRUE),
                         
                         Distant.Failure.from.MDT=anno_oncoprint_barplot(anno_df1$Distant.Failure.from.MDT, ylim = c(0, max(Distant.Failure.from.MDT, na.rm = TRUE)), axis = TRUE),
                         Time.to.Distant.Failure.from.oligomet=anno_oncoprint_barplot(anno_df1$Time.to.Distant.Failure.from.oligomet,ylim = c(0, max(anno_df1$Time.to.Distant.Failure.from.oligomet, na.rm = TRUE)), axis = F),
                         Post.MDT.Local.Failure..time.from.MDT.to.local.failure = anno_oncoprint_barplot(anno_df1$Post.MDT.Local.Failure..time.from.MDT.to.local.failure., ylim = c(0, max(anno_df1$Post.MDT.Local.Failure..time.from.MDT.to.local.failure., na.rm = TRUE)), axis = TRUE),
                         alive=anno_oncoprint_barplot(alive,ylim = c(0, max(alive, na.rm = TRUE)), axis = F),
                        
                         age = anno_points(age, ylim = c(0, max(age, na.rm = TRUE))),
                         col = list(alive = c("1" = "red", "0" = "blue","NA"="grey")),
                         Imaging.Type = c("1"="red","2"="blue"),Number.of.mets=c("1"="red","2"="blue","3"="green","4"="pink"),
                         Risk.Group=c("FI"="red","H"="blue","Metastatic"="green","Node+"="pink","UI"="yellow","VH"="cyan"),
                         Primary.Gleason=c("3"="red","4"="blue","5"="green"),
                         annotation_height = unit(c(1, 1, 1,1), "cm"))
 
 
 anno_df1=anno_df1[order(anno_df1$total_perMB,decreasing=T),]
 
 rownames(anno_df1)=anno_df1$Tumor_Sample_Barcode
 
 
 setwd("~/Desktop")
 
 amp = ifelse(apply(onco_mats, 1, function(x) sum(grepl("Missense_Mutation", x))/length(x) > 0.1), "high Missense_Mutation events", "low Missense_Mutation events")
 amp = factor(amp, levels = c("low Missense_Mutation events", "high Missense_Mutation events"))
 
 gene1=data.frame(maf_2@gene.summary)
 
 gene_plot1=gene1[which(gene1$MutatedSamples>110*.7),]
 
 mat_gs1=onco_mats[which(rownames(onco_mats)%in%gene_plot$Hugo_Symbol),]
 
 onco_mats= mat_gs1[order(match(colnames( mat_gs1),anno_df1$Tumor_Sample_Barcode))]
 #onco_mat.sub1=onco_mat.sub[order(match(colnames(onco_mat.sub),anno_df.sort$Tumor_Sample_Barcode))]
 

 

 
 gene_to_plot=read.table("P_Tran_ref/matation_to_plot.txt",sep="\t",header=T)
 
 mat_gs=onco_mats[which(rownames(onco_mats)%in%gene_to_plot$name),]
 


 ###tcga
 
 tcga=read.table("/Users/ysong/Downloads/prostate_dkfz_2018/onco_matrix.txt",header=T,sep="\t")
 
 
 mat_gs=tcga[which(rownames(tcga)%in%gene_to_plot$name),]
 
 pdf("somatic_tcga.v2.pdf",paper="special",  pointsize=10,  width=69.7/2.54,height=35.9/2.54 * 1,pagecentre=T)
 p0=oncoPrint(  mat_gs,col=col,alter_fun = alter_fun,pct_gp = gpar(fontsize = 8))
 draw(p0,heatmap_legend_side="bottom", annotation_legend_side="bottom")
 dev.off()
 
 
 
 #file:///Users/ysong/Downloads/ComplexHeatmap-supplementary1-4/supplS1_TCGA_OncoPrint/supplS1_OncoPrint.html
 pdf("germ4_metastatic_gene.pdf",paper="special",  pointsize=10,  width=39.7/2.54,height=40.9/2.54 * 1,pagecentre=T)
 #column_title = "Mutate in more than 70% of samples "
# p1=oncoPrint( mat_gs1,col=col,
 #            alter_fun = alter_fun, get_type = function(x) strsplit(x, "\t") [[1]],column_order = sample_order,
  #           remove_empty_columns = TRUE, remove_empty_rows = TRUE,show_column_names = TRUE,column_names_gp = gpar(cex=0.5, font=1, col= "black"),
   #          column_title = column_title,bottom_annotation = ha1,pct_gp = gpar(fontsize = 4))
 #,  top_annotation=ha.top)
 
# draw(p1,heatmap_legend_side="bottom", annotation_legend_side="bottom")
 column_title = "cancer mutant marker genes "
 p=oncoPrint(onco_mat.sub1,col=col,
             alter_fun = alter_fun, get_type = function(x) strsplit(x, "\t") [[1]],column_order = sample_order,
             remove_empty_columns = TRUE, remove_empty_rows = TRUE,show_column_names = TRUE,column_names_gp = gpar(cex=0.5, font=1, col= "black"),
             column_title = column_title,bottom_annotation = ha1,pct_gp = gpar(fontsize = 4))
 #,  top_annotation=ha.top)
 
 draw(p,heatmap_legend_side="bottom", annotation_legend_side="bottom")
 dev.off()
 
 
 
 
 
