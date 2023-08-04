

tempus.old=readRDS("/Users/ysong/Desktop/P_Tran_ref/analysis_result/maftools_post_analysis/maf_to_use/maf_previous_cancer_normal_pair_607.rds")


write.table(tempus.old@clinical.data,"tempus_pair_meta.txt",sep="\t",quote=F)

tempus=readRDS("/Users/ysong/Desktop/combine_tempus_foundation/tempus.filter.209.rds")
tempus_gene=read.table("/Users/ysong/Desktop/P_Tran_ref/analysis_result/freebay_plot/820_correct/tab_file/tempus_marker_clean.txt",sep="\t",header=T)
meta.t=tempus@clinical.data


setwd("~/Desktop")
write.table(meta.t,"meta_tempus.txt",sep="\t",quote=F)




tempus.sub=subsetMaf(tempus,genes=tempus_gene$value,dropLevels = F)
tempus.sub@clinical.data=tempus@clinical.data
var.classes=c("Frame_Shift_Del", 
              "Frame_Shift_Ins",
              "In_Frame_Del", 
              "In_Frame_Ins",
              "Missense_Mutation",
              "Nonsense_Mutation",
              "Nonstop_Mutation", 
              "Splice_Site" )
tempus.sub1=subsetMaf(maf = tempus.sub, includeSyn = F, query = "Variant_Classification%in%var.classes")
tempus.sub1@clinical.data=tempus@clinical.data

saveRDS(tempus.sub1,"tempus_panel_w_meta.212.rds")

write.table(tempus.sub1@data,"tempus_main.matrix.txt",sep="\t",row.names=F)


tempus_gene=read.table("/Users/ysong/Desktop/P_Tran_ref/analysis_result/freebay_plot/820_correct/tab_file/tempus_marker_clean.txt",sep="\t",header=T)

foundation_gene=read.table("/Users/ysong/Desktop/P_Tran_ref/foundation/foundation_gene.tab",sep="\t",header=F)

tempus.found=merge(tempus_gene,foundation_gene,by.y="V1",by.x="value")

pathogen.cat=c("Pathogenic","Pathogenic/Likely_pathogenic","Pathogenic/Likely_pathogenic,_other","Pathogenic/Likely_pathogenic,_risk_factor","Likely_pathogenic","Conflicting_interpretations_of_pathogenicity")
#pathogen.cat=c("Pathogenic","Pathogenic/Likely_pathogenic","Pathogenic/Likely_pathogenic,_other","Pathogenic/Likely_pathogenic,_risk_factor","Likely_pathogenic")



tempus.sub1_patho=subsetMaf(tempus.sub1, query=  "ClinVar_VCF_CLNSIG %in% pathogen.cat",dropLevels =F)
tempus.sub1_patho@clinical.data=tempus@clinical.data
tempus.sub1_patho.f=subsetMaf(tempus.sub1_patho,genes=tempus.found$value,dropLevels=F)

somatic_all=read.delim("/Users/ysong/Desktop/P_Tran_ref/foundation/combine_somatic_germline/somatic_nonfilter.mod.txt",sep="\t",header=T)
somatic_all.maf=read.maf(somatic_all)


meta.f=somatic_all.maf@clinical.data

somatic_f=subsetMaf(maf = somatic_all.maf, includeSyn = F, query = "Variant_Classification%in%var.classes")

pathogen.cat=c("Pathogenic","Pathogenic/Likely_pathogenic","Pathogenic/Likely_pathogenic,_other","Pathogenic/Likely_pathogenic,_risk_factor","Likely_pathogenic","Conflicting_interpretations_of_pathogenicity")

#pathogen.cat=c("Pathogenic","Pathogenic/Likely_pathogenic","Pathogenic/Likely_pathogenic,_other","Pathogenic/Likely_pathogenic,_risk_factor","Likely_pathogenic")


foundation_patho=subsetMaf(somatic_f, query=  "ClinVar_VCF_CLNSIG %in% pathogen.cat",dropLevels =F)


foundation_patho.f=subsetMaf(foundation_patho,genes=tempus.found$value,dropLevels=F)

setwd("combine_tempus_foundation/")


write.table(foundation_patho.f@data,"foundation_filter.matrix.212.txt",sep='\t',quote=F)
write.table(tempus.sub1_patho.f@data,"tempus_filter.matrix.212.txt",sep='\t',quote=F)


#tempus.found1=merge_mafs(c(tempus.sub1_patho.f,foundation_patho.f),removeDuplicatedVariants = T, verbose = TRUE)

f.matrix=foundation_patho.f@data
f.matrix.f=f.matrix[which(f.matrix$Hugo_Symbol%in%tempus.found$value),]

t.matrix=tempus.sub1_patho.f@data
t.matrix.f=t.matrix[which(t.matrix$Hugo_Symbol%in%tempus.found$value),]

f.matrix.f$group="foundation"
t.matrix.f$group="tempus"




p1=oncoplot(foundation_patho.f,removeNonMutated = F,writeMatrix = T,genes=foundation_patho.f@gene.summary$Hugo_Symbol)


p2=oncoplot(tempus.sub1_patho.f,removeNonMutated = F,writeMatrix = T,genes=tempus.sub1_patho.f@gene.summary$Hugo_Symbol)

#f_matrix= maftools:::createOncoMatrix(foundation_patho,g=foundation_patho@gene.summary$Hugo_Symbol,droplevels=F)




clinvar_db=read.table("/Users/ysong/Desktop/P_Tran_ref/clinvar_db/somatic/conflict_pathogent1.txt",sep="\t",header=F)
clinvar_db$V6=gsub("CLNSIG=","",clinvar_db$V6)
clinvar_db$id=paste(clinvar_db$V1,clinvar_db$V2,sep="_")

f.data=foundation_patho.f@data
t.data=tempus.sub1_patho.f@data

saveRDS(foundation_patho.f,"foundation_pathogene_w_conflic.214.rds")
saveRDS(tempus.sub1_patho.f,"tempus_pathogene_w_conflic.214.rds")


tempus.sub1_patho.f@clinical.data=tempus@variants.per.sample








##############################pathogen_only
 `%notin%` <- Negate(`%in%`)


#query="ClinVar_VCF_CLNSIG %in% pathogen.cat"

term=c("Conflicting_interpretations_of_pathogenicity")

foundation_patho.f.sub=subsetMaf(foundation_patho.f,query="ClinVar_VCF_CLNSIG %notin%term",dropLevels =F)

tempus_patho.f.sub=subsetMaf(tempus.sub1_patho.f,query="ClinVar_VCF_CLNSIG %notin%term",dropLevels =F)

setwd("/Users/ysong/Desktop/combine_tempus_foundation/212/pathogen_only/")
pdf("foundation.pathogen_only.pdf")
oncoplot(foundation_patho.f.sub,removeNonMutated = F,writeMatrix = T,genes=foundation_patho.f.sub@gene.summary$Hugo_Symbol)
dev.off()


pdf("tempus.pathogen_only.pdf")
oncoplot(tempus_patho.f.sub,removeNonMutated = F,writeMatrix = T,genes=tempus_patho.f.sub@gene.summary$Hugo_Symbol)
dev.off()


f1=read.table("/Users/ysong/Desktop/combine_tempus_foundation/foundation.onco_matrix.txt",sep="\t",header=T,row.names=1,check.names = F)

t1=read.table("/Users/ysong/Desktop/combine_tempus_foundation/tempus.onco_matrix.txt",sep="\t",header=T,row.names=1,check.names=F)

ft1=merge(f1,t1,all=T,by=0)

ft1[is.na(ft1)]<-""

write.table(ft1,"found_tempus_comb_all.txt",sep="\t",row.names=F,quote=F)
rownames(ft1)=ft1$Row.name

ft1=read.table("found_tempus_comb_all.txt",sep="\t",header=T,row.names=1,check.names=F)
setwd("~/Desktop/combine_tempus_foundation/212/")
input_matrix=ft1
meta=read.table("meta_213.txt",sep="\t",header=T)

meta.sub=meta


#meta.sub=meta[which(meta$Time..2.synchronous..1.metachronous.!="NA"),]

imput.sub=input_matrix[,which(colnames(input_matrix)%in%meta.sub$id)]

library(reshape2)
imput.sub.m=melt(as.matrix(imput.sub))

table(imput.sub.m$value)



input.sub=imput.sub


meta.met=meta.sub[which(meta.sub$Time..2.synchronous..1.metachronous.==2),]

meta.syn=meta.sub[which(meta.sub$Time..2.synchronous..1.metachronous.==1),]

matrix.met=input.sub[,which(colnames(input.sub)%in%meta.met$id)]

matrix.syn=input.sub[,which(colnames(input.sub)%in%meta.syn$id)]





col = c("Missense_Mutation" = "blue", "In_Frame_Ins" = "red", "Splice_Site" = "#008000","In_Frame_Del"="purple","Frame_Shift_Ins"="yellow",
        "Frame_Shift_Del"="grey", "Nonsense_Mutation"="brown", "Multi_Hit"="orange")




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
              gp = gpar(fill = col["Multi_Hit"], col = NA))}
)

matrix.met[is.na(matrix.met)]<-"" 

matrix.syn[is.na(matrix.syn)]<-""

matrix.met[matrix.met==0]<-""

matrix.syn[matrix.syn==0]<-""

pdf("~/Desktop/somatic_tempus.meta_vs_syn_pathogenic_only.pdf",20,20)
t01=oncoPrint(matrix.met, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="Metachronous(82)")
t02=oncoPrint(matrix.syn, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="Synchronous(201)")

t01+t02
dev.off()



#########################################################

f.data.sub1=f.data[which(f.data$ClinVar_VCF_CLNSIG=="Conflicting_interpretations_of_pathogenicity"),]
t.data.sub1=t.data[which(t.data$ClinVar_VCF_CLNSIG=="Conflicting_interpretations_of_pathogenicity"),]

f.data.sub1$new_id=paste(f.data.sub1$Chromosome,f.data.sub1$Start_Position,sep="_")
t.data.sub1$new_id=paste(t.data.sub1$Chromosome,t.data.sub1$Start_Position,sep="_")


f.data.sub1.ann=merge(data.frame(f.data.sub1),data.frame(clinvar_db),by.x="new_id",by.y="id")
f.data.sub1.ann.f=f.data.sub1.ann[grep(paste("pathogenic","Pathogenic",sep="|"),f.data.sub1.ann$V7),]
f.data.sub1.ann.f1=f.data.sub1.ann.f[,c(1:23,34)]

t.data.sub1.ann=merge(data.frame(t.data.sub1),data.frame(clinvar_db),by.x="new_id",by.y="id")

t.data.sub1.ann.f=t.data.sub1.ann[grep(paste("pathogenic","Pathogenic",sep="|"),t.data.sub1.ann$V7),]
t.data.sub1.ann.f1=t.data.sub1.ann.f[,which(colnames(t.data.sub1.ann.f)%in%colnames(f.data.sub1.ann.f1))]
t.data.sub1.ann.f2=t.data.sub1.ann.f1[,order(match(colnames(t.data.sub1.ann.f1),colnames(f.data.sub1.ann.f1)))]




f.data.sub2=f.data[which(f.data$ClinVar_VCF_CLNSIG!="Conflicting_interpretations_of_pathogenicity"),]
t.data.sub2=t.data[which(t.data$ClinVar_VCF_CLNSIG!="Conflicting_interpretations_of_pathogenicity"),]

t.data.sub2$V7=t.data.sub2$ClinVar_VCF_CLNSIG
f.data.sub2$V7=f.data.sub2$ClinVar_VCF_CLNSIG



f.data.sub2.f2=f.data.sub2[,c(1:22,20)]
colnames(f.data.sub2.f2)[23]="V7"



t.data.sub2.f2=t.data.sub2[,c(1:2,5:7,16,9:13,80:84,14,40,42,107,119,57)]

t.data.sub2.f2$V7=t.data.sub2.f2$ClinVar_VCF_CLNSIG






tempus.new.all=rbind(t.data.sub1.ann.f2[-1],t.data.sub2.f2)


foundation.new.all=rbind(f.data.sub1.ann.f1[-1],f.data.sub2.f2)

foundation.new.all$center="Foundation"
tempus.new.all$center="Tempus"

all=rbind(foundation.new.all,tempus.new.all)

allmaf=read.maf(all)

p0=oncoplot(allmaf,genes=allmaf@gene.summary$Hugo_Symbol,writeMatrix=T,removeNonMutated = F)


sample=unique(somatic_f@data$Tumor_Sample_Barcode)
sample1=(tempus@variants.per.sample$Tumor_Sample_Barcode)
sample2=c(sample,sample1)

sample2t=data.frame(sample2)


matrix=read.table("onco_matrix.txt",sep="\t",header=T,row.names=1,check.names=F)
matrix.t=t(matrix)

matrix.t1=merge(sample2t,matrix.t,by.x="sample2",by.y=0,all.x=T)

matrix.t1[is.na(matrix.t1)]<-""
rownames(matrix.t1)=matrix.t1$sample2
colnames(matrix.t1)=gsub("-T","",colnames(matrix.t1))

write.table(t(matrix.t1[-1]),"foundation_tempus_pathogen_comb.212.txt",sep="\t",quote=F)

meta=read.delim("/Users/ysong/Downloads/01_Sample_Metadata/Clinical_overall.txt",sep="\t",header=T)

id=data.frame(matrix.t1$sample2)
names(id)="id"

matrix2=merge(id,meta,by.x="id",by.y="matrix_id",all.x=T)
write.table(matrix2,"~/Desktop/meta_213.txt",sep="\t",quote=F)

#########pathogene+conflict

setwd("~/Desktop/combine_tempus_foundation/212/")
input_matrix=read.table("foundation_tempus_pathogen_comb.212.txt",sep="\t",header=T,row.names=1,check.names=F)
meta=read.table("meta_213.txt",sep="\t",header=T)

meta1=data.frame(meta[,1])
names(meta1)="id"


matrix1=tempus.sub1_patho.f@data

matrix2=foundation_patho.f@data

matrix1f=matrix1[,c(1,   2  , 5:7,  16,   9:13,  80:84,  14,  40,  42, 107, 119,  57, 233:235)]
matrix2f=matrix2[,-24]




matrix.all=rbind(matrix1f,matrix2f)

matrix.all.m=merge(matrix.all,meta1,by.y="id",by.x="Tumor_Sample_Barcode",all.y=T)


all.maf=read.maf(matrix.all.m)
sample=all.maf@clinical.data

sample1=merge(sample,unique(meta[-c(203,206,83),-10]),by.x="Tumor_Sample_Barcode",by.y="id")


sample1f=setDT(sample1)
all.maf@clinical.data=sample1f

saveRDS(all.maf,"tempus_foundation_w_conflicP.rds")


setwd("/Users/ysong/Desktop/combine_tempus_foundation/212")
all.maf=read_rds("/Users/ysong/Desktop/combine_tempus_foundation/212/tempus_foundation_w_conflicP.rds")


OncogenicPathways(maf = all.maf)



pdf("summary.pdf")
plotmafSummary(
  all.maf,
  rmOutlier = TRUE,
  dashboard = TRUE,
  titvRaw = TRUE,
  log_scale = FALSE,
  addStat = NULL,
  showBarcodes = FALSE,
  fs = 1,
  textSize = 0.8,
  color = NULL,
  titleSize = c(1, 0.8),
  titvColor = NULL,
  top = 10
)
dev.off()





library("MesKit")




 all.maf@clinical.data$crpc=as.factor(as.character(all.maf@clinical.data$crpc))
 all.maf@clinical.data$Time..2.synchronous..1.metachronous.=as.factor(as.character(all.maf@clinical.data$Time..2.synchronous..1.metachronous.))


 vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
pdf("spop.mutation.pdf")
oncoplot(all.maf,draw_titv = T,removeNonMutated = F,genes=c("SPOP","BRCA1"),clinicalFeatures=c("crpc","Time..2.synchronous..1.metachronous."),bgCol = "white")
dev.off()



laml.mut.excl = somaticInteractions(maf = all.maf, top = 10)

co=data.frame(laml.mut.excl)


co.f=co[which(co$pValue<0.01&co$oddsRatio>5),]

co.f1=co[which(co$pValue<0.1&co$Event=="Mutually_Exclusive"),]

 pdf("onco_strip.coocurrenc.pdf")
 oncostrip(all.maf,genes =unique(c(co.f$gene1,co.f$gene2)),removeNonMutated = F, showTumorSampleBarcodes = FALSE)
 oncostrip(all.maf,genes =unique(c(co.f1$gene1,co.f1$gene2)),removeNonMutated = F, showTumorSampleBarcodes = FALSE)
 
 
 dev.off()




pdf("co-occurance1.pdf",20,20)
plot1=somaticInteractions(maf = all.maf, top = 50, pvalue = c(0.05,0.01))
dev.off()


#laml.sig = oncodrive(maf = all.maf, AACol = 'Protein_Change', minMut = 3, pvalMethod = 'zscore')

all.maf.sub=subsetMaf(all.maf,query="Hugo_Symbol%in%all.maf@gene.summary$Hugo_Symbol")
laml.sig = oncodrive(maf = all.maf.sub, AACol = 'Protein_Change', minMut = 3, pvalMethod = 'zscore')
pdf("drive.pdf")
plotOncodrive(res = laml.sig, fdrCutOff = 0.05, useFraction = TRUE, labelSize = 0.3)
ggplot(data = laml.sig, aes(x =fract_muts_in_clusters, y = -log10(fdr), label = Hugo_Symbol))+geom_point(color = ifelse(-log10(laml.sig$pval) >2,  "red","grey50"))+ggrepel::geom_text_repel()
dev.off()


p=oncoplot(maf = all.maf, pathways = "auto", gene_mar = 8, fontSize = 0.6)


##laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = laml.titv)



Pathway  N n_affected_genes fraction_affected Mutated_samples Fraction_mutated_samples
1:                   Genome_integrity 14                7         0.5000000             105              0.362068966
2:            Wnt/B-catenin_signaling  8               NA                NA              53              0.182758621
3:          Chromatin_SWI/SNF_complex  8               NA                NA              43              0.148275862
4: Protein_homeostasis/ubiquitination 15               NA                NA              31              0.106896552
5:                    Other_signaling 28                8         0.2857143              29              0.100000000
6:                     PI3K_signaling  9                6         0.6666667              27              0.093103448
7:                     MAPK_signaling  9                5         0.5555556              25              0.086206897
8:                      TOR_signaling  3                3         1.0000000              20              0.068965517
9:                     TGFB_signaling  7                3         0.4285714              19              0.065517241
10:                      RTK_signaling 16                7         0.4375000              15              0.051724138
11:                         Cell_cycle  8                3         0.3750000              15              0.051724138
12:               Transcription_factor 39                5         0.1282051               8              0.027586207
13:                    Chromatin_other 14                2         0.1428571               6              0.020689655
14:                    NOTCH_signaling  1                1         1.0000000               5              0.017241379
15:                           Splicing  6                1         0.1666667               5              0.017241379
16:        Chromatin_histone_modifiers 15                5         0.3333333               3              0.010344828
17:                         Metabolism  2                2         1.0000000               2              0.006896552
18:                     NFKB_signaling  2                1         0.5000000               1              0.003448276



s1=OncogenicPathways(maf = all.maf)


pathway=c("Cell_Cycle", "Hippo"  ,    "MYC" ,       "NOTCH" ,     "NRF2" ,      "PI3K"  ,     "RTK-RAS" ,   "TGF-Beta"  , "TP53"  ,     "WNT"  )

length(pathway)



pdf("pathway.pdf")
i=1
PlotOncogenicPathways(maf = all.maf, pathways = pathway[i], fullPathway = TRUE)
i=2
PlotOncogenicPathways(maf = all.maf, pathways = pathway[i], fullPathway = TRUE)

i=3
PlotOncogenicPathways(maf = all.maf, pathways = pathway[i], fullPathway = TRUE)


i=4
PlotOncogenicPathways(maf = all.maf, pathways = pathway[i], fullPathway = TRUE)


#i=5
#PlotOncogenicPathways(maf = all.maf, pathways = pathway[i], fullPathway = TRUE)


i=6
PlotOncogenicPathways(maf = all.maf, pathways = pathway[i], fullPathway = TRUE)


i=7
PlotOncogenicPathways(maf = all.maf, pathways = pathway[i], fullPathway = TRUE)


i=8
PlotOncogenicPathways(maf = all.maf, pathways = pathway[i], fullPathway = TRUE)

i=9
PlotOncogenicPathways(maf = all.maf, pathways = pathway[i], fullPathway = TRUE)
dev.off()






#######################################################################
setwd("~/Desktop/combine_tempus_foundation/212/")
input_matrix=read.table("foundation_tempus_pathogen_comb.212.txt",sep="\t",header=T,row.names=1,check.names=F)
imput.sub=input_matrix[,which(colnames(input_matrix)%in%meta.sub$id)]

library(reshape2)
imput.sub.m=melt(as.matrix(imput.sub))

table(imput.sub.m$value)



input.sub=imput.sub

meta.met=meta.sub[which(meta.sub$Time..2.synchronous..1.metachronous.==1),]

meta.syn=meta.sub[which(meta.sub$Time..2.synchronous..1.metachronous.==2),]

matrix.met=input.sub[,which(colnames(input.sub)%in%meta.met$id)]

matrix.syn=input.sub[,which(colnames(input.sub)%in%meta.syn$id)]





col = c("Missense_Mutation" = "blue", "In_Frame_Ins" = "red", "Splice_Site" = "#008000","In_Frame_Del"="purple","Frame_Shift_Ins"="yellow",
        "Frame_Shift_Del"="grey", "Nonsense_Mutation"="brown", "Multi_Hit"="orange")




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
              gp = gpar(fill = col["Multi_Hit"], col = NA))}
)

matrix.met[is.na(matrix.met)]<-""

matrix.syn[is.na(matrix.syn)]<-""

pdf("~/Desktop/somatic_tempus.meta_vs_syn_pathogenic_related.pdf",20,20)
t01=oncoPrint(matrix.met, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="Metachronous(201(109F+92T))")
t02=oncoPrint(matrix.syn, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="Synchronous(82(71F+11T))")

t01+t02
dev.off()


meta_use=all.maf@clinical.data
meta.met=meta_use[which(meta_use$Time..2.synchronous..1.metachronous.==1),]

meta.met.s=meta.met[order(match(meta.met$id,colnames(matrix.met))),]

anno_oncoprint_barplot(type = NULL, which = c("column", "row"),
                       bar_width = 0.6, ylim = NULL, show_fraction = FALSE, axis = TRUE,
                       axis_param = if(which == "column") default_axis_param("column") else list(side = "top", labels_rot = 0),
                       width = NULL, height = NULL, border = FALSE)





input_matrix=read.table("foundation_tempus_pathogen_comb.212.txt",sep="\t",header=T,row.names=1,check.names=F)



meta.ncrpc=meta.sub[which(meta.sub$crpc==0),]

meta.crpc=meta.sub[which(meta.sub$crpc==1),]

matrix.ncrpc=input_matrix[,which(colnames(input_matrix)%in%meta.ncrpc$id)]

matrix.crpc=input_matrix[,which(colnames(input_matrix)%in%meta.crpc$id)]





col = c("Missense_Mutation" = "blue", "In_Frame_Ins" = "red", "Splice_Site" = "#008000","In_Frame_Del"="purple","Frame_Shift_Ins"="yellow",
        "Frame_Shift_Del"="grey", "Nonsense_Mutation"="brown", "Multi_Hit"="orange")



matrix.ncrpc[is.na(matrix.ncrpc)]<-""

matrix.crpc[is.na(matrix.crpc)]<-""

pdf("~/Desktop/somatic_tempus.ncrpc_vs_crpc_pathogenic_related.pdf",20,20)
t01=oncoPrint(matrix.ncrpc, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="non-crpc(184(95F+86T))")
t02=oncoPrint(matrix.crpc, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="crpc(94(77F+17F))")

t01+t02
dev.off()








input.sub[is.na(input.sub)]<-""


library(survival)
library(lubridate)
library(ggsurvfit)
library(gtsummary)
#library(tidycmprsk)
#library(condsurv)

##https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html
write.table(input.sub,"input_numeric_matrix.txt",sep="\t",quote=F)


setwd("/Users/ysong/Desktop/combine_tempus_foundation/212/")

setwd("~/Desktop/combine_tempus_foundation/212/")
input_matrix=read.table("foundation_tempus_pathogen_comb.212.txt",sep="\t",header=T,row.names=1,check.names=F)



meta=read.table("meta_213.txt",sep="\t",header=T)
meta.sub=meta[which(meta$Time..2.synchronous..1.metachronous.!="NA"),]

imput.sub=input_matrix[,which(colnames(input_matrix)%in%meta.sub$id)]


pw_rt = data.frame(
  gene = c('ABL1', 'EGFR', 'ERBB2', 'EPHA3',"EPHA2",'MAP3K1','MAP2K7','LIN7A','RASA1','ERBB3', 'ERBB4', 'PDGFRA', 'PDGFRB', 'MET', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FLT3', 'ALK', 'RET', 'ROS1', 'KIT', 'IGF1R', 'NTRK1', 'NTRK2', 'NTRK3', 'SOS1', 'GRB2', 'PTPN11', 'KRAS', 'HRAS', 'NRAS', 'RIT1', 'ARAF', 'BRAF', 'RAF1', 'RAC1', 'MAP2K1', 'MAP2K2', 'MAPK1', 'NF1', 'RASA1', 'CBL', 'ERRFI1', 'CBLB', 'CBLC', 'INSR', 'INSRR', 'IRS1', 'SOS2', 'SHC1', 'SHC2', 'SHC3', 'SHC4', 'RASGRP1', 'RASGRP2', 'RASGRP3', 'RASGRP4', 'RAPGEF1', 'RAPGEF2', 'RASGRF1', 'RASGRF2', 'FNTA', 'FNTB', 'RCE1', 'ICMT', 'MRAS', 'PLXNB1', 'MAPK3', 'ARHGAP35', 'RASA2', 'RASA3', 'RASAL1', 'RASAL2', 'RASAL3', 'SPRED1', 'SPRED2', 'SPRED3', 'DAB2IP', 'SHOC2', 'PPP1CA', 'SCRIB', 'PIN1', 'KSR1', 'KSR2', 'PEBP1', 'ERF', 'PEA15', 'JAK2', 'IRS2','EIF4EBP1', 'AKT1', 'AKT2', 'AKT3', 'AKT1S1', 'DEPDC5', 'DEPTOR', 'INPP4B', 'MAPKAP1', 'MLST8', 'MTOR', 'NPRL2', 'NPRL3', 'PDK1', 'PIK3CA', 'PIK3CB', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PPP2R1A', 'PTEN', 'RHEB', 'RICTOR', 'RPTOR', 'RPS6', 'RPS6KB1', 'STK11', 'TSC1', 'TSC2'),
  pw = c('RTK/PIK3'),
  stringsAsFactors = F
)

pw_cc = data.frame(
  gene = c('CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CDK2', 'CDK4', 'CDK6', 'RB1', 'E2F1', 'E2F3'),
  pw = c('CellCycle'),
  stringsAsFactors = F
)


pw_wnt = data.frame(
  gene = c('CHD8', 'LEF1', 'LGR4', 'LGR5', 'LRP5', 'LRP6', 'LZTR1', 'NDP', 'PORCN', 'RSPO1', 'SFRP1', 'SFRP2', 'SFRP4', 'SFRP5', 'SOST', 'TCF7L1', 'TLE1', 'TLE2', 'TLE3', 'TLE4', 'WIF1', 'ZNRF3', 'CTNNB1', 'DVL1', 'DVL2', 'DVL3', 'FRAT1', 'FRAT2', 'FZD1', 'FZD10', 'FZD2', 'FZD3', 'FZD4', 'FZD5', 'FZD6', 'FZD7', 'FZD8', 'FZD9', 'WNT1', 'WNT10A', 'WNT10B', 'WNT11', 'WNT16', 'WNT2', 'WNT3A', 'WNT4', 'WNT5A', 'WNT5B', 'WNT6', 'WNT7A', 'WNT7B', 'WNT8A', 'WNT8B', 'WNT9A', 'WNT9B', 'AMER1', 'APC', 'AXIN1', 'AXIN2', 'DKK1', 'DKK2', 'DKK3', 'DKK4', 'GSK3B', 'RNF43', 'TCF7', 'TCF7L2', 'CHD4'),
  pw = c('WNT'),
  stringsAsFactors = F
)

pw_tp = data.frame(
  gene = c('TP53', 'MDM2', 'MDM4', 'ATM', 'CHEK2', 'RPS6KA3'),
  pw = c('TP53'),
  stringsAsFactors = F
)

pw_hippo = data.frame(
  gene = c('STK4', 'STK3', 'SAV1', 'LATS1', 'LATS2', 'MOB1A', 'MOB1B', 'YAP1', 'WWTR1', 'TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'PTPN14', 'NF2', 'WWC1', 'TAOK1', 'TAOK2', 'TAOK3', 'CRB1', 'CRB2', 'CRB3', 'LLGL1', 'LLGL2', 'HMCN1', 'SCRIB', 'HIPK2', 'FAT1', 'FAT2', 'FAT3', 'FAT4', 'DCHS1', 'DCHS2', 'CSNK1E', 'CSNK1D', 'AJUBA', 'LIMD1', 'WTIP'),
  pw = c('HIPPO'),
  stringsAsFactors = F
)

pw_tgf = data.frame(
  gene = c('TGFBR1', 'TGFBR2', 'ACVR2A', 'ACVR1B', 'SMAD2', 'SMAD3', 'SMAD4'),
  pw = c('TGF-Beta'),
  stringsAsFactors = F
)

pw_chre = data.frame(
  gene = c('ARID1A', 'ARID2','ARID1B','SMARCA4','PBRM1','KDM6A','MBD6','HIST1H3B'),
  pw = c('ChrRemod'),
  stringsAsFactors = F
)


pw_myc = data.frame(
  gene = c('MAX', 'MGA', 'MLX', 'MLXIP', 'MLXIPL', 'MNT', 'MXD1', 'MXD3', 'MXD4', 'MXI1', 'MYC', 'MYCL', 'MYCN'),
  pw = c('MYC'),
  stringsAsFactors = F
)
pw_notch = data.frame(
  gene = c('ARRDC1', 'CNTN6', 'CREBBP', 'EP300', 'HES1', 'HES2', 'HES3', 'HES4', 'HES5', 'HEY1', 'HEY2', 'HEYL', 'KAT2B', 'KDM5A', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NOV', 'NRARP', 'PSEN2', 'LFNG', 'ITCH', 'NCSTN', 'SPEN', 'JAG1', 'APH1A', 'FBXW7', 'FHL1', 'THBS2', 'HDAC2', 'MFAP2', 'CUL1', 'RFNG', 'NCOR1', 'NCOR2', 'MFAP5', 'HDAC1', 'NUMB', 'JAG2', 'MAML3', 'MFNG', 'CIR1', 'CNTN1', 'MAML1', 'MAML2', 'NUMBL', 'PSEN1', 'PSENEN', 'RBPJ', 'RBPJL', 'RBX1', 'SAP30', 'SKP1', 'SNW1', 'CTBP1', 'CTBP2', 'ADAM10', 'APH1B', 'ADAM17', 'DLK1', 'DLL1', 'DLL3', 'DLL4', 'DNER', 'DTX1', 'DTX2', 'DTX3', 'DTX3L', 'DTX4', 'EGFL7'),
  pw = c('NOTCH'),
  stringsAsFactors = F
)

pw_gi = data.frame(
  gene = c('GATA6', 'GATA4','MUC6'),
  pw = c('GI'),
  stringsAsFactors = F
)

pw = rbind(pw_rt, pw_cc,pw_wnt,pw_tp,pw_hippo,pw_tgf,pw_chre,pw_myc,pw_notch,pw_gi)

pw_col = RColorBrewer::brewer.pal(n = 11, name = 'Set3')[1:11]
names(pw_col) = unique(pw$pw)




ranks <-merge(pw,matrix.met,by.x="gene",by.y=0)

ranks1 <- merge(pw,matrix.syn,by.x="gene",by.y=0)




library("dplyr")
library("tidyr")
rank01=ranks 
rownames(rank01)=ranks$gene

rank01[,-c(1:2)] <- as.integer(rank01[,-c(1:2)] !="")

rank01[rank01 == 0]<-""

write.table(rank01,"path_met.txt",sep="\t",quote=F)






rank01m%>%count(pw,variable)


rank02=ranks1
rownames(rank02)=ranks1$gene

rank02[,-c(1:2)] <- as.integer(rank02[,-c(1:2)] !="")

rank02[rank02 == 0]<-""

write.table(rank02,"path_syn.txt",sep="\t",quote=F)


pathway010=read.table("path_met.path.txt",sep="\t",header=T,row.names=1,check.names=F)
pathway020=read.table("path_syn.pathway2.txt",sep="\t",header=T,row.names=1,check.names=F)
pathway010[is.na(pathway010)]<-""
pathway020[is.na(pathway020)]<-""

pathway010[pathway010==0]<-""
pathway020[pathway020==0]<-""



col0 = c("1" = "red","2" = "red","3" = "red","4"="red","5"="red","6"="red")

alter_fun0 = list(
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
  "2" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col0["2"], col = NA))
    
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

path=data.frame(OncogenicPathways(maf = all.maf))
pathway010.ann=merge(path[,1:3],pathway010,by.x="Pathway",by.y=0)

rownames(pathway010.ann)=paste(pathway010.ann$Pathway,"(",pathway010.ann$n_affected_genes,"/",pathway010.ann$N,")",sep="")

pathway020.ann=merge(path[,1:3],pathway020,by.x="Pathway",by.y=0)
rownames(pathway020.ann)=paste(pathway020.ann$Pathway,"(",pathway020.ann$n_affected_genes,"/",pathway020.ann$N,")",sep="")


pdf("~/Desktop/pathway.215.summary.pdf",20,6)
t01=oncoPrint(pathway010.ann[-c(1:3)], remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun0, col = col0,show_heatmap_legend =F,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="Synchronous(82(71F+11T))")
t02=oncoPrint(pathway020.ann[-c(1:3)], remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun0, col = col0,show_heatmap_legend =F,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="Metachronous(201(109F+92T))")

t01+t02
dev.off()





de.split=split(pw, f=pw$pw)

ranks=""
ranks1=""


for(i in 1:length(de.split)){

  
  
  
  #fgseaRes <- fgseaMultilevel(pathways=c8, stats=ranks,scoreType = "pos")
  filename <- paste(names(de.split)[i], ".pathway.pdf", sep="")




t01=oncoPrint(ranks, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="Synchronous(82(71F+11T))")
t02=oncoPrint(ranks1, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="Metachronous(201(109F+92T))")

p=t01+t02

ggsave(p,filename=paste("onco",names(de.split)[i],".pathway.pdf",sep=""))}








####survival_analysis
numeric_matrix=read.table("input_numeric_matrix.txt",sep="\t",header=T,row.names=1,check.names=F)

numeric_matrix[is.na(numeric_matrix)]<-""


library(dplyr)
#Code
#First make variables in standard format







library(survival);library(survminer)




setwd("/Users/ysong/Desktop/combine_tempus_foundation/212/survival_plot_216")

library(survminer)
library(survival)
# vector with the variables to run through
genes <- rownames(numeric_matrix)



for(i in 1:108){

  tp53=data.frame(t(numeric_matrix[grep(genes[i],rownames(numeric_matrix)),]))
  
  tp53.ann=merge(meta.sub[-c(11:37)],tp53,by.x="id",by.y=0)
  
  
  

  
  fit <- survfit(as.formula(paste0("Surv(met_survival, status) ~", genes[i])),
                 data = tp53.ann)
  
  pdf(paste0(genes[i], "_Survival_w_vs_wo_mutation.pdf"),20,10)
  p1=   ggsurvplot(fit, conf.int=F, pval=TRUE, risk.table=TRUE, 
               legend.labs=c(paste(genes[i],"_WO_mutation",sep="_"), paste(genes[i],"_W_mutation",sep="_")), legend.title="Mutation",  
               palette=c("dodgerblue2", "orchid2","grey"), 
               title = "synchronous/metachronous", subtitle = "Based on Kaplan-Meier estimates",
               
               font.title = c(10, "bold", "darkblue"),
               font.subtitle = c(10, "bold.italic", "purple"),
               font.caption = c(10, "plain", "orange"), 
               risk.table.height=0.2,ggtheme = theme_minimal())
  
  
  
  tp53.ann.m=tp53.ann[which(tp53.ann$Time..2.synchronous..1.metachronous.==1),]
  tp53.ann.s=tp53.ann[which(tp53.ann$Time..2.synchronous..1.metachronous.==2),]
  
  fit1m <- survfit(as.formula(paste0("Surv(met_survival, status) ~", genes[i])),
                 data =  tp53.ann.m)
  fit1s <- survfit(as.formula(paste0("Surv(met_survival, status) ~", genes[i])),
                   data =  tp53.ann.s)
  

  p2=   ggsurvplot(fit1m, conf.int=F, pval=TRUE, risk.table=TRUE, 
                   legend.labs=c(paste(genes[i],"metachronous_WO_mutation",sep="_"), paste(genes[i],"_metachronous_W_mutation",sep="_")), 
                   legend.title="Mutation",  
  palette=c("dodgerblue2", "orchid2","grey"), 
  title = "split by synchronous/metachronous", subtitle = "Based on Kaplan-Meier estimates",
  
  font.title = c(10, "bold", "darkblue"),
  font.subtitle = c(10, "bold.italic", "purple"),
  font.caption = c(10, "plain", "orange"), 
  risk.table.height=0.2,ggtheme = theme_minimal())

  #source("https://raw.githubusercontent.com/BingxinS/survminer-fix/master/ggsurvplot_facet_risktable.R")
  p3=   ggsurvplot(fit1s, conf.int=F, pval=TRUE, risk.table=TRUE, 
                   legend.labs=c(paste(genes[i],"synchronous_WO_mutation",sep="_"), paste(genes[i],"_synchronous_W_mutation",sep="_")), 
                   legend.title="Mutation",  
                   palette=c("dodgerblue2", "orchid2","grey"), 
                   title = "split by synchronous/metachronous", subtitle = "Based on Kaplan-Meier estimates",
                   
                   font.title = c(10, "bold", "darkblue"),
                   font.subtitle = c(10, "bold.italic", "purple"),
                   font.caption = c(10, "plain", "orange"), 
                   risk.table.height=0.2,ggtheme = theme_minimal())
  
  print(p1)
 
    
  require("survminer")
  splots <- list()
  splots[[1]] <- p2
  splots[[2]] <- p3
  
  # Arrange multiple ggsurvplots and print the output
  arrange_ggsurvplots(splots, print = TRUE,
                      ncol = 2, nrow = 1, risk.table.height = 0.4)
  
  dev.off()
}

dev.off()

