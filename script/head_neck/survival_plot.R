
library(maftools)


soma=readRDS("soma.filter.rds")

pdf("tmb.pdf")
tmb=tmb(soma)
dev.off()

meta=read.table("P_Tran_ref/head_neck/20230511_GENOMIC DATA.txt",sep="\t",header=T)

meta0=soma@clinical.data
foo <- data.frame(do.call('rbind', strsplit(as.character(meta0$Tumor_Sample_Barcode),'_',fixed=TRUE)))
meta0$id=foo$X1


meta1=meta[which(meta$DI.Code%in%meta0$id),]
colnames(meta1)[1]=colnames(meta0)[1]

meta1$name=meta0$Tumor_Sample_Barcode
meta1d=meta1[,-c(16,20,23,16,29,32,35,38)]

meta1d=meta1d[,c(32,1:31)]

colnames(meta1d)[1:2]=c("Tumor_Sample_Barcode","id")


tmb.meta=merge(tmb, meta1d,by="Tumor_Sample_Barcode")

tmb.meta$Tissue.Label..1.upfront..2.recurrence.=as.factor(tmb.meta$Tissue.Label..1.upfront..2.recurrence.)
tmb.meta$PFS.Event=as.factor(tmb.meta$PFS.Event)
tmb.meta$DMFS.Event=as.factor(tmb.meta$DMFS.Event)
tmb.meta$P16.Status.Upfront..OPX.only.=as.factor(tmb.meta$P16.Status.Upfront..OPX.only.)
library(ggpubr)
pdf("tmb_tissue_label.pdf")
ggplot(tmb.meta,aes(x=Tissue.Label..1.upfront..2.recurrence.,y=total_perMB))+geom_boxplot()+  geom_jitter(shape=16, position=position_jitter(0.2))+ stat_compare_means()
ggplot(tmb.meta,aes(x=PFS.Event,y=total_perMB))+geom_boxplot()+   geom_jitter(shape=16, position=position_jitter(0.2))+ stat_compare_means()
ggplot(tmb.meta,aes(x=DMFS.Event,y=total_perMB))+geom_boxplot()+  geom_jitter(shape=16, position=position_jitter(0.2))+ stat_compare_means()
ggplot(tmb.meta,aes(x=P16.Status.Upfront..OPX.only.,y=total_perMB))+geom_boxplot()+  geom_jitter(shape=16, position=position_jitter(0.2))+ stat_compare_means()

dev.off()

library(data.table)
meta.f2=setDT(meta1d)

soma@clinical.data=meta.f2


library(survminer)
library(survival)
# vector with the variables to run through
genes <-soma@gene.summary$Hugo_Symbol
soma_matrix=data.frame(mutCountMatrix(maf = soma,removeNonMutated = F))

soma_matrix.f=soma_matrix %>% mutate_if(is.numeric, ~1 * (. > 0))
numeric_matrix=soma_matrix.f
meta=data.frame(soma@clinical.data)
meta$Time.between.initial.RT.and.recurrence..mo.=as.numeric(meta$Time.between.initial.RT.and.recurrence..mo.)

meta$Time.to.DM..death..or.last.FU=as.numeric(meta$Time.to.DM..death..or.last.FU)
meta$Time.to.local.failure.progression=as.numeric(meta$Time.to.local.failure.progression)
meta$OS.Time.to.death.or.last.FU=as.numeric(meta$OS.Time.to.death.or.last.FU)
meta$Age=as.numeric(meta$Age)
meta$Re.RT.Dose=as.numeric(meta$Re.RT.Dose)
meta$GTV.Volume=as.numeric(meta$GTV.Volume)
meta$CTV.Volume=as.numeric(meta$CTV.Volume)



  
coxph(as.formula(paste0("Surv(PFS.Time.to.any.progression,PFS.Event)~",genes[i],"+PFS.Event")),data=tp53.ann)

for(i in 1:nrow(numeric_matrix)){
  
  tp53=data.frame(t(numeric_matrix[grep(genes[i],rownames(numeric_matrix)),]))
  
  tp53.ann=merge(meta,tp53,by.x="Tumor_Sample_Barcode",by.y=0)
  tp53.ann1=tp53.ann[which(tp53.ann$Concurrent.chemotherapy==1),]
  
  fit <- survfit(as.formula(paste0("Surv(PFS.Time.to.any.progression,PFS.Event)~", as.character(genes[i]))),
                 data = tp53.ann1)

  
  pdf(paste0(genes[i], "_PFS.Survival_w_vs_wo_mutation.pdf"))
  p1=   ggsurvplot(fit, conf.int=F, pval=TRUE, risk.table=TRUE, 
                   legend.labs=c(paste(genes[i],"_WO_mutation",sep="_"), paste(genes[i],"_W_mutation",sep="_")), legend.title="Mutation",  
                   palette=c("dodgerblue2", "orchid2","grey"), 
                   title = "Status", subtitle = "Based on Kaplan-Meier estimates",
                   
                   font.title = c(10, "bold", "darkblue"),
                   font.subtitle = c(10, "bold.italic", "purple"),
                   font.caption = c(10, "plain", "orange"), 
                   risk.table.height=0.2,ggtheme = theme_minimal())
  
  
  print(p1)
  
  dev.off()}


for(i in 1:nrow(numeric_matrix)){
  
  tp53=data.frame(t(numeric_matrix[grep(genes[i],rownames(numeric_matrix)),]))
  
  tp53.ann=merge(meta,tp53,by.x="Tumor_Sample_Barcode",by.y=0)
  tp53.ann1=tp53.ann[which(tp53.ann$Concurrent.chemotherapy==0),]
  
  fit <- survfit(as.formula(paste0("Surv(PFS.Time.to.any.progression,PFS.Event)~", as.character(genes[i]))),
                 data = tp53.ann1)
  
  
  pdf(paste0(genes[i], "_PFS.Survival_w_vs_wo_mutation.pdf"))
  p1=   ggsurvplot(fit, conf.int=F, pval=TRUE, risk.table=TRUE, 
                   legend.labs=c(paste(genes[i],"_WO_mutation",sep="_"), paste(genes[i],"_W_mutation",sep="_")), legend.title="Mutation",  
                   palette=c("dodgerblue2", "orchid2","grey"), 
                   title = "Status", subtitle = "Based on Kaplan-Meier estimates",
                   
                   font.title = c(10, "bold", "darkblue"),
                   font.subtitle = c(10, "bold.italic", "purple"),
                   font.caption = c(10, "plain", "orange"), 
                   risk.table.height=0.2,ggtheme = theme_minimal())
  
  
  print(p1)
  
  dev.off()}


prog_geneset = survGroup(maf = soma, genes=soma@gene.summary$Hugo_Symbol, geneSetSize = 2, time = "PFS.Time.to.any.progression", Status = "PFS.Event", verbose = FALSE)
pdf('surv.geneset2.pdf')
mafSurvGroup(maf = soma, geneSet = c("KMT2D","ARID1A"), time = "PFS.Time.to.any.progression", Status = "PFS.Event")
dev.off()

soma_matrix.surv=numeric_matrix[which(rownames(numeric_matrix)%in%unique(c(comp1f$Hugo_Symbol,comp2f$Hugo_Symbol))),]
tp53.ann=merge(meta,t(soma_matrix.surv),by.x="Tumor_Sample_Barcode",by.y=0)


genes=unique(c(comp1f$Hugo_Symbol,comp2f$Hugo_Symbol))

for(i in 1:nrow(soma_matrix.surv)){
  
  tp53=data.frame(t(numeric_matrix[grep(genes[i],rownames(numeric_matrix)),]))
  
  tp53.ann=merge(meta,tp53,by.x="Tumor_Sample_Barcode",by.y=0)
  tp53.ann1=tp53.ann[which(tp53.ann$PFS.Event==1),]
  
  fit <- survfit(as.formula(paste0("Surv(PFS.Time.to.any.progression,PFS.Event)~", as.character(genes[i]))),
                 data = tp53.ann1)
  
  
  pdf(paste0(genes[i], "_PFS.Survival_w_vs_wo_mutation.pdf"))
  p1=   ggsurvplot(fit, conf.int=F, pval=TRUE, risk.table=TRUE, 
                   legend.labs=c(paste(genes[i],"_WO_mutation",sep="_"), paste(genes[i],"_W_mutation",sep="_")), legend.title="Mutation",  
                   palette=c("dodgerblue2", "orchid2","grey"), 
                   title = "Status", subtitle = "Based on Kaplan-Meier estimates",
                   
                   font.title = c(10, "bold", "darkblue"),
                   font.subtitle = c(10, "bold.italic", "purple"),
                   font.caption = c(10, "plain", "orange"), 
                   risk.table.height=0.2,ggtheme = theme_minimal())
  
  
  print(p1)
  
  dev.off()}




for(i in 1:nrow(numeric_matrix)){
  
  tp53=data.frame(t(numeric_matrix[grep(genes[i],rownames(numeric_matrix)),]))
  
  tp53.ann=merge(meta,tp53,by.x="Tumor_Sample_Barcode",by.y=0)
  tp53.ann1=tp53.ann[which(tp53.ann$Salvage.Surgery==0),]
  
  fit <- survfit(as.formula(paste0("Surv(PFS.Time.to.any.progression,PFS.Event)~", as.character(genes[i]))),
                 data = tp53.ann1)
  
  
  pdf(paste0(genes[i], "_PFS.Survival_w_vs_wo_mutation.pdf"))
  p1=   ggsurvplot(fit, conf.int=F, pval=TRUE, risk.table=TRUE, 
                   legend.labs=c(paste(genes[i],"_WO_mutation",sep="_"), paste(genes[i],"_W_mutation",sep="_")), legend.title="Mutation",  
                   palette=c("dodgerblue2", "orchid2","grey"), 
                   title = "Status", subtitle = "Based on Kaplan-Meier estimates",
                   
                   font.title = c(10, "bold", "darkblue"),
                   font.subtitle = c(10, "bold.italic", "purple"),
                   font.caption = c(10, "plain", "orange"), 
                   risk.table.height=0.2,ggtheme = theme_minimal())
  
  
  print(p1)
  
  dev.off()}

laml.sig = oncodrive(maf = soma, AACol = 'Amino_acids', minMut = 3, pvalMethod = 'zscore')

plotOncodrive(res = laml.sig, fdrCutOff = 0.1, useFraction = T)


prog_geneset = survGroup(maf = soma, genes = soma@gene.summary$Hugo_Symbol, geneSetSize = 1, time = "PFS.Time.to.any.progression", Status = "PFS.Event", verbose = FALSE)

prog.f1=prog_geneset[which(prog_geneset$P_value<0.1),]

pdf("survival_gene_set1.PFS.pdf")

for (i in 1:nrow(prog.f1)){
  mafSurvGroup(maf = soma, geneSet = prog.f1$Gene_combination[i],time = "PFS.Time.to.any.progression", Status = "PFS.Event")
}
dev.off()


pdf("rnx1.pdf")
oncoplot(soma,genes=c("ZNRF3","RNX1"))
dev.off()


####survival plot of PFS samples with mutation odd ratio >2 

meta=read.table("20230511_GENOMIC DATA.txt",sep="\t",header=T)

meta0=soma@clinical.data
foo <- data.frame(do.call('rbind', strsplit(as.character(meta0$Tumor_Sample_Barcode),'_',fixed=TRUE)))
meta0$id=foo$X1
library(dplyr)

meta1=meta[which(meta$DI.Code%in%meta0$id),]
colnames(meta1)[1]=colnames(meta0)[1]

meta1$name=meta0$Tumor_Sample_Barcode
meta1d=meta1[,-c(16,20,23,16,29,32,35,38)]

meta1d=meta1d[,c(32,1:31)]

colnames(meta1d)[1:2]=c("Tumor_Sample_Barcode","id")

library(data.table)
meta.f2=setDT(meta1d)

soma@clinical.data=meta.f2


soma_matrix=data.frame(mutCountMatrix(maf = soma,removeNonMutated = F))

soma_matrix.f=soma_matrix %>% mutate_if(is.numeric, ~1 * (. > 0))

meta=data.frame(soma@clinical.data)
meta=meta[order(meta$Age),]

#gene.use=soma@gene.summary$Hugo_Symbol[which(soma@gene.summary$Hugo_Symbol%in%unique(pw$gene))]

meta.plot=data.frame(soma@clinical.data)



sub1=subsetMaf(soma,clinQuery="PFS.Event==1")

sub2=subsetMaf(soma,clinQuery="PFS.Event==0")


pt.vs.rt <- mafCompare(m1 = sub1, m2 = sub2, m1Name = 'PFS.Event==1', m2Name = 'PFS.Event==0', minMut = 0)

comp1=pt.vs.rt$results[which(pt.vs.rt$results$or=="Inf"),]
comp1f=comp1[which(comp1$`PFS.Event==1`>2),]

comp2=pt.vs.rt$results[which(pt.vs.rt$results$or>=2),]
comp2f=comp2[which(comp2$`PFS.Event==1`>2),]

pt.vs.rt <- mafCompare(m2 = sub1, m1 = sub2, m1Name = 'PFS.Event==0', m2Name = 'PFS.Event==1', minMut = 2)
comp3=pt.vs.rt$results[which(pt.vs.rt$results$or=="Inf"),]
comp3f=comp1[which(comp3$`PFS.Event==1`>2),]


write.table(pt.vs.rt$results,"mafcomp_stat_pfs.txt",sep="\t",quote=F)

library(survminer)
library(survival)
# vector with the variables to run through
genes <- rownames(soma_matrix.f)



col = c("Missense_Mutation" = "blue", "In_Frame_Ins" = "red", "Splice_Site" = "pink","In_Frame_Del"="purple","Frame_Shift_Ins"="yellow",
        "Frame_Shift_Del"="green", "Nonsense_Mutation"="brown", "Multi_Hit"="orange")




soma_matrix=data.frame(mutCountMatrix(maf = soma,removeNonMutated = F))

soma_matrix.f=soma_matrix %>% mutate_if(is.numeric, ~1 * (. > 0))

write.table(soma_matrix.f,"soma.f_numeric.matrix",sep='\t',quote=F)


numeric_matrix=soma_matrix.f


soma_matrix.surv=soma_matrix.f[which(rownames(soma_matrix.f)%in%unique(c(comp1f$Hugo_Symbol,comp2f$Hugo_Symbol))),]
sum1=data.frame(colSums(soma_matrix.surv))
names(sum1)="total_mutation"

#sum1=sum1 %>% mutate_if(is.numeric, ~1 * (. > 0))

tp53.ann0=merge(meta,sum1,by.x="Tumor_Sample_Barcode",by.y=0)
tp53.ann=merge(meta,t(soma_matrix.surv),by.x="Tumor_Sample_Barcode",by.y=0)


#tp53.ann0$PFS.Event=as.factor(tp53.ann0$PFS.Event)

tp53.ann0$PFS.Event=as.numeric(tp53.ann0$PFS.Event)



pw_col = RColorBrewer::brewer.pal(n = 12, name = 'Set3')[1:12]

sfita <- survfit(as.formula(paste0("Surv(PFS.Time.to.any.progression,PFS.Event) ~", "total_mutation")), data=tp53.ann0[which(tp53.ann0$PFS.Event==1),])


sfit0 <- survfit(as.formula(paste0("Surv(PFS.Time.to.any.progression,PFS.Event) ~", "total_mutation")), data=tp53.ann0)




tp53.ann0$level <- rep(NA, nrow(tp53.ann0))

tp53.ann0[tp53.ann0$total_mutation<=2, ][, "level"] <- "less mutation"
tp53.ann0[tp53.ann0$total_mutation>2, ][, "level"] <- "more mutation"


pw_col = RColorBrewer::brewer.pal(n = 12, name = 'Paired')[1:12]

sfit0 <- survfit(as.formula(paste0("Surv(PFS.Time.to.any.progression,PFS.Event) ~", "total_mutation")), data=tp53.ann0)


sfit <- survfit(as.formula(paste0("Surv(PFS.Time.to.any.progression,PFS.Event) ~", "level")), data=tp53.ann0)


pdf("~/Desktop/total.surv.pfs.pdf")

p1=   ggsurvplot(sfit, conf.int=F, pval=TRUE, risk.table=TRUE, 
                  legend.labs=c(paste("less_mutation","<=2",sep="_"), paste("more_mutation",">2",sep="_")), legend.title="Mutation",  
                 palette=c("blue","red"), 
                 title = "Status", subtitle = "Based on Kaplan-Meier estimates",
                 
                 font.title = c(10, "bold", "darkblue"),
                 font.subtitle = c(10, "bold.italic", "purple"),
                 font.caption = c(10, "plain", "orange"), 
                 risk.table.height=0.2,ggtheme = theme_minimal())


print(p1)


dev.off()

