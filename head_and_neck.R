library(maftools)
setwd("~/Desktop/head_neck/")

soma=readRDS("soma.filter.rds")

col = c("Missense_Mutation" = "blue", "In_Frame_Ins" = "red", "Splice_Site" = "#008000","In_Frame_Del"="purple","Frame_Shift_Ins"="yellow",
        "Frame_Shift_Del"="cyan", "Nonsense_Mutation"="brown", "Multi_Hit"="orange","Translation_Start_Site"="pink")




pdf("all_maf.pdf",10,30)
oncoplot(soma,removeNonMutated = F,writeMatrix = T,top=50, clinicalFeatures ="Local.Failure.Progression")
dev.off()



rate1 = clinicalEnrichment(maf = soma, clinicalFeature = 'Recurrence.Location')
rate1$pairwise_comparision[fdr < 0.05]

rate1$groupwise_comparision[p_value < 0.05]


soma_matrix=data.frame(mutCountMatrix(maf = soma,removeNonMutated = F))

soma_matrix.f=soma_matrix %>% mutate_if(is.numeric, ~1 * (. > 0))

write.table(soma_matrix.f,"soma_numeric.matrix",sep='\t',quote=F)


meta=data.frame(soma@clinical.data)

meta.sub1=meta[-c(4,10,14,16:18,20,22,24,26,28)]
meta.sub1.m=melt(meta[-c(4,10,14,16:18,20,22,24,26,28)],id.vars="Tumor_Sample_Barcode")



pdf("clinical.pdf")
ggplot(meta.sub1.m,aes(x=variable,y=value))+geom_boxplot()+facet_wrap(~variable)
dev.off()

meta.sub2=meta[c(1,4,10,14,16:18,20,22,24,26,28)]

df=meta.sub1




gapminder=meta.sub2

library("vtable")
gapminder$Tumor_Sample_Barcode=as.factor(gapminder$Tumor_Sample_Barcode)

st(gapminder,group = 'Tumor_Sample_Barcode', group.long = TRUE, factor.counts = T)


sumtable(gapminder)

vartable <- vtable(gapminder,out='return')

rownames(gapminder)=gapminder$ID

df=gapminder[-1]

#install.packages("gtsummary")

library(gtsummary)

X<-split(df, df$index)




meta.sub2.m=melt(meta[c(10,14,16:18,20,22,24,26,28)],id.vars="Tumor_Sample_Barcode")

pdf("clinical.pdf")
ggplot(meta.m,aes(x=variable,y=value))+geom_bar(aes(fill=variable), position = "dodge", stat="identity")
dev.off()


###progression################################################

sub1=subsetMaf(soma,clinQuery="Locoregional.Failure.Progression==0")
sub2=subsetMaf(soma,clinQuery="Locoregional.Failure.Progression==1")

mafCompare(sub1,sub2)


pt.vs.rt <- mafCompare(m1 = sub1, m2 = sub2, m2Name = 'Locoregional.Failure.Progression==0', m1Name = 'Locoregional.Failure.Progression==1', minMut = 0)


comp2=pt.vs.rt$results[which(pt.vs.rt$results$or>=2),]
comp2f=comp2[which(comp2$`Locoregional.Failure.Progression==1`>2),]


##compare local vs progression

col = c("Missense_Mutation" = "blue", "In_Frame_Ins" = "red", "Splice_Site" = "#008000","In_Frame_Del"="purple","Frame_Shift_Ins"="yellow",
        "Frame_Shift_Del"="cyan", "Nonsense_Mutation"="brown", "Multi_Hit"="orange","Translation_Start_Site"="pink")






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
  Translation_Start_Site= function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["Translation_Start_Site"], col = NA))},
  
  Multi_Hit= function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["Multi_Hit"], col = NA))}
)




meta.lf=  meta[which(meta$Locoregional.Failure.Progression ==0),]

meta.pr=  meta[which(meta$Locoregional.Failure.Progression  ==1),]


input_matrix=read.table("soma_onco_matrix.txt",sep="\t",header=T,row.names=1,check.names=F)

lf=input_matrix[which(colnames(input_matrix)%in%meta.lf$Tumor_Sample_Barcode)]
pr=input_matrix[which(colnames(input_matrix)%in%meta.pr$Tumor_Sample_Barcode)]


lf.s=lf[,order(match(colnames(lf),meta.lf$Tumor_Sample_Barcode))]
pr.s=pr[,order(match(colnames(pr),meta.pr$Tumor_Sample_Barcode))]

pdf("~/Desktop/somatic_Locoregional Failure.pdf",20,10)
t01=oncoPrint(lf.s[which(rownames(lf.s)%in%comp2f$Hugo_Symbol),],remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              #bottom_annotation = ha2,top_annotation=ha01,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="No L/R failure (21)")


t02=oncoPrint(pr.s[which(rownames(pr.s)%in%comp2f$Hugo_Symbol),], remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              #bottom_annotation = ha1,top_annotation=ha02,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="Local or regional failure (death/last FU) (15)")

t01+t02

dev.off()

######pathway

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



input_matrix.path=merge(pw,input_matrix,by.x="gene",by.y=0)

#input_matrix.path[input_matrix.path== ''] <- NA
#input_matrix.path$count=rowSums(is.na(input_matrix.path[-c(1:2)]))

#input_matrix.path1=input_matrix.path[which(input_matrix.path$count!=163),]

#input_matrix.path=input_matrix.path1[-ncol(input_matrix.path1)]


input_matrix.path.s=unique(input_matrix.path[order(input_matrix.path$pw),])

rownames(input_matrix.path.s)=input_matrix.path.s$gene


input_matrix.path.s[is.na(input_matrix.path.s)]<-""

#colnames(input_matrix.path.s)=gsub("-T","",colnames(input_matrix.path.s))

input_matrix.path.pr=input_matrix.path.s[which(colnames(input_matrix.path.s)%in%meta.pr$Tumor_Sample_Barcode)]
input_matrix.path.pr$pw=input_matrix.path.s$pw

input_matrix.path2.lf=input_matrix.path.s[which(colnames(input_matrix.path.s)%in%meta.lf$Tumor_Sample_Barcode)]
input_matrix.path2.lf$pw=input_matrix.path.s$pw

pdf("~/Desktop/path.Locoregional_Failure.pdf",20,20)

t01=oncoPrint(input_matrix.path2.lf[,-ncol(input_matrix.path2.lf)], remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,row_order =input_matrix.path.pr$gene,row_split=input_matrix.path.pr$pw, 
              #row_title_rot = switch(row_title_side[1], "left" = 0),
              
              row_title_gp = gpar(fontsize = 10) ,show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 5),column_title="No L/R failure (21)")
#t01s=oncoPrint(matrix.ncrpc[which(rownames(matrix.ncrpc)%in%unique(c(f4$gen1,f4$gene2,f5$gen1,f5$gene2))),], remove_empty_columns = FALSE, remove_empty_rows = T,
#             alter_fun = alter_fun, col = col,
#            show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 5),column_title="non-crpc(184(95F+86T))")

t01a=oncoPrint(input_matrix.path.pr[,-ncol(input_matrix.path.pr)], remove_empty_columns = FALSE, remove_empty_rows = F,
               alter_fun = alter_fun, col = col,row_order =input_matrix.path2.lf$gene,row_split=input_matrix.path2.lf$pw, 
               #row_title_rot = switch(row_title_side[1], "left" = 0),
               
               row_title_gp = gpar(fontsize = 10) ,show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 5),column_title="Local or regional failure (death/last FU) (15)")
#t01s=oncoPrint(matrix.ncrpc[which(rownames(matrix.ncrpc)%in%unique(c(f4$gen1,f4$gene2,f5$gen1,f5$gene2))),], remove_empty_columns = FALSE, remove_empty_rows = T,
t01+t01a

dev.off()


#input_matrix.path1=soma_matrix.f

input_matrix.path1=merge(pw,soma_matrix.f,by.x="gene",by.y=0)


df1=input_matrix.path1[-c(1,2)]
#df1a=replace(df1, is.na(df1), 0)
#df2a=replace(df1a, !is.na(df1), 1)

#df3a=df2a %>% mutate_if(is.character,as.numeric)


df1$pw=input_matrix.path1$pw

path.agg=aggregate(df1[,-ncol(df1)], by=list(pathway=df1$pw), FUN=sum)
rownames(path.agg)=path.agg$pathway

#colnames(path.agg)=gsub("-T","",colnames(path.agg))

pathway010.ann=path.agg[which(colnames(path.agg)%in%colnames(input_matrix.path.pr))]
pathway010.ann1=path.agg[which(colnames(path.agg)%in%colnames(input_matrix.path2.lf))]




library(dplyr)
pathway010.ann[pathway010.ann == 0] <- ""
pathway010.ann1[pathway010.ann1 == 0] <- ""

pathway010.ann[] <- lapply(pathway010.ann, function(x) as.numeric(as.character(x)))
pathway010.ann1[] <- lapply(pathway010.ann1, function(x) as.numeric(as.character(x)))
pathway010.ann=pathway010.ann %>% mutate_if(is.numeric, ~1 * (. != 0))
pathway010.ann1=pathway010.ann1 %>% mutate_if(is.numeric, ~1 * (. != 0))

pathway010.ann[is.na(pathway010.ann)]<-""
pathway010.ann1[is.na(pathway010.ann1)]<-""




pdf("~/Desktop/pathway.locoregional.failure.agg..v3.pdf",15,5)
t01=oncoPrint(pathway010.ann1, remove_empty_columns = FALSE, remove_empty_rows = F,top_annotation = NULL,
              alter_fun = alter_fun0, col = col0,show_heatmap_legend =F,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 5),column_title="No L/R failure (21)")
t01a=oncoPrint(pathway010.ann, remove_empty_columns = FALSE, remove_empty_rows = F,top_annotation = NULL,
               alter_fun = alter_fun0, col = col0,show_heatmap_legend =F,
               show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 5),column_title="Local or regional failure (death/last FU) (15)")


t01+t01a


dev.off()






####survival

prog_geneset = survGroup(maf = soma, genes = soma@gene.summary$Hugo_Symbol, geneSetSize = 1, time = "OS.Time.to.death.or.last.FU", Status = "OS.Event", verbose = FALSE)

prog.f1=prog_geneset[which(prog_geneset$P_value<0.1),]

pdf("~/Desktop/survival_gene_set1.stand.pdf")

for (i in 1:nrow(prog.f1)){
  mafSurvGroup(maf = soma, geneSet = prog.f1$Gene_combination[i],time = "OS.Time.to.death.or.last.FU", Status = "OS.Event")
}
dev.off()






library(survminer)
library(survival)
# vector with the variables to run through
genes <- rownames(soma_matrix.f)

numeric_matrix=soma_matrix.f

meta$Time.between.initial.RT.and.recurrence..mo.=as.numeric(meta$Time.between.initial.RT.and.recurrence..mo.)

meta$Time.to.DM..death..or.last.FU=as.numeric(meta$Time.to.DM..death..or.last.FU)
meta$Time.to.local.failure.progression=as.numeric(meta$Time.to.local.failure.progression)
meta$OS.Time.to.death.or.last.FU=as.numeric(meta$OS.Time.to.death.or.last.FU)
meta$Age=as.numeric(meta$Age)
meta$Re.RT.Dose=as.numeric(meta$Re.RT.Dose)
meta$GTV.Volume=as.numeric(meta$GTV.Volume)
meta$CTV.Volume=as.numeric(meta$CTV.Volume)


for(i in 1:nrow(numeric_matrix)){
  
  tp53=data.frame(t(numeric_matrix[grep(genes[i],rownames(numeric_matrix)),]))
  
  tp53.ann=merge(meta,tp53,by.x="Tumor_Sample_Barcode",by.y=0)
  
  fit <- survfit(as.formula(paste0("Surv(OS.Time.to.death.or.last.FU) ~", genes[i])),
                 data = tp53.ann)
  
  pdf(paste0(genes[i], "_Survival_w_vs_wo_mutation.pdf"))
  p1=   ggsurvplot(fit, conf.int=F, pval=TRUE, risk.table=TRUE, 
                   legend.labs=c(paste(genes[i],"_WO_mutation",sep="_"), paste(genes[i],"_W_mutation",sep="_")), legend.title="Mutation",  
                   palette=c("dodgerblue2", "orchid2","grey"), 
                   title = "Status", subtitle = "Based on Kaplan-Meier estimates",
                   
                   font.title = c(10, "bold", "darkblue"),
                   font.subtitle = c(10, "bold.italic", "purple"),
                   font.caption = c(10, "plain", "orange"), 
                   risk.table.height=0.2,ggtheme = theme_minimal())
  

  print(p1)
  
  

  dev.off()
}

