library(maftools);library(data.table)
setwd("combine_tempus_foundation/212")
all.maf=readRDS("pattern_failure/tempus_foundation_w_conflicP.417.s289.rds")
write.table(all.maf@data,"all_raw_maf_matrix.txt",sep='\t',quote=F)



meta=read.delim("pattern_failure/meta_417_corrected.v1.filter.txt",sep="\t",header=T)
meta.s=meta[order(meta$Tumor_Sample_Barcode),]
all.maf@clinical.data=setDT(meta.s)




setwd("pattern_failure/result_417")

test=mutCountMatrix(all.maf,removeNonMutated = F)

pdf("all.maf.pdf")
oncoplot(all.maf,removeNonMutated = F,writeMatrix = T,genes=all.maf@gene.summary$Hugo_Symbol)
dev.off()

met.maf=subsetMaf(all.maf,genes=all.maf@gene.summary$Hugo_Symbol,clinQuery = "Time..2.synchronous..1.metachronous==1")

meta=read.delim("pattern_failure/meta_417_corrected.v1.filter.txt",sep="\t",header=T)
meta.s=meta[order(meta$Tumor_Sample_Barcode),]

input_matrix=read.table("pattern_failure/result_417/onco_matrix.txt",sep="\t",row.names=1,check.names=F,header=T)


meta.metachro=meta.s[which(meta.s$Time..2.synchronous..1.metachronous==1),]





oligo=meta.metachro[which(meta.metachro$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.==1),]

poly=meta.metachro[which(meta.metachro$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.==2),]
nopro=meta.metachro[which(meta.metachro$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.==3),]

nopro.maf=subsetMaf(maf = met.maf, tsb=  nopro$id)
oligo.maf=subsetMaf(maf = met.maf, tsb=  oligo$id)
poly.maf=subsetMaf(maf = met.maf, tsb=  poly$id)

stat1=mafCompare(nopro.maf,oligo.maf,minMut=3)
stat2=mafCompare(nopro.maf,poly.maf,minMut=3)
stat3=mafCompare(poly.maf,oligo.maf,minMut=3)

input_matrix1=data.frame(t(input_matrix))

colnames(input_matrix)=gsub("-T","",colnames(input_matrix))
library(dplyr)


matrix.nopro=input_matrix[,which(colnames(input_matrix)%in%nopro$Tumor_Sample_Barcode)]

matrix.oligo=input_matrix[,which(colnames(input_matrix)%in%oligo$Tumor_Sample_Barcode)]

matrix.poly=input_matrix[,which(colnames(input_matrix)%in%poly$Tumor_Sample_Barcode)]

input_matrix.pf=cbind(matrix.nopro,matrix.oligo,matrix.poly)



matrix.nopro[is.na(matrix.nopro)]<-""

matrix.oligo[is.na(matrix.oligo)]<-""

matrix.poly[is.na(matrix.poly)]<-""



pdf("somatic_metachronous.Paterns.of.Failure.417.pdf",20,20)


t01=oncoPrint(matrix.nopro, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="no.progression.at.last.fu(48)")
#t01s=oncoPrint(matrix.ncrpc[which(rownames(matrix.ncrpc)%in%unique(c(f4$gen1,f4$gene2,f5$gen1,f5$gene2))),], remove_empty_columns = FALSE, remove_empty_rows = T,
#             alter_fun = alter_fun, col = col,
#            show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="non-crpc(184(95F+86T))")


t02=oncoPrint(matrix.poly, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="polyprogressor(38)")
#t02s=oncoPrint(matrix.crpc[which(rownames(matrix.crpc)%in%unique(c(f5$gen1,f5$gene2,f4$gen1,f4$gene2))),], remove_empty_columns = FALSE, remove_empty_rows = T,
#              alter_fun = alter_fun, col = col,
#             show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="crpc(94(77F+17F))")

t03=oncoPrint(matrix.oligo, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="oligoprogressor(74)")
t01+t03+t02

dev.off()





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





input_matrix.path=merge(pw,input_matrix.pf,by.x="gene",by.y=0)

input_matrix.path[input_matrix.path== ''] <- NA


input_matrix.path$count=rowSums(is.na(input_matrix.path[-c(1:2)]))

input_matrix.path1=input_matrix.path[which(input_matrix.path$count!=160),]

input_matrix.path=input_matrix.path1[-ncol(input_matrix.path1)]


input_matrix.path.s=input_matrix.path[order(input_matrix.path$pw),]

rownames(input_matrix.path.s)=input_matrix.path.s$gene


input_matrix.path.s[is.na(input_matrix.path.s)]<-""

input_matrix.path.nopro=input_matrix.path.s[which(colnames(input_matrix.path.s)%in%colnames(matrix.nopro))]
input_matrix.path.nopro$pw=input_matrix.path.s$pw

input_matrix.path.poly=input_matrix.path.s[which(colnames(input_matrix.path.s)%in%colnames(matrix.poly))]
input_matrix.path.poly$pw=input_matrix.path.s$pw

input_matrix.path.oligo=input_matrix.path.s[which(colnames(input_matrix.path.s)%in%colnames(matrix.oligo))]
input_matrix.path.oligo$pw=input_matrix.path.s$pw








pdf("path.new.pattern.failure.417.pdf",20,10)

t01=oncoPrint(input_matrix.path.nopro[-ncol(input_matrix.path.nopro)], remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,row_order =input_matrix.path.pelvic$gene,row_split=input_matrix.path.nopro$pw, 
              #row_title_rot = switch(row_title_side[1], "left" = 0),
              
              row_title_gp = gpar(fontsize = 6) ,show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="no.progression.at.last.fu(48)")
#t01s=oncoPrint(matrix.ncrpc[which(rownames(matrix.ncrpc)%in%unique(c(f4$gen1,f4$gene2,f5$gen1,f5$gene2))),], remove_empty_columns = FALSE, remove_empty_rows = T,
#             alter_fun = alter_fun, col = col,
#            show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="non-crpc(184(95F+86T))")


t02=oncoPrint(input_matrix.path.poly[-ncol(input_matrix.path.poly)], remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,row_order =input_matrix.path.node$gene,row_split=input_matrix.path,poly$pw, row_title_gp = gpar(fontsize = 6) ,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="polyprogressor(38)")
#t01s=oncoPrint(matrix.ncrpc[which(rownames(matrix.ncrpc)%in%unique(c(f4$gen1,f4$gene2,f5$gen1,f5$gene2))),], remove_empty_columns = FALSE, remove_empty_rows = T,
#             alter_fun = alter_fun, col = col,
#            show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="non-crpc(184(95F+86T))")

t03=oncoPrint(input_matrix.path.oligo[-ncol(input_matrix.path.oligo)], remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,row_order =input_matrix.path.bone$gene,row_split=input_matrix.path.oligo$pw, row_title_gp = gpar(fontsize = 6) ,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="oligoprogressor(74)")
#t01s=oncoPrint(matrix.ncrpc[which(rownames(matrix.ncrpc)%in%unique(c(f4$gen1,f4$gene2,f5$gen1,f5$gene2))),], remove_empty_columns = FALSE, remove_empty_rows = T,
#             alter_fun = alter_fun, col = col,
#            show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="non-crpc(184(95F+86T))")




t01+t03+t02
#t01s
#t02s
dev.off()






df1=input_matrix.path1[-c(1,2,ncol(input_matrix.path1))]
df1a=replace(df1, is.na(df1), 0)
df2a=replace(df1a, !is.na(df1), 1)

df3a=df2a %>% mutate_if(is.character,as.numeric)


df3a$pw=input_matrix.path1$pw

path.agg=aggregate(df3a[,-ncol(df3a)], by=list(pathway=df3a$pw), FUN=sum)
rownames(path.agg)=path.agg$pathway


pathway010.ann=path.agg[which(colnames(path.agg)%in%colnames(matrix.nopro))]
pathway020.ann=path.agg[which(colnames(path.agg)%in%colnames(matrix.poly))]
pathway030.ann=path.agg[which(colnames(path.agg)%in%colnames(matrix.oligo))]

pathway010.ann[pathway010.ann == 0] <- ""
pathway020.ann[pathway020.ann == 0] <- ""
pathway030.ann[pathway030.ann == 0] <- ""


pathway010.ann1=pathway010.ann %>% mutate_if(is.numeric, ~1 * (. != 0))
pathway020.ann1=pathway020.ann %>% mutate_if(is.numeric, ~1 * (. != 0))
pathway030.ann1=pathway030.ann %>% mutate_if(is.numeric, ~1 * (. != 0))


pdf("pathway.pattern.failure.agg.417.pdf",20,6)
t01=oncoPrint(pathway010.ann1, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun0, col = col0,show_heatmap_legend =F,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="no.progression.at.last.fu(48)")
t02=oncoPrint(pathway020.ann1, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun0, col = col0,show_heatmap_legend =F,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="polyprogressor(38)")
t03=oncoPrint(pathway030.ann1, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun0, col = col0,show_heatmap_legend =F,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="oligoprogressor(74)")

t01+t03+t02

dev.off()






###############################################################



ranks <-merge(pw,matrix.nopro,by.x="gene",by.y=0)

ranks1 <- merge(pw,matrix.oligo,by.x="gene",by.y=0)


ranks2 <- merge(pw,matrix.poly,by.x="gene",by.y=0)




library("dplyr")
library("tidyr")
rank01=ranks 
rownames(rank01)=ranks$gene

rank01[,-c(1:2)] <- as.integer(rank01[,-c(1:2)] !="")

rank01[rank01 == 0]<-""
rank01[is.na(rank01)]<-""

write.table(rank01,"path_noprog.txt",sep="\t",quote=F)

setwd("~/Desktop")

rank01[-c(1:2)] %>% mutate_if(is.character,as.numeric)
rank01[is.na(rank01)]<-""


rank02=ranks1
rownames(rank02)=ranks1$gene

rank02[,-c(1:2)] <- as.integer(rank02[,-c(1:2)] !="")

rank02[rank02 == 0]<-""

rank02[-c(1:2)] %>% mutate_if(is.character,as.numeric)
rank02[is.na(rank02)]<-""

write.table(rank02,"path_oligo.txt",sep="\t",quote=F)



rank03=ranks2
rownames(rank03)=ranks2$gene

rank03[,-c(1:2)] <- as.integer(rank03[,-c(1:2)] !="")

rank03[rank03 == 0]<-""
rank03[-c(1:2)] %>% mutate_if(is.character,as.numeric)
rank03[is.na(rank03)]<-""

write.table(rank03,"path_poly.txt",sep="\t",quote=F)




library(dplyr)
library(janitor)
library(tidyr)
library(purrr)

dat %>% 
  group_by(pw) %>%
  group_modify(~ bind_rows(., summarize(., across(where(is.numeric), sum)))) %>%
  ungroup %>%
  mutate(category = coalesce(category, paste("Total", size)),
         size = if_else(startsWith(category, "Total"), "", size)) %>%
  data.frame


pathway010=read.table("path_noprog.m.txt",sep="\t",header=T,row.names=2,check.names=F)
pathway020=read.table("path_oligo.m.txt",sep="\t",header=T,row.names=2,check.names=F)

pathway030=read.table("path_poly.m.txt",sep="\t",header=T,row.names=2,check.names=F)




pathway010[is.na(pathway010)]<-""
pathway020[is.na(pathway020)]<-""
pathway030[is.na(pathway030)]<-""
pathway010[pathway010==0]<-""
pathway020[pathway020==0]<-""
pathway030[pathway030==0]<-""


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
pdf("path.pdf")
path=data.frame(OncogenicPathways(maf = all.maf))
dev.off()

pathway010.ann=merge(path[,1:3],pathway010,by.x="Pathway",by.y=0)
rownames(pathway010.ann)=pathway010.ann$Pathway

#rownames(pathway010.ann)=paste(pathway010.ann$Pathway,"(",pathway010.ann$n_affected_genes,"/",pathway010.ann$N,")",sep="")

pathway020.ann=merge(path[,1:3],pathway020,by.x="Pathway",by.y=0)
rownames(pathway020.ann)=pathway020.ann$Pathway

#rownames(pathway020.ann)=paste(pathway020.ann$Pathway,"(",pathway020.ann$n_affected_genes,"/",pathway020.ann$N,")",sep="")

pathway030.ann=merge(path[,1:3],pathway030,by.x="Pathway",by.y=0)
#rownames(pathway030.ann)=paste(pathway030.ann$Pathway,"(",pathway030.ann$n_affected_genes,"/",pathway030.ann$N,")",sep="")
rownames(pathway030.ann)=pathway030.ann$Pathway

pdf("pathway.new.pattern.pdf",20,6)
t01=oncoPrint(pathway010.ann[-c(1:4)], remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun0, col = col0,show_heatmap_legend =F,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="metachronous.noprogression(46)")
t02=oncoPrint(pathway020.ann[-c(1:4)], remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun0, col = col0,show_heatmap_legend =F,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="metachronous.OLIGO(73))")
t03=oncoPrint(pathway030.ann[-c(1:4)], remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun0, col = col0,show_heatmap_legend =F,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="metachronous.POLY(36)))")

t01+t02+t03
dev.off()








#########################################







meta.sub=meta[which(meta$location.of.failure.1.Pelvic.Node.2.Distant.Node.3.Bone.4.visceral!="#N/A"),]


l1=meta.sub[which(meta.sub$location.of.failure.1.Pelvic.Node.2.Distant.Node.3.Bone.4.visceral==1),]

l2=meta.sub[which(meta.sub$location.of.failure.1.Pelvic.Node.2.Distant.Node.3.Bone.4.visceral==2),]
l3=meta.sub[which(meta.sub$location.of.failure.1.Pelvic.Node.2.Distant.Node.3.Bone.4.visceral==3),]
l4=meta.sub[which(meta.sub$location.of.failure.1.Pelvic.Node.2.Distant.Node.3.Bone.4.visceral==4),]


matrix.l1=all.matrix[which(all.matrix$Tumor_Sample_Barcode%in%l1$id),]


matrix.l2=all.matrix[which(all.matrix$Tumor_Sample_Barcode%in%l2$id),]

matrix.l3=all.matrix[which(all.matrix$Tumor_Sample_Barcode%in%l3$id),]


matrix.l4=all.matrix[which(all.matrix$Tumor_Sample_Barcode%in%l4$id),]


maf.l1=read.maf(matrix.l1)
maf.l2=read.maf(matrix.l2)

maf.l3=read.maf(matrix.l3)
maf.l4=read.maf(matrix.l4)


compare_result=mafCompare(maf.l1,maf.l2,minMut = 3,m1Name = 'pelvic.node',m2Name = 'Distant.Node')

compare_result1=mafCompare(maf.l1,maf.l3,minMut = 3,m1Name = 'pelvic.node',m2Name = 'Bone')


compare_result2=mafCompare(maf.l1,maf.l4,minMut = 3,m1Name = 'pelvic.node',m2Name = 'visceral')




col = c("Missense_Mutation" = "blue", "In_Frame_Ins" = "red", "Splice_Site" = "pink","In_Frame_Del"="purple","Frame_Shift_Ins"="yellow",
        "Frame_Shift_Del"="green", "Nonsense_Mutation"="brown", "Multi_Hit"="orange")




input_matrix=read.table("foundation_tempus_pathogen_comb.212.txt",sep="\t",header=T,row.names=1,check.names=F)
input_matrix[is.na(input_matrix)]<-""

gene=all.maf@gene.summary

gene.f=gene[which(gene$total>=3),]

input_matrix.f=input_matrix[which(rownames(input_matrix)%in%gene.f$Hugo_Symbol),]


matrix.l1=input_matrix[,which(colnames(input_matrix)%in%l1$id)]

matrix.l2=input_matrix[,which(colnames(input_matrix)%in%l2$id)]

matrix.l3=input_matrix[,which(colnames(input_matrix)%in%l3$id)]


matrix.l4=input_matrix[,which(colnames(input_matrix)%in%l4$id)]



matrix.l1[is.na(matrix.l1)]<-""

matrix.l2[is.na(matrix.l2)]<-""

matrix.l3[is.na(matrix.l3)]<-""

matrix.l4[is.na(matrix.l4)]<-""


pdf("somatic_tempus.fail_location_pathogenic.pdf",20,20)


t01=oncoPrint(matrix.l1, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="Pelvic.Node(18)")
#t01s=oncoPrint(matrix.ncrpc[which(rownames(matrix.ncrpc)%in%unique(c(f4$gen1,f4$gene2,f5$gen1,f5$gene2))),], remove_empty_columns = FALSE, remove_empty_rows = T,
#             alter_fun = alter_fun, col = col,
#            show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="non-crpc(184(95F+86T))")


t02=oncoPrint(matrix.l2, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="Distant.Node(30)")
#t02s=oncoPrint(matrix.crpc[which(rownames(matrix.crpc)%in%unique(c(f5$gen1,f5$gene2,f4$gen1,f4$gene2))),], remove_empty_columns = FALSE, remove_empty_rows = T,
#              alter_fun = alter_fun, col = col,
#             show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="crpc(94(77F+17F))")

t03=oncoPrint(matrix.l3, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="Bone(44)")
#t02s=oncoPrint(matrix.crpc[which(rownames(matrix.crpc)%in%unique(c(f5$gen1,f5$gene2,f4$gen1,f4$gene2))),], remove_empty_columns = FALSE, remove_empty_rows = T,
#              alter_fun = alter_fun, col = col,
#             show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="crpc(94(77F+17F))")

t04=oncoPrint(matrix.l4, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="visceral(14)")
#t02s=oncoPrint(matrix.crpc[which(rownames(matrix.crpc)%in%unique(c(f5$gen1,f5$gene2,f4$gen1,f4$gene2))),], remove_empty_columns = FALSE, remove_empty_rows = T,
#              alter_fun = alter_fun, col = col,
#             show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="crpc(94(77F+17F))")


t01+t02+t03+t04
#t01s
#t02s
dev.off()


ltc=meta.sub[which(meta.sub$Deek.Paterns.of.Failure..1.LTC..2.olgio..3.poly.0.blank.==1),]

oligo=meta.sub[which(meta.sub$Deek.Paterns.of.Failure..1.LTC..2.olgio..3.poly.0.blank.==2),]
poly=meta.sub[which(meta.sub$Deek.Paterns.of.Failure..1.LTC..2.olgio..3.poly.0.blank.==3),]




matrix.ltc=input_matrix[,which(colnames(input_matrix)%in%ltc$id)]

matrix.oligo=input_matrix[,which(colnames(input_matrix)%in%oligo$id)]

matrix.poly=input_matrix[,which(colnames(input_matrix)%in%poly$id)]



matrix.ltc[is.na(matrix.ltc)]<-""

matrix.oligo[is.na(matrix.oligo)]<-""

matrix.poly[is.na(matrix.poly)]<-""



pdf("somatic_Deek.Paterns.of.Failure.pdf",20,20)


t01=oncoPrint(matrix.ltc, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="LTC(153)")
#t01s=oncoPrint(matrix.ncrpc[which(rownames(matrix.ncrpc)%in%unique(c(f4$gen1,f4$gene2,f5$gen1,f5$gene2))),], remove_empty_columns = FALSE, remove_empty_rows = T,
#             alter_fun = alter_fun, col = col,
#            show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="non-crpc(184(95F+86T))")


t02=oncoPrint(matrix.oligo, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="OLIGO(34)")
#t02s=oncoPrint(matrix.crpc[which(rownames(matrix.crpc)%in%unique(c(f5$gen1,f5$gene2,f4$gen1,f4$gene2))),], remove_empty_columns = FALSE, remove_empty_rows = T,
#              alter_fun = alter_fun, col = col,
#             show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="crpc(94(77F+17F))")

t03=oncoPrint(matrix.poly, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="ploy(33)")
t01+t02+t03

dev.off()


