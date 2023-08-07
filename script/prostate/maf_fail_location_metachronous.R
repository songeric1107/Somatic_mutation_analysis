library(maftools);library(data.table)
setwd("combine_tempus_foundation/212")


all.maf=readRDS("pattern_failure/tempus_foundation_w_conflicP.417.s289.rds")
write.table(all.maf@data,"all_raw_maf_matrix.txt",sep='\t',quote=F)



meta=read.delim("meta_417_corrected.v1.filter.txt",sep="\t",header=T)
meta.s=meta[order(meta$Tumor_Sample_Barcode),]
all.maf@clinical.data=setDT(meta.s)


input_matrix=read.table("onco_matrix.txt",sep="\t",row.names=1,check.names=F,header=T)


meta.metachro=meta.s[which(meta.s$Time..2.synchronous..1.metachronous==1),]
rownames(meta.metachro)=meta.metachro$Tumor_Sample_Barcode       

gene.interest=input_matrix[which(rownames(input_matrix)%in%c("SPOP","APC","CTNNB1","RNF43")),]


`%notin%` <- Negate(`%in%`)

meta.metachro$Pelvic.Node.fail.location=gsub("1","Pelvic.Node",meta.metachro$Pelvic.Node.fail.location)
meta.metachro$Distant.Node=gsub("1","Distant.Node",meta.metachro$Distant.Node)
meta.metachro$Bone=gsub("1","Bone",meta.metachro$Bone)
meta.metachro$Visceral=gsub("1","Visceral",meta.metachro$Visceral)




pel.meta=meta.metachro[which(meta.metachro$Pelvic.Node.fail.location%in%c(0,"Pelvic.Node")),]
dis.meta=meta.metachro[which(meta.metachro$Distant.Node%in%c(0,"Distant.Node")),]
bone.meta=meta.metachro[which(meta.metachro$Bone%in%c(0,"Bone")),]
vis.meta=meta.metachro[which(meta.metachro$Visceral%in%c(0,"Visceral")),]

pel.meta1=meta.metachro[which(meta.metachro$Pelvic.Node.fail.location%in%c(0)),]
pel.meta2=meta.metachro[which(meta.metachro$Pelvic.Node.fail.location%in%c("Pelvic.Node")),]
dis.meta1=meta.metachro[which(meta.metachro$Distant.Node%in%c(0)),]
dis.meta2=meta.metachro[which(meta.metachro$Distant.Node%in%c("Distant.Node")),]

bone.meta1=meta.metachro[which(meta.metachro$Bone%in%c(0)),]
bone.meta2=meta.metachro[which(meta.metachro$Bone%in%c("Bone")),]
vis.meta1=meta.metachro[which(meta.metachro$Visceral%in%c(0)),]
vis.meta2=meta.metachro[which(meta.metachro$Visceral%in%c("Visceral")),]

meta.maf=subsetMaf(all.maf,genes=all.maf@gene.summary$Hugo_Symbol,clinQuery = "Time..2.synchronous..1.metachronous==1")


pel.maf=subsetMaf(maf = meta.maf, tsb=  as.character(pel.meta$Tumor_Sample_Barcode),mafObj=T)
dis.maf=subsetMaf(maf = meta.maf, tsb=  as.character(dis.meta$Tumor_Sample_Barcode),mafObj=T)
bone.maf=subsetMaf(maf = meta.maf, tsb=  as.character(bone.meta$Tumor_Sample_Barcode),mafObj=T)
vis.maf=subsetMaf(maf = meta.maf, tsb=  as.character(vis.meta$Tumor_Sample_Barcode),mafObj=T)


pel.maf1=subsetMaf(maf = meta.maf, tsb=  pel.meta1$Tumor_Sample_Barcode,mafObj=T)
dis.maf1=subsetMaf(maf = meta.maf, tsb=  dis.meta1$Tumor_Sample_Barcode,mafObj=T)
bone.maf1=subsetMaf(maf = meta.maf, tsb= bone.meta1$Tumor_Sample_Barcode,mafObj=T)
vis.maf1=subsetMaf(maf = meta.maf, tsb=  vis.meta1$Tumor_Sample_Barcode,mafObj=T)


pel.maf2=subsetMaf(maf = meta.maf, tsb=  pel.meta2$Tumor_Sample_Barcode,mafObj=T)
dis.maf2=subsetMaf(maf = meta.maf, tsb= dis.meta2$Tumor_Sample_Barcode,mafObj=T)
bone.maf2=subsetMaf(maf = meta.maf, tsb=  bone.meta2$Tumor_Sample_Barcode,mafObj=T)
vis.maf2=subsetMaf(maf = meta.maf, tsb= vis.meta2$Tumor_Sample_Barcode,mafObj=T)

stat1=mafCompare(pel.maf1,pel.maf2,minMut=3)
stat2=mafCompare(dis.maf1,dis.maf2,minMut=3)
stat3=mafCompare(bone.maf1,bone.maf2,minMut=3)
stat4=mafCompare(vis.maf1,vis.maf2,minMut=3)


input_matrix.location=input_matrix[colnames(input_matrix)%in%unique(c(rownames(pel.meta),rownames(dis.meta),rownames(bone.meta),rownames(vis.meta)))]

matrix.pel=input_matrix[,which(colnames(input_matrix)%in%pel.meta$Tumor_Sample_Barcode)]

matrix.dist=input_matrix[,which(colnames(input_matrix)%in%dis.meta$Tumor_Sample_Barcode)]

matrix.bone=input_matrix[,which(colnames(input_matrix)%in%bone.meta$Tumor_Sample_Barcode)]

matrix.vis=input_matrix[,which(colnames(input_matrix)%in%vis.meta$Tumor_Sample_Barcode)]




pel.sub1=matrix.pel[,which(colnames(matrix.pel)%in%rownames(pel.meta[which(pel.meta$Pelvic.Node.fail.location==0),]))]
pel.sub2=matrix.pel[,which(colnames(matrix.pel)%in%rownames(pel.meta[which(pel.meta$Pelvic.Node.fail.location=="Pelvic.Node"),]))]

dist.sub1=matrix.dist[,which(colnames(matrix.dist)%in%rownames(dis.meta[which(dis.meta$Distant.Node==0),]))]
dist.sub2=matrix.dist[,which(colnames(matrix.dist)%in%rownames(dis.meta[which(dis.meta$Distant.Node=="Distant.Node"),]))]



bone.sub1=matrix.bone[,which(colnames(matrix.bone)%in%rownames(bone.meta[which(bone.meta$Bone==0),]))]
bone.sub2=matrix.bone[,which(colnames(matrix.bone)%in%rownames(bone.meta[which(bone.meta$Bone=="Bone"),]))]

vis.sub1=matrix.vis[,which(colnames(matrix.vis)%in%rownames(vis.meta[which(vis.meta$Visceral==0),]))]
vis.sub2=matrix.vis[,which(colnames(matrix.vis)%in%rownames(vis.meta[which(vis.meta$Visceral=="Visceral"),]))]







col = c("Missense_Mutation" = "blue", "In_Frame_Ins" = "red", "Splice_Site" = "pink","In_Frame_Del"="purple","Frame_Shift_Ins"="yellow",
        "Frame_Shift_Del"="green", "Nonsense_Mutation"="brown", "Multi_Hit"="orange")
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



pdf("~/Desktop/location.failure.oncoprint.409.v1.pdf",20,20)


t00=oncoPrint(pel.sub1, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize =5),column_title="Location_of_failure: no Pelvic Node (137)")
#t01s=oncoPrint(matrix.ncrpc[which(rownames(matrix.ncrpc)%in%unique(c(f4$gen1,f4$gene2,f5$gen1,f5$gene2))),], remove_empty_columns = FALSE, remove_empty_rows = T,
#             alter_fun = alter_fun, col = col,
#            show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 5),column_title="non-crpc(184(95F+86T))")
t01=oncoPrint(pel.sub2, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize =5),column_title="Location_of_failure:Pelvic Node (40)")
#t01

t02=oncoPrint(dist.sub1, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize =5),column_title="Location_of_failure:Distant Node no  (137)")
#t01s=oncoPrint(matrix.ncrpc[which(rownames(matrix.ncrpc)%in%unique(c(f4$gen1,f4$gene2,f5$gen1,f5$gene2))),], remove_empty_columns = FALSE, remove_empty_rows = T,
#             alter_fun = alter_fun, col = col,
#            show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 5),column_title="non-crpc(184(95F+86T))")

t03=oncoPrint(dist.sub2, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize =5),column_title="Location_of_failure:Distant Node yes (40)")


t04=oncoPrint(bone.sub1, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize =5),column_title="Location_of_failure:Bone no (129)")
#t01s=oncoPrint(matrix.ncrpc[which(rownames(matrix.ncrpc)%in%unique(c(f4$gen1,f4$gene2,f5$gen1,f5$gene2))),], remove_empty_columns = FALSE, remove_empty_rows = T,
#             alter_fun = alter_fun, col = col,
#            show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 5),column_title="non-crpc(184(95F+86T))")

t05=oncoPrint(bone.sub2, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize =5),column_title="Location_of_failure:Bone yes (50)")

t06=oncoPrint(vis.sub1, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize =5),column_title="Location_of_failure:Visceral no (163)")
#t01s=oncoPrint(matrix.ncrpc[which(rownames(matrix.ncrpc)%in%unique(c(f4$gen1,f4$gene2,f5$gen1,f5$gene2))),], remove_empty_columns = FALSE, remove_empty_rows = T,
#             alter_fun = alter_fun, col = col,
#            show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 5),column_title="non-crpc(184(95F+86T))")

t07=oncoPrint(vis.sub2, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize =5),column_title="Location_of_failure:Visceral yes (14)")
#t01s=oncoPrint(matrix.ncrpc[which(rownames(matrix.ncrpc)%in%unique(c(f4$gen1,f4$gene2,f5$gen1,f5$gene2))),], remove_empty_columns = FALSE, remove_empty_rows = T,
#             alter_fun = alter_fun, col = col,
#          

t00+t01
t02+t03
t04+t05
t06+t07
#t01s
#t02s
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




input_matrix.path=merge(pw,input_matrix.location,by.x="gene",by.y=0)

input_matrix.path[input_matrix.path== ''] <- NA
input_matrix.path$count=rowSums(is.na(input_matrix.path[-c(1:2)]))
input_matrix.path1=input_matrix.path[which(input_matrix.path$count!=179),]

input_matrix.path=input_matrix.path1[-ncol(input_matrix.path1)]


input_matrix.path.s=input_matrix.path[order(input_matrix.path$pw),]

rownames(input_matrix.path.s)=input_matrix.path.s$gene


input_matrix.path.s[is.na(input_matrix.path.s)]<-""

input_matrix.path.pelvic=input_matrix.path.s[which(colnames(input_matrix.path.s)%in%colnames(pel.sub1))]
input_matrix.path.pelvic$pw=input_matrix.path.s$pw

input_matrix.path2.pelvic=input_matrix.path.s[which(colnames(input_matrix.path.s)%in%colnames(pel.sub2))]
input_matrix.path2.pelvic$pw=input_matrix.path.s$pw


input_matrix.path.node=input_matrix.path.s[which(colnames(input_matrix.path.s)%in%colnames(dist.sub1))]
input_matrix.path.node$pw=input_matrix.path.s$pw

input_matrix.path2.node=input_matrix.path.s[which(colnames(input_matrix.path.s)%in%colnames(dist.sub2))]
input_matrix.path2.node$pw=input_matrix.path.s$pw




input_matrix.path.bone=input_matrix.path.s[which(colnames(input_matrix.path.s)%in%colnames(bone.sub1))]
input_matrix.path.bone$pw=input_matrix.path.s$pw

input_matrix.path2.bone=input_matrix.path.s[which(colnames(input_matrix.path.s)%in%colnames(bone.sub2))]
input_matrix.path2.bone$pw=input_matrix.path.s$pw



input_matrix.path.vis=input_matrix.path.s[which(colnames(input_matrix.path.s)%in%colnames(vis.sub1))]
input_matrix.path.vis$pw=input_matrix.path.s$pw

input_matrix.path2.vis=input_matrix.path.s[which(colnames(input_matrix.path.s)%in%colnames(vis.sub2))]
input_matrix.path2.vis$pw=input_matrix.path.s$pw





pdf("pathway_gene_location_failure.pdf",10,10)

t01=oncoPrint(input_matrix.path.pelvic[-ncol(input_matrix.path.pelvic)], remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,row_order =input_matrix.path.pelvic$gene,row_split=input_matrix.path.pelvic$pw, 
              #row_title_rot = switch(row_title_side[1], "left" = 0),
              
              row_title_gp = gpar(fontsize = 6) ,show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),
              column_names_gp = gpar(fontsize = 5),column_title="Location_of_failure: no Pelvic Node (137)")
#t01s=oncoPrint(matrix.ncrpc[which(rownames(matrix.ncrpc)%in%unique(c(f4$gen1,f4$gene2,f5$gen1,f5$gene2))),], remove_empty_columns = FALSE, remove_empty_rows = T,
#             alter_fun = alter_fun, col = col,
#            show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 5),column_title="non-crpc(184(95F+86T))")

t01a=oncoPrint(input_matrix.path2.pelvic[-ncol(input_matrix.path2.pelvic)], remove_empty_columns = FALSE, remove_empty_rows = F,
               alter_fun = alter_fun, col = col,row_order =input_matrix.path2.pelvic$gene,row_split=input_matrix.path2.pelvic$pw, 
               #row_title_rot = switch(row_title_side[1], "left" = 0),
               
               row_title_gp = gpar(fontsize = 6) ,show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 5),column_title="Location_of_failure: Pelvic Nod (22)")
#t01s=oncoPrint(matrix.ncrpc[which(rownames(matrix.ncrpc)%in%unique(c(f4$gen1,f4$gene2,f5$gen1,f5$gene2))),], remove_empty_columns = FALSE, remove_empty_rows = T,
#             alter_fun = alter_fun, col = col,
#            show_co

t02=oncoPrint(input_matrix.path.node[-ncol(input_matrix.path.node)], remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,row_order =input_matrix.path.node$gene,row_split=input_matrix.path.node$pw, row_title_gp = gpar(fontsize = 6) ,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 5),column_title="Location_of_failure:Distant Node no  (137)")
#t01s=oncoPrint(matrix.ncrpc[which(rownames(matrix.ncrpc)%in%unique(c(f4$gen1,f4$gene2,f5$gen1,f5$gene2))),], remove_empty_columns = FALSE, remove_empty_rows = T,
t02a=oncoPrint(input_matrix.path2.node[-ncol(input_matrix.path2.node)], remove_empty_columns = FALSE, remove_empty_rows = F,
               alter_fun = alter_fun, col = col,row_order =input_matrix.path2.node$gene,row_split=input_matrix.path2.node$pw, row_title_gp = gpar(fontsize = 6) ,
               show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 5),column_title="Location_of_failure: Distance Node (22)")
#t01s=oncoPrint(matrix.ncrpc[which(rownames(matrix.ncrpc)%in%unique(c(f4$gen1,f4$gene2,f5$gen1,f5$gene2))),], remove_empty_columns = FALSE, remove_empty_rows = T,
#             alter_fun = alter_fun, col = col,
#            show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 5),column_title="non-crpc(184(95F+86T))")

t03=oncoPrint(input_matrix.path.bone[-ncol(input_matrix.path.bone)], remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,row_order =input_matrix.path.bone$gene,row_split=input_matrix.path.bone$pw, row_title_gp = gpar(fontsize = 6) ,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize =5),column_title="Location_of_failure:Bone no (129)")
#t01s=oncoPrint(matrix.ncrpc[which(rownames(matrix.ncrpc)%in%unique(c(f4$gen1,f4$gene2,f5$gen1,f5$gene2))),], remove_empty_columns = FALSE, remove_empty_rows = T,
#             alter_fun = alter_fun, col = col,
#            show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 5),column_title="non-crpc(184(95F+86T))")
t03a=oncoPrint(input_matrix.path2.bone[-ncol(input_matrix.path2.bone)], remove_empty_columns = FALSE, remove_empty_rows = F,
               alter_fun = alter_fun, col = col,row_order =input_matrix.path2.bone$gene,row_split=input_matrix.path2.bone$pw, row_title_gp = gpar(fontsize = 6) ,
               show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize =5),column_title="Location_of_failure:Bone yes (50)")

#t01s=oncoPrint(matrix.ncrpc[which(rownames(matrix.ncrpc)%in%unique(c(f4$gen1,f4$gene2,f5$gen1,f5$gene2))),], remove_empty_columns = FALSE, remove_empty_rows = T,
#             alter_fun = alter_fun, col = col,
#            show_column_n


t04=oncoPrint(input_matrix.path.vis[-ncol(input_matrix.path.vis)], remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,row_order =input_matrix.path.vis$gene,row_split=input_matrix.path.vis$pw, row_title_gp = gpar(fontsize = 6) ,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 5),column_title="Location_of_failure:Visceral no (163)")
#t01s=oncoPrint(matrix.ncrpc[which(rownames(matrix.ncrpc)%in%unique(c(f4$gen1,f4$gene2,f5$gen1,f5$gene2))),], remove_empty_columns = FALSE, remove_empty_rows = T,
#             alter_fun = alter_fun, col = col,
#            show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 5),column_title="non-crpc(184(95F+86T))")
t04a=oncoPrint(input_matrix.path2.vis[-ncol(input_matrix.path2.vis)], remove_empty_columns = FALSE, remove_empty_rows = F,
               alter_fun = alter_fun, col = col,row_order =input_matrix.path2.vis$gene,row_split=input_matrix.path2.vis$pw, row_title_gp = gpar(fontsize = 6) ,
               show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 5),column_title="Location_of_failure:Visceral yes (14)")
#t01s=oncoPrint(matrix.ncrpc[which(rownames(matrix.ncrpc)%in%unique(c(f4$gen1,f4$gene2,f5$gen1,f5$gene2))),], remove_empty_columns = FALSE, remove_empty_rows = T,
#             alter_fun = alter_fun, col = col,
#            show_column_names = T


t01+t01a
t02+t02a
t03+t03a
t04+t04a
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


pathway010.ann=path.agg[which(colnames(path.agg)%in%colnames(input_matrix.path.pelvic))]
pathway010.ann1=path.agg[which(colnames(path.agg)%in%colnames(input_matrix.path2.pelvic))]

pathway020.ann=path.agg[which(colnames(path.agg)%in%colnames(input_matrix.path.node))]
pathway020.ann1=path.agg[which(colnames(path.agg)%in%colnames(input_matrix.path2.node))]

pathway030.ann=path.agg[which(colnames(path.agg)%in%colnames(input_matrix.path.bone))]
pathway030.ann1=path.agg[which(colnames(path.agg)%in%colnames(input_matrix.path2.bone))]

pathway040.ann=path.agg[which(colnames(path.agg)%in%colnames(input_matrix.path.vis))]
pathway040.ann1=path.agg[which(colnames(path.agg)%in%colnames(input_matrix.path2.vis))]




pathway010.ann[pathway010.ann == 0] <- ""
pathway010.ann1[pathway010.ann1 == 0] <- ""

pathway020.ann[pathway020.ann == 0] <- ""

pathway020.ann1[pathway020.ann1 == 0] <- ""



pathway030.ann[pathway030.ann == 0] <- ""
pathway030.ann1[pathway030.ann1 == 0] <- ""

pathway040.ann[pathway040.ann == 0] <- ""
pathway040.ann1[pathway040.ann1 == 0] <- ""

pathway010.ann.1=pathway010.ann %>% mutate_if(is.numeric, ~1 * (. != 0))
pathway020.ann.1=pathway020.ann %>% mutate_if(is.numeric, ~1 * (. != 0))
pathway030.ann.1=pathway030.ann %>% mutate_if(is.numeric, ~1 * (. != 0))
pathway040.ann.1=pathway040.ann %>% mutate_if(is.numeric, ~1 * (. != 0))


pathway010.ann1.1=pathway010.ann1 %>% mutate_if(is.numeric, ~1 * (. != 0))
pathway020.ann1.1=pathway020.ann1 %>% mutate_if(is.numeric, ~1 * (. != 0))
pathway030.ann1.1=pathway030.ann1 %>% mutate_if(is.numeric, ~1 * (. != 0))
pathway040.ann1.1=pathway040.ann1 %>% mutate_if(is.numeric, ~1 * (. != 0))


pdf("pathway.location.failure.agg..v2.pdf",15,5)
t01=oncoPrint(pathway010.ann.1, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun0, col = col0,show_heatmap_legend =F,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 5),column_title="Location_of_failure:Pelvic Node no (137)")
t01a=oncoPrint(pathway010.ann1.1, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun0, col = col0,show_heatmap_legend =F,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 5),column_title="Location_of_failure:Pelvic Node (40)")


t02=oncoPrint(pathway020.ann.1, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun0, col = col0,show_heatmap_legend =F,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 5),column_title="Location_of_failure:Distant Node no  (137)")
t02a=oncoPrint(pathway020.ann1.1, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun0, col = col0,show_heatmap_legend =F,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 5),column_title="Location_of_failure:Distant Node yes (40)")

t03=oncoPrint(pathway030.ann.1, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun0, col = col0,show_heatmap_legend =F,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 5),column_title="Location_of_failure:Bone no (129)")

t03a=oncoPrint(pathway030.ann1.1, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun0, col = col0,show_heatmap_legend =F,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 5),column_title="Location_of_failure:Bone yes (50)")



t04=oncoPrint(pathway040.ann.1, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun0, col = col0,show_heatmap_legend =F,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 5),column_title="Location_of_failure:Visceral no (163)")

t04a=oncoPrint(pathway040.ann1.1, remove_empty_columns = FALSE, remove_empty_rows = F,
               alter_fun = alter_fun0, col = col0,show_heatmap_legend =F,
               show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 5),column_title="Location_of_failure:Visceral yes (14)")



t01+t01a
t02+t02a
t03+t03a
t04+t04a

dev.off()

