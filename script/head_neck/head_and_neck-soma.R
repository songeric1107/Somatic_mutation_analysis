###vcf2maf
#for i in *.vcf;do echo perl /usr/local/packages/vcf2maf/vcf2maf.pl --input-vcf $i --output-maf $i.maf --vep-path /usr/local/packages/miniconda3/envs/vep/bin --ref-fasta /local/db/vep/Homo_sapiens.GRCh37.dna.toplevel.fa.gz --ncbi-build GRCh37;done>step1_vcf_maf.pbs
#merge all maf together
#https://github.com/songeric1107/Somatic_mutation_analysis/blob/main/script/prostate/merge_maf.R


library(maftools)

soma=readRDS("/Volumes/projects-t3/PTRAN/Projects_starting_Jan2024/ysong/HN_analysis_2024/maf/maf/")

meta=read.delim("/Volumes/projects-t3/PTRAN/Projects_starting_Jan2024/ysong/HN_analysis_2024/maf/maf/hn_meta.update_ysong.txt",sep="\t",header=T)

sample=soma@clinical.data
sample.sub=sample[which(sample$Tumor_Sample_Barcode%in%meta$xT),]


meta=read.delim("/Volumes/projects-t3/PTRAN/Projects_starting_Jan2024/ysong/HN_analysis_2024/maf/maf/meta_onco.txt",sep="\t",header=T)


# Define a function that applies clinicalEnrichment to each column

pdf("summary.pdf")
plotmafSummary(soma)
dev.off()


##compare local vs progression

col = c("Missense_Mutation" = "blue", "In_Frame_Ins" = "red", "Splice_Site" = "#008000","In_Frame_Del"="purple","Frame_Shift_Ins"="yellow",
        "Frame_Shift_Del"="cyan", "Nonsense_Mutation"="brown", "Multi_Hit"="orange","Translation_Start_Site"="pink")





#prepare the topbar

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),
              gp = gpar(fill = "white", col = NA))
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

meta=read.delim("/Volumes/projects-t3/PTRAN/Projects_starting_Jan2024/ysong/HN_analysis_2024/maf/maf/meta_onco.txt",sep="\t",header=T)





meta$P16.Status..0.neg..1.pos..2.unknown.=as.factor(meta$P16.Status..0.neg..1.pos..2.unknown.)

################################
input_matrix=read.table("/Volumes/projects-t3/PTRAN/Projects_starting_Jan2024/ysong/HN_analysis_2024/maf/maf/onco_matrix.txt",sep="\t",header=T,row.names=1,check.names=F)
input_matrix[input_matrix==0]<-""

input_matrix=input_matrix[,order(match(colnames(input_matrix),meta$Tumor_Sample_Barcode))]


set.seed(10)
library(ComplexHeatmap)

# Function to generate color mappings automatically
column_colors <- lapply(df, function(col) {
  if(is.factor(col) || is.character(col)) {
    # For categorical columns
    factor_levels <- unique(col)
    colors <- rainbow(length(factor_levels))  # Assign different colors
    return(setNames(colors, factor_levels))
  } else {
    # For continuous data, NULL means default color scale
    return(NULL)
  }
})

# Create HeatmapAnnotation using the df argument
ha <- HeatmapAnnotation(df = meta[-c(1:3)], col = column_colors)


pdf("all.pdf",30,30)

p1=oncoPrint(input_matrix, remove_empty_columns = FALSE, remove_empty_rows = F,
          alter_fun = alter_fun, col = col,
          bottom_annotation=ha,

          #column_order = meta.lf$Tumor_Sample_Barcode,
          show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6))

draw(p1,          heatmap_legend_side = "bottom",      # Place heatmap legend at the bottom
     annotation_legend_side = "bottom",)

dev.off()
meta=meta[which(meta$P16.Status..0.neg..1.pos..2.unknown.!=2),]

meta.p16=meta[which(meta$P16.Status..0.neg..1.pos..2.unknown.=="1"),]
meta.p16.n=meta[which(meta$P16.Status..0.neg..1.pos..2.unknown.!="1"),]


p16=input_matrix[which(colnames(input_matrix)%in%meta.p16$Tumor_Sample_Barcode)]
p16.n=input_matrix[which(colnames(input_matrix)%in%meta.p16.n$Tumor_Sample_Barcode)]

p16.s=p16[,order(match(colnames(p16),meta.p16$Tumor_Sample_Barcode))]
p16n.s=p16.n[,order(match(colnames(p16.n),meta.p16.n$Tumor_Sample_Barcode))]

library(ComplexHeatmap)

set.seed(10)
ha11 = HeatmapAnnotation(p16 = meta.p16$P16.Status..0.neg..1.pos..2.unknown.,col=list(p16=c("1"="red","0"="grey","2"="grey")))
ha12 = HeatmapAnnotation(p16 = meta.p16.n$P16.Status..0.neg..1.pos..2.unknown.,col=list(p16=c("1"="red","0"="grey","2"="grey")))


pdf("~/Desktop/somatic_p16.pdf",20,10)
t01=oncoPrint(p16.s, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              bottom_annotation=ha11,
              #column_order = meta.lf$Tumor_Sample_Barcode,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="P16 (Y) (42)")


t02=oncoPrint(p16n.s, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              bottom_annotation = ha12,
             # column_order = meta.pr$Tumor_Sample_Barcode,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="P16 (N) (52)")

t01+t02

dev.off()





########upfront


meta=read.delim("/Volumes/projects-t3/PTRAN/Projects_starting_Jan2024/ysong/HN_analysis_2024/maf/maf/meta_onco.txt",sep="\t",header=T)


meta=meta[which(meta$Upfront.resecton...0.no..1.yes.!=2),]

meta.p16=meta[which(meta$Upfront.resecton...0.no..1.yes.=="1"),]
meta.p16.n=meta[which(meta$Upfront.resecton...0.no..1.yes.!="1"),]


p16=input_matrix[which(colnames(input_matrix)%in%meta.p16$Tumor_Sample_Barcode)]
p16.n=input_matrix[which(colnames(input_matrix)%in%meta.p16.n$Tumor_Sample_Barcode)]

p16.s=p16[,order(match(colnames(p16),meta.p16$Tumor_Sample_Barcode))]
p16n.s=p16.n[,order(match(colnames(p16.n),meta.p16.n$Tumor_Sample_Barcode))]

library(ComplexHeatmap)

set.seed(10)
ha11 = HeatmapAnnotation(p16 = meta.p16$Upfront.resecton...0.no..1.yes.,col=list(p16=c("1"="red","0"="grey","2"="grey")))
ha12 = HeatmapAnnotation(p16 = meta.p16.n$Upfront.resecton...0.no..1.yes.,col=list(p16=c("1"="red","0"="grey","2"="grey")))


pdf("~/Desktop/somatic_upfron.resction.pdf",20,10)
t01=oncoPrint(p16.s, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              bottom_annotation=ha11,
              #column_order = meta.lf$Tumor_Sample_Barcode,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="Upfront.resecton (Y) (50)")


t02=oncoPrint(p16n.s, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              bottom_annotation = ha12,
              # column_order = meta.pr$Tumor_Sample_Barcode,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="Upfront.resecton (N) (47)")

t01+t02

dev.off()



########local.advance


meta=read.delim("/Volumes/projects-t3/PTRAN/Projects_starting_Jan2024/ysong/HN_analysis_2024/maf/maf/meta_onco.txt",sep="\t",header=T)




meta.sub1=meta[which(meta$X1.locally.advanced..2.local.recurrence.3.de.novo.metastatic.definitive.4.recurrent.metastatic.definitive.intent=="1"),]
meta.sub2=meta[which(meta$X1.locally.advanced..2.local.recurrence.3.de.novo.metastatic.definitive.4.recurrent.metastatic.definitive.intent=="2"),]
meta.sub3=meta[which(meta$X1.locally.advanced..2.local.recurrence.3.de.novo.metastatic.definitive.4.recurrent.metastatic.definitive.intent%in%c("1 2","1 3","2/4-node","3","4")),]



m1=input_matrix[which(colnames(input_matrix)%in%meta.sub1$Tumor_Sample_Barcode)]
m2=input_matrix[which(colnames(input_matrix)%in%meta.sub2$Tumor_Sample_Barcode)]
m3=input_matrix[which(colnames(input_matrix)%in%meta.sub3$Tumor_Sample_Barcode)]



m1.s=m1[,order(match(colnames(m1),m1$Tumor_Sample_Barcode))]
m2.s=m2[,order(match(colnames(m2),m2$Tumor_Sample_Barcode))]
m3.s=m3[,order(match(colnames(m3),m3$Tumor_Sample_Barcode))]

library(ComplexHeatmap)

set.seed(10)
ha11 = HeatmapAnnotation(local_advance = meta.sub1$X1.locally.advanced..2.local.recurrence.3.de.novo.metastatic.definitive.4.recurrent.metastatic.definitive.intent,col=list(local_advance=c("1"="red","0"="grey","2"="green")),
                         annotation_legend_param = list(local_advance=list(nrow=1) ))
  #nmf_cluster = list(nrow=1),



ha12 = HeatmapAnnotation(local_advance = meta.sub2$X1.locally.advanced..2.local.recurrence.3.de.novo.metastatic.definitive.4.recurrent.metastatic.definitive.intent,col=list(local_advance=c("1"="red","0"="grey","2"="green")),
                         annotation_legend_param = list(local_advance=list(nrow=1) ))
ha13 = HeatmapAnnotation(local_advance = meta.sub3$X1.locally.advanced..2.local.recurrence.3.de.novo.metastatic.definitive.4.recurrent.metastatic.definitive.intent, annotation_legend_param = list(local_advance=list(nrow=1) ))

pdf("~/Desktop/somatic_local_advance.pdf",20,15)
t01=oncoPrint(m1.s, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              bottom_annotation=ha11,
              #column_order = meta.lf$Tumor_Sample_Barcode,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="local advance (Y) (80)")


t02=oncoPrint(m2.s, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              bottom_annotation = ha12,
              # column_order = meta.pr$Tumor_Sample_Barcode,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="local.recurrence (Y) (12)")


t03=oncoPrint(m3.s, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              bottom_annotation = ha13,
              # column_order = meta.pr$Tumor_Sample_Barcode,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="others (5)")


draw(t01+t02+t03,          heatmap_legend_side = "bottom",      # Place heatmap legend at the bottom
     annotation_legend_side = "bottom",)

dev.off()






####T..0.Tis.1.T1.2.T2.3.T3.4.T4.5.T0.6.Unk...7.recurrent.



meta=soma@clinical.data


meta=meta[order(meta$T..0.Tis.1.T1.2.T2.3.T3.4.T4.5.T0.6.Unk...7.recurrent.),]
input_matrix=input_matrix[,order(match(colnames(input_matrix),meta$Tumor_Sample_Barcode))]
# Split the data frame based on the 'Group' factor
split_df <- split(meta, meta$T..0.Tis.1.T1.2.T2.3.T3.4.T4.5.T0.6.Unk...7.recurrent.)

# Check the split data frames





m1=input_matrix[which(colnames(input_matrix)%in%split_df[[1]]$Tumor_Sample_Barcode)]
m2=input_matrix[which(colnames(input_matrix)%in%split_df[[2]]$Tumor_Sample_Barcode)]
m3=input_matrix[which(colnames(input_matrix)%in%split_df[[3]]$Tumor_Sample_Barcode)]
m4=input_matrix[which(colnames(input_matrix)%in%split_df[[4]]$Tumor_Sample_Barcode)]
m5=input_matrix[which(colnames(input_matrix)%in%split_df[[5]]$Tumor_Sample_Barcode)]
m6=input_matrix[which(colnames(input_matrix)%in%c(split_df[[6]]$Tumor_Sample_Barcode,split_df[[7]]$Tumor_Sample_Barcode,split_df[[8]]$Tumor_Sample_Barcode))]





library(ComplexHeatmap)

set.seed(10)
ha11 = HeatmapAnnotation(local_advance = split_df$`0`$T..0.Tis.1.T1.2.T2.3.T3.4.T4.5.T0.6.Unk...7.recurrent.,
                         annotation_legend_param = list(local_advance=list(nrow=1) ))
#nmf_cluster = list(nrow=1),



ha12 = HeatmapAnnotation(local_advance = split_df$`1`$T..0.Tis.1.T1.2.T2.3.T3.4.T4.5.T0.6.Unk...7.recurrent.,
                         annotation_legend_param = list(local_advance=list(nrow=1) ))
ha13 = HeatmapAnnotation(local_advance = split_df$`2`$T..0.Tis.1.T1.2.T2.3.T3.4.T4.5.T0.6.Unk...7.recurrent.,
                         annotation_legend_param = list(local_advance=list(nrow=1) ))

ha14 = HeatmapAnnotation(local_advance = split_df$`3`$T..0.Tis.1.T1.2.T2.3.T3.4.T4.5.T0.6.Unk...7.recurrent.,
                         annotation_legend_param = list(local_advance=list(nrow=1) ))


ha15 = HeatmapAnnotation(local_advance = split_df$`4`$T..0.Tis.1.T1.2.T2.3.T3.4.T4.5.T0.6.Unk...7.recurrent.,
                         annotation_legend_param = list(local_advance=list(nrow=1) ))


meta.other=rbind(split_df[[6]],split_df[[7]],split_df[[8]])
ha16 = HeatmapAnnotation(local_advance = meta.other$T..0.Tis.1.T1.2.T2.3.T3.4.T4.5.T0.6.Unk...7.recurrent.,
                         annotation_legend_param = list(local_advance=list(nrow=1) ))



pdf("~/Desktop/somatic_T.pdf",50,15)
t01=oncoPrint(m1, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              bottom_annotation=ha11,
              #column_order = meta.lf$Tumor_Sample_Barcode,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="T0 (4)")


t02=oncoPrint(m2, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              bottom_annotation = ha12,
              # column_order = meta.pr$Tumor_Sample_Barcode,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="T1 (14)")


t03=oncoPrint(m3, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              bottom_annotation = ha13,
              # column_order = meta.pr$Tumor_Sample_Barcode,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="T2 (29)")

t04=oncoPrint(m4, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              bottom_annotation=ha14,
              #column_order = meta.lf$Tumor_Sample_Barcode,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="T3 (26)")


t05=oncoPrint(m5, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              bottom_annotation = ha15,
              # column_order = meta.pr$Tumor_Sample_Barcode,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="T4 (23)")


t06=oncoPrint(m6, remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              bottom_annotation = ha16,
              # column_order = meta.pr$Tumor_Sample_Barcode,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="T>4 (12)")


draw(t01+t02+t03+t04+t05+t06,          heatmap_legend_side = "bottom",      # Place heatmap legend at the bottom
     annotation_legend_side = "bottom",)

dev.off()

##############################################################










#Y <- sapply(lf.s, function(x) sum(is.na(x) | x == ""))

pdf("~/Desktop/somatic_PFS_label.pdf",20,10)
t01=oncoPrint(lf.s[which(rownames(lf.s)%in%unique(c(comp1f$Hugo_Symbol,comp2f$Hugo_Symbol))),], remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              bottom_annotation=ha11,
              column_order = meta.lf$Tumor_Sample_Barcode,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="PFS_0_alive_noProg (8)")


t02=oncoPrint(pr.s[which(rownames(pr.s)%in%unique(c(comp1f$Hugo_Symbol,comp2f$Hugo_Symbol))),], remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              bottom_annotation = ha12,
              column_order = meta.pr$Tumor_Sample_Barcode,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="PFS_1_death_wProg (28)")

t01+t02

dev.off()

ha21 = HeatmapAnnotation(PFS.Time = meta.lf$PFS.Time.to.any.progression)
ha22 = HeatmapAnnotation(PFS.Time=meta.pr$PFS.Time.to.any.progression)

pdf("~/Desktop/somatic_PFS_label2.pdf",20,10)
t01=oncoPrint(lf.s[which(rownames(lf.s)%in%unique(c(comp1f$Hugo_Symbol,comp2f$Hugo_Symbol))),], remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              bottom_annotation=ha21,
              # column_order = meta.lf$Tumor_Sample_Barcode,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="PFS_0_alive_noProg (8)")


t02=oncoPrint(pr.s[which(rownames(pr.s)%in%unique(c(comp1f$Hugo_Symbol,comp2f$Hugo_Symbol))),], remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              bottom_annotation = ha22,
              #  column_order = meta.pr$Tumor_Sample_Barcode,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="PFS_1_death_wProg (28)")

t01+t02

dev.off()


meta.lf=  meta[which(meta$PFS.Event ==0),]
meta.lf=meta.lf[order(meta.lf$OS.Time.to.death.or.last.FU),]

meta.pr=  meta[which(meta$PFS.Event  ==1),]
meta.pr=meta.pr[order(meta.pr$OS.Time.to.death.or.last.FU),]

################################
input_matrix=read.table("soma.f.onco_matrix.txt",sep="\t",header=T,row.names=1,check.names=F)


input_matrix=input_matrix[,order(match(colnames(input_matrix),meta1$Tumor_Sample_Barcode))]


lf=input_matrix[which(colnames(input_matrix)%in%meta.lf$Tumor_Sample_Barcode)]
pr=input_matrix[which(colnames(input_matrix)%in%meta.pr$Tumor_Sample_Barcode)]


lf.s=lf[,order(match(colnames(lf),meta.lf$Tumor_Sample_Barcode))]
pr.s=pr[,order(match(colnames(pr),meta.pr$Tumor_Sample_Barcode))]

ha11 = HeatmapAnnotation(OS.Time = anno_points(meta.lf$OS.Time.to.death.or.last.FU))
ha12 = HeatmapAnnotation(OS.Time=anno_points(meta.pr$OS.Time.to.death.or.last.FU))

#Y <- sapply(lf.s, function(x) sum(is.na(x) | x == ""))

pdf("~/Desktop/somatic_PFS_label.os.pdf",20,10)
t01=oncoPrint(lf.s[which(rownames(lf.s)%in%unique(c(comp1f$Hugo_Symbol,comp2f$Hugo_Symbol))),], remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              bottom_annotation=ha11,
              column_order = meta.lf$Tumor_Sample_Barcode,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="PFS_0_alive_noProg (8)")


t02=oncoPrint(pr.s[which(rownames(pr.s)%in%unique(c(comp1f$Hugo_Symbol,comp2f$Hugo_Symbol))),], remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,
              bottom_annotation = ha12,
              column_order = meta.pr$Tumor_Sample_Barcode,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 6),column_title="PFS_1_death_wProg (28)")

t01+t02

dev.off()



file1=lf.s[which(rownames(lf.s)%in%unique(c(comp1f$Hugo_Symbol,comp2f$Hugo_Symbol))),]
file2=pr.s[which(rownames(pr.s)%in%unique(c(comp1f$Hugo_Symbol,comp2f$Hugo_Symbol))),]


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

input_matrix.path.pr=input_matrix.path.s[,which(colnames(input_matrix.path.s)%in%meta.pr$Tumor_Sample_Barcode)]
input_matrix.path.pr$pw=input_matrix.path.s$pw

input_matrix.path2.lf=input_matrix.path.s[,which(colnames(input_matrix.path.s)%in%meta.lf$Tumor_Sample_Barcode)]
input_matrix.path2.lf$pw=input_matrix.path.s$pw


input_matrix.path.prs=input_matrix.path.pr[,order(match(colnames(input_matrix.path.pr),meta.pr$Tumor_Sample_Barcode)),]
input_matrix.path2.lfs=input_matrix.path2.lf[order(match(colnames(input_matrix.path2.lf),meta.lf$Tumor_Sample_Barcode))]



ha1=HeatmapAnnotation(OS=meta.lf$OS.Time.to.death.or.last.FU)

ha2=HeatmapAnnotation(OS=meta.pr$OS.Time.to.death.or.last.FU)

pdf("~/Desktop/path.PFS.pdf",20,20)

t01=oncoPrint(input_matrix.path2.lfs[,-ncol(input_matrix.path2.lfs)], remove_empty_columns = FALSE, remove_empty_rows = F,
              alter_fun = alter_fun, col = col,row_split=input_matrix.path2.lf$pw,bottom_annotation = ha1,
              #row_title_rot = switch(row_title_side[1], "left" = 0),

              row_title_gp = gpar(fontsize = 10) ,show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),
              column_names_gp = gpar(fontsize = 5),column_title="PFS_0_alive_noProgression (8)")
#t01s=oncoPrint(matrix.ncrpc[which(rownames(matrix.ncrpc)%in%unique(c(f4$gen1,f4$gene2,f5$gen1,f5$gene2))),], remove_empty_columns = FALSE, remove_empty_rows = T,
#             alter_fun = alter_fun, col = col,
#            show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 5),column_title="non-crpc(184(95F+86T))")

t01a=oncoPrint(input_matrix.path.prs[,-ncol(input_matrix.path.prs)], remove_empty_columns = FALSE, remove_empty_rows = F,
               alter_fun = alter_fun, col = col,row_split=input_matrix.path.pr$pw ,bottom_annotation = ha2,
               #row_title_rot = switch(row_title_side[1], "left" = 0),

               row_title_gp = gpar(fontsize = 10) ,show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),
               column_names_gp = gpar(fontsize = 5),column_title="PFS_1_death_w_progression (28)")
#t01s=oncoPrint(matrix.ncrpc[which(rownames(matrix.ncrpc)%in%unique(c(f4$gen1,f4$gene2,f5$gen1,f5$gene2))),], remove_empty_columns = FALSE, remove_empty_rows = T,
t01+t01a

dev.off()


#input_matrix.path1=soma_matrix.f

soma_matrix=data.frame(mutCountMatrix(maf = soma,removeNonMutated = F))

soma_matrix.f=soma_matrix %>% mutate_if(is.numeric, ~1 * (. > 0))





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




pdf("~/Desktop/pathway.PFS_label.agg..v3.pdf",15,5)
t01=oncoPrint(pathway010.ann1, remove_empty_columns = FALSE, remove_empty_rows = F,top_annotation = NULL,
              alter_fun = alter_fun0, col = col0,show_heatmap_legend =F,
              show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 5),column_title="PFS_0_alive_noProg (8)")
t01a=oncoPrint(pathway010.ann, remove_empty_columns = FALSE, remove_empty_rows = F,top_annotation = NULL,
               alter_fun = alter_fun0, col = col0,show_heatmap_legend =F,
               show_column_names = TRUE,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 5),column_title="PFS_1_death_wProg (28)")


t01+t01a


dev.off()






#####################################

meta$Tissue.Label..1.upfront..2.recurrence.=as.character(meta$Tissue.Label..1.upfront..2.recurrence.)
meta$Gender=as.factor(meta$Gender)

pdf("age_PFS_lable.pdf")
ggplot(meta,aes(x=Tissue.Label..1.upfront..2.recurrence.,y=Age))+geom_boxplot()+geom_point(aes(fill=factor(Gender)))

ggplot(meta,aes(x=Tissue.Label..1.upfront..2.recurrence.,y=Age,fill=Gender))+geom_boxplot()+facet_wrap(~Gender)+geom_point()

dev.off()
meta$Recurrence.Location=as.character(meta$Recurrence.Location)


ggplot(meta,aes(x=Gender,y=Age))+geom_boxplot()+geom_point()


pdf("age_PFS")
ggplot(meta,aes(x=meta$Tissue.Label..1.upfront..2.recurrence.,y=Age))+geom_boxplot()+geom_point(aes(fill=factor(Gender)))

ggplot(meta,aes(x=meta$Tissue.Label..1.upfront..2.recurrence.,y=Gender))+geom_boxplot()+geom_point()

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


soma_matrix.surv=soma_matrix.f[which(rownames(soma_matrix.f)%in%unique(c(comp1f$Hugo_Symbol,comp2f$Hugo_Symbol))),]
sum1=data.frame(colSums(soma_matrix.surv))
names(sum1)="total_mutation"

#sum1=sum1 %>% mutate_if(is.numeric, ~1 * (. > 0))



tp53.ann0=merge(meta,sum1,by.x="Tumor_Sample_Barcode",by.y=0)
tp53.ann=merge(meta,t(soma_matrix.surv),by.x="Tumor_Sample_Barcode",by.y=0)


tp53.ann0$PFS.Event=as.factor(tp53.ann0$PFS.Event)
pdf("mutation.pdf")
ggplot(tp53.ann0,aes(x=PFS.Event,y=total_mutation))+geom_violin()+geom_jitter(height = 0, width = 0.1)
dev.off()


pw_col = RColorBrewer::brewer.pal(n = 12, name = 'Set3')[1:12]

sfita <- survfit(as.formula(paste0("Surv(PFS.Time.to.any.progression,PFS.Event) ~", "total_mutation")), data=tp53.ann0[which(tp53.ann0$PFS.Event==1),])




pdf("~/Desktop/total.surv.pfs.pdf",20,15)
p1=   ggsurvplot(sfit0, conf.int=F, pval=TRUE, risk.table=TRUE,
                 # legend.labs=c(paste("total_mutation","_WO_mutation",sep="_"), paste("total_mutation","_W_mutation",sep="_")), legend.title="Mutation",
                 palette=pw_col,
                 title = "Status", subtitle = "Based on Kaplan-Meier estimates",

                 font.title = c(10, "bold", "darkblue"),
                 font.subtitle = c(10, "bold.italic", "purple"),
                 font.caption = c(10, "plain", "orange"),
                 risk.table.height=0.2,ggtheme = theme_minimal())


print(p1)
dev.off()


tp53.ann0$level <- rep(NA, nrow(tp53.ann0))

tp53.ann0[tp53.ann0$total_mutation<=2, ][, "level"] <- "less mutation"
tp53.ann0[tp53.ann0$total_mutation>2, ][, "level"] <- "more mutation"


pw_col = RColorBrewer::brewer.pal(n = 12, name = 'Paired')[1:12]

sfit0 <- survfit(as.formula(paste0("Surv(PFS.Time.to.any.progression,PFS.Event) ~", "total_mutation")), data=tp53.ann0)


sfit <- survfit(as.formula(paste0("Surv(PFS.Time.to.any.progression,PFS.Event) ~", "level")), data=tp53.ann0)


pdf("~/Desktop/total.surv.pfs.pdf",20,15)
p1=   ggsurvplot(sfit0, conf.int=F, pval=TRUE, risk.table=TRUE,
                # legend.labs=c(paste("total_mutation","_WO_mutation",sep="_"), paste("total_mutation","_W_mutation",sep="_")), legend.title="Mutation",
                 palette=pw_col,
                 title = "Status", subtitle = "Based on Kaplan-Meier estimates",

                 font.title = c(10, "bold", "darkblue"),
                 font.subtitle = c(10, "bold.italic", "purple"),
                 font.caption = c(10, "plain", "orange"),
                 risk.table.height=0.2,ggtheme = theme_minimal())


print(p1)

p1=   ggsurvplot(sfit, conf.int=F, pval=TRUE, risk.table=TRUE,
                 # legend.labs=c(paste("total_mutation","_WO_mutation",sep="_"), paste("total_mutation","_W_mutation",sep="_")), legend.title="Mutation",
                 palette=c("red","blue"),
                 title = "Status", subtitle = "Based on Kaplan-Meier estimates",

                 font.title = c(10, "bold", "darkblue"),
                 font.subtitle = c(10, "bold.italic", "purple"),
                 font.caption = c(10, "plain", "orange"),
                 risk.table.height=0.2,ggtheme = theme_minimal())


print(p1)


dev.off()

pdf("path.pdf")
path1=OncogenicPathways(maf =sub2)
path0=OncogenicPathways(maf = sub1)
dev.off()

path1=data.frame(path1)
path0=data.frame(path0)
names(path1)[2:6]=paste(names(path1)[2:6],"PFS_0",sep=".")
names(path0)[2:6]=paste(names(path0)[2:6],"PFS_1",sep=".")


sub1.2=merge(data.frame(path1),data.frame(path0),by="Pathway")
sub1.2$fraction_PFS_0_OTHER=1-sub1.2$fraction_affected.PFS_0
sub1.2$fraction_PFS_1_OTHER=1-sub1.2$fraction_affected.PFS_1

sub1.2$ratio=paste(sub1.2$n_affected_genes.PFS_0,sub1.2$N.PFS_0,sep="/")
sub1.2$ratio1=paste(sub1.2$n_affected_genes.PFS_1,sub1.2$N.PFS_1,sep="/")
sub1.2$id1=paste(sub1.2$Pathway,"(",sub1.2$ratio,")",sep="")
sub1.2$id2=paste(sub1.2$Pathway,"(",sub1.2$ratio1,")",sep="")


sub1.2.m=melt(sub1.2[c(17,18,4,9,12:13)])

sub1.2.m$variable=factor(sub1.2.m$variable,levels=c("fraction_affected.PFS_0","fraction_PFS_0_OTHER","fraction_affected.PFS_1","fraction_PFS_1_OTHER"))
a1=sub1.2.m[grep("PFS_0",sub1.2.m$variable),]
a2=sub1.2.m[grep("PFS_1",sub1.2.m$variable),]

a1$variable <- factor(a1$variable, levels=c("fraction_PFS_0_OTHER" , "fraction_affected.PFS_0"))

a2$variable <- factor(a2$variable, levels=c("fraction_PFS_1_OTHER" , "fraction_affected.PFS_1"))



pdf("fraction.path.pfs.pdf")
p1=ggplot(a1[-2],aes(x=id1,y=value,fill= variable))+geom_bar( position = "stack", stat="identity",show.legend = FALSE)+
  coord_flip()+scale_fill_manual(values=c("grey","red"))+xlab("pathway")+ylab("")

p2=ggplot(a2[-1],aes(x=id2,y=value,fill=variable))+
  geom_bar( position = "stack", stat="identity",show.legend = FALSE)+
  coord_flip()+scale_fill_manual(values=c("grey","red"))+xlab("pathway")+ylab("")
p1+p2


dev.off()

sig = oncodrive(maf = soma, AACol = pchange, bgEstimat=F,minMut = 2, pvalMethod = 'zscore')

pdf("~/Desktop/all_drive_mutation.pdf")
plotOncodrive(res = sig, fdrCutOff = 0.05, useFraction = TRUE, labelSize = 0.5,bubbleSize=1)

laml.sig=sig
laml.sig$level <- rep(NA, nrow(laml.sig))
laml.sig1=data.frame(laml.sig)
laml.sig1[laml.sig1$fdr<0.05, ][, "level"] <- "sig"
laml.sig1[laml.sig1$fdr>0.05, ][, "level"] <- "nonsig"

ggplot(data = laml.sig1, aes(x = muts_in_clusters, y = -log10(fdr), label = Hugo_Symbol,color=level))+
geom_point()+geom_hline(yintercept=1.3, linetype="dashed",      color = "red", size=1)+
ggrepel::geom_text_repel(max.overlaps=100)

dev.off()


library(ActiveDriverWGS)


input=soma@data[,c(6:8,12,13,1)]

input$Source_MAF=gsub(".soma.freebayes.vcf.hg38-filtered.vep","",input$Source_MAF)

this_genome = BSgenome.Hsapiens.UCSC.hg38::Hsapiens
library(GenomicRanges)

# format_muts


gene=soma@data[,c(6:8,2)]
colnames(gene)=c("chr",   "start",     "end",       "id")
results = ActiveDriverWGS(mutations = input,elements=unique(gene))


fit <- survfit(as.formula(paste0("Surv(PFS.Time.to.any.progression,PFS.Event) ~", "TP53")),  data = tp53.ann)

pdf("sur.pfs.621.pdf",15,10)
library(survminer)

genes=rownames(soma_matrix.surv)


for(i in 1:length(soma_matrix.surv)){



  fit <- survfit(as.formula(paste0("Surv(PFS.Time.to.any.progression,PFS.Event)  ~", as.character(genes[i]))),
                 data = tp53.ann[which(tp53.ann$PFS.Event==1),])

  pdf(paste0(genes[i], "_Survival_w_vs_wo_mutation.pfs1.pdf"))
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


sfita <- survfit(as.formula(paste0("Surv(PFS.Time.to.any.progression,PFS.Event) ~", "total_mutation")), data=tp53.ann0[which(tp53.ann0$PFS.Event==1),])


pdf("sur.pfs.621.pdf",15,10)
library(survminer)

genes=rownames(soma_matrix.surv)


for(i in 1:length(soma_matrix.surv)){



  fit <- survfit(as.formula(paste0("Surv(PFS.Time.to.any.progression,PFS.Event)  ~", genes[i])),
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

