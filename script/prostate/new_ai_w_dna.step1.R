setwd("/local/projects-t3/PTRAN/Projects_starting_Jan2025/dnaseq/ysong/IGS_pipeline_output/maf/new1/")
path=getwd()
file.names <- dir(path, pattern =".maf")

#meta.soma=read.table("meta.soma.txt",sep="\t",header=T)
#maf_3=maftools:::merge_mafs(file.names,
 #                           vc_nonSyn=c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation"))



list.mafs <- list.files(path, pattern = "\\.maf$", full.names = TRUE)

list.all.maf.files <- lapply(list.mafs, function(file_path) {
  # Extract sample name from filename
  sample_name <- tools::file_path_sans_ext(basename(file_path))
  
  # Read MAF file
  maf_df <- read.delim(file_path, sep = "\t", header = TRUE, fill = TRUE, comment.char = "#")
  
  # Replace Tumor_Sample_Barcode with sample_name
  if ("Tumor_Sample_Barcode" %in% colnames(maf_df)) {
    maf_df$Tumor_Sample_Barcode <- sample_name
  }
  
  return(maf_df)
})




var.classes=c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation")

var.classes1 <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation","Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation","3'Flank", "3'UTR", "5'Flank",
                "5'UTR","Intron","RNA","Silent","Splice_Region", "DE_NOVO_START_IN_FRAME", "DE_NOVO_START_OUT_FRAME", "START_CODON_SNP")


merged_mafs_TN_1000G <- maftools::merge_mafs(list.all.maf.files, vc_nonSyn=var.classes)


merged_mafs_all <- maftools::merge_mafs(list.all.maf.files, vc_nonSyn=var.classes1)

maf1.nonsyn=merged_mafs_TN_1000G
maf1.syn=merged_mafs_all
___________________________________________________________________________

setwd("/local/projects-t3/PTRAN/Projects_starting_Jan2025/dnaseq/ysong/IGS_pipeline_output/maf/new0/")
path=getwd()
file.names <- dir(path, pattern =".maf")

#meta.soma=read.table("meta.soma.txt",sep="\t",header=T)
#maf_3=maftools:::merge_mafs(file.names,
 #                           vc_nonSyn=c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation"))



list.mafs <- list.files(path, pattern = "\\.maf$", full.names = TRUE)

list.all.maf.files <- lapply(list.mafs, function(file_path) {
  # Extract sample name from filename
  sample_name <- tools::file_path_sans_ext(basename(file_path))
  
  # Read MAF file
  maf_df <- read.delim(file_path, sep = "\t", header = TRUE, fill = TRUE, comment.char = "#")
  
  # Replace Tumor_Sample_Barcode with sample_name
  if ("Tumor_Sample_Barcode" %in% colnames(maf_df)) {
    maf_df$Tumor_Sample_Barcode <- sample_name
  }
  
  return(maf_df)
})




var.classes=c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation")

var.classes1 <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation","Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation","3'Flank", "3'UTR", "5'Flank",
                "5'UTR","Intron","RNA","Silent","Splice_Region", "DE_NOVO_START_IN_FRAME", "DE_NOVO_START_OUT_FRAME", "START_CODON_SNP")


merged_mafs_TN_1000G <- maftools::merge_mafs(list.all.maf.files, vc_nonSyn=var.classes)


merged_mafs_all <- maftools::merge_mafs(list.all.maf.files, vc_nonSyn=var.classes1)

maf2.nonsyn=merged_mafs_TN_1000G
maf2.syn=merged_mafs_all



________________________________________________________________________________________________


setwd("/local/projects-t3/PTRAN/Projects_starting_Jan2025/dnaseq/ysong/IGS_pipeline_output/maf/new2/")
path=getwd()
file.names <- dir(path, pattern =".maf")

#meta.soma=read.table("meta.soma.txt",sep="\t",header=T)
#maf_3=maftools:::merge_mafs(file.names,
 #                           vc_nonSyn=c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation"))



list.mafs <- list.files(path, pattern = "\\.maf$", full.names = TRUE)

list.all.maf.files <- lapply(list.mafs, function(file_path) {
  # Extract sample name from filename
  sample_name <- tools::file_path_sans_ext(basename(file_path))
  
  # Read MAF file
  maf_df <- read.delim(file_path, sep = "\t", header = TRUE, fill = TRUE, comment.char = "#")
  
  # Replace Tumor_Sample_Barcode with sample_name
  if ("Tumor_Sample_Barcode" %in% colnames(maf_df)) {
    maf_df$Tumor_Sample_Barcode <- sample_name
  }
  
  return(maf_df)
})




var.classes=c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation")

var.classes1 <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation","Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation","3'Flank", "3'UTR", "5'Flank",
                "5'UTR","Intron","RNA","Silent","Splice_Region", "DE_NOVO_START_IN_FRAME", "DE_NOVO_START_OUT_FRAME", "START_CODON_SNP")


merged_mafs_TN_1000G <- maftools::merge_mafs(list.all.maf.files, vc_nonSyn=var.classes)


merged_mafs_all <- maftools::merge_mafs(list.all.maf.files, vc_nonSyn=var.classes1)

maf3.nonsyn=merged_mafs_TN_1000G
maf3.syn=merged_mafs_all


________________________________________________________________


library(maftools)

maf.all.nonsy <- maftools::merge_mafs(maf = list(maf1.nonsyn,maf2.nonsyn,maf3.nonsyn))
saveRDS(maf.all.nonsy,"/Volumes/projects-t3/PTRAN/Projects_starting_Jan2025/dnaseq/ysong/IGS_pipeline_output/maf/all.nonsyn.rds")

maf.all.sy=readRDS("/local/projects-t3/PTRAN/Projects_starting_Jan2025/dnaseq/ysong/IGS_pipeline_output/maf/all.include.syn.rds")


merged_mafs_TN_1000G_sub=maf.all.nonsy@data
merged_mafs_TN_1000G_sub$tot_depth_tumor <- merged_mafs_TN_1000G_sub$t_alt_count + merged_mafs_TN_1000G_sub$t_ref_count
merged_mafs_TN_1000G_sub$tot_depth_normal <- merged_mafs_TN_1000G_sub$n_alt_count + merged_mafs_TN_1000G_sub$n_ref_count
merged_mafs_TN_1000G_sub$DNA.Mutations <- paste(merged_mafs_TN_1000G_sub$Hugo_Symbol,merged_mafs_TN_1000G_sub$Protein_Change)


merged_mafs_TN_1000G_sub <- merged_mafs_TN_1000G_sub %>% filter(tumor_f < 0.40 & tot_depth_tumor>=10)
library(maftools)
sub1=read.maf(merged_mafs_TN_1000G_sub)
clin=unique( sub1@data$ClinVar_VCF_CLNSIG)
sub1f=subsetMaf(sub1, query=  "ClinVar_VCF_CLNSIG %in% clin[c(3,4,6,7,11,12,13)]",dropLevels =F)
#d <- merge_mafs(lapply(Sys.glob("mafs/Patient*.maf"), read.maf))
saveRDS(sub1f,"/Volumes/projects-t3/PTRAN/Projects_starting_Jan2025/dnaseq/ysong/IGS_pipeline_output/maf/all.patho.w.conflic.patho.rds")



tempus.gene=read.table("~/Desktop/P_Tran_ref/tempus/tempus_gene_update.version.txt",sep="\t",header=F)
foundation.gene=read.table("/Users/ysong/Desktop/P_Tran_ref/foundation/maf_validation/foundation_gene_to_plot.txt",,sep="\t",header=F)

sub1f.sub <- subsetMaf(maf = sub1f, genes = tempus.gene$V1,dropLevels=F)

saveRDS(sub1f.sub,"/Volumes/projects-t3/PTRAN/Projects_starting_Jan2025/dnaseq/ysong/IGS_pipeline_output/maf/all.patho.w.conflic.patho.update.rds")
pdf("/local/projects-t3/PTRAN/Projects_starting_Jan2025/dnaseq/ysong/IGS_pipeline_output/maf/plot.patho.pdf",10,20)
oncoplot(maf = sub1f.sub,genes = unique(sub1f.sub@gene.summary$Hugo_Symbol),writeMatrix = T)
dev.off()

mutation.new=read.delim("/Volumes/projects-t3/PTRAN/Projects_starting_Jan2025/dnaseq/ysong/IGS_pipeline_output/maf/new_onco_matrix.txt",sep="\t",header=T)


meta=read.delim("//Users/ysong/Desktop/P_Tran_ref/score/to_use/score_111/image_predict_score/Complete Artera Database_v5.mod.txt",sep="\t",header=T)
meta$Tumor_Sample_Barcode=gsub("ORD-1108071-01-T","ORD-1108071-01",meta$Tumor_Sample_Barcode)

meta.therapy=meta[order(meta$Tumor_Sample_Barcode),]    

mutation=read.table("/Users/ysong/Desktop/P_Tran_ref/score/to_use/score_Phil_update_v2/onco_matrix.txt",sep="\t",check.names=F,header=T,row.names=1)    
mutation.sub=mutation[which(colnames(mutation)%in%meta.therapy$Tumor_Sample_Barcode)]
mutation.sub=mutation.sub[order(colnames(mutation.sub))]

mutation.all=merge(mutation.sub,mutation.new,by=0,all=T)

write.table(mutation.all,"~/Desktop/all_mutation.combine.710.txt",sep="\t",quote=F)

meta=read.delim("/local/projects-t3/PTRAN/Projects_starting_Jan2024/ysong/clinical_from_phil/clinical_2024/Complete_Artera_Database_v5.mod.txt",sep="\t",header=T)
meta$Tumor_Sample_Barcode=gsub("ORD-1108071-01-T","ORD-1108071-01",meta$Tumor_Sample_Barcode)

meta.therapy=meta[order(meta$Tumor_Sample_Barcode),]    

