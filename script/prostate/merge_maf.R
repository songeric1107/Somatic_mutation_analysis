library(maftools)
#setwd("/Users/ysong/Desktop/P_Tran_ref/maf/soma/")
path=getwd()
clinical_tmb=read.table("/Users/ysong/Desktop/P_Tran_ref/onco_508/soma/meta_no97.txt",sep="\t",header=T,check.names = F)
file.names <- dir(path, pattern =".maf.mod.txt")


maf_3=maftools:::merge_mafs(file.names,clinicalData=clinical_tmb,
                            vc_nonSyn=c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation","3'Flank",
"3'UTR", "5'Flank","5'UTR","Intron","RNA","Silent","IGR","Splice_Region"))
saveRDS(maf_3,"~/Desktop/maf_all_.var.515.rds")


pdf("soma3.pdf")
setwd("projects-t3/PTRAN/Projects_starting_Jan2022/dnaseq/analysis/DNA/gatk/ysong_test/analysis_foundation_tempus/tempus_IGS_MAF/")
path=getwd()
file.names <- dir(path, pattern =".maf")

meta.soma=read.table("meta.soma.txt",sep="\t",header=T)
maf_3=maftools:::merge_mafs(file.names,
                            vc_nonSyn=c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation"))



list.mafs= list.files(path, pattern=".maf", full.names=T)
list.all.maf.files <- lapply(list.mafs, function(i){
  read.delim(i, sep = "\t", header = TRUE, fill = TRUE, comment.char = "#")
})

var.classes=c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation")

#var.classes <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation","Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation","3'Flank", "3'UTR", "5'Flank",
 #                "5'UTR","Intron","RNA","Silent","Splice_Region", "DE_NOVO_START_IN_FRAME", "DE_NOVO_START_OUT_FRAME", "START_CODON_SNP")


merged_mafs_TN_1000G <- maftools::merge_mafs(list.all.maf.files, vc_nonSyn=var.classes)

clin=unique(merged_mafs_TN_1000G@data$ClinVar_VCF_CLNSIG)
sub=subsetMaf(merged_mafs_TN_1000G, query=  "ClinVar_VCF_CLNSIG %in% clin[c(2,4,5,8,11,13,16)]",dropLevels =F)

saveRDS(sub,"non_synonous_tempus.patho.like-patho.patho-conflict.gatk.rds")
#d <- merge_mafs(lapply(Sys.glob("mafs/Patient*.maf"), read.maf))

# Load sample information
#c <- read.table(file="sample-information.tsv", sep="\t", header=T)  
#setDT

# Combine MAF and sample info
#d@clinical.data <- c
merged_mafs_TN_1000G_df <- merged_mafs_TN_1000G@data
merged_mafs_TN_1000G_sub<-merged_mafs_TN_1000G_df[,c("Hugo_Symbol","Entrez_Gene_Id", "Chromosome", "Start_Position", "End_Position", "Tumor_Sample_Barcode", "Variant_Classification", "Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2", "tumor_f","t_alt_count","t_ref_count","n_alt_count","n_ref_count","dbSNP_RS","cDNA_Change","Protein_Change","ClinVar_VCF_CLNSIG","ClinVar_VCF_ID","COSMIC_overlapping_mutations")]



merged_mafs_TN_1000G_sub$tot_depth_tumor <- merged_mafs_TN_1000G_sub$t_alt_count + merged_mafs_TN_1000G_sub$t_ref_count
merged_mafs_TN_1000G_sub$tot_depth_normal <- merged_mafs_TN_1000G_sub$n_alt_count + merged_mafs_TN_1000G_sub$n_ref_count
merged_mafs_TN_1000G_sub$DNA.Mutations <- paste(merged_mafs_TN_1000G_sub$Hugo_Symbol,merged_mafs_TN_1000G_sub$Protein_Change)


merged_mafs_TN_1000G_sub <- merged_mafs_TN_1000G_sub %>% filter(tumor_f < 0.40 & tot_depth_tumor>=10)
