library(maftools)
#setwd("/Users/ysong/Desktop/P_Tran_ref/maf/soma/")
path=getwd()
clinical_tmb=read.table("/Users/ysong/Desktop/P_Tran_ref/onco_508/soma/meta_no97.txt",sep="\t",header=T,check.names = F)
#file.names <- dir(path, pattern =".maf.mod.txt")

setwd("/Volumes/projects-t3/PTRAN/Projects_starting_Jan2022/dnaseq/analysis/DNA/gatk/ysong_test/analysis_foundation_tempus/tempus_IGS_MAF/")
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

merged_mafs_TN_1000G_df <- merged_mafs_TN_1000G@data
merged_mafs_TN_1000G_sub<-merged_mafs_TN_1000G_df[,c("Hugo_Symbol","Entrez_Gene_Id", "Chromosome", "Start_Position", "End_Position", "Tumor_Sample_Barcode", "Variant_Classification", "Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2", "tumor_f","t_alt_count","t_ref_count","n_alt_count","n_ref_count","dbSNP_RS","cDNA_Change","Protein_Change","ClinVar_VCF_CLNSIG","ClinVar_VCF_ID","COSMIC_overlapping_mutations")]



merged_mafs_TN_1000G_sub$tot_depth_tumor <- merged_mafs_TN_1000G_sub$t_alt_count + merged_mafs_TN_1000G_sub$t_ref_count
merged_mafs_TN_1000G_sub$tot_depth_normal <- merged_mafs_TN_1000G_sub$n_alt_count + merged_mafs_TN_1000G_sub$n_ref_count
merged_mafs_TN_1000G_sub$DNA.Mutations <- paste(merged_mafs_TN_1000G_sub$Hugo_Symbol,merged_mafs_TN_1000G_sub$Protein_Change)


merged_mafs_TN_1000G_sub <- merged_mafs_TN_1000G_sub %>% filter(tumor_f < 0.40 & tot_depth_tumor>=10)

sub1=read.maf(merged_mafs_TN_1000G_sub)
sub1=subsetMaf(sub1, query=  "ClinVar_VCF_CLNSIG %in% clin[c(2,4,5,8,11,13,16)]",dropLevels =F)
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
saveRDS(merged_maf_TN_1000G_sub,"/local/projects-t3/PTRAN/Projects_starting_Jan2022/dnaseq/analysis/DNA/gatk/ysong_test/analysis_foundation_tempus/tempus_IGS_MAF.rds")


#tempus=readRDS("/Users/ysong/Desktop/combine_tempus_foundation/tempus.filter.209.rds")



### foundation


path="/local/projects-t3/PTRAN/Projects_starting_Jan2022/dnaseq/analysis/DNA/gatk/ysong_test/analysis_foundation_tempus/maf/"

list.mafs= list.files(path, pattern=".maf", full.names=T)
list.all.maf.files <- lapply(list.mafs, function(i){
  read.delim(i, sep = "\t", header = TRUE, fill = TRUE, comment.char = "#")
})

#var.classes=c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation","3'Flank",
 #           "3'UTR", "5'Flank","5'UTR","Intron","RNA","Silent","IGR","Splice_Region","DE_NOVO_START_IN_FRAME", "DE_NOVO_START_OUT_FRAME",
  #          "START_CODON_SNP")

merged_mafs_TN_1000G <- maftools::merge_mafs(list.all.maf.files , vc_nonSyn=var.classes)

library(data.table)
meta=read.table("maf_finish_163.txt",sep="\t",header=T)

names(meta)[1:3]=c("sample_ID1","sample_ID2","Done")
meta=meta[1:69]
library(data.table)
meta.sub1=setDT(meta)
merged_mafs_TN_1000G@clinical.data= meta.sub1

merged_mafs_TN_1000G_df <- merged_mafs_TN_1000G@data
merged_mafs_TN_1000G_sub<-merged_mafs_TN_1000G_df[,c("Hugo_Symbol","Entrez_Gene_Id", "Chromosome", "Start_Position", "End_Position", "Tumor_Sample_Barcode", "Variant_Classification", "Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2", "tumor_f","t_alt_count","t_ref_count","n_alt_count","n_ref_count","dbSNP_RS","cDNA_Change","Protein_Change","ClinVar_VCF_CLNSIG","ClinVar_VCF_ID","COSMIC_overlapping_mutations")]

#Added columns
merged_mafs_TN_1000G_sub$tot_depth_tumor <- merged_mafs_TN_1000G_sub$t_alt_count + merged_mafs_TN_1000G_sub$t_ref_count
merged_mafs_TN_1000G_sub$tot_depth_normal <- merged_mafs_TN_1000G_sub$n_alt_count + merged_mafs_TN_1000G_sub$n_ref_count
merged_mafs_TN_1000G_sub$DNA.Mutations <- paste(merged_mafs_TN_1000G_sub$Hugo_Symbol,merged_mafs_TN_1000G_sub$Protein_Change)


merged_mafs_TN_1000G_sub <- merged_mafs_TN_1000G_sub %>% filter(tumor_f < 0.40 & tot_depth_tumor>=10)
saveRDS(merged_maf_TN_1000G_sub,"/local/projects-t3/PTRAN/Projects_starting_Jan2022/dnaseq/analysis/DNA/gatk/ysong_test/analysis_foundation_tempus/foundation_IGS_MAF.rds")





##########not use#####





library(dplyr)
found=read.table("foundation_gene.tab",header=F,sep="\t")

merged_mafs_TN_1000G_sub2 <- merged_mafs_TN_1000G_sub %>% 
  mutate(type = case_when(
    merged_mafs_TN_1000G_sub$Hugo_Symbol %in% found$V1 ~ "Foundation", 
    !merged_mafs_TN_1000G_sub$Hugo_Symbol %in% found$V1 ~ "Non-Foundation"))


sub.maf <- merged_mafs_TN_1000G_sub2 %>% filter(type %in% "Foundation" )

sub_filtered <- sub.maf %>% filter(tumor_f < 0.40 & tot_depth_tumor>=10)
saveRDS(sub

tempus=readRDS("/Users/ysong/Desktop/combine_tempus_foundation/tempus.filter.209.rds")
tempus_gene=read.table("/Users/ysong/Desktop/P_Tran_ref/analysis_result/freebay_plot/820_correct/tab_file/tempus_marker_clean.txt",sep="\t",header=T)
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
