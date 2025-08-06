library(maftools)
meta=read.delim("/Volumes/projects-t3/PTRAN/Projects_starting_Jan2025/dnaseq/ysong/IGS_pipeline_output/maf/meta_w_ai_711.txt",header=T,sep="\t")
foundation=readRDS("/Volumes/projects-t3/PTRAN/Projects_starting_Jan2022/dnaseq/analysis/DNA/gatk/ysong_test/analysis_foundation_tempus/foundation_IGS_MAF_correctd.rds")
tempus=readRDS("/Volumes/projects-t3/PTRAN/Projects_starting_Jan2022/dnaseq/analysis/DNA/gatk/ysong_test/analysis_foundation_tempus/tempus_IGS_MAF/RDS/non_synonous_tempus.gatk.rds")

merged_mafs_TN_1000G_sub=foundation@data
merged_mafs_TN_1000G_sub$tot_depth_tumor <- merged_mafs_TN_1000G_sub$t_alt_count + merged_mafs_TN_1000G_sub$t_ref_count
merged_mafs_TN_1000G_sub$tot_depth_normal <- merged_mafs_TN_1000G_sub$n_alt_count + merged_mafs_TN_1000G_sub$n_ref_count
merged_mafs_TN_1000G_sub$DNA.Mutations <- paste(merged_mafs_TN_1000G_sub$Hugo_Symbol,merged_mafs_TN_1000G_sub$Protein_Change)

library(dplyr)
merged_mafs_TN_1000G_sub <- merged_mafs_TN_1000G_sub %>% filter(tumor_f < 0.40 & tot_depth_tumor>=10)
library(maftools)
tempus.gene=read.table("~/Desktop/P_Tran_ref/tempus/tempus_gene_update.version.txt",sep="\t",header=F)
foundation.gene=read.table("/Users/ysong/Desktop/P_Tran_ref/foundation/maf_validation/foundation_gene_to_plot.txt",,sep="\t",header=F)
share_gene=merge(foundation.gene,tempus.gene,by="V1")


input_found=merged_mafs_TN_1000G_sub[which(merged_mafs_TN_1000G_sub$Hugo_Symbol%in%share_gene$V1),]


foundation.nonsy.f=read.maf(input_found)
saveRDS(foundation.nonsy.f,"~/Desktop/foundation.nonsyn.share.gene.rds")


merged_mafs_TN_1000G_sub=tempus@data
merged_mafs_TN_1000G_sub$tot_depth_tumor <- merged_mafs_TN_1000G_sub$t_alt_count + merged_mafs_TN_1000G_sub$t_ref_count
merged_mafs_TN_1000G_sub$tot_depth_normal <- merged_mafs_TN_1000G_sub$n_alt_count + merged_mafs_TN_1000G_sub$n_ref_count
merged_mafs_TN_1000G_sub$DNA.Mutations <- paste(merged_mafs_TN_1000G_sub$Hugo_Symbol,merged_mafs_TN_1000G_sub$Protein_Change)

library(dplyr)
merged_mafs_TN_1000G_sub <- merged_mafs_TN_1000G_sub %>% filter(tumor_f < 0.40 & tot_depth_tumor>=10)

library(dplyr)

# List of all unique barcodes (samples) before filtering
all_samples <- unique(merged_mafs_TN_1000G_sub$Tumor_Sample_Barcode)

# Filter for shared genes
shared_gene_maf <- merged_mafs_TN_1000G_sub %>%
  filter(Hugo_Symbol %in% share_gene$V1)

# Add missing samples back in with NA rows
missing_samples <- setdiff(all_samples, unique(shared_gene_maf$Tumor_Sample_Barcode))

# Create NA rows for missing samples
na_rows <- data.frame(Tumor_Sample_Barcode = missing_samples,
                      Hugo_Symbol = NA,
                      stringsAsFactors = FALSE)

# Bind real data with NA placeholders
input_tempus <- bind_rows(shared_gene_maf, na_rows)
tempus.nonsy.f=read.maf(input_tempus)
saveRDS(tempus.nonsy.f,"~/Desktop/tempus.nonsyn.w.share.gene.rds")

new=readRDS("/Volumes/projects-t3/PTRAN/Projects_starting_Jan2025/dnaseq/ysong/IGS_pipeline_output/maf/all.nonsyn.rds")

merged_mafs_TN_1000G_sub=new@data
merged_mafs_TN_1000G_sub$tot_depth_tumor <- merged_mafs_TN_1000G_sub$t_alt_count + merged_mafs_TN_1000G_sub$t_ref_count
merged_mafs_TN_1000G_sub$tot_depth_normal <- merged_mafs_TN_1000G_sub$n_alt_count + merged_mafs_TN_1000G_sub$n_ref_count
merged_mafs_TN_1000G_sub$DNA.Mutations <- paste(merged_mafs_TN_1000G_sub$Hugo_Symbol,merged_mafs_TN_1000G_sub$Protein_Change)

library(dplyr)
merged_mafs_TN_1000G_sub <- merged_mafs_TN_1000G_sub %>% filter(tumor_f < 0.40 & tot_depth_tumor>=10)

library(dplyr)

# List of all unique barcodes (samples) before filtering
library(dplyr)

# List of all unique barcodes (samples) before filtering
all_samples <- unique(merged_mafs_TN_1000G_sub$Tumor_Sample_Barcode)

# Filter for shared genes
shared_gene_maf <- merged_mafs_TN_1000G_sub %>%
  filter(Hugo_Symbol %in% share_gene$V1)

# Add missing samples back in with NA rows
missing_samples <- setdiff(all_samples, unique(shared_gene_maf$Tumor_Sample_Barcode))

# Create NA rows for missing samples
na_rows <- data.frame(Tumor_Sample_Barcode = missing_samples,
                      Hugo_Symbol = NA,
                      stringsAsFactors = FALSE)

# Bind real data with NA placeholders
input_new <- bind_rows(shared_gene_maf, na_rows)


new.nonsy.f=read.maf(input_new)

saveRDS(new.nonsy.f,"~/Desktop/new.nonsyn.share.gene.rds")

all.non=merge_mafs(maf=list(new.nonsy.f,tempus.nonsy.f,foundation.nonsy.f))
saveRDS(all.non,"~/Desktop/all.nonsyn.share.gene.805.rds")

pdf("plot.all.nonsy.pdf",10,20)
oncoplot(maf = all.non,genes = unique(all.non@gene.summary$Hugo_Symbol),writeMatrix = T)
dev.off()
