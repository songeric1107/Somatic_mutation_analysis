all.maf=readRDS("/local/projects-t3/PTRAN/Projects_starting_Jan2025/dnaseq/ysong/tempus_foundation_w_conflicP.417.s289.rds")
#all.maf.sub=subsetMaf(all.maf,query="Tumor_Sample_Barcode%in%meta1f$Tumor_Sample_Barcode",dropLevels = F)

data.all=all.maf@data

id=data.frame(meta1f$Tumor_Sample_Barcode)
colnames(id)="Tumor_Sample_Barcode"
data.all.ann=merge(id,data.all,by="Tumor_Sample_Barcode",all.x=T)
batch1=read.maf(data.all.ann,clinicalData = meta1f[-c(1:2,19:39)])

#saveRDS(batch1,"batch1.rds")


batch2=readRDS("/local/projects-t3/PTRAN/Projects_starting_Jan2025/dnaseq/ysong/IGS_pipeline_output/maf/all.patho.w.conflic.patho.update.rds")
meta2=read.table("/local/projects-t3/PTRAN/Projects_starting_Jan2025/dnaseq/ysong/IGS_pipeline_output/maf/raw_batch1/batch2_w_ai.txt",sep="\t",header=T)

data2=batch2@data

id=data.frame(meta2$Tumor_Sample_Barcode)
colnames(id)="Tumor_Sample_Barcode"
data2a=merge(id,data2,by="Tumor_Sample_Barcode",all.x=T)

batch2=read.maf(data2a,clinicalData = meta2)



foundation_gene=read.table("/local/projects-t3/PTRAN/Projects_starting_Jan2025/dnaseq/ysong/IGS_pipeline_output/maf/raw_batch1/Download FoundationOne_CDx_324_genes.txt",sep="\t",header=F)
tempus_gene=read.table("/local/projects-t3/PTRAN/Projects_starting_Jan2025/dnaseq/ysong/IGS_pipeline_output/maf/raw_batch1/tempus_gene_update.version.txt",sep="\t",header=F)
share_gene=merge(foundation_gene,tempus_gene,by="V1")

batch1.sub=subsetMaf(batch1,genes=share_gene$V1,dropLevels = F)
saveRDS(batch1.sub,"/local/projects-t3/PTRAN/Projects_starting_Jan2025/dnaseq/ysong/IGS_pipeline_output/maf/raw_batch1/batch1.share.gene.rds")
batch2.sub=subsetMaf(batch2,genes=share_gene$V1,dropLevels = F)
saveRDS(batch2.sub,"/local/projects-t3/PTRAN/Projects_starting_Jan2025/dnaseq/ysong/IGS_pipeline_output/maf/raw_batch1/batch2.share.gene.rds")
library(maftools)

# Step 1: Merge MAFs
batch.all <- merge_mafs(list(batch1.sub, batch2.sub))

# Step 2: Merge clinical data manually (preserve all sample barcodes)
m1=data.frame(batch1.sub@clinical.data[,c(1:3)])
m2=data.frame(batch2.sub@clinical.data[,c(1:3)])
pdf("batch1.pdf")
oncoplot(batch1.sub,genes=batch1.sub@gene.summary$Hugo_Symbol,removeNonMutated = F, writeMatrix = T)
dev.off()

pdf("batch2.pdf")
oncoplot(batch2.sub,genes=batch2.sub@gene.summary$Hugo_Symbol,removeNonMutated = F, writeMatrix = T)
dev.off()


colnames(m2)=colnames(m1)
m12=rbind(m1,m2)

# Step 3: Overwrite clinical data in the merged object
batch.all@clinical.data <- as.data.table(m12)
id.all=data.frame(m12$Tumor_Sample_Barcode)
colnames(id.all)="Tumor_Sample_Barcode"
mall=data.frame(batch.all@data)
mut_data <- batch.all@data

all_samples <- unique(m12$Tumor_Sample_Barcode)
mutated_samples <- unique(mut_data$Tumor_Sample_Barcode)
nonmutated_samples <- setdiff(all_samples, mutated_samples)

# Create placeholder rows (minimal valid structure)
dummy_rows <- mut_data[1, , drop = FALSE]
dummy_rows <- dummy_rows[rep(1, length(nonmutated_samples)), ]
dummy_rows$Tumor_Sample_Barcode <- nonmutated_samples

# Fill other required columns with safe defaults
cols_to_na <- setdiff(names(dummy_rows), "Tumor_Sample_Barcode")
dummy_rows[, (cols_to_na) := lapply(.SD, function(x) NA), .SDcols = cols_to_na]

dummy_rows$Variant_Classification <- "Silent"
dummy_rows$Hugo_Symbol <- "NA"



# Combine real + dummy mutation rows
merged_mut_data <- rbind(mut_data, dummy_rows)




# Determine columns to blank
cols_to_na <- setdiff(names(dummy_rows), "Tumor_Sample_Barcode")

# Set NA for those columns
dummy_rows[, (cols_to_na) := lapply(.SD, function(x) NA), .SDcols = cols_to_na]

# Add minimal valid fields
dummy_rows$Variant_Classification <- "Silent"
dummy_rows$Hugo_Symbol <- "NA"

merged_mut_data <- rbind(mut_data, dummy_rows, fill = TRUE)
all <- read.maf(maf = merged_mut_data, clinicalData = m12)

saveRDS(all,"/local/projects-t3/PTRAN/Projects_starting_Jan2025/dnaseq/ysong/IGS_pipeline_output/maf/raw_batch1/combine_all_ai.711.rds")
