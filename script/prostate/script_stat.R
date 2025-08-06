
maf_data=readRDS("/local/projects-t3/PTRAN/Projects_starting_Jan2025/dnaseq/ysong/IGS_pipeline_output/maf/raw_batch1/combine_all_ai.711.rds")


# Count mutations per gene
gene_summary <- maf_data@data %>%
  group_by(Hugo_Symbol) %>%
  summarise(MutatedSamples = n_distinct(Tumor_Sample_Barcode)) %>%
  arrange(desc(MutatedSamples))

# Get the top 10 mutated genes


# Plot the top mutated genes
sample_load_plot=ggplot(gene_summary, aes(x = reorder(Hugo_Symbol, MutatedSamples), y = MutatedSamples)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top Mutated Genes", x = "Genes", y = "Number of Mutated Samples") +
  theme_minimal()



# Count variant classifications
variant_summary <- maf_data@data %>%
  group_by(Variant_Classification) %>%
  summarise(Count = n()) %>%
  arrange(desc(Count))

# Plot the variant classification distribution
variant_class=ggplot(variant_summary, aes(x = reorder(Variant_Classification, Count), y = Count, fill = Variant_Classification)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Variant Classification Distribution", x = "variant_class", y = "Count") +
  theme_minimal() +
  theme(legend.position = "none")



# Count mutations per sample
sample_summary <- maf_data@data%>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(TotalMutations = n()) %>%
  arrange(desc(TotalMutations))

# Plot the mutation load across samples
count_mutation=mutation_load=ggplot(sample_summary, aes(x = reorder(Tumor_Sample_Barcode, TotalMutations), y = TotalMutations)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Mutation Load", x = "Samples", y = "Total Mutations") +
  theme_minimal()






# Define transitions and transversions
transitions <- c("A>G", "G>A", "C>T", "T>C")
transversions <- c("A>C", "A>T", "G>C", "G>T", "C>A", "C>G", "T>A", "T>G")

# Create mutation types
# Access the mutation data
mutation_data <- maf_data@data

# Add a new column for mutation type
mutation_data$MutationType <- paste(mutation_data$Reference_Allele, mutation_data$Tumor_Seq_Allele2, sep = ">")

library(dplyr)
titv_summary <- mutation_data %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(
    Transitions = sum(MutationType %in% transitions),
    Transversions = sum(MutationType %in% transversions),
    TiTvRatio = ifelse(Transversions == 0, NA, Transitions / Transversions)  # Avoid division by zero
  )

# Plot Ti/Tv Ratio
library(ggplot2)
titv_plot=ggplot(titv_summary, aes(x = Tumor_Sample_Barcode, y = TiTvRatio)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Transition/Transversion Ratio", x = "Sample", y = "Ti/Tv Ratio") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))




# Access mutation data from the MAF object
mutation_data <- maf_data@data

# Summarize counts of Variant_Type
variant_type_stats <- mutation_data %>%
  group_by(Variant_Type) %>%
  summarise(Count = n()) %>%
  arrange(desc(Count))

# Print the summary
print(variant_type_stats)

# Visualize the statistics
library(ggplot2)
variant_type=ggplot(variant_type_stats, aes(x = reorder(Variant_Type, Count), y = Count, fill = Variant_Type)) +
  geom_bar(stat = "identity") +
  labs(title = "Distribution of Variant Types", x = "Variant Type", y = "Count") +
  theme_minimal() + coord_flip() +
  theme(legend.position = "none")



# Create SNV_Class if not present
mutation_data$SNV_Class <- paste(mutation_data$Reference_Allele, mutation_data$Tumor_Seq_Allele2, sep = ">")

# Summarize counts of SNV_Class
snv_class_stats <- mutation_data %>%
  filter(Variant_Type == "SNP") %>%  # Focus only on SNPs
  group_by(SNV_Class) %>%
  summarise(Count = n()) %>%
  arrange(desc(Count))

# Print the summary
print(snv_class_stats)

# Visualize the statistics
snv_class=ggplot(snv_class_stats, aes(x = reorder(SNV_Class, Count), y = Count, fill = SNV_Class)) +
  geom_bar(stat = "identity") +
  labs(title = "Distribution of SNV Classes", x = "SNV Class", y = "Count") +
  theme_minimal() +coord_flip()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")


#install.packages("patchwork")
library(patchwork)

# Combine plots



# Check for open devices
print(dev.list())


while (!is.null(dev.list())) {
  dev.off()
}



# Save the grid to PDF
pdf(paste(index.name,"mutation.summary.pdf",sep="."), width = 12, height = 6)

# Arrange ggplots in a grid
grid_plot <- grid.arrange(variant_class, snv_class, variant_type, count_mutation, ncol = 2)
# Close any lingering devices



# Plot 1: grid.draw (gridExtra plot)
#grid.draw(grid_plot)

dev.off()  # Start another new page

while (!is.null(dev.list())) {
  dev.off()
}




#tumor mutation burden based on all mutation
#tmb.maf=data.frame(tmb(maf_data,captureSize =2.6))



# Close any lingering devices
while (!is.null(dev.list())) {
  dev.off()
}
mutation_stat=data.frame(maf_data@summary)[,1:2]
#mutation_stat <- rbind(mutation_stat, c("total_TMB_perMB", tmb.maf$total_perMB) ) # Use mean TMB for simplicity

write.table(mutation_stat,"mutation_summary.txt",sep="\t",quote=F,row.names=F)
write.table(maf_data@data,"mutation.txt",sep="\t",quote=F,row.names=F)

write.table(maf_data@gene.summary,"mutation_gene_stat.txt",sep="\t",quote=F,row.names=F)
