#conda activate r_cibersort
#https://bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/oncoplots.html


library(maftools)
all.syn=readRDS("/local/projects-t3/PTRAN/Projects_starting_Jan2025/dnaseq/ysong/IGS_pipeline_output/maf/all.include.syn.rds")

path="/Volumes/projects-t3/PTRAN/Projects_starting_Jan2025/dnaseq/ysong/IGS_pipeline_output/maf/raw_batch1/combine_all_ai.711.rds"

all.p=readRDS(path)

pdf("pathway.pdf")
pws = pathways(maf = all.p, plotType = 'treemap')
plotPathways(maf = all.p, pathlist = pws)
dev.off()


mutation_matrix <- mutCountMatrix(maf = all.p, includeSyn = FALSE,removeNonMutated = F) 
write.table(mutation_matrix,"/local/projects-t3/PTRAN/Projects_starting_Jan2025/dnaseq/ysong/IGS_pipeline_output/maf/raw_batch1/numeric.mat.txt",sep="\t",quote=F)

mutation_rate <- rowSums(mutation_matrix > 0, na.rm = TRUE) / ncol(mutation_matrix)
cat("Mutation rates:\n")
print(mutation_rate)

# Filter genes with mutation rate > 0.1% (0.001)
num.matrix.f <- mutation_matrix[mutation_rate > 0.01, ]



pdf("summary.pdf")
plotmafSummary(all)
oncoplot(all,genes=all@gene.summary$Hugo_Symbol,removeNonMutated = F,writeMatrix = T)
dev.off()

library(maftools)
meta=read.delim("/local/projects-t3/PTRAN/Projects_starting_Jan2025/dnaseq/ysong/IGS_pipeline_output/maf/meta_w_ai_711.txt",header=T,sep="\t")

input=read.table("/local/projects-t3/PTRAN/Projects_starting_Jan2025/dnaseq/ysong/IGS_pipeline_output/maf/onco_matrix.all.patho.txt",sep="\t",header=T,row.names=1,check.names=F)



col = c("Missense_Mutation" = "blue", "In_Frame_Ins" = "red", "Splice_Site" = "#008000","In_Frame_Del"="purple","Frame_Shift_Ins"="yellow",
        "Frame_Shift_Del"="grey", "Nonsense_Mutation"="brown", "Multi_Hit"="orange")


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
  Multi_Hit= function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["Multi_Hit"], col = NA))}
)

meta.all=meta
meta.all=meta.all[order(meta.all$ARTERA.Score),]
mutation2f=input[which(colnames(input)%in%meta.all$ID)]


mutation.s=mutation2f[order(match(colnames(mutation2f),meta.all$ID))]        
df=mutation.s
df_clean <- df[!apply(df, 1, function(row) all(is.na(row) | row == "")), ]




ha = HeatmapAnnotation(batch = meta.all$batch, group=meta.all$met.syn.1.2.,score = anno_points(meta.all$ARTERA.Score, ylim = c(0, 1), axis = TRUE),
                       col = list(group=c("#N/A"="grey","1"="cyan","2"="orange"),batch = c("batch1" = "green", "batch2(PMH)" = "blue","batch3(Oligopelvis)"="yellow","batch4 (GETUG)"="orange","MatthewD"="red")),
                       #annotation_height = unit(c(5, 5, 15), "mm"),
                       annotation_legend_param = list(batch = list(title = "batch")))


library(dplyr)
library(ComplexHeatmap)

# Ensure sample IDs in df_clean match meta.all
library(dplyr)

# Ensure sample IDs are in the metadata and match column names of df_clean
meta.all$sample_id <- colnames(df_clean)

# Sort by group (met.syn.1.2.) and ARTERA.Score (increasing order)
meta.reorder <- meta.all %>%
  dplyr::select(sample_id, met.syn.1.2., ARTERA.Score) %>%
  dplyr::group_by(met.syn.1.2.) %>%
  dplyr::arrange(ARTERA.Score, .by_group = TRUE)

# Apply the new column order to your matrix and metadata
ordered_samples <- meta.reorder$sample_id
df_clean <- df_clean[, ordered_samples]
meta.all <- meta.all[match(ordered_samples, meta.all$sample_id), ]

# Create heatmap annotation
ha <- HeatmapAnnotation(
  batch = meta.all$batch,
  group = meta.all$met.syn.1.2.,
  score = anno_points(
    meta.all$ARTERA.Score,
    ylim = c(0, 1),
    axis = TRUE,
    gp = gpar(fontsize = 8)  # Adjust font size for readability
  ),
  col = list(
    group = c("#N/A" = "grey", "1" = "cyan", "2" = "orange"),
    batch = c(
      "batch1" = "green",
      "batch2(PMH)" = "blue",
      "batch3(Oligopelvis)" = "yellow",
      "batch4 (GETUG)" = "orange",
      "MatthewD" = "red"
    )
  ),
  annotation_legend_param = list(batch = list(title = "batch"))
)

# Generate heatmap
pdf("~/Desktop/mutation_all.score.pdf", width = 10, height = 20)
oncoPrint(
  df_clean,
  col = col,  # Ensure col is defined (e.g., colors for alterations)
  alter_fun = alter_fun,  # Ensure alter_fun is defined (e.g., alteration functions)
  column_split = meta.reorder$met.syn.1.2.,  # Use meta.reorder for consistency
  column_order = ordered_samples,  # Explicitly set column order
  remove_empty_columns = FALSE,
  remove_empty_rows = FALSE,
  column_title = "",
  top_annotation = ha
)
dev.off()


  

# -------------------------------------------------------------------------

library(ComplexHeatmap)
library(circlize)

library(dplyr)
library(ComplexHeatmap)
library(grid)  # For grid.rect and gpar

df.clean.f=df_clean[which(rownames(df_clean)%in%rownames(num.matrix.f)),]


meta.all$sample_id <- colnames(df.clean.f)
if (any(!meta.all$sample_id %in% colnames(df.clean.f))) stop("Sample IDs in meta.all do not match df_clean columns")

# Preprocess df_clean to handle NA and empty strings
df.clean.f[is.na(df.clean.f)] <- ""  # Replace NA with empty string
df.clean.f[df.clean.f == ""] <- ""  # Ensure consistency (optional)


# Define colors for each alteration type
col <- c(
  "Missense_Mutation" = "blue",
  "Nonsense_Mutation" = "red",
  "Frame_Shift_Del" = "purple",
  "Splice_Site" = "orange",
  "Frame_Shift_Ins" = "green",
  "Multi_Hit" = "pink"
)

# Define graphic functions for each alteration type
alter_fun <- list(
  Missense_Mutation = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["Missense_Mutation"], col = NA)),
  Nonsense_Mutation = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["Nonsense_Mutation"], col = NA)),
  Frame_Shift_Del = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.45, gp = gpar(fill = col["Frame_Shift_Del"], col = NA)),
  Splice_Site = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["Splice_Site"], col = NA)),
  Frame_Shift_Ins = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.45, gp = gpar(fill = col["Frame_Shift_Ins"], col = NA)),
  Multi_Hit = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["Multi_Hit"], col = NA))
)

# Check alteration types in df_clean
alter_types <- unique(unlist(df.clean.f))
if (any(!alter_types %in% c(names(col), ""))) stop("Undefined alteration types in df_clean: ", paste(setdiff(alter_types, c(names(col), "")), collapse = ", "))

# Data integrity checks
if (!is.numeric(meta.all$ARTERA.Score)) stop("ARTERA.Score is not numeric")
if (any(is.na(meta.all$ARTERA.Score))) stop("NA values found in ARTERA.Score")

# Split data by met.syn.1.2. and sort ARTERA.Score in ascending order
meta.split <- meta.all %>%
  dplyr::select(sample_id, met.syn.1.2., ARTERA.Score, batch) %>%
  dplyr::group_by(met.syn.1.2.) %>%
  dplyr::arrange(ARTERA.Score, .by_group = TRUE) %>%
  dplyr::ungroup() %>%
  split(.$met.syn.1.2.)

# Verify expected groups
expected_groups <- c("#N/A", "1", "2")
#if (!all(expected_groups %in% names(meta.split))) stop("Missing expected met.syn.1.2. groups: ", paste(setdiff(expected_groups, names(meta.split)), collapse = ", "))

# Verify sorting within each group
print("Verifying ARTERA.Score sorting within each met.syn.1.2. group:")
lapply(names(meta.split), function(group) {
  scores <- meta.split[[group]]$ARTERA.Score
  is_sorted <- all(diff(scores) >= 0)
  cat(sprintf("Group %s: is_sorted = %s\n", group, is_sorted))
  return(data.frame(group = group, is_sorted = is_sorted))
})

# Create a list to store heatmaps
onco_plots <- list()

# Generate oncoPrint for each met.syn.1.2. group
for (g in names(meta.split)) {
  # Subset data for the current group
  meta.sub <- meta.split[[g]]
  sample_ids <- meta.sub$sample_id
  df_sub <- df.clean.f[, sample_ids, drop = FALSE]
  
  # Verify subset data
  if (any(!sample_ids %in% colnames(df.clean.f))) stop(sprintf("Sample IDs in group %s do not match df_clean columns", g))
  
  # Create heatmap annotation
  ha <- HeatmapAnnotation(
    batch = meta.sub$batch,
    score = anno_points(
      meta.sub$ARTERA.Score,
      ylim = c(0, 1),
      axis = TRUE,
      gp = gpar(fontsize = 18)  # Increased from 12 to 14 for axis labels
    ),
    col = list(
      batch = c(
        "batch1" = "green",
        "batch2(PMH)" = "blue",
        "batch3(Oligopelvis)" = "yellow",
        "batch4 (GETUG)" = "orange",
        "MatthewD" = "red"
      )
    ),
    annotation_legend_param = list(batch = list(title = "batch")),
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontsize = 16)  # Added to increase annotation name label size
  )
  # Create oncoPrint
  onco_plots[[g]] <- oncoPrint(
    df_sub,
    col = col,row_names_gp = gpar(fontsize = 14),
    alter_fun = alter_fun,
    column_order = sample_ids,
    remove_empty_columns = FALSE,
    remove_empty_rows = FALSE,
    column_title = paste("Group", g),
    top_annotation = ha
  )
}

# Draw heatmaps side by side in a single PDF
pdf("~/Desktop/mutation_all.score_split.filter.pdf", width = 30, height = 15)
ht_list <- onco_plots[["#N/A"]] + onco_plots[["1"]] + onco_plots[["2"]]
draw(ht_list, ht_gap = unit(1, "cm"), heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# -------------------------------------------------------------------------


# Ensure sample IDs in df_clean match meta.all
meta.all$sample_id <- colnames(df_clean)
if (any(!meta.all$sample_id %in% colnames(df_clean))) stop("Sample IDs in meta.all do not match df_clean columns")

# Preprocess df_clean to handle NA and empty strings
df_clean[is.na(df_clean)] <- ""  # Replace NA with empty string
df_clean[df_clean == ""] <- ""  # Ensure consistency (optional)


# Define colors for each alteration type
col <- c(
  "Missense_Mutation" = "blue",
  "Nonsense_Mutation" = "red",
  "Frame_Shift_Del" = "purple",
  "Splice_Site" = "orange",
  "Frame_Shift_Ins" = "green",
  "Multi_Hit" = "pink"
)

# Define graphic functions for each alteration type
alter_fun <- list(
  Missense_Mutation = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["Missense_Mutation"], col = NA)),
  Nonsense_Mutation = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["Nonsense_Mutation"], col = NA)),
  Frame_Shift_Del = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.45, gp = gpar(fill = col["Frame_Shift_Del"], col = NA)),
  Splice_Site = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["Splice_Site"], col = NA)),
  Frame_Shift_Ins = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.45, gp = gpar(fill = col["Frame_Shift_Ins"], col = NA)),
  Multi_Hit = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["Multi_Hit"], col = NA))
)

# Check alteration types in df_clean
alter_types <- unique(unlist(df_clean))
if (any(!alter_types %in% c(names(col), ""))) stop("Undefined alteration types in df_clean: ", paste(setdiff(alter_types, c(names(col), "")), collapse = ", "))

# Data integrity checks
if (!is.numeric(meta.all$ARTERA.Score)) stop("ARTERA.Score is not numeric")
if (any(is.na(meta.all$ARTERA.Score))) stop("NA values found in ARTERA.Score")

# Split data by met.syn.1.2. and sort ARTERA.Score in ascending order
meta.split <- meta.all %>%
  dplyr::select(sample_id, met.syn.1.2., ARTERA.Score, batch) %>%
  dplyr::group_by(met.syn.1.2.) %>%
  dplyr::arrange(ARTERA.Score, .by_group = TRUE) %>%
  dplyr::ungroup() %>%
  split(.$met.syn.1.2.)

# Verify expected groups
expected_groups <- c("#N/A", "1", "2")
#if (!all(expected_groups %in% names(meta.split))) stop("Missing expected met.syn.1.2. groups: ", paste(setdiff(expected_groups, names(meta.split)), collapse = ", "))

# Verify sorting within each group
print("Verifying ARTERA.Score sorting within each met.syn.1.2. group:")
lapply(names(meta.split), function(group) {
  scores <- meta.split[[group]]$ARTERA.Score
  is_sorted <- all(diff(scores) >= 0)
  cat(sprintf("Group %s: is_sorted = %s\n", group, is_sorted))
  return(data.frame(group = group, is_sorted = is_sorted))
})

# Create a list to store heatmaps
onco_plots <- list()

# Generate oncoPrint for each met.syn.1.2. group
for (g in names(meta.split)) {
  # Subset data for the current group
  meta.sub <- meta.split[[g]]
  sample_ids <- meta.sub$sample_id
  df_sub <- df_clean[, sample_ids, drop = FALSE]
  
  # Verify subset data
  if (any(!sample_ids %in% colnames(df_clean))) stop(sprintf("Sample IDs in group %s do not match df_clean columns", g))
  
  # Create heatmap annotation
  ha <- HeatmapAnnotation(
    batch = meta.sub$batch,
    score = anno_points(
      meta.sub$ARTERA.Score,
      ylim = c(0, 1),
      axis = TRUE,
      gp = gpar(fontsize = 8)
    ),
    col = list(
      batch = c(
        "batch1" = "green",
        "batch2(PMH)" = "blue",
        "batch3(Oligopelvis)" = "yellow",
        "batch4 (GETUG)" = "orange",
        "MatthewD" = "red"
      )
    ),
    annotation_legend_param = list(batch = list(title = "batch")),
    show_annotation_name = TRUE
  )
  
  # Create oncoPrint
  onco_plots[[g]] <- oncoPrint(
    df_sub,
    col = col,
    alter_fun = alter_fun,
    column_order = sample_ids,
    remove_empty_columns = FALSE,
    remove_empty_rows = FALSE,
    column_title = paste("Group", g),
    top_annotation = ha
  )
}

# Draw heatmaps side by side in a single PDF
pdf("~/Desktop/mutation_all.score_split.pdf", width = 30, height = 20)
ht_list <- onco_plots[["#N/A"]] + onco_plots[["1"]] + onco_plots[["2"]]
draw(ht_list, ht_gap = unit(1, "cm"), heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()


# pathway -----------------------------------------------------------------


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

pw_ddr=data.frame(gene= c("BRCA1", "BRCA2", "ATM", "RAD51", "PALB2", "CHEK2",  # HRR
                                 "PRKDC", "XRCC6", "XRCC5", "LIG4",                  # NHEJ
                                 "PARP1", "MUTYH", "OGG1",                          # BER
                                 "MLH1", "MSH2", "MSH6", "PMS2",                    # MMR
                                 "ERCC1", "XPC", "ERCC2",                           # NER
                                 "FANCG", "FANCA", "FANCD2",                        # FA
                                 "TP53", "CHEK1", "ATR", "MDM2", "WEE1", "H2AFX"),
                  pw = c('DDR'),
                  stringsAsFactors = F )
pw = rbind(pw_rt, pw_cc,pw_wnt,pw_tp,pw_hippo,pw_tgf,pw_chre,pw_myc,pw_notch,pw_gi,pw_ddr)


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


#input_matrix.n=read.table("~/Desktop/numeric.mat.txt",sep="\t",header=T,check.names=F,row.names=1)


input_matrix.path=merge(pw,df.clean.f,by.x="gene",by.y=0)

input_matrix.path[input_matrix.path== ''] <- NA
input_matrix.path$count=rowSums(is.na(input_matrix.path[-c(1:2)]))

input_matrix.path1=input_matrix.path[which(input_matrix.path$count!=223),]

input_matrix.path=input_matrix.path1[-ncol(input_matrix.path1)]


input_matrix.path.s=input_matrix.path[order(input_matrix.path$pw),]

rownames(input_matrix.path.s)=make.unique(input_matrix.path.s$gene)


input_matrix.path.s[is.na(input_matrix.path.s)]<-""

df_clean=input_matrix.path.s




# Remove gene and pw columns for oncoPrint
df_clean = input_matrix.path.s[, !(colnames(input_matrix.path.s) %in% c("gene", "pw"))]

# Create row annotation for pathways
row_anno = rowAnnotation(
  Pathway = input_matrix.path.s$pw,
  col = list(Pathway = pw_col),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 16)
)

# Generate oncoPrint for each met.syn.1.2. group
onco_plots <- list()

for (g in names(meta.split)) {
  # Subset data for the current group
  meta.sub <- meta.split[[g]]
  sample_ids <- meta.sub$sample_id
  df_sub <- input_matrix.path.s[, sample_ids, drop = FALSE]
  df_sub$pw=input_matrix.path.s$pw
  # Verify subset data
  if (any(!sample_ids %in% colnames(df_clean))) stop(sprintf("Sample IDs in group %s do not match df_clean columns", g))
  
  # Create heatmap annotation with increased label sizes
  ha <- HeatmapAnnotation(
    batch = meta.sub$batch,
    score = anno_points(
      meta.sub$ARTERA.Score,
      ylim = c(0, 1),
      axis = TRUE,
      gp = gpar(fontsize = 14)  # Increased axis label size
    ),
    col = list(
      batch = c(
        "batch1" = "green",
        "batch2(PMH)" = "blue",
        "batch3(Oligopelvis)" = "yellow",
        "batch4 (GETUG)" = "orange",
        "MatthewD" = "red"
      )
    ),
    annotation_legend_param = list(batch = list(title = "batch")),
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontsize = 16)  # Increased annotation name label size
  )
  
  # Create oncoPrint
  onco_plots[[g]] <- oncoPrint(
    df_sub[-ncol(df_sub)],
    col = col,
    alter_fun = alter_fun,row_split=df_sub$pw, 
    column_order = sample_ids,row_names_gp = gpar(fontsize = 14),
    remove_empty_columns = FALSE,
    remove_empty_rows = FALSE,
    column_title = paste("Group", g),
    top_annotation = ha,pct_gp=gpar(size=14),
    right_annotation = row_anno  # Add pathway annotation
  )
}

# Draw heatmaps side by side in a single PDF
pdf("~/Desktop/mutation_all.pathway.pdf", width = 30, height = 15)
ht_list <- onco_plots[["#N/A"]] + onco_plots[["1"]] + onco_plots[["2"]]
draw(ht_list, ht_gap = unit(1, "cm"), heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
