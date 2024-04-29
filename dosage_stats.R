#Load packages
library(edgeR)  # v 3.42.4
library(ggplot2)  # v 3.4.4
library(reshape2)  # v 1.4.4
library(ggpubr)  # v 0.6.0
library(tidyverse) # v 2.0.0

# Define function to process each stage
process_stage <- function(stage_name, pattern, tissue) {
  
  br <- counts[, grepl(pattern, colnames(counts))]
  br <- cbind(counts$Gene_name, br)
  colnames(br)[1] <- "Gene_name"
  
  # Identify columns containing "F" or "M"
  female_cols <- grep("F", colnames(br), value = TRUE)
  male_cols <- grep("M", colnames(br), value = TRUE)
  
  # Normalize and calculate RPKM
  d <- DGEList(counts = br, genes = rownames(br))
  d$genes$Length <- as.numeric(scaffold_gene_size$exonic.gene.sizes)
  keep <- rowSums(cpm(d) > 0.5) >= 2
  d <- d[keep, , keep.lib.sizes = FALSE]
  d <- calcNormFactors(d)
  br_rpkm <- rpkm(d)
  
  # Compute average RPKM for males and females
  average_M <- rowMeans(br_rpkm[, male_cols])
  average_F <- rowMeans(br_rpkm[, female_cols])
  df <- data.frame(average_F, average_M, d$genes$Gene_name)
  colnames(df)[3] <- "Gene_name"
  df$m_f_rpkm <- log2(average_M / average_F)
  df$stage <- stage_name
  df$tissue <- tissue
  
  return(df)
}

# Read data
counts <- read.csv("/Users/jdjordje/Documents/GitHub/Dosage_compensation/dosage_counts.csv", sep = "\t")
scaffold_gene_size <- read.csv("/Users/jdjordje/Documents/GitHub/Dosage_compensation/scaffold_gene_size_temp.txt", sep = "")

# Rename columns
for (col in 1:ncol(counts)) {
  colnames(counts)[col] <- sub("_extR.*", "", colnames(counts)[col])
}

# Subset data
counts <- subset(counts, Gene_name != "Tps_037574")
names<-as.character(counts$Gene_name)
scaffold_gene_size<-scaffold_gene_size[scaffold_gene_size$group_name %in% names,]


###### TISSUES ###### 
# Process each stage and tissue 
#brain
dosage_ad_br <- process_stage("adult", "_B_Ad_18\\.", "brain")
dosage_ha_br <- process_stage("hatchling", "_B_Ha_", "brain")
dosage_j3_br <- process_stage("juv3", "_B_J3_", "brain")
dosage_j2_br <- process_stage("juv2", "_B_J2_", "brain")
dosage_j4_br <- process_stage("juv4", "_B_J4_", "brain")

#antenna
dosage_ad_ant <- process_stage("adult", "_A_Ad_18\\.", "antenna")
dosage_ha_ant <- process_stage("hatchling", "_A_Ha_", "antenna")
dosage_j3_ant <- process_stage("juv3", "_A_J3_", "antenna")
dosage_j2_ant <- process_stage("juv2", "_A_J2_", "antenna")
dosage_j4_ant <- process_stage("juv4", "_A_J4_", "antenna")

#gut
dosage_ad_gut <- process_stage("adult", "_Gu_Ad_18\\.", "gut")
dosage_ha_gut <- process_stage("hatchling", "_Gu_Ha_", "gut")
dosage_j3_gut <- process_stage("juv3", "_Gu_J3_", "gut")
dosage_j2_gut <- process_stage("juv2", "_Gu_J2_", "gut")
dosage_j4_gut <- process_stage("juv4", "_Gu_J4_", "gut")

#leg
dosage_ha_leg <- process_stage("hatchling", "_Lg_Ha_", "leg")
dosage_j3_leg <- process_stage("juv3", "_Lg_J3_", "leg")
dosage_j2_leg <- process_stage("juv2", "_Lg_J2_", "leg")
dosage_j4_leg <- process_stage("juv4", "_Lg_J4_", "leg")

#gonad
dosage_ad_gonad <- process_stage("adult", "_Go_Ad_", "gonad")
dosage_ha_gonad <- process_stage("hatchling", "_Go_Ha_", "gonad")
dosage_j4_gonad <- process_stage("juv4", "_Go_J4_", "gonad")

# Combine results for all tissues
dosage_all <- rbind(dosage_ad_br, dosage_ha_br, dosage_j3_br, dosage_j2_br, dosage_j4_br,
                    dosage_ad_ant, dosage_ha_ant, dosage_j3_ant, dosage_j2_ant, dosage_j4_ant,
                    dosage_ad_gut, dosage_ha_gut, dosage_j3_gut, dosage_j2_gut, dosage_j4_gut,
                    dosage_ha_leg, dosage_j3_leg, dosage_j2_leg, dosage_j4_leg, dosage_ad_gonad, 
                    dosage_ha_gonad, dosage_j4_gonad)


# Add chromosome information
dosage_all$scaff <- scaffold_gene_size$seqnames[match(dosage_all$Gene_name, scaffold_gene_size$group_name)]
dosage_all$chrom <- ifelse(dosage_all$scaff == "Tps_LRv5b_scf3", "X_chrom", "Auto")

# Rename stages
dosage_all$stage <- ifelse(dosage_all$stage == "hatchling", "N1",
                           ifelse(dosage_all$stage == "juv2", "N2",
                                  ifelse(dosage_all$stage == "juv3", "N3",
                                         ifelse(dosage_all$stage == "juv4", "N4",
                                                ifelse(dosage_all$stage == "adult", "A", dosage_all$stage)))))



# List of tissues
tissues <- c("brain", "antenna", "gut", "leg", "gonad")

# Perform statistical analysis for each tissue 
stats_ratio <- map_dfr(tissues, ~ {
  tissue_data <- filter(dosage_all, .data$tissue == .x)
  result <- compare_means(m_f_rpkm ~ chrom, group.by = "stage", data = tissue_data, p.adjust.method = "BH")
  mutate(result, tissue = .x)
})


# Get plot for average RPKM for each sex and each type of chromosome (X and autosomes)
# Melt the data excluding "m_f_rpkm" and "scaff"
melted_data<- reshape2::melt(dosage_all[, !names(dosage_all) %in% c("m_f_rpkm", "scaff")], 
                                  id.vars = c("Gene_name", "stage", "tissue", "chrom"))

# Create the average_value_chrom column by combining and chrom and variable
melted_data$sex_chrom<- as.factor(paste(melted_data$chrom, melted_data$variable, sep = "_"))


# Perform statistical analysis for each tissue 
stats_averageRPKM <- map_dfr(tissues, ~ {
  tissue_data <- filter(melted_data, .data$tissue == .x)
  result <- compare_means(value ~ sex_chrom, group.by = "stage", data = tissue_data, p.adjust.method = "BH")
  mutate(result, tissue = .x)
})

# Save results as a table
write.table(stats_ratio, file = "stats_dosage_ratio.txt", sep = "\t", row.names = FALSE)
write.table(stats_averageRPKM, file = "stats_dosage_average_rpkm.txt", sep = "\t", row.names = FALSE)
