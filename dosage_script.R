#Load packages
library(edgeR)  # Version 3.42.4
library(ggplot2)  # Version 3.4.4
library(reshape2)  # Version 1.4.4
library(ggpubr)  # Version 0.6.0
library(cowplot)  # Version 1.1.1

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
counts <- read.csv("/Users/jdjordje/Documents/GitHub/Dosage_compensation/Data/dosage_counts.csv",sep = "\t")
scaffold_gene_size <- read.csv("/Users/jdjordje/Documents/GitHub/Dosage_compensation/Data/scaffold_gene_size_temp.txt", sep = "")

# Rename columns
for (col in 1:ncol(counts)) {
  colnames(counts)[col] <- sub("_extR.*", "", colnames(counts)[col])
}

# Subset data
counts <- subset(counts, Gene_name != "Tps_037574")
names<-as.character(counts$Gene_name)
names
scaffold_gene_size<-scaffold_gene_size[scaffold_gene_size$group_name %in% names,]

###### SOMATIC TISSUES ###### 
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



# Combine results for all tissues
dosage_all <- rbind(dosage_ad_br, dosage_ha_br, dosage_j3_br, dosage_j2_br, dosage_j4_br,
                    dosage_ad_ant, dosage_ha_ant, dosage_j3_ant, dosage_j2_ant, dosage_j4_ant,
                    dosage_ad_gut, dosage_ha_gut, dosage_j3_gut, dosage_j2_gut, dosage_j4_gut,
                    dosage_ha_leg, dosage_j3_leg, dosage_j2_leg, dosage_j4_leg)


# Add chromosome information
dosage_all$scaff <- scaffold_gene_size$seqnames[match(dosage_all$Gene_name, scaffold_gene_size$group_name)]
dosage_all$chrom <- ifelse(dosage_all$scaff == "Tps_LRv5b_scf3", "X_chrom", "Auto")

# Rename stages
dosage_all$stage <- ifelse(dosage_all$stage == "hatchling", "N1",
                           ifelse(dosage_all$stage == "juv2", "N2",
                                  ifelse(dosage_all$stage == "juv3", "N3",
                                         ifelse(dosage_all$stage == "juv4", "N4",
                                                ifelse(dosage_all$stage == "adult", "A", dosage_all$stage)))))

dosage_all$tissue <- ifelse(dosage_all$tissue == "antenna", "Antenna",
                           ifelse(dosage_all$tissue == "brain", "Brain",
                                  ifelse(dosage_all$tissue == "gut", "Gut",
                                         ifelse(dosage_all$tissue == "leg", "Leg",
                                                ifelse(dosage_all$tissue == "gonad", "Gonad", dosage_all$tissue)))))
#relevel
dosage_all$stage <-factor(dosage_all$stage,
                   levels = c('N1','N2','N3','N4','A'),ordered = TRUE)

# Visualization
p <- ggplot(dosage_all, aes(x = stage, y = m_f_rpkm, fill = chrom)) + 
  geom_boxplot(outlier.shape = NA) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),  
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_hline(yintercept = -1, linetype = "dashed", color = "black") +
  ylim(-4, 4) +
  ylab("log2(Male RPKM/Female RPKM)") + 
  xlab("Developmental stage")+
  facet_wrap(~tissue, scales = "free_y", ncol = 1) +
  #ggtitle("Expression Variation Across Tissues") +
  scale_fill_manual(values = c("gray", "orange"),
                    labels=c('Autosome', 'X'))+
  theme_bw(base_size = 16) +
  theme(legend.title = element_blank())

print(p)

# Plot for average RPKM for each sex and each type of chromosome (X and autosomes)
# Melt the data excluding "m_f_rpkm" and "scaff"
melted_data <- reshape2::melt(dosage_all[, !names(dosage_all) %in% c("m_f_rpkm", "scaff")], 
                              id.vars = c("Gene_name", "stage", "tissue", "chrom"))

# Create the average_value_chrom column by combining and chrom and variable
melted_data$average_value_chrom <- as.factor(paste(melted_data$chrom, melted_data$variable, sep = "_"))

# Plot of Average RPKM
somatic_avr <- ggplot(melted_data, aes(x = stage, y = value, fill = average_value_chrom)) + 
  geom_boxplot(outlier.shape = NA) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),  
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  ylab("Average RPKM")+ 
  xlab("Developmental stage")+
  ylim(0, 80)+ # Modify the y-axis limits as needed
  scale_fill_manual(values=c("gray","white","red","blue"),
                    labels=c('Autosome (F)', 'Autosome (M)', 'X (F)', 'X (M)'))+
  facet_wrap(~tissue, ncol = 1) +
  theme_bw(base_size = 16)+
  theme(legend.title = element_blank())

# Print the plot
print(somatic_avr)

####combine plots (ratio and average RPKM) 
ggarrange(p, somatic_avr,
          labels = c("A", "B"),
          ncol = 2)


#################### GONAD TISSUE #################################

dosage_ad_gonad <- process_stage("adult", "_Go_Ad_", "gonad")
dosage_ha_gonad <- process_stage("hatchling", "_Go_Ha_", "gonad")
dosage_j4_gonad <- process_stage("juv4", "_Go_J4_", "gonad")

# Combine results for all tissues
dosage_gon <- rbind(dosage_ad_gonad, dosage_ha_gonad, dosage_j4_gonad)

# Add chromosome information
dosage_gon$scaff <- scaffold_gene_size$seqnames[match(dosage_gon$Gene_name, scaffold_gene_size$group_name)]
dosage_gon$chrom <- ifelse(dosage_gon$scaff == "Tps_LRv5b_scf3", "X_chrom", "Auto")

# Rename stages
dosage_gon$stage <- ifelse(dosage_gon$stage == "hatchling", "N1",
                           ifelse(dosage_gon$stage == "juv2", "N2",
                                  ifelse(dosage_gon$stage == "juv3", "N3",
                                         ifelse(dosage_gon$stage == "juv4", "N4",
                                                ifelse(dosage_gon$stage == "adult", "A", dosage_gon$stage)))))

#relevel
dosage_gon$stage <-factor(dosage_gon$stage,
                          levels = c('N1','N2','N3','N4','A'),ordered = TRUE)

# Visualization
gon_ratio <- ggplot(dosage_gon, aes(x = stage, y = m_f_rpkm, fill = chrom)) + 
  geom_boxplot(outlier.shape = NA) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),  
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_hline(yintercept = -1, linetype = "dashed", color = "black") +
  ylim(-4, 4) +
  ylab("log2(Male RPKM/ Female RPKM)") + 
  xlab("Developmental stage")+
  #facet_wrap(~tissue, scales = "free_y", ncol = 1) +
  ggtitle("Gonads") +
  scale_fill_manual(values = c("gray", "orange"),
                    labels=c('Autosome','X'))+
  theme_bw(base_size = 16)+
  theme(legend.title = element_blank())

print(gon_ratio)


# Get plot for average RPKM for each sex and each type of chromosome (X and autosomes)
# Melt the data excluding "m_f_rpkm" and "scaff"
melted_data_gon <- reshape2::melt(dosage_gon[, !names(dosage_gon) %in% c("m_f_rpkm", "scaff")], 
                              id.vars = c("Gene_name", "stage", "tissue", "chrom"))

# Create the average_value_chrom column by combining and chrom and variable
melted_data_gon$sex_chrom<- as.factor(paste(melted_data_gon$chrom, melted_data_gon$variable, sep = "_"))

# Plot of Average RPKM
gonad_avr <- ggplot(melted_data_gon, aes(x = stage, y = value, fill = sex_chrom)) + 
  geom_boxplot(outlier.shape = NA) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),  
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  ylab("Average RPKM")+ 
  xlab("Developmental stage")+
  ylim(0, 80)+ # Modify the y-axis limits as needed
  scale_fill_manual(values=c("gray","white","red","blue"),
                    labels=c('Autosome (F)', 'Autosome (M)', 'X (F)', 'X (M)'))+
  #facet_wrap(~tissue, ncol = 1) +
  theme_bw(base_size = 16) +
  theme(legend.title = element_blank())

# Print the plot
print(gonad_avr)

# Extract only X chromosome data
dosage_x <- melted_data_gon[melted_data_gon$chrom == "X_chrom", ]

# Modify stage names
dosage_x$stage[dosage_x$stage == "adult"] <- "A"
dosage_x$stage[dosage_x$stage == "hatch"] <- "N1"

levels(dosage_x$sex_chrom)

# Remove unused levels from factors
dosage_x$stage <- droplevels(dosage_x$stage)
dosage_x$sex_chrom <- droplevels(dosage_x$sex_chrom)

levels(dosage_x$sex_chrom)

# Convert factor to character
dosage_x$sex_chrom <- as.character(dosage_x$sex_chrom)

# Modify sex_chrom names
dosage_x$sex_chrom[dosage_x$sex_chrom == "X_chrom_average_F"] <- "Females"
dosage_x$sex_chrom[dosage_x$sex_chrom == "X_chrom_average_M"] <- "Males"

# Convert back to factor
dosage_x$sex_chrom <- factor(dosage_x$sex_chrom)

# Reorder stage factor
dosage_x$stage <- factor(dosage_x$stage, 
                         levels = c('N1','N4','A'), ordered = TRUE)

# Density plot
legend_title <- "X chr"
gon_dens <- ggplot(dosage_x, aes(x = log2(value), fill = sex_chrom)) + 
  geom_density(alpha = 0.5) +   
  scale_fill_manual(legend_title, values = c("red", "blue")) +
  xlab("Log2(RPKM)") +
  ylab("Density") +
  facet_wrap(~stage, nrow = 1) +
  theme_bw(base_size = 16) +
  theme(legend.title = element_blank())

gon_dens

## combine all plots for gonads
ggdraw() +
  draw_plot(gon_ratio, x = 0, y = .5, width = .5, height = .5) +
  draw_plot(gonad_avr, x = .5, y = .5, width = .5, height = .5) +
  draw_plot(gon_dens, x = 0, y = 0, width = 1, height = 0.5) +
  draw_plot_label(label = c("A)", "B)", "C)"), size = 16,
                  x = c(0, 0.5, 0), y = c(1, 1, 0.5))
