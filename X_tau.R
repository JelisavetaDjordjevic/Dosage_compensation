# Load necessary libraries
library(edgeR)
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)

# Read data
counts <- read.csv("/Users/jdjordje/Documents/GitHub/Dosage_compensation/Data/dosage_counts.csv", sep = "\t")
scaffold_gene_size <- read.csv("/Users/jdjordje/Documents/GitHub/Dosage_compensation/Data/scaffold_gene_size_temp.txt", sep = "")

# Rename each column
for (col in 1:ncol(counts)) {
  colnames(counts)[col] <- sub("_extR.*", "", colnames(counts)[col])
}

# Remove specific gene
counts <- subset(counts, Gene_name != "Tps_037574")

# Get names
names <- as.character(counts$Gene_name)

# Subset scaffold_gene_size based on names
scaffold_gene_size <- scaffold_gene_size[scaffold_gene_size$group_name %in% names,]

# Extract sample information
samples <- colnames(counts[,2:ncol(counts)])
sp_sex_tissue_stage <- sapply(strsplit(samples, split='_ind', fixed=TRUE), function(x) (x[1]))
sp <- sapply(strsplit(samples, split='_', fixed=TRUE), function(x) (x[1]))
sex <- sapply(strsplit(samples, split='_', fixed=TRUE), function(x) (x[2]))
tissue <- sapply(strsplit(samples, split='_', fixed=TRUE), function(x) (x[3]))
stage <- sapply(strsplit(samples, split='_', fixed=TRUE), function(x) (x[4]))
sex_tissue_stage <- paste(sex, tissue, stage, sep="_")

# Normalization
mycounts1 <- as.numeric(scaffold_gene_size$exonic.gene.sizes)
tmycounts1 <- t(mycounts1)
rownames(counts)<-counts$Gene_name
d <- DGEList(counts = counts[,2:ncol(counts)], genes = rownames(counts), group = sex_tissue_stage)
d$genes$Length <- c(tmycounts1)
keep_a <- rowSums(cpm(d) > 0.5) >= 1
d <- d[keep_a, , keep.lib.sizes = FALSE]
d <- calcNormFactors(d)
d_rpkm <- rpkm(d)
df1 <- as.data.frame(d_rpkm)
df1 <- cbind(row.names(d_rpkm), d_rpkm)
colnames(df1)[1] <- "Gene_name"

# Convert df1 to a data frame
df1 <- as.data.frame(df1)

# Gather data into long format
newcounts <- df1 %>%
  gather(var, rpkm, Tps_F_A_J6_ind59:Tps_M_A_Ha_ind4) %>%
  separate(var, c('sp', 'sex', 'tissue', 'stage', 'individual'), sep = '_')

# Add sex_tissue_stage column
newcounts$sex_tissue_stage <- paste(newcounts$sex, newcounts$tissue, newcounts$stage, sep = "_")

# Convert rpkm to numeric
newcounts$rpkm <- as.numeric(as.character(newcounts$rpkm))
levels(as.factor(newcounts$tissue))

# Summarize data
data <- newcounts %>% 
  group_by(Gene_name, sex, stage, tissue) %>% 
  summarise(rpkm = median(rpkm))

# Spread data into wide format
df <- spread(data, key = tissue, value = rpkm)
df <- df[, 1:7]

# Select relevant stages
yes <- df[df$stage == "Ha" | df$stage == "Ad" | df$stage == "J3" | df$stage == "J4",]

# Define Tau function
fTau <- function(x) {
  if(all(!is.na(x))) {
    if(min(x, na.rm=TRUE) >= 0) {
      if(max(x) != 0) {
        x <- (1 - (x / max(x)))
        res <- sum(x, na.rm=TRUE)
        res <- res / (length(x) - 1)
      } else {
        res <- 0
      }
    } else {
      res <- NA
    }
  } else {
    res <- NA
  }
  return(res)
}

# Separate tables by sex
#females
fem <- yes[yes$sex == "F",]
f_df <- fem[, c(4, 5, 6, 7)]
rownames(f_df) <- paste(fem$Gene_name, fem$stage, sep="_")

#males
mal <- yes[yes$sex == "M",]
m_df <- mal[, c(4, 5, 6, 7)]
rownames(m_df) <- paste(mal$Gene_name, mal$stage, sep="_")

# Calculate Tau for females and males
Tau_f <- apply(f_df, 1, fTau)
Tau_m <- apply(m_df, 1, fTau)

# Bind Tau results
Tau <- cbind(Tau_f, Tau_m)
Tau_fin <- Tau
names <- rownames(Tau)
rownames(Tau_fin) <- NULL
Tau_fin <- cbind(names, Tau_fin, fem$Gene_name)
Tau_fin <- as.data.frame(Tau_fin)

# Use "scaffold_gene_size.txt" for assigning genes to scaffolds
Tau_fin$seq_name <- scaffold_gene_size$seqnames[match(Tau_fin$V4, scaffold_gene_size$group_name)]

# Filter based on scaffold names
Tau_fin <- Tau_fin[Tau_fin$seq_name %in% paste0("Tps_LRv5b_scf", 1:13),]

# Convert seq_name to ordered factor
Tau_fin$seq_name <- factor(Tau_fin$seq_name, levels = paste0("Tps_LRv5b_scf", 1:13), ordered = TRUE)

# Melt the data frame
melted <- melt(Tau_fin, id.vars = c("names", "V4", "seq_name"), variable.name = "variable", value.name = "Tau_value")

# Print the melted data frame
colnames(Tau_fin)

# Separate Tau_fin based on variable
sub_f <- subset(melted, variable == "Tau_f")
melted$Tau_value <- as.numeric(melted$Tau_value)

# Create ggplot for visualization
scaff <- ggplot(melted, aes(x = variable, y = Tau_value, fill = seq_name)) + 
  geom_boxplot() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),  
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  scale_fill_manual(values = c("gray", "grey", "orange", "grey", "grey", "gray", "gray", "grey", "grey", "gray", "grey", "grey"))

scaff
# Statistical test
test_scaff <- compare_means(Tau_value ~ seq_name, group.by = "variable", data = melted, p.adjust.method = "BH")

# Add another factor to color differently boxes
melted$chrom <- rep("auto", nrow(melted))
melted$chrom[melted$seq_name == "Tps_LRv5b_scf3"] <- "x_chrom"
melted$fact <- paste(melted$chrom, melted$variable, sep = "_")

melted$variable<-as.character(melted$variable)
melted$variable[melted$variable=="Tau_f"]<-"F"
melted$variable[melted$variable=="Tau_m"]<- "M"

# Create ggplot for comparison
X_vs_A <- ggplot(melted, aes(x = variable, y = Tau_value, fill = fact)) + 
  geom_boxplot() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),  
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  xlab("Sex") +
  ylab("Tissue specificity [Tau]") +
  scale_fill_manual(values = c("gray", "white", "red", "blue"),
                    labels = c('Autosomes (F)', 'Autosomes (M)', 'X (F)', 'X (M)')) +
  theme_bw(base_size = 24) +
  theme(legend.title = element_blank())

X_vs_A

# Statistical test
test_x_vs_a <- compare_means(Tau_value ~ fact, group.by = "variable", data = melted, p.adjust.method = "BH")

########### repeating Tau for only somatic tissues (= excluding gonads) #########

f_df<-f_df[,c(1,2,4)]
rownames(f_df) <- paste(fem$Gene_name, fem$stage, sep="_")
m_df<-m_df[,c(1,2,4)]
rownames(m_df) <- paste(mal$Gene_name, mal$stage, sep="_")
# Calculate Tau for females and males
Tau_f <- apply(f_df, 1, fTau)
Tau_m <- apply(m_df, 1, fTau)

# Bind Tau results
Tau <- cbind(Tau_f, Tau_m)
Tau_fin <- Tau
names <- rownames(Tau)
rownames(Tau_fin) <- NULL
Tau_fin <- cbind(names, Tau_fin, fem$Gene_name)
Tau_fin <- as.data.frame(Tau_fin)

# Use "scaffold_gene_size.txt" for assigning genes to scaffolds
Tau_fin$seq_name <- scaffold_gene_size$seqnames[match(Tau_fin$V4, scaffold_gene_size$group_name)]

# Filter based on scaffold names
Tau_fin <- Tau_fin[Tau_fin$seq_name %in% paste0("Tps_LRv5b_scf", 1:13),]

# Convert seq_name to ordered factor
Tau_fin$seq_name <- factor(Tau_fin$seq_name, levels = paste0("Tps_LRv5b_scf", 1:13), ordered = TRUE)

# Melt the data frame
melted <- melt(Tau_fin, id.vars = c("names", "V4", "seq_name"), variable.name = "variable", value.name = "Tau_value")

# Print the melted data frame
colnames(Tau_fin)

# Separate Tau_fin based on variable
sub_f <- subset(melted, variable == "Tau_f")
melted$Tau_value <- as.numeric(melted$Tau_value)

# Create ggplot for visualization
scaff <- ggplot(melted, aes(x = variable, y = Tau_value, fill = seq_name)) + 
  geom_boxplot() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),  
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  scale_fill_manual(values = c("gray", "grey", "orange", "grey", "grey", "gray", "gray", "grey", "grey", "gray", "grey", "grey"))+
  ylab("Tissue specificity [Tau]")+
  theme_bw(base_size = 20) +
  theme(legend.title = element_blank())


scaff
# Statistical test
test_scaff <- compare_means(Tau_value ~ seq_name, group.by = "variable", data = melted, p.adjust.method = "BH")

# Add another factor to color differently boxes
melted$chrom <- rep("auto", nrow(melted))
melted$chrom[melted$seq_name == "Tps_LRv5b_scf3"] <- "x_chrom"
melted$fact <- paste(melted$chrom, melted$variable, sep = "_")

melted$variable<-as.character(melted$variable)
melted$variable[melted$variable=="Tau_f"]<-"F"
melted$variable[melted$variable=="Tau_m"]<- "M"

# Create ggplot for comparison
X_vs_A <- ggplot(melted, aes(x = variable, y = Tau_value, fill = fact)) + 
  geom_boxplot() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),  
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  xlab("Sex") +
  ylab("Tissue specificity [Tau]") +
  scale_fill_manual(values = c("gray", "white", "red", "blue"),
                    labels = c('Autosomes (F)', 'Autosomes (M)', 'X (F)', 'X (M)')) +
  theme_bw(base_size = 20) +
  theme(legend.title = element_blank())

X_vs_A

####combine plots 
ggarrange(X_vs_A, scaff,
          labels = c("A", "B"),
          ncol = 2)
###separate by stage

# Statistical test
test_x_vs_a <- compare_means(Tau_value ~ fact, group.by = "variable", data = melted, p.adjust.method = "BH")

# Save results as a table
write.table(test_scaff, file = "stats_tau_scaff.txt", sep = "\t", row.names = FALSE)
write.table(test_x_vs_a, file = "stats_tau_no_GO.txt", sep = "\t", row.names = FALSE)
getwd()
