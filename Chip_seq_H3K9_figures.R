library(ggplot2) #v 3.42.4

#Function to read data and calculate corrected coverage
read_and_correct_coverage <- function(input_file, h3k9_file, lib_size_input, lib_size_h3k9) {
  input <- read.delim(input_file, header=FALSE)
  h3k9 <- read.delim(h3k9_file, header=FALSE)
  
  h3k9$corrected_coverage <- h3k9$V4 / lib_size_h3k9
  input$corrected_coverage <- input$V4 / lib_size_input
  
  ratio <- log2(h3k9$corrected_coverage / input$corrected_coverage)
  
  df <- input[,1:3]
  df$ratio <- ratio
  colnames(df) <- c("scaff", "start", "end", "log2ratio")
  
  # Filter for specific scaffolds
  df <- df[df$scaff %in% paste0("Tps_LRv5b_scf", 1:12), ]
  
  return(df)
}

# Function to plot data
plot_data <- function(df) {
  df$scaff_number <- as.numeric(gsub("^.*scf", "", df$scaff))
  
  df2 <- df[order(df$scaff_number, df$start),]
  df2$cumul_start <- seq(0, by = 30000, length.out = nrow(df2))
  df2$chrom <- ifelse(df2$scaff_number == 3, "sex-chrom", "autosome")
  df3 <- df2[df2$scaff_number == 3,]

#modify y limits    
  ggplot(df2, aes(x = cumul_start, y = log2ratio, color = as.factor(scaff_number))) +
    geom_point(size = 0.5) +
    scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(df2$scaff_number))) / 2)) +
    scale_y_continuous(limits = c(-1, 1)) +
    geom_point(data = df3, aes(x = cumul_start, y = log2ratio), colour = "blue", size = 0.5) +
    ylab("log2(H3K9/Input)") +
    xlab("Genomic coordinates") +
    geom_hline(yintercept = c(-0.5, 0, 0.5), linetype = c("dashed", "solid", "dashed"), color = "black") +
    theme(
      axis.ticks.x = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "none"
    )
}

# Data frame that contains scaffold information start and end of each 30K window and covearge
# Testes data
testes_input <- "~/Desktop/Tps_testes_input1_GW_coverage_30K_windows_DR.txt"
testes_h3k9 <- "~/Desktop/Tps_testes_h3k9A_GW_coverage_30K_windows_DR.txt"
lib_size_testes_h3k9 <- 57363039 
lib_size_testes_input <- 53495659

# Plot testes data
testes_df <- read_and_correct_coverage(testes_input, testes_h3k9, lib_size_testes_input, lib_size_testes_h3k9)
plot_data(testes_df)

# Males somatic data
input_males_somatic <- "~/Desktop/Tps_male_organs_input3_GW_coverage_30W_DR.txt"
males_somatic_h3k9 <- "~/Desktop/Tps_male_organs_h3k9A_GW_coverage_30W_DR.txt"
lib_size_male_somatic_H3K9 <- 62833968
lib_size_male_somatic_input <- 78355731

# Plot males somatic data
males_somatic_df <- read_and_correct_coverage(input_males_somatic, males_somatic_h3k9, lib_size_male_somatic_input, lib_size_male_somatic_H3K9)
plot_data(males_somatic_df)
