library(edgeR)
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(tidyr)

newcounts <- read.delim("~/Documents/GitHub/Dosage_compensation/Data/newcounts")

get_factors <- function(sub) {
  fac<-colnames(sub)
  print(fac)
  rows <- sapply(strsplit(fac, split='_', fixed=TRUE), function(x) (x[c(1,2,4)]))
  fact <- paste(rows[1,],  rows[2,], rows[3,], sep='_')
  return(fact)
}

#code to do analysis 
analysis <- function(subset, element, contrast) { 
  print("---------------------------------")
  print("element")
  print(element)
  print("contrast")
  print(contrast)
  #subset is a variable with  rownames=genenames 
  subset <- as.data.frame(subset)
  subset <- subset  %>%  tidyr::unite("variables", sp:individual)
  subset <-  tidyr::spread(  subset, variables, names)
  print (head(subset))
  rownames(subset) <- subset$Gene_name
  subset <- as.matrix(subset[,-1])
  sp_sex_stage<-get_factors(subset) # what you subsetted by
  print("sp_sex_stage")
  print(sp_sex_stage)
  d <- DGEList(counts=subset , genes =rownames(subset),  group = sp_sex_stage)
  #normalization and filtering of the genes with low expression
  d$samples
  keep_a <- rowSums(cpm(d[,grep("F", colnames(d))])>0.5)>=(ncol(d[,grep("F", colnames(d))])-1) | rowSums(cpm(d[,grep("M", colnames(d))])>0.5) >=(ncol(d[,grep("M", colnames(d))])-1)
  head(keep_a)
  d<- d[keep_a, , keep.lib.sizes=FALSE]
  filtered<-as.data.frame(table(keep_a))
  fil<-filtered$Freq[filtered$keep_a=="FALSE"]
  kept<-filtered$Freq[filtered$keep_a=="TRUE"]
  fil
  kept
  d<- calcNormFactors(d)
  d$samples
  # writing the result
  #write.table(br_rpkm, file = "br_adult_rpkm.txt", dec = ".", sep = "\t", quote = FALSE)  
  #design a model 
  design<- model.matrix(~ 0 + sp_sex_stage, data= d$samples)
  colnames(design)<- levels(d$samples$group)
  print("design")
  print(design)
  design
  #Dispersion
  d <- estimateDisp(d, design, robust=TRUE)
  summary(d$trended.dispersion)
  sqrt(d$common.dispersion) #this is BCV value
  plot_bcv<-plotBCV(d)
  
  #Fitting a GLM model 
  fit <- glmQLFit(d, design, robust=TRUE)
  head(fit$coefficients)
  plotQLDisp(fit)
  summary(fit$df.prior)
  
  #between m and F
  #my_contrast<-get_contrast(sp_sex_stage)
  
  astr=paste(contrast, collapse=",")
  prestr="makeContrasts("
  poststr=",levels=design)"
  commandstr=paste(prestr,astr,poststr,sep="")
  FM<-eval(parse(text=commandstr))
  
  #FM<- makeContrasts(cont, levels=design)
  #FM<- makeContrasts(Tps_F_Ad-Tps_M_Ad, levels=design)
  res_FM<- glmQLFTest(fit, contrast=FM)
  is.de_FM <- decideTestsDGE(res_FM, p.value=0.05)
  summary(is.de_FM)
  SBG<-as.data.frame(summary(is.de_FM))
  colnames(SBG)<-c("Dirrection", "Contrast", "Number_of_G")
  
  total_sbg<-(SBG$Number_of_G[[1]])+SBG$Number_of_G[[3]]
  total_sbg
  total_genes<-SBG$Number_of_G[[1]]+SBG$Number_of_G[[3]]+SBG$Number_of_G[[2]]
  total_sbg_perc<-(total_sbg/total_genes)*100
  down<-SBG$Number_of_G[[1]]
  down
  up<-SBG$Number_of_G[[3]]
  up
  down_perc<-(SBG$Number_of_G[[1]]/(SBG$Number_of_G[[1]]+SBG$Number_of_G[[2]]+SBG$Number_of_G[[3]]))*100
  down_perc
  up_perc<-(SBG$Number_of_G[[3]]/(SBG$Number_of_G[[1]]+SBG$Number_of_G[[2]]+SBG$Number_of_G[[3]]))*100
  up_perc
  
  cols<-colnames(subset)[1] 
  cols
  ok<-strsplit(cols, "_")[[1]]
  tissue<-ok[3]
  print("Tissue")
  print(tissue)
  stage<-ok[4]
  print("Stage")
  print(stage)
  topTags(res_FM)
  message_filtering<-sprintf("In the %s tissue, within the %s stage, %d of gene were filtered out, while %d were kept: ", tissue, stage, fil, kept)
  message <- sprintf("In the %s tissue, within the %s stage, between the %s, there are %d sex-biased genes out of %d genes in total, in percentages %f,with %d up and %d  down genes, or %f percentage and %f percentage of expressed genes", tissue, stage, contrast, total_sbg, total_genes,total_sbg_perc, up, down, up_perc, down_perc)
  
  
  cat(message_filtering,message, "\n")
  
  top_fm<- topTags(res_FM, n= 15000)
  top_fm<- top_fm$table
  #get a table with FDR and gene names for the TopGO= GO terms
  #hat_top_go<- hat_top_fm[c(1,6)]
  
  ###lets add SB factor 
  top_fm$SB[top_fm$logFC>0 & top_fm$FDR<0.05]<- "FB"
  top_fm$SB[top_fm$logFC>0 & top_fm$FDR>0.05]<- "UB"
  top_fm$SB[top_fm$logFC<0 & top_fm$FDR<0.05]<- "MB"
  top_fm$SB[top_fm$logFC<0 & top_fm$FDR>0.05]<- "UB"
  factor_counts <- table(top_fm$SB)
  print(factor_counts)
  #read just_genes file to assign X vs Auto
  just_genes<-read.csv("/Users/jdjordje/Documents/Jeli/job/Tps_Tdi_development/Tps_LRv5b_mtDNAv350_v2.1/Just_genes.txt", sep="")
  top_fm$chrom<- just_genes$seqnames[match(top_fm$genes,just_genes$gene_id)]
  top_fm$chrom[top_fm$chrom=="Tps_LRv5b_scf3"]<-"X"
  top_fm$chrom[!top_fm$chrom=="X"]<-"Auto"
  top_fm$chrom
  #get distribution of sex biased genes on X and auto
  xvsa<-table(top_fm$SB,top_fm$chrom)
  print("Sex- bias on X vs A")
  print(xvsa)
  xvsa<-as.data.frame(xvsa)
  xvsa
  colnames(xvsa)<-c("SB", "Chrom", "Freq")
  xvsa<-reshape(xvsa, idvar = "SB", timevar = "Chrom", direction = "wide")
  xvsa$tissue<-rep(tissue, nrow(xvsa))
  xvsa$stage<-rep(stage, nrow(xvsa))
  #return(list(top_table=top_fm,prop_onX=xvsa,message = message))
  res<-xvsa
  return(res)
}

# Your vector of strings

my_vector<- c("Tps_Go_Ad", "Tps_Go_Ha", "Tps_Go_J4")


# Example list of contrasts
contrast_list <- c("Tps_F_Ad-Tps_M_Ad", "Tps_F_Ha-Tps_M_Ha", "Tps_F_J4-Tps_M_J4")

# List to store results
#results_list <- list()
df<-data.frame()

# Function to perform analysis
perform_analysis <- function(element, contrast) {
  # Split the element based on "_"
  elements <- strsplit(element, "_")[[1]]
  
  # Assign the elements to variables
  Species <- elements[1]
  Tissue <- elements[2]
  Stage <- elements[3]
  
  subset <- newcounts %>% filter(sp == Species, tissue == Tissue, stage == Stage)
  #subset <- newcounts %>% filter(sp == "Tps", tissue == "Go", stage == "Ad")
  #contrast<-"Tps_F_Ad-Tps_M_Ad"
  
  # Call the analysis function with the assigned values and current contrast
  result <- analysis(subset, element, contrast)
  
  # Store the result in the list
  print("=================")
  print(element)
  print(contrast)
  #results_list[[paste0(element, "_", contrast)]]  <- result
  #    row<-paste0(element, "_", contrast)
  df<-rbind(df,result)
}

# Apply the analysis function to each element and each contrast
number_analysis<-length(my_vector)

for (i in 1:number_analysis) {
  perform_analysis(my_vector[i], contrast_list[i])
}


