# Dosage compensation in *Timema Popense*

Repository for the collected scripts used in study: 

*DATA*

Contains 
1. Timema popense read counts- **dosage_counts.csv**
2. Gene length and scaffold information- **scaffold_gene_size_temp.txt**
3. Chip-seq coverage for Input- **Tps_testes_input1_GW_coverage_30K_windows_DR.txt**
4. Chip-seq coverage for H3K9- **Tps_testes_h3k9A_GW_coverage_30K_windows_DR.txt**

*SCRIPTS* 

*To analyze Dosage compensation*

Calculate RPKM values, ratio of male to female expression for X chromosome and autosomes and generate figures

   1. **dosage_script.R** 
   
To perform statistical tests 

   2. **dosage_stats.R**

To calculate specificity of gene expression (Tau index) 

   3. **X_tau.R**

To plot the coverage across scaffolds for H3K9 chip seq data

   4. **Chip_seq_H3K9_figures.R**

To calculate sex-bias gene expression in Gonad tissue
   5. **Sex_bias_gonads.R**
   

