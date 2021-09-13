library(dplyr)

setwd("/omics/groups/OE0540/internal/users/ruehle/rnavelocity/vqtls/data/velocyto/")


top_iPSC = read.delim("./velocity/deterministic/celltypes/iPSC/top_qtl_results_all_FDR0.05.txt", sep ="", as.is = T)
top_mes = read.delim("./velocity/deterministic/celltypes/mesendoderm/top_qtl_results_all_FDR0.05.txt", sep="", as.is = T)
top_pre = read.delim("./velocity/deterministic/celltypes/predefendoderm/top_qtl_results_all_FDR0.05.txt", sep="", as.is = T)
top_def = read.delim("./velocity/deterministic/celltypes/defendoderm/top_qtl_results_all_FDR0.05.txt", sep="", as.is = T)

full_iPSC = read.delim("./velocity/deterministic/celltypes/iPSC/qtl_results_all.txt", sep="", as.is = T)
full_mes = read.delim("./velocity/deterministic/celltypes/mesendoderm/qtl_results_all.txt", sep="", as.is = T)
full_pre = read.delim("./velocity/deterministic/celltypes/predefendoderm/qtl_results_all.txt", sep="", as.is = T)
full_def = read.delim("./velocity/deterministic/celltypes/defendoderm/qtl_results_all.txt", sep="", as.is = T)

#create df, with all significant SNPs (FDR > 0.05) in iPSC 
iPSC_df = select(top_iPSC, snp_id, feature_id)

#for those SNPs in iPSC_df look at p-value, and sign of SNP in other stages => decide if specific to one, two or three stages. 
for (gene_index in 1:nrow(top_iPSC)){
  iPSC_df[gene_index, 3] = "Stage specific"
  snp_mes = filter(full_mes, snp_id == top_iPSC[gene_index, 1], feature_id == top_iPSC[gene_index, 19])
  snp_pre = filter(full_pre, snp_id == top_iPSC[gene_index, 1], feature_id == top_iPSC[gene_index, 19])
  snp_def = filter(full_def, snp_id == top_iPSC[gene_index, 1], feature_id == top_iPSC[gene_index, 19])
  
  if (sign(top_iPSC[gene_index, 3]) == sign(snp_mes$beta) && snp_mes$p_value < 0.05 && nrow(snp_mes) == 1){
    iPSC_df[gene_index, 3] = "Two stages"
    if (sign(top_iPSC[gene_index, 3]) == sign(snp_pre$beta) && snp_pre$p_value < 0.05 && nrow(snp_pre) == 1) {
      iPSC_df[gene_index, 3] = "Three stages"
      if (sign(top_iPSC[gene_index, 3]) == sign(snp_def$beta) && snp_def$p_value < 0.05 && nrow(snp_def) == 1) {
        iPSC_df[gene_index, 3] = "Four stages"
      }
    }
  }
  else {
    if (sign(top_iPSC[gene_index, 3]) == sign(snp_pre$beta) && snp_pre$p_value < 0.05 && nrow(snp_pre) == 1){
      iPSC_df[gene_index, 3] = "Two stages"
      if (sign(top_iPSC[gene_index, 3]) == sign(snp_def$beta) && snp_def$p_value < 0.05 && nrow(snp_def) ==1) {
        iPSC_df[gene_index, 3] = "Three stages"
      }
    }
    else if (sign(top_iPSC[gene_index, 3]) == sign(snp_def$beta) && snp_def$p_value < 0.05 && nrow(snp_def) ==1) {
      iPSC_df[gene_index, 3] = "Two stages"
    }
  }
}


#counting mesendoderm, same way as for iPSC QTL

mes_df = select(top_mes, snp_id, feature_id)

for (gene_index in 1:nrow(top_mes)){
  mes_df[gene_index, 3] = "Stage specific"
  snp_iPSC = filter(full_iPSC, snp_id == top_mes[gene_index, 1], feature_id == top_mes[gene_index, 19])
  snp_pre = filter(full_pre, snp_id == top_mes[gene_index, 1], feature_id == top_mes[gene_index, 19])
  snp_def = filter(full_def, snp_id == top_mes[gene_index, 1], feature_id == top_mes[gene_index, 19])
  
  if (sign(top_mes[gene_index, 3]) == sign(snp_iPSC$beta) && snp_iPSC$p_value < 0.05 && nrow(snp_iPSC) == 1){
    mes_df[gene_index, 3] = "Two stages"
    if (sign(top_mes[gene_index, 3]) == sign(snp_pre$beta) && snp_pre$p_value < 0.05 && nrow(snp_pre) == 1) {
      mes_df[gene_index, 3] = "Three stages"
      if (sign(top_mes[gene_index, 3]) == sign(snp_def$beta) && snp_def$p_value < 0.05 && nrow(snp_def) == 1) {
        mes_df[gene_index, 3] = "Four stages"
      }
    }
  }
  else {
    if (sign(top_mes[gene_index, 3]) == sign(snp_pre$beta) && snp_pre$p_value < 0.05 && nrow(snp_pre) == 1){
      mes_df[gene_index, 3] = "Two stages"
      if (sign(top_mes[gene_index, 3]) == sign(snp_def$beta) && snp_def$p_value < 0.05 && nrow(snp_def) ==1) {
        mes_df[gene_index, 3] = "Three stages"
      }
    }
    else if (sign(top_mes[gene_index, 3]) == sign(snp_def$beta) && snp_def$p_value < 0.05 && nrow(snp_def) ==1) {
      mes_df[gene_index, 3] = "Two stages"
    }
  }
}

#predef counting

pre_df = select(top_pre, snp_id, feature_id)

for (gene_index in 1:nrow(top_pre)){
  pre_df[gene_index, 3] = "Stage specific"
  snp_iPSC = filter(full_iPSC, snp_id == top_pre[gene_index, 1], feature_id == top_pre[gene_index, 19])
  snp_mes = filter(full_mes, snp_id == top_pre[gene_index, 1], feature_id == top_pre[gene_index, 19])
  snp_def = filter(full_def, snp_id == top_pre[gene_index, 1], feature_id == top_pre[gene_index, 19])
  
  if (sign(top_pre[gene_index, 3]) == sign(snp_iPSC$beta) && snp_iPSC$p_value < 0.05 && nrow(snp_iPSC) == 1){
    pre_df[gene_index, 3] = "Two stages"
    if (sign(top_pre[gene_index, 3]) == sign(snp_mes$beta) && snp_mes$p_value < 0.05 && nrow(snp_mes) == 1) {
      pre_df[gene_index, 3] = "Three stages"
      if (sign(top_pre[gene_index, 3]) == sign(snp_def$beta) && snp_def$p_value < 0.05 && nrow(snp_def) == 1) {
        pre_df[gene_index, 3] = "Four stages"
      }
    }
  }
  else {
    if (sign(top_pre[gene_index, 3]) == sign(snp_mes$beta) && snp_mes$p_value < 0.05 && nrow(snp_mes) == 1){
      pre_df[gene_index, 3] = "Two stages"
      if (sign(top_pre[gene_index, 3]) == sign(snp_def$beta) && snp_def$p_value < 0.05 && nrow(snp_def) ==1) {
        pre_df[gene_index, 3] = "Three stages"
      }
    }
    else if (sign(top_pre[gene_index, 3]) == sign(snp_def$beta) && snp_def$p_value < 0.05 && nrow(snp_def) ==1) {
      pre_df[gene_index, 3] = "Two stages"
    }
  }
}

def_df = select(top_def, snp_id, feature_id)

for (gene_index in 1:nrow(top_def)){
  def_df[gene_index, 3] = "Stage specific"
  snp_iPSC = filter(full_iPSC, snp_id == top_def[gene_index, 1], feature_id == top_def[gene_index, 19])
  snp_mes = filter(full_mes, snp_id == top_def[gene_index, 1], feature_id == top_def[gene_index, 19])
  snp_pre = filter(full_pre, snp_id == top_def[gene_index, 1], feature_id == top_def[gene_index, 19])
  
  if (sign(top_def[gene_index, 3]) == sign(snp_iPSC$beta) && snp_iPSC$p_value < 0.05 && nrow(snp_iPSC) == 1){
    def_df[gene_index, 3] = "Two stages"
    if (sign(top_def[gene_index, 3]) == sign(snp_mes$beta) && snp_mes$p_value < 0.05 && nrow(snp_mes) == 1) {
      def_df[gene_index, 3] = "Three stages"
      if (sign(top_def[gene_index, 3]) == sign(snp_pre$beta) && snp_pre$p_value < 0.05 && nrow(snp_pre) == 1) {
        def_df[gene_index, 3] = "Four stages"
      }
    }
  }
  else {
    if (sign(top_def[gene_index, 3]) == sign(snp_mes$beta) && snp_mes$p_value < 0.05 && nrow(snp_mes) == 1){
      def_df[gene_index, 3] = "Two stages"
      if (sign(top_def[gene_index, 3]) == sign(snp_pre$beta) && snp_pre$p_value < 0.05 && nrow(snp_pre) ==1) {
        def_df[gene_index, 3] = "Three stages"
      }
    }
    else if (sign(top_def[gene_index, 3]) == sign(snp_pre$beta) && snp_pre$p_value < 0.05 && nrow(snp_pre) ==1) {
      def_df[gene_index, 3] = "Two stages"
    }
  }
}

# Write the files, to then do the barplots 

#write.csv(iPSC_df,"/omics/groups/OE0540/internal/users/ruehle/rnavelocity/vqtls/downstream/stage_specificity/det_velocity/top_iPSC_qtls.csv", row.names = FALSE)
#write.csv(mes_df,"/omics/groups/OE0540/internal/users/ruehle/rnavelocity/vqtls/downstream/stage_specificity/det_velocity/top_mes_qtls.csv", row.names = FALSE)
#write.csv(pre_df,"/omics/groups/OE0540/internal/users/ruehle/rnavelocity/vqtls/downstream/stage_specificity/det_velocity/top_pre_qtls.csv", row.names = FALSE)
#write.csv(def_df,"/omics/groups/OE0540/internal/users/ruehle/rnavelocity/vqtls/downstream/stage_specificity/det_velocity/top_def_qtls.csv", row.names = FALSE)


setwd("/omics/groups/OE0540/internal/users/ruehle/rnavelocity/vqtls/downstream/stage_specificity/unspliced")

iPSC_df = read.csv("./top_iPSC_qtls.csv")
mes_df = read.csv("./top_mes_qtls.csv")
pre_df = read.csv("./top_pre_qtls.csv")
def_df = read.csv("./top_def_qtls.csv")


iPSC_unspliced_counts = as.vector(table(iPSC_df$V3))
mes_unspliced_counts = as.vector(table(mes_df$V3))
pre_unspliced_counts = as.vector(table(pre_df$V3))
def_unspliced_counts = as.vector(table(def_df$V3))

rel_iPSC_expression_counts = iPSC_unspliced_counts/sum(iPSC_unspliced_counts)
rel_mes_expression_counts = mes_unspliced_counts/sum(mes_unspliced_counts)
rel_pre_expression_counts = pre_unspliced_counts/sum(pre_unspliced_counts)  
rel_def_expression_counts = def_unspliced_counts/sum(def_unspliced_counts)  

# create a dataset
type = c(rep("iPSC", 4), rep("Mesendoderm", 4), rep("Predefendoderm", 4), rep("Defendoderm", 4))
condition = rep(c("Four-stages" , "Stage-specific" , "Three-stages", "Two-stages"), 4)

#value = c(iPSC_expression_counts, mes_expression_counts, pre_expression_counts, def_expression_counts)
value = c(rel_iPSC_expression_counts, rel_mes_expression_counts, rel_pre_expression_counts, rel_def_expression_counts)
def_data <- data.frame(type,condition,value)
uns_data = data.frame(type, condition, value)

uns_data$condition = factor(def_data$condition, levels = c("Stage-specific" , "Two-stages","Three-stages",  "Four-stages"), ordered = TRUE)
uns_data$type = factor(def_data$type, levels = c("iPSC", "Mesendoderm", "Predefendoderm", "Defendoderm"), ordered = TRUE)
def_data$type = factor(def_data$type, levels = c("iPSC", "Mesendoderm", "Predefendoderm", "Defendoderm"), ordered = TRUE)

# Stacked
ggplot(uns_data, aes(fill=condition, y=value, x=type)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette="Reds") +  theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(title="Stage-specificity of u-eQTL",
       x ="Cell type", y = "Relative number of u-eQTL")


ggplot(ts[order(ts$y, decreasing = T),],
       aes(z, x, fill=factor(y, levels=c("blue","white" )))) + 
  geom_bar(stat = "identity")
