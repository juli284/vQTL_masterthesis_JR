library(dplyr)

setwd("/omics/groups/OE0540/internal/users/ruehle/rnavelocity/vqtls/data/velocyto/")

#run this seperately for every developmental stage (iPSC, mes, pre, def)

top_iPSC_eQTL = read.delim("./expression/celltypes/defendoderm/top_qtl_results_all_FDR0.05.txt", sep ="", as.is = T)
top_iPSC_ueQTL = read.delim("./expression/celltypes/defendoderm/top_qtl_results_all_FDR0.05.txt", sep="", as.is = T)
top_iPSC_vQTL = read.delim("./expression/celltypes/defendoderm/top_qtl_results_all_FDR0.05.txt", sep="", as.is = T)

full_iPSC_exp = read.delim("./expression/celltypes/defendoderm/qtl_results_all.txt", sep = "", as.is = T)
full_iPSC_un = read.delim("./unspliced/celltypes/defendoderm/qtl_results_all.txt")
full_iPSC_velo = read.delim("./velocity/deterministic/celltypes/defendoderm/qtl_results_all.txt", sep="", as.is = T)


eQTL_df = select(top_iPSC_eQTL, snp_id, feature_id)

#counting overlap ueQTL/vQTL with eQTL 
for (gene_index in 1:nrow(top_iPSC_eQTL)){
  eQTL_df[gene_index, 3] = "Type-specific"
  snp_un = filter(full_iPSC_un, snp_id == top_iPSC_eQTL[gene_index, 1], feature_id == top_iPSC_eQTL[gene_index, 19])
  snp_velo = filter(full_iPSC_velo, snp_id == top_iPSC_eQTL[gene_index, 1], feature_id == top_iPSC_eQTL[gene_index, 19])
  
  if (sign(top_iPSC_eQTL[gene_index, 3]) == sign(snp_un$beta) && snp_un$p_value < 0.05 && nrow(snp_un) == 1){
    eQTL_df[gene_index, 3] = "Two types"
    eQTL_df[gene_index,4] = "eQTL & u-eQTL"
    if (sign(top_iPSC_eQTL[gene_index, 3]) == sign(snp_velo$beta) && snp_velo$p_value < 0.05 && nrow(snp_velo) == 1) {
      eQTL_df[gene_index, 3] = "Three types"
    }
  }
  else if (sign(top_iPSC_eQTL[gene_index, 3]) == sign(snp_velo$beta) && snp_velo$p_value < 0.05 && nrow(snp_velo) ==1) {
    eQTL_df[gene_index, 3] = "Two types"
    eQTL_df[gene_index, 4] = "eQTL & vQTL"
  }
}


#counting overlap eQTL/vQTL with u-eQTL

ueQTL_df = select(top_iPSC_ueQTL, snp_id, feature_id)

for (gene_index in 1:nrow(top_iPSC_ueQTL)){
  ueQTL_df[gene_index, 3] = "Type-specific"
  snp_exp = filter(full_iPSC_exp, snp_id == top_iPSC_ueQTL[gene_index, 1], feature_id == top_iPSC_ueQTL[gene_index, 19])
  snp_velo = filter(full_iPSC_velo, snp_id == top_iPSC_ueQTL[gene_index, 1], feature_id == top_iPSC_ueQTL[gene_index, 19])
  
  if (sign(top_iPSC_ueQTL[gene_index, 3]) == sign(snp_velo$beta) && snp_velo$p_value < 0.05 && nrow(snp_velo) == 1){
    ueQTL_df[gene_index, 3] = "Two types"
    ueQTL_df[gene_index, 4] = "u-eQTL & vQTL"
    if (sign(top_iPSC_ueQTL[gene_index, 3]) == sign(snp_exp$beta) && snp_exp$p_value < 0.05 && nrow(snp_exp) == 1) {
      ueQTL_df[gene_index, 3] = "Three types"
    }
  }
  else if (sign(top_iPSC_ueQTL[gene_index, 3]) == sign(snp_exp$beta) && snp_exp$p_value < 0.05 && nrow(snp_exp) ==1) {
    ueQTL_df[gene_index, 3] = "Two types"
    ueQTL_df[gene_index, 4] = "u-eQTL & eQTL"
  }
}

# counting overlap eQTL/ueQTL with vQTL 

veQTL_df = select(top_iPSC_vQTL, snp_id, feature_id)

for (gene_index in 1:nrow(top_iPSC_vQTL)){
  veQTL_df[gene_index, 3] = "Type-specific"
  snp_exp = filter(full_iPSC_exp, snp_id == top_iPSC_vQTL[gene_index, 1], feature_id == top_iPSC_vQTL[gene_index, 19])
  snp_un = filter(full_iPSC_un, snp_id == top_iPSC_vQTL[gene_index, 1], feature_id == top_iPSC_vQTL[gene_index, 19])
  
  if (sign(top_iPSC_vQTL[gene_index, 3]) == sign(snp_un$beta) && snp_un$p_value < 0.05 && nrow(snp_un) == 1){
    veQTL_df[gene_index, 3] = "Two types"
    veQTL_df[gene_index, 4] = "vQTL & u-eQTL"
    if (sign(top_iPSC_vQTL[gene_index, 3]) == sign(snp_exp$beta) && snp_exp$p_value < 0.05 && nrow(snp_exp) == 1) {
      veQTL_df[gene_index, 3] = "Three types"
    }
  }
  else if (sign(top_iPSC_vQTL[gene_index, 3]) == sign(snp_exp$beta) && snp_exp$p_value < 0.05 && nrow(snp_exp) ==1) {
    veQTL_df[gene_index, 3] = "Two types"
    veQTL_df[gene_index,4] = "vQTL & eQTL"
  }
}


# Hier erst mal Stopp 
#write.csv(eQTL_df,"/omics/groups/OE0540/internal/users/ruehle/rnavelocity/vqtls/downstream/type_specificity/def/top_exp_qtls.csv", row.names = FALSE)
#write.csv(ueQTL_df,"/omics/groups/OE0540/internal/users/ruehle/rnavelocity/vqtls/downstream/type_specificity/def/top_un_qtls.csv", row.names = FALSE)
#write.csv(vQTL_df,"/omics/groups/OE0540/internal/users/ruehle/rnavelocity/vqtls/downstream/type_specificity/def/top_velo_qtls.csv", row.names = FALSE)




setwd("/omics/groups/OE0540/internal/users/ruehle/rnavelocity/vqtls/downstream/type_specificity/iPSC")

eQTL_df = read.csv("./top_exp_qtls.csv")
ueQTL_df = read.csv("./top_unspliced_qtls.csv")
veQTL_df = read.csv("./top_velo_qtls.csv")

iPSC_unspliced_counts = as.vector(table(eQTL_df$V3))
mes_unspliced_counts = as.vector(table(ueQTL_df$V3))
pre_unspliced_counts = as.vector(table(veQTL_df$V3))

rel_iPSC_expression_counts = iPSC_unspliced_counts/sum(iPSC_unspliced_counts)
rel_mes_expression_counts = mes_unspliced_counts/sum(mes_unspliced_counts)
rel_pre_expression_counts = pre_unspliced_counts/sum(pre_unspliced_counts)  

# create a dataset
type = c(rep("eQTL", 3), rep("u-eQTL", 3), rep("vQTL", 3))
condition = rep(c("Three types", "Two types", "Type-specific"), 3)

#value = c(iPSC_expression_counts, mes_expression_counts, pre_expression_counts)
value = c(rel_iPSC_expression_counts, rel_mes_expression_counts, rel_pre_expression_counts)
def_data <- data.frame(type,condition,value)
uns_data = data.frame(type, condition, value)

uns_data$condition = factor(def_data$condition, levels = c("Type-specific", "Two types", "Three types"), ordered = TRUE)
uns_data$type = factor(def_data$type, levels = c("eQTL", "u-eQTL", "vQTL"), ordered = TRUE)
def_data$type = factor(def_data$type, levels = c("eQTL", "u-eQTL", "vQTL"), ordered = TRUE)

# Stacked
ggplot(uns_data, aes(fill=condition, y=value, x=type)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette="Reds") +  theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(title="Type-specificity of QTL",
       x ="iPSC", y = "Relative number of QTL")


ggplot(ts[order(ts$y, decreasing = T),],
       aes(z, x, fill=factor(y, levels=c("blue","white" )))) + 
  geom_bar(stat = "identity")