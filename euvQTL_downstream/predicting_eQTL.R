library(dplyr)

setwd("/omics/groups/OE0540/internal/users/ruehle/rnavelocity/vqtls/data/velocyto/")

#here we check, if vQTL occur as eQTL in later developmental stages 
#for vQTL in iPSC, we check if eQTL in mes, pre and/or def
#for vQTL in mes, we check if eQTL in pre and/or def 
#for vQTL in pre, we check if eQTL in def 

top_iPSC = read.delim("./velocity/deterministic/celltypes/iPSC/top_qtl_results_all_FDR0.05.txt", sep ="", as.is = T)
top_mes = read.delim("./velocity/deterministic/celltypes/mesendoderm/top_qtl_results_all_FDR0.05.txt", sep="", as.is = T)
top_pre = read.delim("./velocity/deterministic/celltypes/predefendoderm/top_qtl_results_all_FDR0.05.txt", sep="", as.is = T)
top_def = read.delim("./velocity/deterministic/celltypes/defendoderm/top_qtl_results_all_FDR0.05.txt", sep="", as.is = T)

full_iPSC = read.delim("./expression/celltypes/iPSC/qtl_results_all.txt", sep="", as.is = T)
full_mes = read.delim("./expression/celltypes/mesendoderm/qtl_results_all.txt", sep="", as.is = T)
full_pre = read.delim("./expression/celltypes/predefendoderm/qtl_results_all.txt", sep="", as.is = T)
full_def = read.delim("./expression/celltypes/defendoderm/qtl_results_all.txt", sep="", as.is = T)

#counting iPSC vQTL
genes_df = dplyr::select(top_iPSC, snp_id, feature_id)

for (gene_index in 1:nrow(top_iPSC)){
  snp_mes = filter(full_mes, snp_id == top_iPSC[gene_index, 1], feature_id == top_iPSC[gene_index, 19])
  snp_pre = filter(full_pre, snp_id == top_iPSC[gene_index, 1], feature_id == top_iPSC[gene_index, 19])
  snp_def = filter(full_def, snp_id == top_iPSC[gene_index, 1], feature_id == top_iPSC[gene_index, 19])
  
  if (sign(top_iPSC[gene_index, 3]) == sign(snp_mes$beta) && snp_mes$p_value < 0.05 && nrow(snp_mes) == 1){
    genes_df[gene_index, 3] = "mes-eQTL"
    if (sign(top_iPSC[gene_index, 3]) == sign(snp_pre$beta) && snp_pre$p_value < 0.05 && nrow(snp_pre) == 1) {
      genes_df[gene_index, 4] = "pre-eQTL"
      if (sign(top_iPSC[gene_index, 3]) == sign(snp_def$beta) && snp_def$p_value < 0.05 && nrow(snp_def) == 1) {
        genes_df[gene_index, 5] = "def-eQTL"
      }
    }
  }
  else {
    if (sign(top_iPSC[gene_index, 3]) == sign(snp_pre$beta) && snp_pre$p_value < 0.05 && nrow(snp_pre) == 1){
      genes_df[gene_index, 4] = "pre-eQTL"
      if (sign(top_iPSC[gene_index, 3]) == sign(snp_def$beta) && snp_def$p_value < 0.05 && nrow(snp_def) ==1) {
        genes_df[gene_index, 5] = "def-eQTL"
      }
    }
    else if (sign(top_iPSC[gene_index, 3]) == sign(snp_def$beta) && snp_def$p_value < 0.05 && nrow(snp_def) ==1) {
      genes_df[gene_index, 5] = "def-eQTL"
    }
  }
}


#counting mesendoderm vQTL
mes_df = dplyr::select(top_mes, snp_id, feature_id)

for (gene_index in 1:nrow(top_mes)){
  #mes_df[gene_index, 3] = "no eQTL"
  
  snp_pre = filter(full_pre, snp_id == top_mes[gene_index, 1], feature_id == top_mes[gene_index, 19])
  snp_def = filter(full_def, snp_id == top_mes[gene_index, 1], feature_id == top_mes[gene_index, 19])
  
  if (sign(top_mes[gene_index, 3]) == sign(snp_pre$beta) && snp_pre$p_value < 0.05 && nrow(snp_pre) == 1){
    mes_df[gene_index, 3] = "pre-eQTL"
    if (sign(top_mes[gene_index, 3]) == sign(snp_def$beta) && snp_def$p_value < 0.05 && nrow(snp_def) == 1) {
      mes_df[gene_index, 4] = "def-eQTL"
    }
  }
  else {
    if (sign(top_mes[gene_index, 3]) == sign(snp_def$beta) && snp_def$p_value < 0.05 && nrow(snp_def) == 1){
      mes_df[gene_index, 4] = "def-eQTL"
    }
  }
}

#predef counting

pre_df = dplyr::select(top_pre, snp_id, feature_id)

for (gene_index in 1:nrow(top_pre)){
  #pre_df[gene_index, 3] = "Stage spec
  snp_def = filter(full_def, snp_id == top_pre[gene_index, 1], feature_id == top_pre[gene_index, 19])
  
  if (sign(top_pre[gene_index, 3]) == sign(snp_def$beta) && snp_def$p_value < 0.05 && nrow(snp_def) == 1){
    pre_df[gene_index, 3] = "def-eQTL"
  }
}

#check if iPSC-vQTL also iPSC eQTL, mes vQTL also mes eQTL,...

for (gene_index in 1:nrow(top_iPSC)){
  #genes_df[gene_index, 3] = "no eQTL"
  snp_iPSC = filter(full_iPSC, snp_id == top_iPSC[gene_index, 1], feature_id == top_iPSC[gene_index, 19])
  if (sign(top_iPSC[gene_index, 3]) == sign(snp_iPSC$beta) && snp_iPSC$p_value < 0.05 && nrow(snp_iPSC) == 1){
    genes_df[gene_index, 6] = "iPSC-eQTL"}
}

for (gene_index in 1:nrow(top_mes)){
  #genes_df[gene_index, 3] = "no eQTL"
  snp_mes = filter(full_mes, snp_id == top_mes[gene_index, 1], feature_id == top_mes[gene_index, 19])
  if (sign(top_mes[gene_index, 3]) == sign(snp_mes$beta) && snp_mes$p_value < 0.05 && nrow(snp_mes) == 1){
    mes_df[gene_index, 5] = "mes-eQTL"}
}

for (gene_index in 1:nrow(top_pre)){
  #genes_df[gene_index, 3] = "no eQTL"
  snp_pre = filter(full_pre, snp_id == top_pre[gene_index, 1], feature_id == top_pre[gene_index, 19])
  if (sign(top_pre[gene_index, 3]) == sign(snp_pre$beta) && snp_pre$p_value < 0.05 && nrow(snp_pre) == 1){
    pre_df[gene_index, 4] = "pre-eQTL"}
}


#counting later stages 
iPSC_not_later = genes_df[(is.na(genes_df$V3) & is.na(genes_df$V4) & is.na(genes_df$V5)), ]
number_iPSC_later = nrow(genes_df) - nrow(iPSC_not_later)

#counting later stages 
mes_not_later = mes_df[(is.na(mes_df$V3) & is.na(mes_df$V4)), ]
number_mes_later = nrow(mes_df) - nrow(mes_not_later)

pre_not_later = pre_df[is.na(pre_df$V3), ]
number_pre_later = nrow(pre_df) - nrow(pre_not_later)


not_eQTL_and_vQTL_iPSC = nrow(genes_df[is.na(genes_df$V6),])
not_eQTL_and_vQTL_mes = nrow(mes_df[is.na(mes_df$V5),])
not_eQTL_and_vQTL_pre = nrow(pre_df[is.na(pre_df$V4),])

eQTL_and_vQTL_iPSC = nrow(genes_df) - not_eQTL_and_vQTL_iPSC
eQTL_and_vQTL_mes = nrow(mes_df) - not_eQTL_and_vQTL_mes
eQTL_and_vQTL_pre = nrow(pre_df) - not_eQTL_and_vQTL_pre




