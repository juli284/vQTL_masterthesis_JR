library(UpSetR)
library(ggplot2)

DES = read.csv("/omics/groups/OE0540/internal/users/ruehle/rnavelocity/velocyto/scvelo_output/use/2000DEpercelltype_velocities_neighbors_pca/list_of_DEs_per_celltype.txt", sep = "\t", header = TRUE)
VES = read.csv("/omics/groups/OE0540/internal/users/ruehle/rnavelocity/velocyto/scvelo_output/use/2000DEpercelltype_velocities_neighbors_pca/list_of_velocity_genes.txt", sep = "\t", header = TRUE)

setwd("/omics/groups/OE0540/internal/users/ruehle/rnavelocity/vqtls/data/velocyto/")

e_iPSC = read.delim("./expression/celltypes/iPSC/top_qtl_results_all_FDR0.05.txt", sep="", as.is = T)
e_mes = read.delim("./expression/celltypes/mesendoderm/top_qtl_results_all_FDR0.05.txt", sep="", as.is = T)
e_pre = read.delim("./expression/celltypes/predefendoderm/top_qtl_results_all_FDR0.05.txt", sep="", as.is = T)
e_def = read.delim("./expression/celltypes/defendoderm/top_qtl_results_all_FDR0.05.txt", sep="", as.is = T)

u_iPSC = read.delim("./unspliced/celltypes/iPSC/top_qtl_results_all_FDR0.05.txt", sep="", as.is = T)
u_mes = read.delim("./unspliced/celltypes/mesendoderm/top_qtl_results_all_FDR0.05.txt", sep="", as.is = T)
u_pre = read.delim("./unspliced/celltypes/predefendoderm/top_qtl_results_all_FDR0.05.txt", sep="", as.is = T)
u_def = read.delim("./unspliced/celltypes/defendoderm/top_qtl_results_all_FDR0.05.txt", sep="", as.is = T)

v_iPSC = read.delim("./velocity/deterministic/celltypes/iPSC/top_qtl_results_all_FDR0.05.txt", sep="", as.is = T)
v_mes = read.delim("./velocity/deterministic/celltypes/mesendoderm/top_qtl_results_all_FDR0.05.txt", sep="", as.is = T)
v_pre = read.delim("./velocity/deterministic/celltypes/predefendoderm/top_qtl_results_all_FDR0.05.txt", sep="", as.is = T)
v_def = read.delim("./velocity/deterministic/celltypes/defendoderm/top_qtl_results_all_FDR0.05.txt", sep="", as.is = T)

par(mfrow = c(3,1))

listInput <- list(one = c(1, 2, 3, 5, 7, 8, 11, 12, 13), two = c(1, 2, 4, 5, 10), three = c(1, 5, 6, 7, 8, 9, 10, 12, 13))

e_genes_listi = list(eGenes_iPSC = as.vector(e_iPSC$feature_id), eGenes_mes = as.vector(e_mes$feature_id), eGenes_pre = as.vector(e_pre$feature_id), eGenes_def = as.vector(e_def$feature_id))
v_genes_listi = list(vGenes_iPSC = as.vector(v_iPSC$feature_id), vGenes_mes = as.vector(v_mes$feature_id), vGenes_pre = as.vector(v_pre$feature_id), vGenes_def = as.vector(v_def$feature_id))
u_genes_listi = list(uGenes_iPSC = as.vector(u_iPSC$feature_id), uGenes_mes = as.vector(u_mes$feature_id), uGenes_pre = as.vector(u_pre$feature_id), uGenes_def = as.vector(u_def$feature_id))

e_and_u_listi = list(eGenes_iPSC = as.vector(e_iPSC$feature_id), eGenes_mes = as.vector(e_mes$feature_id), eGenes_pre = as.vector(e_pre$feature_id), eGenes_def = as.vector(e_def$feature_id), uGenes_iPSC = as.vector(u_iPSC$feature_id), uGenes_mes = as.vector(u_mes$feature_id), uGenes_pre = as.vector(u_pre$feature_id), uGenes_def = as.vector(u_def$feature_id))

listi = list(eGenes_iPSC = as.vector(e_iPSC$feature_id), eGenes_mes = as.vector(e_mes$feature_id), vGenes_iPSC = as.vector(v_pre$feature_id), vGenes_mes = as.vector(v_mes$feature_id))

v_genes_listi = list(eGenes_def = as.vector(e_def$feature_id), vGenes_iPSC = as.vector(v_iPSC$feature_id), vGenes_mes = as.vector(v_mes$feature_id), vGenes_pre = as.vector(v_pre$feature_id), vGenes_def = as.vector(v_def$feature_id))


upset(fromList(e_genes_listi), order.by = "freq")
upset(fromList(v_genes_listi), order.by = "freq")
upset(fromList(u_genes_listi), order.by = "freq", group.by = "set")
upset(fromList(e_and_u_listi), order.by = "freq")


library(ggplot2)

# create a dataset
genes <- c(e_iPSC$feature_id, e_mes$feature_id, u_iPSC$feature_id, u_mes$feature_id)
stage <- c(rep("iPSC", length(e_iPSC$feature_id)), rep("mesendoderm", length(e_mes$feature_id)), rep("iPSC", length(u_iPSC$feature_id)), rep("mesendoderm", length(u_mes$feature_id)))
QTL <- c(rep("expression", length(e_iPSC$feature_id)), rep("expression", length(e_mes$feature_id)), rep("unspliced", length(u_iPSC$feature_id)), rep("unspliced", length(u_mes$feature_id)))

specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
value <- c(rep(1, length(genes)))
data <- data.frame(stage,QTL,value)



# Grouped

p <- data %>%
  ggplot( aes(x=value, fill=QTL, x = stage)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_ipsum() +
  labs(fill="")

ggplot(data, aes(fill=QTL, y=values, x=stage)) + 
  geom_bar(position="dodge", stat="identity")


df <- data.frame(QTL=c("expression", "unspliced", "velocity", "expression", "unspliced", "velocity", "expression", "unspliced", "velocity", "expression", "unspliced", "velocity"),
                 counts=c(1474, 790, 770, 1237, 475, 199, 984, 616, 367, 947, 651, 375),
                 stage = c("iPSC", "iPSC", "iPSC", "mesendoderm", "mesendoderm", "mesendoderm", "predefendoderm","predefendoderm", "predefendoderm", "defendoderm", "defendoderm", "defendoderm"))

df$stage <- factor(df$stage,                                    # Change ordering manually
                  levels = c("iPSC", "mesendoderm", "predefendoderm", "defendoderm"))

ggplot(data=df, aes(fill=qtl, y=counts, x = stage)) +
geom_bar(stat="identity", position = "dodge") +
scale_fill_manual(values=c("#ABE1A8", "#EAA084", "#898FF1")) +
theme_classic() + theme(axis.text.x = element_text(angle = 15, vjust = 1.1, hjust=1), axis.text=element_text(size=10),
                                                                                       axis.title=element_text(size=12,face="bold"),
                                                                                       legend.title = element_text(size=11), 
                                                                                       legend.text = element_text(size=10), 
                                                                                       legend.key.size = unit(0.8, 'cm')) +
labs(title="", x ="Cell type", y = "# genes associated with QTL") +
guides(fill=guide_legend(title="Type of QTL")) 


p + theme(legend.text = element_text(colour="blue", size=10, 
                                     face="bold"))
theme(legend.title = element_text(size=30))

library(VennDiagram)

# Generate 3 sets of 200 words
eGenes = c(e_iPSC$feature_id, e_mes$feature_id, e_pre$feature_id, e_def$feature_id)
uGenes = c(u_iPSC$feature_id, u_mes$feature_id, u_pre$feature_id, u_def$feature_id)
vGenes = c(v_iPSC$feature_id, v_mes$feature_id, v_pre$feature_id, v_def$feature_id)

# Chart
x = list(eGenes, uGenes, vGenes)

venn.diagram(x,
  category.names = c("eGenes" , "uGenes" , "vGenes"),
  filename = './14_venn_diagramm.png',
  output=TRUE
)

# Helper function to display Venn diagram
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}


# Change fill color
display_venn(
  x,
  category.names = c("eGenes" , "uGenes" , "vGenes"),
  fill = c("#ABE1A8", "#EAA084", "#898FF1")
)

#ABE1A8", "#EAA084", "#898FF1"

bla = e_iPSC$feature_id %in% e_mes$feature_id 
neu = e_iPSC$feature_id
neu = neu[bla]

blub = neu %in% e_pre$feature_id
neu = neu[blub]

bli = neu %in% e_def$feature_id
neu = neu[bli]

sum(bla, na.rm = TRUE)


eGenes_iPSC  =  e_mes$feature_id
genes_df = data.frame(row.names = eGenes_iPSC)

for (gene_index in 1:length(eGenes_iPSC)){
  if (eGenes_iPSC[gene_index] %in% v_mes$feature_id) {
    genes_df[gene_index, 1] = "2"
    if (eGenes_iPSC[gene_index] %in% v_pre$feature_id) {
      genes_df[gene_index, 1] = "3"
      if (eGenes_iPSC[gene_index] %in% v_def$feature_id) {
        genes_df[gene_index,1] = "4" 
      }
    }
    else if (eGenes_iPSC[gene_index] %in% v_def$feature_id) {
      genes_df[gene_index, 1] = "3"
    }
  }
  else if (eGenes_iPSC[gene_index] %in% v_pre$feature_id) {
    genes_df[gene_index,1] = "2"
    if (eGenes_iPSC[gene_index] %in% v_def$feature_id) {
    genes_df[gene_index,1] = "3"
    }
  }
  else if (eGenes_iPSC[gene_index] %in% v_def$feature_id) {
    genes_df[gene_index,1] = "2"
  }
  else {
  genes_df[gene_index, 1] = "1"
  }
}

unspliced_counts = as.vector(table(genes_df$V1))
rel_unspliced_counts = unspliced_counts/length(u_iPSC$feature_id)

expression_counts = as.vector(table(genes_df$V1))
rel_spliced_counts = expression_counts/length(e_iPSC$feature_id)

velocity_counts = as.vector(table(genes_df$V1))
rel_velocity_counts = velocity_counts/length(v_iPSC$feature_id)

eGenes_pre  =  v_def$feature_id
genes_df = data.frame(row.names = eGenes_pre)

for (gene_index in 1:length(eGenes_pre)){
  if (eGenes_pre[gene_index] %in% v_iPSC$feature_id) {
    genes_df[gene_index, 1] = "2"
    if (eGenes_pre[gene_index] %in% v_mes$feature_id) {
      genes_df[gene_index, 1] = "3"
      if (eGenes_pre[gene_index] %in% v_pre$feature_id) {
        genes_df[gene_index,1] = "4" 
      }
    }
    else if (eGenes_pre[gene_index] %in% v_pre$feature_id) {
      genes_df[gene_index, 1] = "3"
    }
  }
  else if (eGenes_pre[gene_index] %in% v_mes$feature_id) {
    genes_df[gene_index,1] = "2"
    if (eGenes_pre[gene_index] %in% v_pre$feature_id) {
      genes_df[gene_index,1] = "3"
    }
  }
  else if (eGenes_pre[gene_index] %in% v_pre$feature_id) {
    genes_df[gene_index,1] = "2"
  }
  else {
    genes_df[gene_index, 1] = "1"
  }
}

unspliced_counts = as.vector(table(genes_df$V1))
rel_unspliced_counts = unspliced_counts/length(u_def$feature_id)

expression_counts = as.vector(table(genes_df$V1))
rel_spliced_counts = expression_counts/length(e_def$feature_id)

velocity_counts = as.vector(table(genes_df$V1))
rel_velocity_counts = velocity_counts/length(v_def$feature_id)

# create a dataset
type <- c(rep("expression" , 4) , rep("unspliced" , 4) , rep("velocity" , 4) )
condition <- rep(c("stage-specific" , "two-stage" , "three-stage", "four-stage") , 3)
stage = c(rep("defendoderm", 12))

value <- c(rel_spliced_counts, rel_unspliced_counts, rel_velocity_counts)
def_data <- data.frame(type,condition,value, stage)

par(mfrow = c(2,1))
# Stacked
ggplot(def_data, aes(fill=condition, y=value, x=type)) + 
  geom_bar(position="stack", stat="identity")

