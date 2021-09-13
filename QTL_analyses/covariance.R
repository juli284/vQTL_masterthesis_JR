celltypes = c("iPSC", "mesendoderm", "predefendoderm", "defendoderm")

celltype = "defendoderm"

inpath = "/omics/groups/OE0540/internal/users/ruehle/rnavelocity/vqtls/data/velocyto/velocity/deterministic/celltypes/"

expMat = read.csv(paste0(inpath, celltype, "/input/pheno_file_det_velocity_", celltype, "_DE_genes.txt"), sep = "\t", header = TRUE, row.names = 1)

#t(apply(expMat, 1, function(x) {ifelse(is.na(x), min(x, na.rm = TRUE), x)}))
#try = t(apply(expMat, 1, function(x) {ifelse(is.na(x), min(x, na.rm = TRUE), x)}))
#expMat <- try[!is.infinite(rowSums(try)),]
colnames(expMat) = gsub("\\.", "-", colnames(expMat))

expMat_pcs <- prcomp(t(expMat))
expMat_pcs.df <- expMat_pcs$x
row.names(expMat_pcs.df) = gsub("\\.", "-", row.names(expMat_pcs.df))

write.table(expMat_pcs.df[,1:20], paste0(inpath, celltype, "/input/PCA_Covariates_20.txt") ,quote=FALSE,row.names=TRUE,sep="\t")
#write.table(expMat, paste0(inpath, celltype, "/input/pheno_file_logodds_velocity_", celltype, "_DE_genes.txt"),quote=FALSE,row.names=TRUE,sep="\t")


