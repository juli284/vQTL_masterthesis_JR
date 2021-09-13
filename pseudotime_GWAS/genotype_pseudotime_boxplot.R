## read in data 

sample_mapping <- read.delim("/omics/groups/OE0540/internal/users/ruehle/rnavelocity/vqtls/data/velocyto/pseudotime/harmony/input_day2/sample_mapping_file_pseudotime_day2.txt", header = F)
exprs <-read.delim("/omics/groups/OE0540/internal/users/ruehle/rnavelocity/vqtls/data/velocyto/pseudotime/harmony/input_day2/pheno_file_pseudotime_day2.txt",as.is=T, row.names=1,check.names = F)
cov <-read.delim("/omics/groups/OE0540/internal/users/ruehle/rnavelocity/vqtls/data/velocyto/pseudotime/harmony/input_day2/covariance_pseudotime_day2.txt",as.is=T, row.names=1,check.names = F)
geno <-read.delim("/omics/groups/OE0540/internal/users/ruehle/rnavelocity/genotypes/genome_harmonizer/pseudo_snps.genotypes.txt",as.is=T, row.names=1)
topQTLs <- read.delim("/omics/groups/OE0540/internal/users/ruehle/rnavelocity/vqtls/data/velocyto/pseudotime/harmony/output/day2/qtl_results_all.txt")

mes = topQTLs$snp_id[1:6]
SNPOI = mes[1]
GOI = "pc_pseudotime"

#prepare geno to match sample mapping
colnames(geno) = gsub("\\.", "-", colnames(geno))
geno = geno[c(SNPOI), ]
mat = match(sample_mapping$V1, colnames(geno))

# ## for the validity of the test I drop duplicated donors.
sample_mapping = sample_mapping[which(!duplicated(sample_mapping[,1])),]

#make sure only relevant entries are in the sample mapping file.
#sample_mapping = sample_mapping[-which(sample_mapping[,1]=="zero"),]##Not needed no matches to zero.
sample_mapping = sample_mapping[which(sample_mapping[,1] %in% colnames(geno)),]
sample_mapping = sample_mapping[which(sample_mapping[,2] %in% colnames(exprs)),]

#Filter the actual matrices.
exprs = exprs[,which(colnames(exprs) %in% sample_mapping[,2])]
exprs = exprs[,order(colnames(exprs))]
cov = cov[which(rownames(cov) %in% sample_mapping[,2]),]
cov = cov[order(rownames(cov)),]
all(rownames(cov)==colnames(exprs))

exprsCorrected = residuals(lm(t(exprs)~as.matrix((cov))))

geno = geno[,which(colnames(geno) %in% sample_mapping[,1])]
geno2 = geno[,match(sample_mapping[,1], colnames(geno))]
colnames(geno2) = sample_mapping[,2]
geno2 = geno2[,order(colnames(geno2)),]
all(colnames(exprs)  == colnames(geno2))


boxplot(as.numeric(exprs[1,])~as.factor(geno2[1,]))
boxplot(as.numeric(exprsCorrected)~as.factor(geno2[1,]))

final_mat = data.frame("pseudotime"= as.numeric(exprsCorrected), "variant" = as.factor(geno2[1,]))

#c("#354f5f", "#627682", "#8f9da6")
#c("#999999", "#E69F00", "#56B4E9")

ggplot(final_mat, aes(x=variant, y=pseudotime, fill = variant)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill = "white") + scale_fill_manual(values=c( "#d2d8db", "#8f9da6","#354f5f")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(title="GWAS on pseudotime", x =SNPOI, y = "Pseudotime") #+ scale_y_continuous(breaks=c(0.4,0.5, 0.6, 0.7))
