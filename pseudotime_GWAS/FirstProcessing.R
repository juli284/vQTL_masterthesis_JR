library(qvalue)
setwd("/omics/groups/OE0540/internal/users/ruehle/rnavelocity/vqtls/data/velocyto/velocity/deterministic/celltypes/defendoderm")

VES = read.csv("/omics/groups/OE0540/internal/users/ruehle/rnavelocity/velocyto/scvelo_output/use/2000DEpercelltype_velocities_neighbors_pca/list_of_velocity_genes.txt", sep = "\t", header = TRUE)
VES = VES$Accession
Highly_variable_marc = read.csv("/omics/groups/OE0540/internal/users/ruehle/rnavelocity/vqtls/data/marcs/tested_genes.txt")
subset_marc_12k = read.delim("/omics/groups/OE0540/internal/users/ruehle/rnavelocity/vqtls/data/marcs/12Ksharedgenes.txt")


##Sup table 1 correlations
qtlResultFiles= c("./Mean/kinshipOnly_Neuro/", "./Mean/Neuro_baynorm_kinshipOnly/","./Mean/Neuro_scTransform_kinshipOnly/")
fdr = 0.05

##Initial pass. Full results.
allGenes <- NULL
qtls = NULL
#Gene selection pass.
folder = "./"

subset = VES

#Read top eQTLs only.

topQtlResults <- read.delim(paste(folder,"top_qtl_results_all.txt",sep=""),as.is=T)
topQtlResults["QTL"] <- paste(topQtlResults$snp_id,topQtlResults$feature_id,sep="-")
topQtlResults <- topQtlResults[order(topQtlResults$empirical_feature_p_value,topQtlResults$p_value),]               
#make sure that you don't have a feature more than once.
topQtlResults <- topQtlResults[which(!duplicated(topQtlResults$feature_id)),]

#here subset relevant number of genes 
#topQtlResults = topQtlResults[topQtlResults$feature_id %in% as.vector(subset),]

##Corrected for the number of genes tested for eQTL.
topQtlResults["global_corrected_pValue"] <- qvalue(topQtlResults$empirical_feature_p_value)$qvalues

#Read in all eQTL.
qtlResults <- read.delim(paste(folder,"qtl_results_all.txt",sep=""),as.is=T)
qtlResults["QTL"] <- paste(qtlResults$snp_id,qtlResults$feature_id,sep="-")
qtlResults <- qtlResults[order(qtlResults$empirical_feature_p_value,qtlResults$p_value),]
#make sure QTLs are not more than twice.
qtlResults <- qtlResults[which(!duplicated(qtlResults$QTL)),]


##Write out the relevant sections.
write.table(topQtlResults,paste(folder,"top_qtl_results_all_FDR.txt",sep=""), quote=F, sep="\t",row.names=F)
write.table(topQtlResults[which(topQtlResults$global_corrected_pValue<fdr),],paste(folder,"top_qtl_results_all_FDR",fdr,".txt",sep=""), quote=F, sep="\t",row.names=F)
print(nrow(topQtlResults[which(topQtlResults$global_corrected_pValue<fdr),]))
print(nrow(topQtlResults))

##Minimal P (That we find interesting)
minP = max(topQtlResults$empirical_feature_p_value[which(topQtlResults$global_corrected_pValue<fdr)])
#check what is the Pvalue of interest.
#filter out again subset of genes
#qtlResults <- qtlResults[which(qtlResults$feature_id %in% as.vector(subset)),]

print(nrow(qtlResults[which(qtlResults$empirical_feature_p_value<=minP),]))

#Write out the QTL results under a certain FDR level.

#write.table(qtlResults[which(qtlResults$empirical_feature_p_value<=minP),],paste(folder,"qtl_results_all_FDR",fdr,".txt",sep=""), quote=F, sep="\t",row.names=F)
dim(qtlResults[which(qtlResults$empirical_feature_p_value<=minP),])

               