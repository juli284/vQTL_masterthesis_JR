#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

library(stringr)
library(tictoc)
library(tidyverse)

#read in annotation file 
anno = read.table(file = "/omics/groups/OE0540/internal/users/ruehle/rnavelocity/pre_kallisto/hg38_release101/velocity_files_coll/anno.csv", sep = ",", header = TRUE)

  abundance <- read.table(file = args[1], header = TRUE, sep = '\t')
  #separate intronic and  exonic reads in separate dfs 
  quantintrons <- function(abundance) {
    tic("sleeping")
    cols = colnames(abundance)
    rownames(abundance) = abundance$target_id
    abundance$spliced = grepl("I.*", abundance$target_id)
    blub <- split(abundance, abundance$spliced)
    
    unspliced = data.frame(blub [2])
    
    unspliced = unspliced [,-6]
    colnames(unspliced) = cols
    unspliced$target_id = gsub("[.].*", "", unspliced$target_id)
    
    #counts in unspliced 
    rowlist = list()
    unspliced_counts = data.frame()
    effgeneLength = 0
    tpm = 0
    Countspergene = 0
    geneLength = 0
    
    dfpergene <- split(unspliced, unspliced$target_id)
    
    lapply(dfpergene, function(pergene) {
      geneCount = 0
      pergene = as.data.frame(pergene, col.names = cols)
      for (i in 1:nrow(pergene)) {
        geneCount = geneCount + (pergene$est_counts[i]/pergene$eff_length[i])  
      }
      geneID = pergene[1,1] #ich glaub des geht so net
      geneLength = mean(pergene$length)
      effgeneLength = mean(pergene$eff_length)
      tpm = geneCount/1000000
      Countspergene = geneCount*effgeneLength
      
      rowlist = c(geneID, geneLength, effgeneLength, Countspergene, tpm)
      unspliced_counts <<- rbind(unspliced_counts, rowlist)
    })
    colnames(unspliced_counts) = c("GeneID","AverageLength", "EffectiveLength", "Counts", "TPM")
    unspliced_counts = data.frame(unspliced_counts)
    print("...waking up")
    toc()
    return(unspliced_counts)
  }  
  
  
  quantexons <- function(abundance) {
    #prepare spliced dataframe from whole abundance df
    tic("sleeping")
    
    cols = colnames(abundance)
    rownames(abundance) = abundance$target_id
    abundance$spliced = grepl("I.*", abundance$target_id)
    blub <- split(abundance, abundance$spliced)
    
    spliced = data.frame(blub[1])
    
    spliced = spliced[,-6]
    colnames(spliced) = cols
    spliced$target_id = gsub("[.].*", "", spliced$target_id)
    
    #count in spliced   
    rowlist = list()
    spliced_counts = data.frame()
    tpm = 0
    Countspergene = 0
    geneLength = 0
    
    dfpergene <- split(anno, anno$GENEID)
    
    lapply(dfpergene, function(pergene) {
      geneCount = 0
      new = data.frame()
      tx.list = vector()
      
      pergene = as.data.frame(pergene)
      colnames(pergene) = c("ENTREZID", "GENEID", "GENENAME", "TXID")
      
      #select transcripts belonging to one gene in spliced
      tx.list = as.vector(pergene$TXID)
      new = spliced %>% filter(
        target_id %in% tx.list)
      
      for (i in 1:nrow(new)) {
        geneCount = geneCount + (new$est_counts[i]/new$eff_length[i])  
      }
      geneID = pergene[1,2]
      geneLength = mean(new$length)
      effgeneLength = mean(new$eff_length)
      tpm = geneCount/1000000
      Countspergene = geneCount*effgeneLength
      
      rowlist = c(geneID, geneLength, effgeneLength, Countspergene, tpm)
      spliced_counts <<- rbind(spliced_counts, rowlist)
    })
    colnames(spliced_counts) = c("GeneID","AverageLength", "EffectiveLength", "Counts", "TPM")
    spliced_counts = data.frame(spliced_counts)
    print("...waking up")
    toc()
    return(spliced_counts)
  }  

outfile = args[2]  
outpath = gsub("/spliced.csv", "", outfile)
  
introncount <- quantintrons(abundance)
exoncount <- quantexons(abundance)
#outpath = paste0("/omics/groups/OE0540/internal/users/ruehle/rnavelocity/pre_kallisto/hg38_release101/velocity_files_coll/counts/run21999/", args[2])
write.csv(introncount, file=paste0(outpath,"/unspliced.csv"), row.names = FALSE)
print("wrote unspliced file into")
print(outpath)
write.csv(exoncount, file=paste0(outpath, "/spliced.csv"), row.names = FALSE)
  


