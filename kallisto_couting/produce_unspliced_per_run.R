#wir machen erst mal nur die einfachen counts

library(tidyverse)
library(stringr)

merge_all_cells_from_run <- function(files) {
  #read in all unspliced/spliced files from one run 
  ldf <- lapply(files, read.table, sep = ",", header = TRUE)
  names(ldf) <- cellids
  gene_ids <- ldf[[1]]$GeneID
  
  #check if list of genes per cell identical 
  comp = list()
  comp <- lapply(ldf, function(df) {
    to_compare = df$GeneID
    comp = c(comp, identical(gene_ids, to_compare))
  })
  comp = as.vector(unlist(comp))
  print(all(comp))
  
  if (all(comp)) { #all(comp) checks if all elements in comp vector are TRUE 
    unspliced = sapply(ldf, function(df) {as.numeric(df$Counts)})
  }
  row.names(unspliced) = gene_ids
  colnames(unspliced) = as.vector(cellids)
  return(unspliced)
}
 


#read in all spliced and unspliced.csv files of one run 
runs <- dir("/omics/groups/OE0540/internal/users/ruehle/rnavelocity/pre_kallisto/hg38_release101/velocity_files_coll/counts", recursive=FALSE, full.names=TRUE, pattern="run")

#initialize final dfs
final_spliced = matrix(, nrow = 22793)
final_unspliced = matrix(, nrow = 22793)

lapply(runs, function(run) {
  print(run)
  filenames <- dir(run, recursive=TRUE, full.names=TRUE, pattern="spliced.csv")
  
  cellids <- str_extract(filenames, "....._.#.*/")
  cellids = str_replace(cellids, '/','')
  cellids <<- cellids[c(TRUE, FALSE)]
  
  #subset list of filenames to filenames_spliced and filenames_unspliced
  filenames_spliced = filenames[c(TRUE, FALSE)]
  print(paste0("spliced files read in ", filenames_spliced[1:6]))
  
  filenames_unspliced = filenames[c(FALSE, TRUE)]
  print(paste0("unspliced files read in ", filenames_unspliced[1:6]))
  
  #apply function to both unspliced and spliced 
  intron_counts = merge_all_cells_from_run(filenames_unspliced)
  exon_counts = merge_all_cells_from_run(filenames_spliced)
  
  #subsetting to genes that are shared among both matrices
  shared = intersect(rownames(intron_counts), rownames(exon_counts))
  intron_counts = intron_counts[shared,]
  exon_counts = exon_counts[shared,]
  print(length(shared))
  
  final_spliced <<- cbind(final_spliced, exon_counts)
  final_unspliced <<- cbind(final_unspliced, intron_counts)
})

write.csv(final_spliced, file="/omics/groups/OE0540/internal/users/ruehle/rnavelocity/pre_kallisto/hg38_release101/velocity_files_coll/counts/final_spliced_all_cells.csv", row.names = TRUE)
write.csv(final_unspliced, file="/omics/groups/OE0540/internal/users/ruehle/rnavelocity/pre_kallisto/hg38_release101/velocity_files_coll/counts/final_unspliced_all_cells.csv", row.names = TRUE)

