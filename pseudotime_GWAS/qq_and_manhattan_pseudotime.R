library(qqman)
library(dplyr)
library(ggplot2)
library(plotly)
library(tidyverse) 
library(RColorBrewer) 
library(ggrepel) 
library(kableExtra)

setwd("/omics/groups/OE0540/internal/users/ruehle/rnavelocity/vqtls/data/velocyto/pseudotime/harmony/output/")

FDR = 0.05

day0 = read.table("./day0/qtl_results_all.txt", sep="\t", header=TRUE)
day0['timepoint'] = "day0"
day1 = read.table("./day1/qtl_results_all.txt", sep="\t", header=TRUE)
day1['timepoint'] = "day1"
day2 = read.table("./day2/qtl_results_all.txt", sep="\t", header=TRUE)
day2['timepoint'] = "day2"
day3 = read.table("./day3/qtl_results_all.txt", sep="\t", header=TRUE)
day3['timepoint'] = "day3"


day0_filtered = day0 %>% 
  filter(-log10(p_value)>1) #filtert alle raus die kleiner sind als 0.001

day1_filtered = day1 %>% 
  filter(-log10(p_value)>1)

day2_filtered = day2 %>% 
  filter(-log10(p_value)>3)

day3_filtered = day3 %>%
  filter(-log10(p_value)>1)

#get an overview per day 

par(mfrow=c(2,4))
#manhattan(gwasResults, chr = "CHR", bp = "BP", snp = "SNP", p = "P")
manhattan(day0_filtered, chr = "snp_chromosome", bp = "snp_position", snp = "snp_id", p = "p_value", main ="day0", col = c("grey", "skyblue"), suggestiveline = -log10(5e-6), genomewideline = -log10(5e-8))
manhattan(day0_filtered, chr = "snp_chromosome", bp = "snp_position", snp = "snp_id", p = "empirical_feature_p_value", logp = FALSE, main ="day0", col = c("grey", "skyblue"), genomewideline = FDR)

manhattan(day1_filtered, chr = "snp_chromosome", bp = "snp_position", snp = "snp_id", p = "p_value", main = "day1",col = c("grey", "skyblue"), suggestiveline = -log10(5e-6), genomewideline = -log10(5e-8))
manhattan(day1_filtered, chr = "snp_chromosome", bp = "snp_position", snp = "snp_id", p = "empirical_feature_p_value", logp = FALSE, main = "day1",col = c("grey", "skyblue"), genomewideline = FDR)

manhattan(day2_filtered, chr = "snp_chromosome", bp = "snp_position", snp = "snp_id", p = "p_value", main = "day2", col = c("grey", "skyblue"), suggestiveline = -log10(5e-6), genomewideline = -log10(5e-8))
manhattan(day2_filtered, chr = "snp_chromosome", bp = "snp_position", snp = "snp_id", p = "empirical_feature_p_value", logp = FALSE, main = "day2", col = c("grey", "skyblue"), genomewideline = FDR)

manhattan(day3_filtered, chr = "snp_chromosome", bp = "snp_position", snp = "snp_id", p = "p_value", main = "day3", col = c("grey", "skyblue"), suggestiveline = -log10(5e-6), genomewideline = -log10(5e-8))
manhattan(day3_filtered, chr = "snp_chromosome", bp = "snp_position", snp = "snp_id", p = "empirical_feature_p_value", logp = FALSE, main = "day3", col = c("grey", "skyblue"), genomewideline = FDR)

#make qq plot => colroed according to day 

observed_day0 <- day0$p_value
lobs_day0 <- -(log10(observed_day0))
expected_day0 <- c(1:length(observed_day0)) 
lexp_day0 <- -(log10(expected_day0 / (length(expected_day0)+1)))

observed_day1 <- day1$p_value
lobs_day1 <- -(log10(observed_day1))
expected_day1 <- c(1:length(observed_day1)) 
lexp_day1 <- -(log10(expected_day1 / (length(expected_day1)+1)))

observed_day2 <- day2$p_value
lobs_day2 <- -(log10(observed_day2))
expected_day2 <- c(1:length(observed_day2)) 
lexp_day2 <- -(log10(expected_day2 / (length(expected_day2)+1)))

observed_day3 <- day3$p_value
lobs_day3 <- -(log10(observed_day3))
expected_day3 <- c(1:length(observed_day3)) 
lexp_day3 <- -(log10(expected_day3 / (length(expected_day3)+1)))

#color_map = {'day0': "#440154FF", "day1": "#33638DFF", "day2": "#29AF7FFF", "day3": "#FDE725FF"}

plot.new()
png("/home/j890e/scvelo/figures/qqplot_pseudotime.png")
plot(c(0,8), c(0,8), col="grey", lwd=2, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,7), ylim=c(0,7), las=1, xaxs="i", yaxs="i", bty="l", lty=2)
points(lexp_day0, lobs_day0, pch=23, cex=.3, col="#440154FF") 
points(lexp_day1, lobs_day1, pch=23, cex=.3, col="#33638DFF") 
points(lexp_day2, lobs_day2, pch=23, cex=.3, col="#29AF7FFF") 
points(lexp_day3, lobs_day3, pch=23, cex=.3, col="#FDE725FF") 
dev.off()


par(new=T) 
#manhattan plot of specific region using ggplot-----------------------------------------------------------------------


day3_chr_subset = subset(day3, snp_chromosome == 12)
day0_chr_subset = subset(day0, snp_chromosome == 12)
day1_chr_subset = subset(day1, snp_chromosome == 12)
day2_chr_subset = subset(day2, snp_chromosome == 12)

snp_pos = 83511984

day0_chr_subset_region = filter(day0_chr_subset, between(snp_position, 83300000, 83700000))#4.36e7, 4.37e7))
day1_chr_subset_region = filter(day1_chr_subset, between(snp_position, 83300000, 83700000))#4.36e7, 4.37e7))
day2_chr_subset_region = filter(day2_chr_subset, between(snp_position, 83300000, 83700000))#4.36e7, 4.37e7))
day3_chr_subset_region = filter(day3_chr_subset, between(snp_position, 83300000, 83700000))#4.36e7, 4.37e7))

new = rbind(day0_chr_subset_region, day1_chr_subset_region, day2_chr_subset_region, day3_chr_subset_region)

QTLs = new[,c("snp_chromosome", "snp_position", "snp_id", "p_value", "timepoint")]
colnames(QTLs) = c("CHR", "BP", "SNP", "P", "timepoint")
subset = filter(QTLs, CHR == 12)
df_clean = subset

manh_plot <- function(df, threshold) {
  
  ### 1. Compute the cumulative position of SNP ### 
  plot_data <- df %>%   
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=as.numeric(max(BP))) %>% 
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    dplyr::select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(df_clean, ., by=c("CHR"="CHR")) %>%
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=as.numeric(BP+tot))
  
  ### 2. Generate x-axis ###
  axisdf <- plot_data %>% 
    group_by(CHR) %>% 
    summarize(center=(max(BPcum) + min(BPcum)) / 2 )
  
  ### 3. create plot ###
  plot <- ggplot(plot_data, aes(x=BPcum, y=-log10(P))) + 
    #specify the y and x values
    geom_point( aes(color=as.factor(timepoint)), alpha=0.8, size=0.6) + 
    # create scatterplot colored by chromosome
    scale_color_manual(values = c("#440154", "#173F5F", "#3CAEA3", "#F6D55C")) + 
    # set a colour pattern 
    scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center) + 
    # scale the x-axis
    scale_y_continuous(expand = c(0, 0)) + 
    # remove space between plot area and x axis
    ylim(0,8) +
    theme_light() +
    theme(legend.position="none",
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.line = element_line(color = "black")) +
    xlab("Chromosome") + 
    # add x label
    geom_label_repel( data=plot_data %>% filter(P < threshold), # add annotation value
                      aes(label=SNP), size=3) + # add annotation
    geom_point(data= plot_data %>% filter(P < threshold), # add annotation value
               color="orange", size=2) + # Add highlighted points 
    geom_hline(yintercept = -log10(threshold), linetype="dashed") # threshold line
  
  return(plot) # return the final plot
}

pl <- manh_plot(df_clean, 0)

pdf(file = "/home/j890e/scvelo/figures/GWAS_pseudotime.pdf", width = 5, height = 4)
pl +
  coord_cartesian(xlim = c(83300000, 83700000))
dev.off()