# Manhattan plot for gea

# Load the libraries 
library(qqman)
library(cowplot)
library(ggrepel)
library(dplyr)

# Read GWAS output file created using gemma
gwas_data_bio19 <- read.table("gwas_project_result.assoc.txt", header = TRUE)

# Make the Manhattan plot on the gwas Results dataset
manhattan(gwas_data_bio19, chr="chr", bp="ps", snp="rs", p="p_wald" )

# modify colour of points
manhattan(x = gwas_data_bio19, chr = "chr", bp = "ps", p = "p_wald", snp = "rs", col = c("blue", "red"), logp = TRUE)

# Multiple testing correction in GWAS (Bonferroni correction)
pval_bonf = 0.05/dim(gwas_data_bio19)[[1]]


#adding multiple-testing-corrected p-value line to the plot
manhattan(gwas_data_bio19, chr="chr", bp="ps", snp="rs", p="p_wald",
          suggestiveline = -log10(pval_bonf), genomewideline = FALSE,
          annotatePval = -log10(pval_bonf), col = c("blue4", "red"),
          ylim = c(0, max(-log10(pval_bonf)) + 5),
          main = "GEA Manhattan Plot with Bonferroni correction                                                                                A")


# Read GWAS output file created using gemma
gwas_data_bio19 <- read.table("C:/Users/aless/OneDrive/Desktop/Projekt/gwas_project_result.assoc.txt", header = TRUE)

# Make the Manhattan plot on the gwas Results dataset
manhattan(gwas_data_bio19, chr="chr", bp="ps", snp="rs", p="p_wald" )

# modify colour of points
manhattan(x = gwas_data_bio19, chr = "chr", bp = "ps",
          p = "p_wald", snp = "rs", col = c("blue4", "red"), logp = TRUE)

#Bonferroni correction

# Make the Manhattan plot using FDR-corrected p-values
pvals_bonf <- p.adjust(gwas_data_bio19$p_wald, method = "bonferroni")

# Add FDR-corrected p-value threshold line to the plot
pval_bonf <- 0.05  # Desired FDR threshold, you can adjust this if needed

# Add the adjusted p-values as a new column to the data frame
gwas_data_bio19$pvals_bonf <- pvals_bonf

# Create the Manhattan plot with the new adjusted p-values column
manhattan(gwas_data_bio19, chr = "chr", bp = "ps", snp = "rs", p = "pvals_bonf",
          suggestiveline = -log10(pval_bonf), genomewideline = FALSE,
          annotatePval = -log10(pval_bonf),
          col = c("blue4", "red"),
          main = "GEA Manhattan Plot with FDR correction                                                                                 B")


#FDR correction

# Make the Manhattan plot using FDR-corrected p-values
pvals_fdr <- p.adjust(gwas_data_bio19$p_wald, method = "BH")

# Add FDR-corrected p-value threshold line to the plot
pval_fdr <- 0.05  # Desired FDR threshold

# Add the adjusted p-values as a new column to the data frame
gwas_data_bio19$pvals_fdr <- pvals_fdr

# Create the Manhattan plot with the new adjusted p-values column
manhattan(gwas_data_bio19, chr = "chr", bp = "ps", snp = "rs", p = "pvals_fdr",
          suggestiveline = -log10(pval_fdr), genomewideline = FALSE,
          annotatePval = -log10(pval_fdr),
          col = c("blue4", "red"),
          main = "GEA Manhattan Plot with FDR correction                                                                                 B")


####################################
# GWAS manhatten Plot & SNP annotation
# using ggrepel, dplyr, ggplot

# Prepare the dataset
gwas <- gwas_data_bio19 %>% 
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(ps)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gwas_data_bio19, ., by=c("chr"="chr")) %>%
  # Add a cumulative position of each SNP
  arrange(chr, ps) %>%
  mutate(Chromosome=ps+tot,   # Fixed BPcum assignment
         is_annotate=ifelse(-log10(pvals_fdr)>1.2, "yes", "no")) 

# Prepare X axis
axisdf <- gwas %>% group_by(chr) %>% summarize(center=( max(Chromosome) + min(Chromosome) ) / 2 )
# extract labels column
rs<-gwas$rs

# plot using ggplot

ggplot(gwas, aes(x=Chromosome, y=-log10(pvals_fdr))) +
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("red", "blue"), 22 )) +
  scale_x_continuous(label = axisdf$chr, breaks= axisdf$center) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.8)), breaks = seq(0, 5, 0.5), 
                     limits = c(0,1.5)) +
  geom_label_repel(data=subset(gwas, is_annotate=="yes"), aes(label=rs), size=2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 8)
  )


