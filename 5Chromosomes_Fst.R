# Load required libraries

library(PopGenome) 
library(ggplot2)
library (readr)
library(tibble)
library(vcfR)
library(adegenet)
library(factoextra)
library(FactoMineR)
library(tidyverse)
library(ggrepel)
library(gplots)
library(StAMPP)
library(RColorBrewer)
library(dplyr)
library(VariantAnnotation)
library(tidyr)


# set working directory

setwd("C:/Users/aless/OneDrive/Desktop/Project_pop&qGen/Fst")

# FST calculation between populations 

# Load VCF file
vcf_file <- read.vcfR("group_3_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz")

# Convert VCF to genlight object
genlight_vcf <- vcfR2genlight(vcf_file)

# Read population information from file
pop <- read.table("sample_pop_it_swe.txt", header = TRUE)
str(pop)

# Extract population data
pop1 <- pop$pop
pop2 = as.factor(pop1)
genlight_vcf$pop = pop2

# Convert genlight to stampp object
stampp_vcf = stamppConvert(genlight_vcf, type = "genlight")

# Calculate FST between populations
stamppFst = stamppFst(stampp_vcf, nboots = 100, percent = 95, nclusters = 8)
stamppFst_matrix = as.matrix(stamppFst$Fsts)

# Set diagonal to 0 and add symmetrical upper tri part of matrix
diag(stamppFst_matrix) <- 0
stamppFst_matrix[upper.tri(stamppFst_matrix)]  <- t(stamppFst_matrix)[upper.tri(stamppFst_matrix)]
# Optional: order the names
stamppFst_matrix = stamppFst_matrix[order(row.names(stamppFst_matrix)), order(colnames(stamppFst_matrix))]

# Make an FST heatmap
heatmap(stamppFst_matrix,
        symm = TRUE,
        margins = c(10, 10),
        main = "Genetic Divergence (FST) b/w A. thaliana Pop from Sweden & Italy ")


# Genetic distance calculation between individuals - nei's distance

stamppNeisD = stamppNeisD(stampp_vcf, pop = FALSE)
stamppNeisD_matrix = as.matrix(stamppNeisD)

# Set diagonal to 0 and add symmetrical upper tri part of matrix
diag(stamppNeisD_matrix) <- 0
stamppNeisD_matrix[upper.tri(stamppNeisD_matrix)]  <- t(stamppNeisD_matrix)[upper.tri(stamppNeisD_matrix)]
heatmap(stamppNeisD_matrix)

# add row names
colnames(stamppNeisD_matrix) <- rownames(stamppNeisD_matrix)

# Create a heatmap with symmetric color scale
heatmap(stamppNeisD_matrix,
        symm = TRUE,
        main = "Genetic Divergence (FST) b/w A. thaliana individuals from Sweden & Italy ")



#get chromosomes start and end points
# Read the VCF file

vcf <- readVcf("group_3_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz") 
# Extract chromosome names

chromosomes <- seqlevels(vcf) 
# Initialize an empty data frame to store results

chromosome_ranges <- data.frame(CHROM = character(), Start = numeric(), End = numeric(), stringsAsFactors = FALSE) 
# Iterate over chromosomes 
for (chrom in chromosomes) { 
  # Extract positions for the current chromosome
  
  positions <- start(vcf[seqnames(vcf) == chrom]) 
  # Append results to the data frame
  
  chromosome_ranges <- rbind(chromosome_ranges, data.frame(CHROM = chrom, Start = min(positions), End = max(positions))) 
} 
# Display the results

print(chromosome_ranges)

# Estimate and plot Fst and Tajima'D and Neutrality stats using PopGenome

At_Chr1 <-readVCF("group_3_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz",numcols=89,tid="1", frompos=1373683, topos=30003409, include.unknown =  TRUE)
At_Chr2 <-readVCF("group_3_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz",numcols=89,tid="2", frompos=2332947, topos=11746232, include.unknown =  TRUE)
At_Chr3 <-readVCF("group_3_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz",numcols=89,tid="3", frompos=14093721, topos=20782721, include.unknown =  TRUE)
At_Chr4 <-readVCF("group_3_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz",numcols=89,tid="4", frompos=269237, topos=18328342, include.unknown =  TRUE)
At_Chr5 <-readVCF("group_3_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz",numcols=89,tid="5", frompos=1602625, topos=25995710, include.unknown =  TRUE)

#Examining the variant data
get.sum.data(At_Chr1) 
get.sum.data(At_Chr2) 
get.sum.data(At_Chr3) 
get.sum.data(At_Chr4) 
get.sum.data(At_Chr5) 

At_Chr1@n.biallelic.sites #529
At_Chr2@n.biallelic.sites #222
At_Chr3@n.biallelic.sites #84
At_Chr4@n.biallelic.sites #281
At_Chr5@n.biallelic.sites #403

At_Chr1@n.polyallelic.sites #0
At_Chr2@n.polyallelic.sites #0
At_Chr3@n.polyallelic.sites #0
At_Chr4@n.polyallelic.sites #0
At_Chr5@n.polyallelic.sites #0

#To see what slots in Genome class
#show.slots(At_Chr1)

#To check total number of sites
At_Chr1@n.sites # 28629727
At_Chr2@n.sites # 9413286
At_Chr3@n.sites # 6689001
At_Chr4@n.sites # 18059106
At_Chr5@n.sites # 24393086

#To check starting position and last position of genome class
#At_Chr1@region.names # "1373683 - 30003409"

##Deine populations in your dataset

library (readr)

population_info <- read_delim("sample_pop_it_swe.txt", delim = "\t")

# Get the data for the populations
populations1 <- split(population_info$sample, population_info$pop)
populations1

At_Chr1 <- set.populations(At_Chr1, populations1, diploid = T)
At_Chr1@populations

At_Chr2 <- set.populations(At_Chr2, populations1, diploid = T)
At_Chr2@populations

At_Chr3 <- set.populations(At_Chr3, populations1, diploid = T)
At_Chr3@populations

At_Chr4 <- set.populations(At_Chr4, populations1, diploid = T)
At_Chr4@populations

At_Chr5 <- set.populations(At_Chr5, populations1, diploid = T)
At_Chr5@populations

#Setting up sliding windows

# set chromosome size (= total number of sites)
chr1 <- 28629727
chr2 <- 9413286
chr3 <- 6689002  
chr4 <- 18059106
chr5 <- 24393086

# set window size and window jump
window_size <- 100
window_jump <- 50

#Chr 1
# use seq to find the start points of each window
window_start1 <- seq(from = 1, to = chr1, by = window_jump)
# add the size of the window to each start point 
window_stop1 <- window_start1 + window_size

# no windows start before the end of chromosome 4
sum(window_start1 > chr1)
# 0

# but some window stop positions do occur past the final point
sum(window_stop1 > chr1)
# 2

# remove windows from the start and stop vectors
window_start1 <- window_start1[which(window_stop1 < chr1)]
window_stop1 <- window_stop1[which(window_stop1 < chr1)]

chr1 - window_stop1[length(window_stop1)]
# 26

# save as a data.frame
windows1 <- data.frame(start = window_start1, stop = window_stop1, 
                      mid = window_start1 + (window_stop1-window_start1)/2)


# chr 2
# use seq to find the start points of each window
window_start2 <- seq(from = 1, to = chr2, by = window_jump)
# add the size of the window to each start point 
window_stop2 <- window_start2 + window_size

# no windows start before the end of chromosome 4
sum(window_start2 > chr2)
# 0

# but some window stop positions do occur past the final point
sum(window_stop2 > chr2)
# 2

# remove windows from the start and stop vectors
window_start2 <- window_start2[which(window_stop2 < chr2)]
window_stop2 <- window_stop2[which(window_stop2 < chr2)]

chr2 - window_stop2[length(window_stop2)]
# 35

# save as a data.frame
windows2 <- data.frame(start = window_start2, stop = window_stop2, 
                       mid = window_start2 + (window_stop2-window_start2)/2)


# chr 3 
# use seq to find the start points of each window
window_start3 <- seq(from = 1, to = chr3, by = window_jump)
# add the size of the window to each start point 
window_stop3 <- window_start3 + window_size

# no windows start before the end of chromosome 4
sum(window_start3 > chr3)
# 0

# but some window stop positions do occur past the final point
sum(window_stop3 > chr3)
# 2

# remove windows from the start and stop vectors
window_start3 <- window_start3[which(window_stop3 < chr3)]
window_stop3 <- window_stop3[which(window_stop3 < chr3)]

chr3 - window_stop3[length(window_stop3)]
# 1

# save as a data.frame
windows3 <- data.frame(start = window_start3, stop = window_stop3, 
                       mid = window_start3 + (window_stop3-window_start3)/2)


# chr 4 
# use seq to find the start points of each window
window_start4 <- seq(from = 1, to = chr4, by = window_jump)
# add the size of the window to each start point 
window_stop4 <- window_start4 + window_size

# no windows start before the end of chromosome 4
sum(window_start4 > chr4)
# 0

# but some window stop positions do occur past the final point
sum(window_stop4 > chr4)
# 2

# remove windows from the start and stop vectors
window_start4 <- window_start4[which(window_stop4 < chr4)]
window_stop4 <- window_stop4[which(window_stop4 < chr4)]

chr4 - window_stop4[length(window_stop4)]
# 5

# save as a data.frame
windows4 <- data.frame(start = window_start4, stop = window_stop4, 
                       mid = window_start4 + (window_stop4-window_start4)/2)


# chr 5 
# use seq to find the start points of each window
window_start5 <- seq(from = 1, to = chr5, by = window_jump)
# add the size of the window to each start point 
window_stop5 <- window_start5 + window_size

# no windows start before the end of chromosome 4
sum(window_start5 > chr5)
# 0

# but some window stop positions do occur past the final point
sum(window_stop5 > chr5)
# 2

# remove windows from the start and stop vectors
window_start5 <- window_start5[which(window_stop5 < chr5)]
window_stop5 <- window_stop5[which(window_stop5 < chr5)]

chr5 - window_stop5[length(window_stop5)]
# 35

# save as a data.frame
windows5 <- data.frame(start = window_start5, stop = window_stop5, 
                       mid = window_start5 + (window_stop5-window_start5)/2)


# make sliding window datasets
At_sw1 <- sliding.window.transform(At_Chr1, width = 100, jump = 50, type = 2)
At_sw2 <- sliding.window.transform(At_Chr2, width = 100, jump = 50, type = 2)
At_sw3 <- sliding.window.transform(At_Chr3, width = 100, jump = 50, type = 2)
At_sw4 <- sliding.window.transform(At_Chr4, width = 100, jump = 50, type = 2)
At_sw5 <- sliding.window.transform(At_Chr5, width = 100, jump = 50, type = 2)

# calculate diversity statistics - nd
At_sw1 <- diversity.stats(At_sw1, pi = TRUE)
At_sw2 <- diversity.stats(At_sw2, pi = TRUE)
At_sw3 <- diversity.stats(At_sw3, pi = TRUE)
At_sw4 <- diversity.stats(At_sw4, pi = TRUE)
At_sw5 <- diversity.stats(At_sw5, pi = TRUE)

# calculate diversity statistics - FST
At_sw1 <- F_ST.stats(At_sw1, mode = "nucleotide") 
At_sw2 <- F_ST.stats(At_sw2, mode = "nucleotide") 
At_sw3 <- F_ST.stats(At_sw3, mode = "nucleotide") 
At_sw4 <- F_ST.stats(At_sw4, mode = "nucleotide") 
At_sw5 <- F_ST.stats(At_sw5, mode = "nucleotide") 


#Extracting statistics for visualization

# get the nucleotide diversity data.
# extract nucleotide diversity and correct for window size
nd1 <- At_sw1@nuc.diversity.within/100
nd2 <- At_sw2@nuc.diversity.within/100
nd3 <- At_sw3@nuc.diversity.within/100
nd4 <- At_sw4@nuc.diversity.within/100
nd5 <- At_sw5@nuc.diversity.within/100

#estimates need to be corrected for window size - so we divide them by 100 bp.

# Add the population names to each of estimate
# make population name vector
pops <- c("IT","SWE") 

# set population names
colnames(nd1) <- paste0(pops, "_pi")
colnames(nd2) <- paste0(pops, "_pi")
colnames(nd3) <- paste0(pops, "_pi")
colnames(nd4) <- paste0(pops, "_pi")
colnames(nd5) <- paste0(pops, "_pi")

# extract fst values
fst1 <- t(At_sw1@nuc.F_ST.pairwise)
fst2 <- t(At_sw2@nuc.F_ST.pairwise)
fst3 <- t(At_sw3@nuc.F_ST.pairwise)
fst4 <- t(At_sw4@nuc.F_ST.pairwise)
fst5 <- t(At_sw5@nuc.F_ST.pairwise)


# extract dxy - pairwise absolute nucleotide diversity
dxy1 <- get.diversity(At_sw1, between = T)[[2]]/100
dxy2 <- get.diversity(At_sw2, between = T)[[2]]/100
dxy3 <- get.diversity(At_sw3, between = T)[[2]]/100
dxy4 <- get.diversity(At_sw4, between = T)[[2]]/100
dxy5 <- get.diversity(At_sw5, between = T)[[2]]/100


#As with nucleotide diversity, we also corrected d_XY_ for the window size.

# get column names
x1 <- colnames(fst1)
# does the same thing as above but by indexing the pops vector
x1 <- sub("pop1", pops[1], x1)
x1 <- sub("pop2", pops[2], x1)

# replace forward slash
x1 <- sub("/", "_", x1)
# look at x1 to confirm the replacement has occurred
x1

# Change the column names of our two data matrices, before we put 
#everything together in our final dataset.
colnames(fst1) <- paste0(x1, "_fst1")
colnames(dxy1) <- paste0(x1, "_dxy1")

colnames(fst2) <- paste0(x1, "_fst2")
colnames(dxy2) <- paste0(x1, "_dxy2")

colnames(fst3) <- paste0(x1, "_fst3")
colnames(dxy3) <- paste0(x1, "_dxy3")

colnames(fst4) <- paste0(x1, "_fst4")
colnames(dxy4) <- paste0(x1, "_dxy4")

colnames(fst5) <- paste0(x1, "_fst5")
colnames(dxy5) <- paste0(x1, "_dxy5")

# Combine nd, FST and d_XY_ datasets with our windows information from earlier into a big dataset.

At_data1 <- as_tibble(data.frame(windows1, nd1, fst1, dxy1))
dim(At_data1) #572593     7

At_data2 <- as_tibble(data.frame(windows2, nd2, fst2, dxy2))
At_data3 <- as_tibble(data.frame(windows3, nd3, fst3, dxy3))
At_data4 <- as_tibble(data.frame(windows4, nd4, fst4, dxy4))
At_data5 <- as_tibble(data.frame(windows5, nd5, fst5, dxy5))

#####Visualizing the data - distributions

#For the purposes of this session, we will focus mainly on the difference between Italian and Swedish
#Arabidopsis pop.
#For example, let's say we want to look at mean nucleotide diversity, we can do that like so:

# select nucleotide diversity data and calculate means
dplyr::select(At_data1, contains("pi")) %>% dplyr::summarise_all(mean)
dplyr::select(At_data2, contains("pi")) %>% dplyr::summarise_all(mean)
dplyr::select(At_data3, contains("pi")) %>% dplyr::summarise_all(mean)
dplyr::select(At_data4, contains("pi")) %>% dplyr::summarise_all(mean)
dplyr::select(At_data5, contains("pi")) %>% dplyr::summarise_all(mean)


# To plot this we need to use "gather" on the data

pi_g1 <- At_data1 %>% dplyr::select(contains("pi")) %>% gather(key = "populations", value = "pi")
pi_g2 <- At_data2 %>% dplyr::select(contains("pi")) %>% gather(key = "populations", value = "pi")
pi_g3 <- At_data3 %>% dplyr::select(contains("pi")) %>% gather(key = "populations", value = "pi")
pi_g4 <- At_data4 %>% dplyr::select(contains("pi")) %>% gather(key = "populations", value = "pi")
pi_g5 <- At_data5 %>% dplyr::select(contains("pi")) %>% gather(key = "populations", value = "pi")

# make boxplots
library(ggplot2)
pi_g1$log_pi <- log10(pi_g1$pi)
pi_g2$log_pi <- log10(pi_g2$pi)
pi_g3$log_pi <- log10(pi_g3$pi)
pi_g4$log_pi <- log10(pi_g4$pi)
pi_g5$log_pi <- log10(pi_g5$pi)

library(cowplot)

# Create the individual boxplot objects
a1 <- ggplot(pi_g1, aes(populations, log_pi, fill = populations)) + 
  geom_boxplot(color = "black") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_light() + 
  xlab(NULL) +
  ylab("Log10(pi)") +
  ggtitle("Chromosome 1") +
  annotate("text", x = 1.5, y = max(pi_g1$log_pi), label = paste("p =", format.pval(wilcox_test_result1$p.value, digits = 3)))

a2 <- ggplot(pi_g2, aes(populations, log_pi, fill = populations)) + 
  geom_boxplot(color = "black") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_light() + 
  xlab(NULL) +
  ylab("Log10(pi)") +
  ggtitle("Chromosome 2") +
  annotate("text", x = 1.5, y = max(pi_g2$log_pi), label = paste("p =", format.pval(wilcox_test_result2$p.value, digits = 3)))

a3 <- ggplot(pi_g3, aes(populations, log_pi, fill = populations)) + 
  geom_boxplot(color = "black") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_light() + 
  xlab(NULL) +
  ylab("Log10(pi)") +
  ggtitle("Chromosome 3") +
  annotate("text", x = 1.5, y = max(pi_g3$log_pi), label = paste("p =", format.pval(wilcox_test_result3$p.value, digits = 3)))

a4 <- ggplot(pi_g4, aes(populations, log_pi, fill = populations)) + 
  geom_boxplot(color = "black") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_light() + 
  xlab(NULL) +
  ylab("Log10(pi)") +
  ggtitle("Chromosome 4") +
  annotate("text", x = 1.5, y = max(pi_g4$log_pi), label = paste("p =", format.pval(wilcox_test_result4$p.value, digits = 3)))

a5 <- ggplot(pi_g5, aes(populations, log_pi, fill = populations)) + 
  geom_boxplot(color = "black") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_light() + 
  xlab(NULL) +
  ylab("Log10(pi)") +
  ggtitle("Chromosome 5") +
  annotate("text", x = 1.5, y = max(pi_g5$log_pi), label = paste("p =", format.pval(wilcox_test_result5$p.value, digits = 3)))

# Combine the individual boxplots into a single plot
combined_plots <- plot_grid(a1, a2, a3, a4, a5, labels = "AUTO", nrow = 2)

overall_title <- ggdraw() +
  draw_label("Nucleotide diversity as a function of sample location", size = 20, hjust = 0.5) +
  theme(plot.title = element_text(face = "bold"))

# Combine the overall title and individual boxplots
combined_plots_with_title <- plot_grid(
  overall_title,
  combined_plots,
  ncol = 1,
  rel_heights = c(0.1, 0.9)
)

# Print the combined plot with an overall title
print(combined_plots_with_title)



#Visualizing patterns along the chromosome

#Let's have a look at how FST between Italian and Swedish populations varies along chromosomes.

# chr 1
b1 <- ggplot(At_data1, aes(mid/10^6, IT_SWE_fst1)) +
  geom_line(colour = "red") +
  xlab("Position (Mb)") +
  ylab(expression(italic(F)[ST])) +
  theme_light() +
  ggtitle("Chromosome 1")

# chr 2
b2 <- ggplot(At_data2, aes(mid/10^6, IT_SWE_fst2)) +
  geom_line(colour = "orange") +
  xlab("Position (Mb)") +
  ylab(expression(italic(F)[ST])) +
  theme_light() +
  ggtitle("Chromosome 2")

# chr 3
b3 <- ggplot(At_data3, aes(mid/10^6, IT_SWE_fst3)) +
  geom_line(colour = "yellow2") +
  xlab("Position (Mb)") +
  ylab(expression(italic(F)[ST])) +
  theme_light() +
  ggtitle("Chromosome 3")

# chr 4
b4 <- ggplot(At_data4, aes(mid/10^6, IT_SWE_fst4)) +
  geom_line(colour = "green") +
  xlab("Position (Mb)") +
  ylab(expression(italic(F)[ST])) +
  theme_light() +
  ggtitle("Chromosome 4")

# chr 5
b5 <- ggplot(At_data5, aes(mid/10^6, IT_SWE_fst5)) +
  geom_line(colour = "blue") +
  xlab("Position (Mb)") +
  ylab(expression(italic(F)[ST])) +
  theme_light() +
  ggtitle("Chromosome 5")

combined_plots <- plot_grid(b1, b2, b3, b4, b5, labels = "AUTO", nrow = 2)

overall_title <- ggdraw() +
  draw_label("Variation of FST between Italian and Swedish populations along chromosomes", size = 20, hjust = 0.5) +
  theme(plot.title = element_text(face = "bold"))

# Combine the overall title and individual boxplots
combined_plots_with_title <- plot_grid(
  overall_title,
  combined_plots,
  ncol = 1,
  rel_heights = c(0.1, 0.9)
)

# Print the combined plot with an overall title
print(combined_plots_with_title)


#to plot nd, FST and d_XY_ to examine how they co-vary along the genome. 
#This requires a bit of data manipulation, but is relatively straightforward. We will break it down into steps.
# select data of interest
hs1 <- At_data1 %>% dplyr::select(mid, IT_pi, SWE_pi, IT_SWE_fst1, IT_SWE_dxy1)
hs2 <- At_data2 %>% dplyr::select(mid, IT_pi, SWE_pi, IT_SWE_fst2, IT_SWE_dxy2)
hs3 <- At_data3 %>% dplyr::select(mid, IT_pi, SWE_pi, IT_SWE_fst3, IT_SWE_dxy3)
hs4 <- At_data4 %>% dplyr::select(mid, IT_pi, SWE_pi, IT_SWE_fst4, IT_SWE_dxy4)
hs5 <- At_data5 %>% dplyr::select(mid, IT_pi, SWE_pi, IT_SWE_fst5, IT_SWE_dxy5)

# To set Fst values smaller than zero to zero in the specified columns of a data frame using dplyr and
# the pipe operator %>%, you can use the mutate function along with across

suppressWarnings({
hs1 <- At_data1 %>%
  dplyr::select(mid, IT_pi, SWE_pi, IT_SWE_fst1, IT_SWE_dxy1) %>%
  mutate(across(c(IT_SWE_fst1, IT_SWE_dxy1), ~ ifelse(. < 0, 0, .)))})

suppressWarnings({
hs2 <- At_data2 %>%
  dplyr::select(mid, IT_pi, SWE_pi, IT_SWE_fst2, IT_SWE_dxy2) %>%
  mutate(across(c(IT_SWE_fst2, IT_SWE_dxy2), ~ ifelse(. < 0, 0, .)))})

suppressWarnings({
hs3 <- At_data3 %>%
  dplyr::select(mid, IT_pi, SWE_pi, IT_SWE_fst3, IT_SWE_dxy3) %>%
  mutate(across(c(IT_SWE_fst3, IT_SWE_dxy3), ~ ifelse(. < 0, 0, .)))})

suppressWarnings({
hs4 <- At_data4 %>%
  dplyr::select(mid, IT_pi, SWE_pi, IT_SWE_fst4, IT_SWE_dxy4) %>%
  mutate(across(c(IT_SWE_fst4, IT_SWE_dxy4), ~ ifelse(. < 0, 0, .)))})

suppressWarnings({
hs5 <- At_data5 %>%
  dplyr::select(mid, IT_pi, SWE_pi, IT_SWE_fst5, IT_SWE_dxy5) %>%
  mutate(across(c(IT_SWE_fst5, IT_SWE_dxy5), ~ ifelse(. < 0, 0, .)))})


# use gather to rearrange everything
hs_g1 <- gather(hs1, -mid, key = "stat", value = "value")
hs_g2 <- gather(hs2, -mid, key = "stat", value = "value")
hs_g3 <- gather(hs3, -mid, key = "stat", value = "value")
hs_g4 <- gather(hs4, -mid, key = "stat", value = "value")
hs_g5 <- gather(hs5, -mid, key = "stat", value = "value")







# To take the logarithm of the value variable in your ggplot code, you can use the log10() function
# within the aes() mapping.
hs_g1$log_value <- log10(hs_g1$value)

#However, before we examine our plot in detail, it would also be easier if we rearranged everything 
#so FST came at the top, pi beneath it and then finally, d_XY_. How can we do that? Well we need to
#reorder the stat factor in our hs_g dataset.
# first make a factor
x <- factor(hs_g1$stat)
# then reorder the levels
x <- factor(x, levels(x)[c(3, 1, 4, 2)])
# add to data.frame
hs_g1$stat <- x
#This looks a little complicated, but in the square brackets above we simply rearranged what order
#our facets are displayed. We can replot our figure to demonstrate this:

# construct a plot with facets
a <- ggplot(hs_g1, aes(mid/10^6, value, colour = stat)) + geom_line()
a <- a + facet_grid(stat~., scales = "free_y")
a <- a + xlab("Position (Mb)")
a + theme_light() + theme(legend.position = "none")

#### calculate neutrality statistics####
At_sw1 <- neutrality.stats(At_sw1)

get.neutrality(At_sw1)
#      neurality stats
#pop 1 numeric,5153337
#pop 2 numeric,5153337

#Let's look at the first population [[1]].
get.neutrality(At_sw1)[[1]]

#Let's look at the second population [[2]].
get.neutrality(At_sw1)[[2]]

#extract Tajma's D
td1 <- At_sw1@Tajima.D/100

# set population names
colnames(td1) <- paste0(pops, "_td")


###Delimitate windows on chromosome

# set chromosome start and end position
chr1_start<- 1373683
chr1_end <- 30003409

library(tibble)

#as_tibble: Coerce lists and matrices to data frames
ara_data1 <- as.tibble(data.frame(windows1, td1,nd1)) #tajima, nucleotide diversity pi
nrow(windows)
nrow(nd)
(chr1_end-chr1_start)/50
nrow(ara_data1)
head(ara_data1)
ara_data1 %>% dplyr::select(contains("pi")) %>% summarise_all(mean)

### load selected positions from chromosome e.g., gene 4 5kb upstream and down stream of Defense related genes

bed<-read.table("At_defense_only.bed")
View(bed)

colnames(bed)<-c("chr", "begin","end")
DF1<-vector(length=nrow(ara_data1))
length(DF1)

#ara_data <- as.tibble(data.frame(windows, nd, DF))###if you only want to look at pi
ara_data1 <- as.tibble(data.frame(windows1, nd1, td1, DF1))##if you want to look at tajima D and nucleotide diversity

for (i in 2:nrow(bed)){ara_data1$DF1[which(ara_data1$start>bed$begin[i]&ara_data1$stop<bed$end[i]) ]<-"DF"}##each window that overlaps a DF is tagged
ara_data1$DF1<-as.factor(ara_data1$DF1)
summary(ara_data1)

####
# Italy
###

##Kolmogorov smirnov test - compare the distributions of Pi
sub1<-ara_data1$IT_pi[ara_data1$DF1=="DF"]
sub2<-ara_data1$IT_pi[ara_data1$DF1!="DF"]
ks.test(sub1, sub2)###difference is very significant if windows are small, otherwise not. 

#Draw Density plot "Pi"
plot(density(log(ara_data1$IT_pi)), main="Distribution log Pi")
lines(density(log(ara_data1$IT_pi[ara_data1$DF1=="DF"])), col="red")

##Kolmogorov smirnov test - compare the distributions of Tajima's D
sub1<-ara_data1$IT_td[ara_data1$DF1=="DF"]
sub2<-ara_data1$IT_td[ara_data1$DF1!="DF"]
ks.test(sub1, sub2) 

# Draw Density plots "Tajima's D"
plot(density((ara_data1$IT_td), na.rm=T), main="Distribution Tajima D", ylim = c(0, 230))
lines(density((ara_data1$IT_td[ara_data1$DF1=="DF"]), na.rm = T), col="red")

#p<-ggplot(ara_data1, aes(x=IT_td, fill=DF))
#p+geom_density(alpha=0.4)

# To estimate the lowest value of the x-axis for a density plot
lowest_x <- min(ara_data1$IT_td, na.rm = TRUE)
lowest_x

p <- ggplot(ara_data1, aes(x = IT_td, fill = DF)) +
  geom_density(alpha = 0.4) +
  scale_x_continuous(limits = c(-0.03, max(ara_data1$IT_td, na.rm = TRUE))) +  # Set x-axis limit
  ggtitle("Distribution Tajima D")
# Center the title
p <- p + ggtitle("Distribution Tajima D") +
  theme(plot.title = element_text(hjust = 0.5)) # Adjust the hjust value for centering
  
# plot distribution
p

##Plot along chromosome using ggplot function
sub1<-(ara_data1[ara_data1$DF1=="DF",])
sub2<-ara_data1[ara_data1$DF1!="DF",]
p<-ggplot(sub2, aes(mid,IT_pi))
p+geom_point(size=2)+geom_point(data=sub1, color="red", size=3)

p<-ggplot(sub2, aes(mid,IT_td))
p+geom_point(size=2)+geom_point(data=sub1, color="red", size=3)+ theme_bw()

####
# Sweden
###

##Kolmogorov smirnov test - compare the distributions of Pi
sub1<-ara_data1$SWE_pi[ara_data1$DF1=="DF"]
sub2<-ara_data1$SWE_pi[ara_data1$DF1!="DF"]
ks.test(sub1, sub2)###difference is very significant if windows are small, otherwise not. 

#Draw Density plot "Pi"
plot(density(log(ara_data$SWE_pi)), main="Distribution log Pi")
lines(density(log(ara_data$SWE_pi[ara_data$DF=="DF"])), col="red")

##Kolmogorov smirnov test - compare the distributions of Tajima's D
sub1<-ara_data$SWE_td[ara_data$DF=="DF"]
sub2<-ara_data$SWE_td[ara_data$DF!="DF"]
ks.test(sub1, sub2)

# Draw Density plots "Tajima's D"
plot(density((ara_data$SWE_td), na.rm=T), main="Distribution Tajima D")
lines(density((ara_data$SWE_td[ara_data$DF=="DF"]), na.rm = T), col="red")

p<-ggplot(ara_data, aes(x=SWE_td, fill=DF))
p+geom_density(alpha=0.4)

##
# Base R plot
plot(density(ara_data$SWE_td, na.rm = TRUE), main = "Distribution Tajima D")
lines(density(ara_data$SWE_td[ara_data$DF == "DF"], na.rm = TRUE), col = "red")

# ggplot version

# To estimate the lowest value of the x-axis for a density plot
lowest_x <- min(ara_data$SWE_td, na.rm = TRUE)
lowest_x

p <- ggplot(ara_data, aes(x = SWE_td, fill = DF)) +
  geom_density(alpha = 0.4) +
  scale_x_continuous(limits = c(-0.010, max(ara_data$SWE_td, na.rm = TRUE))) +  # Set x-axis limit
  ggtitle("Distribution Tajima D") + theme_bw()
# plot distribution
p

##Plot along chromosome using ggplot function
sub1<-(ara_data[ara_data$DF=="DF",])
sub2<-ara_data[ara_data$DF!="DF",]
p<-ggplot(sub2, aes(mid,SWE_pi))
p+geom_point(size=2)+geom_point(data=sub1, color="red", size=3) + theme_bw()

p<-ggplot(sub2, aes(mid,SWE_td))
p+geom_point(size=2)+geom_point(data=sub1, color="red", size=3) + theme_bw()



# FST IT/SWE

#as_tibble: Coerce lists and matrices to data frames
ara_data1.2 <- as.tibble(data.frame(windows1, fst1))
nrow(windows)
nrow(fst)
(chr1_end-chr1_start)/50
nrow(ara_data1.2)
head(ara_data1.2)
ara_data1.2 %>% dplyr::select(contains("fst")) %>% summarise_all(mean)






### load selected positions from chromosome -> flowering time genes ####

bed2<-read.table("At_defense_only.bed")
head(bed2)

colnames(bed2)<-c("chr", "begin","end")
DF1.2<-vector(length=nrow(ara_data1.2))
ara_data1.2 <- as.tibble(data.frame(windows1, fst1, DF1.2))

for (i in 2:nrow(bed2)){
  ara_data1.2$DF1.2 <- "all"
} #horrible but works

for (i in 2:nrow(bed2)){
  ara_data1.2$DF1.2[which(ara_data1.2$start>bed2$begin[i]&ara_data1.2$stop<bed2$end[i])]<-"DF2"
}

DF1.2<-vector(length=nrow(ara_data1.2))

ara_data1.2$DF1.2<-as.factor(ara_data1.2$DF1.2)
summary(ara_data1.2)


Defense <- ara_data1.2 %>% filter(DF1.2 == "DF2")
Defense_fst <- Defense %>% filter(IT_SWE_fst1 >= 0)

a <- ggplot(Defense_fst, aes(mid/10^6, IT_SWE_fst1)) + geom_line(colour = "red")
a <- a + xlab("Position (Mb)") + ylab(expression(italic(F)[ST]))
a + theme_light()

# plot with ara_data FST (<=0 values not removed) and FLOWER_fst (all flowering time FSTs also with <=0 values not removed)
tip <- ggplot() + 
  geom_line(data=ara_data1.2, aes(mid/10^6, IT_SWE_fst1), colour = "blue") + 
  geom_line(data=Defense, aes(mid/10^6, IT_SWE_fst1), colour="pink")
tip <- tip + xlab("Position (Mb)") + ylab(expression(italic(F)[ST]))
tip + theme_light()

# remove fst values <=0 
ara_d2 <- ara_data1.2 %>% filter(IT_SWE_fst1 >= 0)
ara_d2

# calculate means
mean_fst <- mean(ara_d2$IT_SWE_fst1)
mean_defense <- mean(Defense_fst$IT_SWE_fst1)

ks.test(ara_d2$IT_SWE_fst1, Defense_fst$IT_SWE_fst1) 
#p-value: 1

#outliers 95% quantile
threshold_95 <- quantile(Defense_fst$IT_SWE_fst1[Defense_fst$DF1.2=="DF2"], 0.975, na.rm = T)
Defense_fst <- Defense_fst %>% mutate(outlier_95 = ifelse(Defense_fst$IT_SWE_fst1 > threshold_95, "outlier", "background"))

#outliers 99% quantile
threshold_99 <- quantile(Defense_fst$IT_SWE_fst1[Defense_fst$DF1.2=="DF2"], 0.995, na.rm = T)
Defense_fst <- Defense_fst %>% mutate(outlier_99 = ifelse(Defense_fst$IT_SWE_fst1 > threshold_99, "outlier", "background"))

# plot with ara data and Defense_fst (all fst values below 0 removed)
top <- ggplot() + 
  geom_point(data=ara_d2, aes(mid/10^6, IT_SWE_fst1), colour = "lightblue") + 
  geom_point(data=Defense_fst, aes(mid/10^6, IT_SWE_fst1), colour="blue") +
  geom_point(data=Defense_fst[Defense_fst$outlier_95 == "outlier",], aes(mid/10^6, IT_SWE_fst1), color="orange") +
  geom_point(data=Defense_fst[Defense_fst$outlier_99 == "outlier",], aes(mid/10^6, IT_SWE_fst1), color="red") +
  geom_hline(yintercept = mean_fst) +
  geom_hline(yintercept = mean_defense, colour="orange")

top <- top + xlab("Position (Mb)") + ylab(expression(italic(F)[ST]))
top + theme_light()

#################################################################################
# Draw Density plots

plot(density((ara_d2$IT_SWE_fst1), na.rm=T), main="Distribution FST", )
lines(density((Defense_fst$IT_SWE_fst1[Defense_fst$DF1.2=="DF2"]), na.rm = T), col="red")

p<-ggplot(ara_d2, aes(x=IT_SWE_fst1, fill=DF1.2))
p+geom_density(alpha=0.4)


# To estimate the lowest value of the x-axis for a density plot
lowest_x <- min(ara_d2$IT_SWE_fst1, na.rm = TRUE)
lowest_x

p <- ggplot(ara_d2, aes(x = IT_SWE_fst1, fill = DF1.2)) +
  geom_density(alpha = 0.4) +
  scale_x_continuous(limits = c(-0.05, max(ara_d2$IT_SWE_fst1, na.rm = TRUE))) +  # Set x-axis limit
  ggtitle("Distribution ESP-SWE FST")
# plot distribution
p

# ggplot version with log scale for x-axis
p <- ggplot(ara_d2, aes(x = IT_SWE_fst1, fill = DF1.2)) +
  geom_density(alpha = 0.4) +
  scale_x_log10() +  # Set log scale for x-axis
  ggtitle("Distribution IT-SWE FST") + 
  theme_bw()

# plot distribution
p

############################################################
### To check where these FST outliers are in the genome ####
###                                                     ####
############################################################

# Outliers 95% quantile whole genome
threshold_95 <- quantile(ara_d2$IT_SWE_fst1, 0.95, na.rm = T)
ara_d2 <- ara_d2 %>% mutate(outlier_95 = ifelse(ara_d2$IT_SWE_fst1 > threshold_95, "outlier", "background"))
ara_d2

# Outliers 99% quantile whole genome
threshold_99 <- quantile(ara_d2$IT_SWE_fst1, 0.99, na.rm = T)
ara_d2 <- ara_d2 %>% mutate(outlier_99 = ifelse(ara_d2$IT_SWE_fst1 > threshold_99, "outlier", "background"))
ara_d2

out <- ara_d2 %>% filter(outlier_95 == "outlier")
out2 <- out %>% filter(DF1.2 == "DF2")
out2 #no 99% outliers but 2 95% outliers which both correspond to sucrose synthase 3

print(out2,n=50)













