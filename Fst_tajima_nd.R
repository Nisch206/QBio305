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
library(cowplot)

#Stampp to calculate FST between populations
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


#calculate genetic distance between individuals - nei's

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


# FST calculation between populations

at.VCF <- read.vcfR("group_3_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz")


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

#To check total number of sites
At_Chr1@n.sites # 28629727
At_Chr2@n.sites # 9413286
At_Chr3@n.sites # 6689001
At_Chr4@n.sites # 18059106
At_Chr5@n.sites # 24393086

##Deine populations in your dataset

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

# set chromosome size
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


# make a sliding window dataset
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
At_data2 <- as_tibble(data.frame(windows2, nd2, fst2, dxy2))
At_data3 <- as_tibble(data.frame(windows3, nd3, fst3, dxy3))
At_data4 <- as_tibble(data.frame(windows4, nd4, fst4, dxy4))
At_data5 <- as_tibble(data.frame(windows5, nd5, fst5, dxy5))

#Visualizing the data - distributions

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

pi_g1$log_pi <- log10(pi_g1$pi)
pi_g2$log_pi <- log10(pi_g2$pi)
pi_g3$log_pi <- log10(pi_g3$pi)
pi_g4$log_pi <- log10(pi_g4$pi)
pi_g5$log_pi <- log10(pi_g5$pi)

wilcox_test_result1 <- wilcox.test(log_pi ~ populations, data = pi_g1)
wilcox_test_result2 <- wilcox.test(log_pi ~ populations, data = pi_g2)
wilcox_test_result3 <- wilcox.test(log_pi ~ populations, data = pi_g3)
wilcox_test_result4 <- wilcox.test(log_pi ~ populations, data = pi_g4)
wilcox_test_result5 <- wilcox.test(log_pi ~ populations, data = pi_g5)


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
  ggtitle("Chromosome 1") +
  coord_cartesian(xlim = c(0, 10))

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
hs_g2$log_value <- log10(hs_g2$value)
hs_g3$log_value <- log10(hs_g3$value)
hs_g4$log_value <- log10(hs_g4$value)
hs_g5$log_value <- log10(hs_g5$value)

# rearrange everything so FST came at the top, pi beneath it and then finally, d_XY_
# first make a factor
x1 <- factor(hs_g1$stat)
x2 <- factor(hs_g2$stat)
x3 <- factor(hs_g3$stat)
x4 <- factor(hs_g4$stat)
x5 <- factor(hs_g5$stat)

# then reorder the levels
x1 <- factor(x1, levels(x1)[c(3, 1, 4, 2)])
x2 <- factor(x2, levels(x2)[c(3, 1, 4, 2)])
x3 <- factor(x3, levels(x3)[c(3, 1, 4, 2)])
x4 <- factor(x4, levels(x4)[c(3, 1, 4, 2)])
x5 <- factor(x5, levels(x5)[c(3, 1, 4, 2)])

# add to data.frame
hs_g1$stat <- x1
hs_g2$stat <- x2
hs_g3$stat <- x3
hs_g4$stat <- x4
hs_g5$stat <- x5

# construct a plot with facets of fst, pi and dxy for each chromosome

a1 <- ggplot(hs_g1, aes(mid/10^6, value, colour = stat)) + geom_line() +
  ggtitle("Fst and Nucleotide diversity on Chromosome 1") +
  facet_grid(stat~., scales = "free_y") +
  xlab("Position (Mb)") +
  theme_light() + theme(legend.position = "none") +
  theme(plot.title = element_text(face = "bold"))
a1

a2 <- ggplot(hs_g2, aes(mid/10^6, value, colour = stat)) + geom_line() +
  ggtitle("Fst and Nucleotide diversity on Chromosome 2") +
  facet_grid(stat~., scales = "free_y") +
  xlab("Position (Mb)") +
  theme_light() + theme(legend.position = "none") +
  theme(plot.title = element_text(face = "bold"))
a2

a3 <- ggplot(hs_g3, aes(mid/10^6, value, colour = stat)) + geom_line() +
  ggtitle("Fst and Nucleotide diversity on Chromosome 3") +
  facet_grid(stat~., scales = "free_y") +
  xlab("Position (Mb)") +
  theme_light() + theme(legend.position = "none") +
  theme(plot.title = element_text(face = "bold"))
a3

a4 <- ggplot(hs_g4, aes(mid/10^6, value, colour = stat)) + geom_line() +
  ggtitle("Fst and Nucleotide diversity on Chromosome 4") +
  facet_grid(stat~., scales = "free_y") +
  xlab("Position (Mb)") +
  theme_light() + theme(legend.position = "none") +
  theme(plot.title = element_text(face = "bold"))
a4

a5 <- ggplot(hs_g5, aes(mid/10^6, value, colour = stat)) + geom_line() +
  ggtitle("Fst and Nucleotide diversity on Chromosome 5") +
  facet_grid(stat~., scales = "free_y") +
  xlab("Position (Mb)") +
  theme_light() + theme(legend.position = "none") +
  theme(plot.title = element_text(face = "bold"))
a5


#### calculate neutrality statistics####
At_sw1 <- neutrality.stats(At_sw1)
At_sw2 <- neutrality.stats(At_sw2)
At_sw3 <- neutrality.stats(At_sw3)
At_sw4 <- neutrality.stats(At_sw4)
At_sw5 <- neutrality.stats(At_sw5)

#extract Tajma's D
td1 <- At_sw1@Tajima.D/100
td2 <- At_sw2@Tajima.D/100
td3 <- At_sw3@Tajima.D/100
td4 <- At_sw4@Tajima.D/100
td5 <- At_sw5@Tajima.D/100

# set population names
colnames(td1) <- paste0(pops, "_td")
colnames(td2) <- paste0(pops, "_td")
colnames(td3) <- paste0(pops, "_td")
colnames(td4) <- paste0(pops, "_td")
colnames(td5) <- paste0(pops, "_td")

#Coerce lists and matrices to data frames
ara_data1 <- as.tibble(data.frame(windows1, td1,nd1))
ara_data2 <- as.tibble(data.frame(windows2, td2,nd2))
ara_data3 <- as.tibble(data.frame(windows3, td3,nd3))
ara_data4 <- as.tibble(data.frame(windows4, td4,nd4))
ara_data5 <- as.tibble(data.frame(windows5, td5,nd5))

# TajimaD

# Select data of interest
hs_td1 <- ara_data1 %>%
  dplyr::select(mid, IT_td, SWE_td)
hs_td2 <- ara_data2 %>%
  dplyr::select(mid, IT_td, SWE_td)
hs_td3 <- ara_data3 %>%
  dplyr::select(mid, IT_td, SWE_td)
hs_td4 <- ara_data4 %>%
  dplyr::select(mid, IT_td, SWE_td)
hs_td5 <- ara_data5 %>%
  dplyr::select(mid, IT_td, SWE_td)

# Extract td values for IT and SWE populations
td_chr1 <- hs_td1 %>%
  summarize(mean_IT_td = mean(na.omit(IT_td)), mean_SWE_td = mean(na.omit(SWE_td)))
td_chr1

td_chr2 <- hs_td2 %>%
  summarize(mean_IT_td = mean(na.omit(IT_td)), mean_SWE_td = mean(na.omit(SWE_td)))
td_chr2

td_chr3 <- hs_td3 %>%
  summarize(mean_IT_td = mean(na.omit(IT_td)), mean_SWE_td = mean(na.omit(SWE_td)))
td_chr3

td_chr4 <- hs_td4 %>%
  summarize(mean_IT_td = mean(na.omit(IT_td)), mean_SWE_td = mean(na.omit(SWE_td)))
td_chr4

td_chr5 <- hs_td5 %>%
  summarize(mean_IT_td = mean(na.omit(IT_td)), mean_SWE_td = mean(na.omit(SWE_td)))
td_chr5

# Use gather to rearrange everything
hs_td_g1 <- gather(hs_td1, -mid, key = "populations", value = "value")
hs_td_g2 <- gather(hs_td2, -mid, key = "populations", value = "value")
hs_td_g3 <- gather(hs_td3, -mid, key = "populations", value = "value")
hs_td_g4 <- gather(hs_td4, -mid, key = "populations", value = "value")
hs_td_g5 <- gather(hs_td5, -mid, key = "populations", value = "value")

# Reorder the levels of the populations factor
hs_td_g1$populations <- factor(hs_td_g1$populations, levels = c("IT_td", "SWE_td"))
hs_td_g2$populations <- factor(hs_td_g2$populations, levels = c("IT_td", "SWE_td"))
hs_td_g3$populations <- factor(hs_td_g3$populations, levels = c("IT_td", "SWE_td"))
hs_td_g4$populations <- factor(hs_td_g4$populations, levels = c("IT_td", "SWE_td"))
hs_td_g5$populations <- factor(hs_td_g5$populations, levels = c("IT_td", "SWE_td"))

# Take the logarithm of the value variable
hs_td_g1$log_value <- log10(hs_td_g1$value)
hs_td_g2$log_value <- log10(hs_td_g2$value)
hs_td_g3$log_value <- log10(hs_td_g3$value)
hs_td_g4$log_value <- log10(hs_td_g4$value)
hs_td_g5$log_value <- log10(hs_td_g5$value)

# Construct plots with facets

a_td1 <- ggplot(hs_td_g1, aes(mid / 10^6, log_value, colour = populations)) +  geom_line() +
  facet_wrap(~populations, scales = "free_y", ncol = 1) +
  xlab("Position (Mb)") +
  theme_light()
a_td1

a_td2 <- ggplot(hs_td_g2, aes(mid / 10^6, log_value, colour = populations)) +  geom_line() +
  facet_wrap(~populations, scales = "free_y", ncol = 1) +
  xlab("Position (Mb)") +
  theme_light()
a_td2

a_td3 <- ggplot(hs_td_g3, aes(mid / 10^6, log_value, colour = populations)) +  geom_line() +
  facet_wrap(~populations, scales = "free_y", ncol = 1) +
  xlab("Position (Mb)") +
  theme_light()
a_td3

a_td4 <- ggplot(hs_td_g4, aes(mid / 10^6, log_value, colour = populations)) +  geom_line() +
  facet_wrap(~populations, scales = "free_y", ncol = 1) +
  xlab("Position (Mb)") +
  theme_light()
a_td4

a_td5 <- ggplot(hs_td_g5, aes(mid / 10^6, log_value, colour = populations)) +  geom_line() +
  facet_wrap(~populations, scales = "free_y", ncol = 1) +
  xlab("Position (Mb)") +
  theme_light()
a_td5

# Boxplots

wilcox_test_result1 <- wilcox.test(log_value ~ populations, data = hs_td_g1 )
wilcox_test_result2 <- wilcox.test(log_value ~ populations, data = hs_td_g2 )
wilcox_test_result3 <- wilcox.test(log_value ~ populations, data = hs_td_g3 )
wilcox_test_result4 <- wilcox.test(log_value ~ populations, data = hs_td_g4 )
wilcox_test_result5 <- wilcox.test(log_value ~ populations, data = hs_td_g5 )


# Create the individual boxplot objects
a1 <- ggplot(hs_td_g1, aes(populations, log_value, fill = populations)) +
  geom_boxplot(color = "black") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_light() +
  xlab(NULL) +
  ylab("Log10(td)") +
  ggtitle("Chromosome 1") +
  annotate("text", x = 2, y = -5.4, label = paste("p =", format.pval(wilcox_test_result1$p.value, digits = 3)))
a1 <- a1 + expand_limits(y = c(-6, max(hs_td_g1$log_value) + 1))


a2 <- ggplot(hs_td_g2, aes(populations, log_value, fill = populations)) +
  geom_boxplot(color = "black") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_light() +
  xlab(NULL) +
  ylab("Log10(td)") +
  ggtitle("Chromosome 2") +
  annotate("text", x = 2, y = -4.5, label = paste("p =", format.pval(wilcox_test_result2$p.value, digits = 3)))
a2 <- a2 + expand_limits(y = c(-5, max(hs_td_g2$log_value) + 1))


a3 <- ggplot(hs_td_g3, aes(populations, log_value, fill = populations)) +
  geom_boxplot(color = "black") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_light() +
  xlab(NULL) +
  ylab("Log10(td)") +
  ggtitle("Chromosome 3") +
  annotate("text", x = 2, y = -3.25, label = paste("p =", format.pval(wilcox_test_result3$p.value, digits = 3)))
a3 <- a3 + expand_limits(y = c(-3.5, max(hs_td_g3$log_value) + 1))


a4 <- ggplot(hs_td_g4, aes(populations, log_value, fill = populations)) +
  geom_boxplot(color = "black") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_light() +
  xlab(NULL) +
  ylab("Log10(td)") +
  ggtitle("Chromosome 4") +
  annotate("text", x = 2, y = -4.5, label = paste("p =", format.pval(wilcox_test_result4$p.value, digits = 3)))
a4 <- a4 + expand_limits(y = c(-5, max(hs_td_g4$log_value) + 1))

a5 <- ggplot(hs_td_g5, aes(populations, log_value,fill = populations)) +
  geom_boxplot(color = "black") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_light() +
  xlab(NULL) +
  ylab("Log10(td)") +
  ggtitle("Chromosome 5") +
  annotate("text", x = 2, y = -5.5, label = paste("p =", format.pval(wilcox_test_result5$p.value, digits = 3)))
a5 <- a5 + expand_limits(y = c(-6, max(hs_td_g5$log_value) + 1))


# Combine the individual boxplots into a single plot
combined_plots <- plot_grid(a1, a2, a3, a4, a5, labels = "AUTO", nrow = 2)

overall_title <- ggdraw() +
  draw_label("Tajima's D as a function of sample location", size = 20, hjust = 0.5) +
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



#  Plotting FST along Chromosome
#  Visualizing FST outliers

# Chr. 1
head(At_data1)

# Provide summary statistics for the IT_SW_fst column
summary(At_data1$IT_SWE_fst1)

#Filter out rows where IT_SW_fst is not equal to zero and create a  histogram using the non-zero values.
non_zero_data1 <- At_data1[At_data1$IT_SWE_fst1 != 0, ]

# Plot histogram of non-zero IT_SW_fst values
ggplot(non_zero_data1, aes(x = IT_SWE_fst1)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black",  alpha = 0.7) +
  labs(title = "Histogram of non-zero IT_SW_fst", x = "IT_SW_fst", y  = "Frequency")

# Filter out zero and negative values
positive_data1 <- At_data1[At_data1$IT_SWE_fst1 > 0, ]

# Display the first few rows of the positive_data
head(positive_data1)


# Plot histogram of positive IT_SW_fst values
ggplot(positive_data1, aes(x = IT_SWE_fst1)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black",  alpha = 0.7) +
  labs(title = "Histogram of positive IT_SW_fst", x = "IT_SW_fst", y  = "Frequency")

# Set threshold for 95% and 99%
threshold_95_1 <- quantile(positive_data1$IT_SWE_fst1, 0.95, na.rm = TRUE)
threshold_99_1 <- quantile(positive_data1$IT_SWE_fst1, 0.99, na.rm = TRUE)

# Plot histogram with threshold marked
ggplot(positive_data1, aes(x = IT_SWE_fst1)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black",  alpha = 0.7) +
  labs(title = "Histogram of positive IT_SW_fst", x = "IT_SW_fst", y  = "Frequency") +
  geom_vline(xintercept = threshold_95_1, colour = "orange", linetype =  "dashed", size = 1)+
  geom_vline(xintercept = threshold_99_1, colour = "red", linetype =  "dashed", size = 1)

# Add outlier columns based on thresholds
positive_data1$outlier_95 <- ifelse(positive_data1$IT_SWE_fst1 >  threshold_95_1, "Outlier", "Non-outlier")
positive_data1$outlier_99 <- ifelse(positive_data1$IT_SWE_fst1 >  threshold_99_1, "Outlier", "Non-outlier")

# Display the first few rows of the updated positive_data
head(positive_data1)

# Create a FST scatter plot along the chromosome
a <- ggplot(positive_data1, aes(mid, IT_SWE_fst1)) + geom_point()
a <- a + xlab("Position (Mb)") + ylab(expression(italic(F)[ST]))
a <- a + theme_light()
a

# Estimate mean FST value
mean_fst_Chr1 <- mean(positive_data1$IT_SWE_fst1)
mean_fst_Chr1

# Filter rows where the 'outlier_95' column has the value "Outlier"
out_95_1 <- filter(positive_data1, outlier_95 == "Outlier")

# Filter rows where the 'outlier_99' column has the value "Outlier"
out_99_1 <- filter(positive_data1, outlier_99 == "Outlier")

# Filter rows where the 'outlier_99' column has the value "Outlier"
out_99_1 <- filter(positive_data1, outlier_99 == "Outlier")

# Print the first 20 rows of the 'out_95' data frame
print(out_95_1, n = 20)


# Print the first 20 rows of the 'out_99' data frame
print(out_99_1, n = 20)

# Adding a chromosome column to the data
out_99_1$chromosome <- 'Chr1'
out_95_1$chromosome <- 'Chr1'

###############################

#Chr. 2

head(At_data2)

# Provide structure information about the dataset
str(At_data2)

# Provide summary statistics for the IT_SW_fst column
summary(At_data2$IT_SWE_fst2)

# Filter out rows where IT_SW_fst is not equal to zero and create a  histogram using the non-zero values.
non_zero_data2 <- At_data2[At_data2$IT_SWE_fst2 != 0, ]

# Filter out rows where IT_SW_fst is not equal to zero and create a  histogram using the non-zero values.
non_zero_data2 <- At_data2[At_data2$IT_SWE_fst2 != 0, ]

# Plot histogram of non-zero IT_SW_fst values
ggplot(non_zero_data2, aes(x = IT_SWE_fst2)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black",  alpha = 0.7) +
  labs(title = "Histogram of non-zero IT_SW_fst", x = "IT_SW_fst", y  = "Frequency")

# Filter out zero and negative values
positive_data2 <- At_data2[At_data2$IT_SWE_fst2 > 0, ]

# Display the first few rows of the positive_data
head(positive_data2)

# Plot histogram of positive IT_SW_fst values
ggplot(positive_data2, aes(x = IT_SWE_fst2)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black",  alpha = 0.7) +
  labs(title = "Histogram of positive IT_SW_fst", x = "IT_SW_fst", y  = "Frequency")

# Set threshold for 95% and 99%
threshold_95_2 <- quantile(positive_data2$IT_SWE_fst2, 0.95, na.rm = TRUE)
threshold_99_2 <- quantile(positive_data2$IT_SWE_fst2, 0.99, na.rm = TRUE)

# Plot histogram with threshold marked
ggplot(positive_data2, aes(x = IT_SWE_fst2)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black",  alpha = 0.7) +
  labs(title = "Histogram of positive IT_SW_fst", x = "IT_SW_fst", y  = "Frequency") +
  geom_vline(xintercept = threshold_95_2, colour = "orange", linetype =  "dashed", size = 1)+
  geom_vline(xintercept = threshold_99_2, colour = "red", linetype =  "dashed", size = 1)

# Add outlier columns based on thresholds
positive_data2$outlier_95 <- ifelse(positive_data2$IT_SWE_fst2 >  threshold_95_2, "Outlier", "Non-outlier")
positive_data2$outlier_99 <- ifelse(positive_data2$IT_SWE_fst2 >  threshold_99_2, "Outlier", "Non-outlier")

# Display the first few rows of the updated positive_data
head(positive_data2)

# Create a FST scatter plot along the chromosome
a <- ggplot(positive_data2, aes(mid, IT_SWE_fst2)) + geom_point()
a <- a + xlab("Position (Mb)") + ylab(expression(italic(F)[ST]))
a <- a + theme_light()
a

# Estimate mean FST value
mean_fst_Chr2 <- mean(positive_data2$IT_SWE_fst2)
mean_fst_Chr2

# Filter rows where the 'outlier_95' column has the value "Outlier"
out_95_2 <- filter(positive_data2, outlier_95 == "Outlier")

# Filter rows where the 'outlier_99' column has the value "Outlier"
out_99_2 <- filter(positive_data2, outlier_99 == "Outlier")

# Print the first 20 rows of the 'out_95' data frame
print(out_95_2, n = 20)

# Print the first 20 rows of the 'out_99' data frame
print(out_99_2, n = 20)

# Adding a chromosome column to the data
out_99_2$chromosome <- 'Chr2'
out_95_2$chromosome <- 'Chr2'

###########################

# Chr.3

# Provide structure information about the dataset
str(At_data3)

summary(At_data3$IT_SWE_fst3)

# Filter out rows where IT_SW_fst is not equal to zero and create a  histogram using the non-zero values.
non_zero_data3 <- At_data3[At_data3$IT_SWE_fst3 != 0, ]

# Plot histogram of non-zero IT_SW_fst values
ggplot(non_zero_data3, aes(x = IT_SWE_fst3)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black",  alpha = 0.7) +
  labs(title = "Histogram of non-zero IT_SW_fst", x = "IT_SW_fst", y  = "Frequency")

# Filter out zero and negative values
positive_data3 <- At_data3[At_data3$IT_SWE_fst3 > 0, ]

# Display the first few rows of the positive_data
head(positive_data3)

# Plot histogram of positive IT_SW_fst values
ggplot(positive_data3, aes(x = IT_SWE_fst3)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black",  alpha = 0.7) +
  labs(title = "Histogram of positive IT_SW_fst", x = "IT_SW_fst", y  = "Frequency")

# Set threshold for 95% and 99%
threshold_95_3 <- quantile(positive_data3$IT_SWE_fst3, 0.95, na.rm = TRUE)
threshold_99_3 <- quantile(positive_data3$IT_SWE_fst3, 0.99, na.rm = TRUE)

# Plot histogram with threshold marked
ggplot(positive_data3, aes(x = IT_SWE_fst3)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black",  alpha = 0.7) +
  labs(title = "Histogram of positive IT_SW_fst", x = "IT_SW_fst", y  = "Frequency") +
  geom_vline(xintercept = threshold_95_3, colour = "orange", linetype =  "dashed", size = 1)+
  geom_vline(xintercept = threshold_99_3, colour = "red", linetype =  "dashed", size = 1)

# Add outlier columns based on thresholds
positive_data3$outlier_95 <- ifelse(positive_data3$IT_SWE_fst3 >  threshold_95_3, "Outlier", "Non-outlier")
positive_data3$outlier_99 <- ifelse(positive_data3$IT_SWE_fst3 >  threshold_99_3, "Outlier", "Non-outlier")

# Display the first few rows of the updated positive_data
head(positive_data3)

# Create a FST scatter plot along the chromosome
a <- ggplot(positive_data3, aes(mid, IT_SWE_fst3)) + geom_point()
a <- a + xlab("Position (Mb)") + ylab(expression(italic(F)[ST]))
a <- a + theme_light()
a

# Estimate mean FST value
mean_fst_Chr3 <- mean(positive_data3$IT_SWE_fst3)
mean_fst_Chr3

# Filter rows where the 'outlier_95' column has the value "Outlier"
out_95_3 <- filter(positive_data3, outlier_95 == "Outlier")

# Filter rows where the 'outlier_99' column has the value "Outlier"
out_99_3 <- filter(positive_data3, outlier_99 == "Outlier")

# Print the first 20 rows of the 'out_95' data frame
print(out_95_3, n = 20)

# Print the first 20 rows of the 'out_99' data frame
print(out_99_3, n = 20)

# Adding a chromosome column to the data
out_99_3$chromosome <- 'Chr3'
out_95_3$chromosome <- 'Chr3'

##########################################

# Chr.4

# Display the first few rows of the dataset
head(At_data4)

# Provide structure information about the dataset
str(At_data4)

# Provide summary statistics for the IT_SW_fst column
summary(At_data4$IT_SWE_fst4)

# Filter out rows where IT_SW_fst is not equal to zero and create a  histogram using the non-zero values.
non_zero_data4 <- At_data4[At_data4$IT_SWE_fst4 != 0, ]

# Plot histogram of non-zero IT_SW_fst values
ggplot(non_zero_data4, aes(x = IT_SWE_fst4)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black",  alpha = 0.7) +
  labs(title = "Histogram of non-zero IT_SW_fst", x = "IT_SW_fst", y  = "Frequency")

# Filter out zero and negative values
positive_data4 <- At_data4[At_data4$IT_SWE_fst4 > 0, ]

# Display the first few rows of the positive_data
head(positive_data4)

# Plot histogram of positive IT_SW_fst values
ggplot(positive_data4, aes(x = IT_SWE_fst4)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black",  alpha = 0.7) +
  labs(title = "Histogram of positive IT_SW_fst", x = "IT_SW_fst", y  = "Frequency")

# Set threshold for 95% and 99%
threshold_95_4 <- quantile(positive_data4$IT_SWE_fst4, 0.95, na.rm = TRUE)
threshold_99_4 <- quantile(positive_data4$IT_SWE_fst4, 0.99, na.rm = TRUE)

# Plot histogram with threshold marked
ggplot(positive_data4, aes(x = IT_SWE_fst4)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black",  alpha = 0.7) +
  labs(title = "Histogram of positive IT_SW_fst", x = "IT_SW_fst", y  = "Frequency") +
  geom_vline(xintercept = threshold_95_4, colour = "orange", linetype =  "dashed", size = 1)+
  geom_vline(xintercept = threshold_99_4, colour = "red", linetype =  "dashed", size = 1)

# Add outlier columns based on thresholds
positive_data4$outlier_95 <- ifelse(positive_data4$IT_SWE_fst4 >  threshold_95_4, "Outlier", "Non-outlier")
positive_data4$outlier_99 <- ifelse(positive_data4$IT_SWE_fst4 >  threshold_99_4, "Outlier", "Non-outlier")

# Display the first few rows of the updated positive_data
head(positive_data4)

# Create a FST scatter plot along the chromosome
a <- ggplot(positive_data4, aes(mid, IT_SWE_fst4)) + geom_point()
a <- a + xlab("Position (Mb)") + ylab(expression(italic(F)[ST]))
a <- a + theme_light()
a

# Estimate mean FST value
mean_fst_Chr4 <- mean(positive_data4$IT_SWE_fst4)
mean_fst_Chr4

# Filter rows where the 'outlier_95' column has the value "Outlier"
out_95_4 <- filter(positive_data4, outlier_95 == "Outlier")

# Filter rows where the 'outlier_99' column has the value "Outlier"
out_99_4 <- filter(positive_data4, outlier_99 == "Outlier")

# Print the first 20 rows of the 'out_95' data frame
print(out_95_4, n = 20)

# Print the first 20 rows of the 'out_99' data frame
print(out_99_4, n = 20)

# Adding a chromosome column to the data
out_99_4$chromosome <- 'Chr4'
out_95_4$chromosome <- 'Chr4'

#################################

# Chr.5

# Provide structure information about the dataset
str(At_data5)

summary(At_data5$IT_SWE_fst5)

# Filter out rows where IT_SW_fst is not equal to zero and create a  histogram using the non-zero values.
non_zero_data5 <- At_data5[At_data5$IT_SWE_fst5 != 0, ]

# Plot histogram of non-zero IT_SW_fst values
ggplot(non_zero_data5, aes(x = IT_SWE_fst5)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black",  alpha = 0.7) +
  labs(title = "Histogram of non-zero IT_SW_fst", x = "IT_SW_fst", y  = "Frequency")

# Filter out zero and negative values
positive_data5 <- At_data5[At_data5$IT_SWE_fst5 > 0, ]

# Display the first few rows of the positive_data
head(positive_data5)

# Plot histogram of positive IT_SW_fst values
ggplot(positive_data5, aes(x = IT_SWE_fst5)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black",  alpha = 0.7) +
  labs(title = "Histogram of positive IT_SW_fst", x = "IT_SW_fst", y  = "Frequency")

# Set threshold for 95% and 99%
threshold_95_5 <- quantile(positive_data5$IT_SWE_fst5, 0.95, na.rm = TRUE)
threshold_99_5 <- quantile(positive_data5$IT_SWE_fst5, 0.99, na.rm = TRUE)

# Plot histogram with threshold marked
ggplot(positive_data5, aes(x = IT_SWE_fst5)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black",  alpha = 0.7) +
  labs(title = "Histogram of positive IT_SW_fst", x = "IT_SW_fst", y  = "Frequency") +
  geom_vline(xintercept = threshold_95_5, colour = "orange", linetype =  "dashed", size = 1)+
  geom_vline(xintercept = threshold_99_5, colour = "red", linetype =  "dashed", size = 1)

# Add outlier columns based on thresholds
positive_data5$outlier_95 <- ifelse(positive_data5$IT_SWE_fst5 >  threshold_95_5, "Outlier", "Non-outlier")
positive_data5$outlier_99 <- ifelse(positive_data5$IT_SWE_fst5 >  threshold_99_5, "Outlier", "Non-outlier")

# Display the first few rows of the updated positive_data
head(positive_data5)

# Create a FST scatter plot along the chromosome
a <- ggplot(positive_data5, aes(mid, IT_SWE_fst5)) + geom_point()
a <- a + xlab("Position (Mb)") + ylab(expression(italic(F)[ST]))
a <- a + theme_light()
a

# Estimate mean FST value
mean_fst_Chr5 <- mean(positive_data5$IT_SWE_fst5)
mean_fst_Chr5

# Filter rows where the 'outlier_95' column has the value "Outlier"
out_95_5 <- filter(positive_data5, outlier_95 == "Outlier")
# Filter rows where the 'outlier_99' column has the value "Outlier"
out_99_5 <- filter(positive_data5, outlier_99 == "Outlier")
# Print the first 20 rows of the 'out_95' data frame
print(out_95_5, n = 20)

# Print the first 20 rows of the 'out_99' data frame
print(out_99_5, n = 20)

# Adding a chromosome column to the "out_99" data
out_99_5$chromosome <- 'Chr5'
out_95_5$chromosome <- 'Chr5'


## All scatterplots
# Create a plot with marked outliers, FST mean line, and customized legend

top1 <- ggplot(positive_data1, aes(mid/10^6, IT_SWE_fst1)) +
  geom_point(aes(colour = ifelse(outlier_99 == "Outlier", "Outlier 99%",
                                 ifelse(outlier_95 == "Outlier", "Outlier 95%", "Non-outlier"))),
             size = 1.5) +
  scale_colour_manual(values = c("Non-outlier" = "black", "Outlier 95%" = "orange", "Outlier 99%" = "red")) +
  geom_hline(yintercept = mean_fst_Chr1, colour = "blue") +
  xlab("Position on Chromosome (Mb)") + ylab(expression(italic(F)[ST])) +
  theme_light() +
  ggtitle("Chromosome 1") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

top2 <- ggplot(positive_data2, aes(mid/10^6, IT_SWE_fst2)) +
  geom_point(aes(colour = ifelse(outlier_99 == "Outlier", "Outlier 99%",
                                 ifelse(outlier_95 == "Outlier", "Outlier 95%", "Non-outlier"))),
             size = 1.5) +
  scale_colour_manual(values = c("Non-outlier" = "black", "Outlier 95%" = "orange", "Outlier 99%" = "red")) +
  geom_hline(yintercept = mean_fst_Chr2, colour = "blue") +
  xlab("Position on Chromosome (Mb)") + ylab(expression(italic(F)[ST])) +
  theme_light() +
  ggtitle("Chromosome 2") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

top3 <- ggplot(positive_data3, aes(mid/10^6, IT_SWE_fst3)) +
  geom_point(aes(colour = ifelse(outlier_99 == "Outlier", "Outlier 99%",
                                 ifelse(outlier_95 == "Outlier", "Outlier 95%", "Non-outlier"))),
             size = 1.5) +
  scale_colour_manual(values = c("Non-outlier" = "black", "Outlier 95%" = "orange", "Outlier 99%" = "red")) +
  geom_hline(yintercept = mean_fst_Chr3, colour = "blue") +
  xlab("Position on Chromosome (Mb)") + ylab(expression(italic(F)[ST])) +
  theme_light() +
  ggtitle("Chromosome 3") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

top4 <- ggplot(positive_data4, aes(mid/10^6, IT_SWE_fst4)) +
  geom_point(aes(colour = ifelse(outlier_99 == "Outlier", "Outlier 99%",
                                 ifelse(outlier_95 == "Outlier", "Outlier 95%", "Non-outlier"))),
             size = 1.5) +
  scale_colour_manual(values = c("Non-outlier" = "black", "Outlier 95%" = "orange", "Outlier 99%" = "red")) +
  geom_hline(yintercept = mean_fst_Chr4, colour = "blue") +
  xlab("Position on Chromosome (Mb)") + ylab(expression(italic(F)[ST])) +
  theme_light() +
  ggtitle("Chromosome 4") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

top5 <- ggplot(positive_data5, aes(mid/10^6, IT_SWE_fst5)) +
  geom_point(aes(colour = ifelse(outlier_99 == "Outlier", "Outlier 99%",
                                 ifelse(outlier_95 == "Outlier", "Outlier 95%", "Non-outlier"))),
             size = 1.5) +
  scale_colour_manual(values = c("Non-outlier" = "black", "Outlier 95%" = "orange", "Outlier 99%" = "red")) +
  geom_hline(yintercept = mean_fst_Chr5, colour = "blue") +
  xlab("Position on Chromosome (Mb)") + ylab(expression(italic(F)[ST])) +
  theme_light() +
  ggtitle("Chromosome 5") +
  guides(color = guide_legend(override.aes = list(size = 3, shape = 16), title = "Outlier Status")) +
  theme(plot.title = element_text(hjust = 0.5))

# Extract the legend
legend <- get_legend(top5)

# Remove the legend from plot5 to prevent skewing
top5 <- top5 + theme(legend.position = "none")

# Combine the individual scatter plots into a single plot
combined_plots <- plot_grid(top1, top2, top3, top4, top5, legend, labels = "AUTO", nrow = 2)

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


###### Asignment to Genes ##########
library(GenomicRanges)

# Reference genome file
gff <-read.delim("/Users/julia/Desktop/3. Semester/QGenetics/GFF_final_with_annotation_2.gff", header = FALSE)

# Add column names
colnames(gff) <- c("model_name", "chromosome", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
# Remove first row that is without data
gff <- gff[-1,]

# Transform the data type of the start and end positions into numeric data
gff$start <- as.numeric(as.character(gff$start))
gff$end <- as.numeric(as.character(gff$end))

# Filter to exclude UTR, intron, and exon types
filtered_annotations <- subset(gff, !(type %in% c("UTR", "intron", "exon", "five_prime_UTR", "three_prime_UTR", "chromosome")))

# Preperation of outlier data
# 95% outlier
# Equalize all column names
colnames(out_95_1)[6] <- "IT_SWE_fst"
colnames(out_95_1)[7] <- "IT_SWE_dxy"

colnames(out_95_2) <- colnames(out_95_1)
colnames(out_95_3) <- colnames(out_95_1)
colnames(out_95_4) <- colnames(out_95_1)
colnames(out_95_5) <- colnames(out_95_1)

# Combine them into one data frame
all_outs_95 <- rbind(out_95_1, out_95_2, out_95_3, out_95_4, out_95_5)
# Remove rows with 99% outliers, since we look at them afterwards
all_outs_95 <- all_outs_95[all_outs_95$outlier_99 != "Outlier", ]

# 99% outlier
# Equalize all column names
colnames(out_99_1)[6] <- "IT_SWE_fst"
colnames(out_99_1)[7] <- "IT_SWE_dxy"

colnames(out_99_2) <- colnames(out_99_1)
colnames(out_99_3) <- colnames(out_99_1)
colnames(out_99_4) <- colnames(out_99_1)
colnames(out_99_5) <- colnames(out_99_1)

# Combine them into one data frame
all_outs <- rbind(out_99_1, out_99_2, out_99_3, out_99_4, out_99_5)

# Search for overlaps
# Create a GRanges object for the annotation file
gene_ranges <- GRanges(
  seqnames = filtered_annotation$chromosome,
  ranges = IRanges(start = filtered_annotation$start, end = filtered_annotation$end)
)


# Create a GRanges object for the outlier regions
# For 95% outliers
outlier_ranges_95 <- GRanges(
  seqnames = all_outs_95$chromosome,
  ranges = IRanges(start = all_outs_95$start, end = all_outs_95$stop)
)

# For 99% outliers
outlier_ranges <- GRanges(
  seqnames = all_outs$chromosome,
  ranges = IRanges(start = all_outs$start, end = all_outs$stop)
)

# Find overlapps of the data files
overlaps_95 <- findOverlaps(outlier_ranges_95, gene_ranges, type = "any", select = "all", ignore.strand = TRUE)
overlaps <- findOverlaps(outlier_ranges, gene_ranges, type = "any", select = "all", ignore.strand = TRUE)

# Extract indices of overlapping genes
overlapping_genes_indices_95 <- subjectHits(overlaps_95)
overlapping_genes_indices <- subjectHits(overlaps)

# Extract information for overlapping genes
overlapping_genes_95 <- filtered_annotation[overlapping_genes_indices_95, ]
overlapping_genes <- filtered_annotation[overlapping_genes_indices, ]
















