
# Load required libraries

library(PopGenome) 
library(vcfR)
library(VariantAnnotation)
library (readr)
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(dunn.test)

# Read the VCF file

FLC <- readVCF("group_3_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz", numcols=89,tid="5", frompos=3173497, topos=3179448,include.unknown = TRUE)

population_info <- read_delim("pop1_sample_pop.txt", delim = "\t") 

# now get the data for the populations
populations <- split(population_info$sample, population_info$pop)

# now set 
FLC <- set.populations(FLC, populations, diploid = T)
##check if it worked
FLC@populations

##To check total number of sites
FLC@n.sites
chr <- 5952

####Setting up sliding windows###

# set window size and window jump
window_size <- 100
window_jump <- 50

# use seq to find the start points of each window
window_start_Chr <- seq(from = 1, to = chr, by = window_jump)
# add the size of the window to each start point 
window_stop_Chr <- window_start_Chr + window_size

# no windows start before the end of chromosome 4
sum(window_start_Chr > chr)
# but some window stop positions do occur past the final point
sum(window_stop_Chr > chr)

# remove windows from the start and stop vectors
window_start_Chr <- window_start_Chr[which(window_stop_Chr < chr)]
window_stop_Chr <- window_stop_Chr[which(window_stop_Chr < chr)]

chr - window_stop_Chr[length(window_stop_Chr)]

# save as a data.frame
windows_Chr <- data.frame(start = window_start_Chr, stop = window_stop_Chr, 
                          mid = window_start_Chr + (window_stop_Chr-window_start_Chr)/2)

# make a sliding window dataset
At_sw_Chr <- sliding.window.transform(FLC, width = 100, jump = 50, type = 2)

# calculate diversity statistics
At_sw_Chr <- diversity.stats(At_sw_Chr, pi = TRUE)


#Next we will calculate FST, which again is very straight forward with a single command.

### calculate diversity statistics
At_sw_Chr <- F_ST.stats(At_sw_Chr, mode = "nucleotide")

#### calculate neutrality statistics####
At_sw_Chr <- neutrality.stats(At_sw_Chr)

####Extracting statistics for visualization####

#extract nucleotide diversity and correct for window size
nd_Chr <- At_sw_Chr@nuc.diversity.within/100


# make population name vector
pops <- c("IT-N", "IT-S","SW-N", "SW-S") 
# set population names
colnames(nd_Chr) <- paste0(pops, "_pi")

# extract fst values
fst_Chr <- t(At_sw_Chr@nuc.F_ST.pairwise)

# extract dxy - pairwise absolute nucleotide diversity
dxy_Chr <- get.diversity(At_sw_Chr, between = T)[[2]]/100

# get column names 
x <- colnames(fst_Chr)
fst_Chr
# Loop through each population and replace the corresponding population name in the column names
for (i in 1:length(pops)) {
  pattern <- paste0("pop", i)
  x <- sub(pattern, pops[i], x)
}

# look at x to confirm the replacement has occurred
print(x)

# replace forward slash
x <- sub("/", "_", x)
# look at x to confirm the replacement has occurred
x

#Make clear these names are for either FST or d_XY_. 
paste0(x, "_fst")
paste0(x, "_dxy")

#Change the column names of our two data matrices
colnames(fst_Chr) <- paste0(x, "_fst")
colnames(dxy_Chr) <- paste0(x, "_dxy")

#extract Tajma's D and set population names
td_Chr <- At_sw_Chr@Tajima.D/100

colnames(td_Chr) <- paste0(pops, "_td")


#Combine the data into a big dataset

At_data_Chr <- as.tibble(data.frame(windows_Chr, td_Chr, nd_Chr, fst_Chr, dxy_Chr))

# select nucleotide diversity data and calculate means
At_data_Chr %>% dplyr::select(contains("pi")) %>% summarise_all(mean)

At_data_Chr

################################

#Nucleotide diversity pi


# select nucleotide diversity data and calculate means
At_data_Chr %>% dplyr::select(contains("pi")) %>% summarise_all(mean)


# To plot this we need to use "gather" on the data

pi_g_Chr <- At_data_Chr %>% dplyr::select(contains("pi")) %>% gather(key = "populations", value = "pi")
pi_g_Chr


pi_g_Chr$log_pi <- log10(pi_g_Chr$pi)

a_pi_Chr <- ggplot(pi_g_Chr, aes(populations, log_pi, fill = populations)) + 
  geom_boxplot(color = "black") +  # Border color
  scale_fill_manual(values = c("red", "orange", "blue", "lightblue")) +  
  theme_light() + 
  xlab(NULL) +
  ggtitle("Nucleotide diversity of FLC gene") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
        legend.text = element_text(size = 14),  
        legend.title = element_text(size = 16)) +
  ylab("Log10(pi)") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) 

a_pi_Chr


# Wilcoxon rank-sum test after filtering the data and two select two populations of your choice
comparison_data_Chr <- pi_g_Chr %>% filter(populations %in% c("IT.S_pi", "SW.N_pi"))

# Perform Wilcoxon rank-sum test
wilcox_test_pi_Chr <- wilcox.test(log_pi ~ populations, data = comparison_data_Chr)

# Print the test result
print(wilcox_test_pi_Chr)

# Kruskal-Wallis test
kruskal_test_pi_Chr <- kruskal.test(log_pi ~ populations, data = pi_g_Chr)

# Print the test result
print(kruskal_test_pi_Chr)

# generate boxplot again
a_pi_Chr

# Add Kruskal-Wallis test p-value to the plot
a_pi_Chr + annotate("text", x = 1.5, y = max(pi_g_Chr$log_pi)-0.1, label = paste("Kruskal-Wallis p =", format.pval(kruskal_test_pi_Chr$p.value, digits = 3))) +
  annotate("text", x = 1.6, y = max(pi_g_Chr$log_pi), label = paste("Wilcox test IT.S SW.N p =", format.pval(wilcox_test_pi_Chr$p.value, digits = 3)))




###
# Compare pi density distributions of populations
###

subset_IT_N_pi <- log(At_data_Chr$IT.N_pi[!is.na(At_data_Chr$IT.N_pi) & !is.nan(log(At_data_Chr$IT.N_pi))])
subset_SW_N_pi <- log(At_data_Chr$SW.N_pi[!is.na(At_data_Chr$SW.N_pi) & !is.nan(log(At_data_Chr$SW.N_pi))])
subset_IT_S_pi <- log(At_data_Chr$IT.S_pi[!is.na(At_data_Chr$IT.S_pi) & !is.nan(log(At_data_Chr$IT.S_pi))])
subset_SW_S_pi <- log(At_data_Chr$SW.S_pi[!is.na(At_data_Chr$SW.S_pi) & !is.nan(log(At_data_Chr$SW.S_pi))])

# Perform Kruskal-Wallis test
pi_data <- list(
  SW.N = subset_SW_N_pi,
  IT.N = subset_IT_N_pi,
  SW.S = subset_SW_S_pi,
  IT.S = subset_IT_S_pi
)

kruskal_pi_dist <- kruskal.test(pi_data)

# Print the Kruskal-Wallis test result
print(kruskal_pi_dist)

if (kruskal_pi_dist$p.value < 0.05) {
  dunn_result <- dunn.test(pi_data)
  
  # Print the post-hoc Dunn's test results
  print(dunn_result)
}

# Plot the density
plot(density(subset_IT_S_pi), col="orange", main = "Distribution log Nucleotide Diversity")
lines(density(subset_SW_S_pi), col="lightblue")
lines(density(subset_SW_N_pi), col="blue")
lines(density(subset_IT_N_pi), col = "red")

# Add legend
legend("topright", legend=c("IT.N", "IT.S", "SW.N", "SW.S"),
       col=c("red", "orange", "blue", "lightblue"), lty=1,
       title="Population")
# Add legend for Kruskal-Wallis test p-value
legend("topleft", 
       legend=paste("Kruskal-Wallis p-value:", format(kruskal_pi_dist$p.value, digits=3)), 
       bty="n", 
       cex=1)
legend("topleft", 
       legend=paste("Wilcox test IT.S SW.N p-value:", format(wilcox_test_pi_Chr$p.value, digits=3)), 
       bty="n", 
       cex=1,  
       y.intersp = 3, 
       yjust = 1.0)


###########################################

# TajimaD (td)


# select nucleotide diversity data and calculate means
At_data_Chr %>% dplyr::select(contains("td_Chr")) %>% summarise_all(mean)


# To plot this we need to use "gather" on the data
td_g_Chr <- At_data_Chr %>% dplyr::select(contains("td")) %>% gather(key = "populations", value = "td")

# Remove rows with missing values
td_g2_Chr <- na.omit(td_g_Chr)
td_g2_Chr

# Take the logarithm 
td_g2_Chr$log_td <- log10(td_g2_Chr$td)

a2_td_g2_Chr <- ggplot(td_g2_Chr, aes(populations, log_td, fill = populations)) + 
  geom_boxplot(color = "black") +  # Border color
  scale_fill_manual(values = c("red", "orange", "blue", "lightblue")) +  
  theme_light() + 
  xlab(NULL) +
  ggtitle("Tajma's D of FLC gene") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
        legend.text = element_text(size = 14),  
        legend.title = element_text(size = 16)) +
  ylab("Log10(Tajima's D)") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) 

a2_td_g2_Chr


# Kruskal-Wallis test
kruskal_test_td_Chr <- kruskal.test(log_td ~ populations, data = td_g2_Chr)

# Print the test result
print(kruskal_test_td_Chr)

# generate boxplot again
a2_td_g2_Chr

# adjust text positions by changing x and y below
a2_td_g2_Chr + 
  annotate("text", x = 0.8, y = -1.6, 
           label = paste("Kruskal-Wallis p =", format.pval(kruskal_test_td_Chr$p.value, digits = 3)),
           size = 5)



###
# Compare "Tajima's D" density distributions of populations
###

subset_IT_td <- log(At_data_Chr$IT_td[!is.na(At_data_Chr$IT_td) & !is.nan(log(At_data_Chr$IT_td))])
subset_SW_td <- log(At_data_Chr$SW_td[!is.na(At_data_Chr$SW_td) & !is.nan(log(At_data_Chr$SW_td))])

# Perform Kruskal-Wallis test
td_data <- list(
  SW = subset_SW_td,
  IT = subset_IT_td,
)

kruskal_td_dist <- kruskal.test(td_data)

# Print the Kruskal-Wallis test result
print(kruskal_td_dist)

if (kruskal_td_dist$p.value < 0.05) {
  dunn_result <- dunn.test(td_data)
  
  # Print the post-hoc Dunn's test results
  print(dunn_result)
}



# modify code change xlim to see complete distribution
plot(density(subset_SW_S_td), col="lightblue", main = "Distribution log TajimaD")
lines(density(subset_IT_S_td), col="orange")
lines(density(subset_SW_N_td), col="blue")
lines(density(subset_IT_N_td), col = "red")

# Add legend
legend("topright", legend=c("IT.N", "IT.S", "SW.N", "SW.S"),
       col=c("red", "orange", "blue", "lightblue"), lty=1,
       title="Population")
# Add legend for Kruskal-Wallis test p-value
legend("topleft", 
       legend=paste("Kruskal-Wallis p-value:", format(kruskal_td_dist$p.value, digits=4)), 
       bty="n", 
       cex=1)

