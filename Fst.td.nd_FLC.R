# Load required libraries

library(PopGenome) 
library(dplyr)
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
library (readr)
library(tidyr)


# set working directory

setwd("C:/Users/aless/OneDrive/Desktop/Project_pop&qGen/Fst")
getwd()

At_Chr <-readVCF("group_3_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz",numcols=89,tid="5", frompos=3173382, topos=3179448, include.unknown =  TRUE)

#To the class of object At_Chr
class(At_Chr)

At_Chr

#Examining the variant data
get.sum.data(At_Chr)

At_Chr@n.biallelic.sites + At_Chr@n.polyallelic.sites
# 22

#To see what slots in Genome class
show.slots(At_Chr)

#To check total number of sites
At_Chr@n.sites
# 6067

#To check starting position and last position of genome class
At_Chr@region.names
# "3173382 - 3179448"

##Deine populations in your dataset

population_info <- read_delim("sample_pop_it_swe.txt", delim = "\t")

# Get the data for the populations
populations <- split(population_info$sample, population_info$pop)
populations
At_Chr <- set.populations(At_Chr, populations, diploid = T)
At_Chr@populations

#Setting up sliding windows

# set chromosome size
chr <- 6067

# set window size and window jump
window_size <- 100
window_jump <- 50

# use seq to find the start points of each window
window_start <- seq(from = 1, to = chr, by = window_jump)
# add the size of the window to each start point 
window_stop <- window_start + window_size

# no windows start before the end of chromosome 4
sum(window_start > chr)
# 0

# but some window stop positions do occur past the final point
sum(window_stop > chr)
# 2

# remove windows from the start and stop vectors
window_start <- window_start[which(window_stop < chr)]
window_stop <- window_stop[which(window_stop < chr)]

chr - window_stop[length(window_stop)]
# 16

# save as a data.frame
windows <- data.frame(start = window_start, stop = window_stop, 
                      mid = window_start + (window_stop-window_start)/2)

# make a sliding window dataset
At_sw <- sliding.window.transform(At_Chr, width = 100, jump = 50, type = 2)

# calculate diversity statistics - nd
At_sw <- diversity.stats(At_sw, pi = TRUE)

# calculate diversity statistics - FST
At_sw <- F_ST.stats(At_sw, mode = "nucleotide") 
#"nucleotide" to specify we want it to be calculated sliding averages
#of nucleotides, rather than using haplotype data


#Extracting statistics for visualization

# get the nucleotide diversity data.
# extract nucleotide diversity and correct for window size
nd <- At_sw@nuc.diversity.within/100

#estimates need to be corrected for window size - so we divide them by 100 bp.

# Add the population names to each of estimate
# make population name vector
pops <- c("IT","SWE") 

# set population names
colnames(nd) <- paste0(pops, "_pi")

# extract fst values
fst <- t(At_sw@nuc.F_ST.pairwise)
# t() to transpose the F_ST matrix so that each column is a pairwise
#comparison and each row is an estimate for a genome window.

# extract dxy - pairwise absolute nucleotide diversity
dxy <- get.diversity(At_sw, between = T)[[2]]/100
#As with nucleotide diversity, we also corrected d_XY_ for the window size.

# get column names 
x <- colnames(fst)
# does the same thing as above but by indexing the pops vector
x <- sub("pop1", pops[1], x)
x <- sub("pop2", pops[2], x)

# replace forward slash
x <- sub("/", "_", x)
# look at x to confirm the replacement has occurred
x

# Change the column names of our two data matrices, before we put 
#everything together in our final dataset.
colnames(fst) <- paste0(x, "_fst")
colnames(dxy) <- paste0(x, "_dxy")

# Combine nd, FST and d_XY_ datasets with our windows information from earlier into a big dataset.

At_data <- as_tibble(data.frame(windows, nd, fst, dxy))
dim(At_data)
#120      7


#####Visualizing the data - distributions#####

#For the purposes of this session, we will focus mainly on the difference between Italian and Swedish
#Arabidopsis pop.
#For example, let's say we want to look at mean nucleotide diversity, we can do that like so:

# select nucleotide diversity data and calculate means
At_data %>% dplyr::select(contains("pi")) %>% summarise_all(mean)
#we used select and contains to select columns from our main dataset that contain 
#pi - i.e. nucleotide diversity columns. We then used summarise_all and mean to calculate
#the mean value for all of the columns we selected.

# To plot this we need to use "gather" on the data

pi_g <- At_data %>% dplyr::select(contains("pi")) %>% gather(key = "populations", value = "pi")

# make a boxplot
a <- ggplot(pi_g, aes(populations, pi)) + geom_boxplot() + theme_light() + xlab(NULL)
a

# Taking the logarithm (in this case, log base 10) of the values can be useful when dealing with data
# that spans several orders of magnitude. This transformation helps in visually emphasizing relative
# differences in the data, especially when there are large variations in scale.

pi_g$log_pi <- log10(pi_g$pi)

a <- ggplot(pi_g, aes(populations, log_pi, fill = populations)) + 
  geom_boxplot(color = "black") +  # Border color
  scale_fill_manual(values = c("red", "blue")) +  # Box fill colors
  theme_light() + 
  xlab(NULL) +
  ylab("Log10(pi)")

a

#This makes it much clearer how nucleotide diversity differs among the populations.

#When comparing two boxplots to determine if they are statistically different, one can perform
# statistical tests such as the t-test or Wilcoxon rank-sum test. 
# Wilcoxon rank-sum test
wilcox_test_result <- wilcox.test(log_pi ~ populations, data = pi_g)

# Print the test result
print(wilcox_test_result)

# Add p-value to the plot
a + annotate("text", x = 1.5, y = max(pi_g$log_pi), label = paste("p =", format.pval(wilcox_test_result$p.value, digits = 3)))

#####Visualizing patterns along the chromosome ####
#Let's have a look at how FST between Italian and Swedish populations varies along chromosomes.
#We can do this very simply with ggplot.

a <- ggplot(At_data, aes(mid/10^6, IT_SWE_fst)) + geom_line(colour = "red")
a <- a + xlab("Position (Mb)") + ylab(expression(italic(F)[ST]))
a + theme_light()


#to plot nd, FST and d_XY_ to examine how they co-vary along the genome. 
#This requires a bit of data manipulation, but is relatively straightforward. We will break it down into steps.
# select data of interest
hs <- At_data %>% dplyr::select(mid, IT_pi, SWE_pi, IT_SWE_fst, IT_SWE_dxy)

# To set Fst values smaller than zero to zero in the specified columns of a data frame using dplyr and
# the pipe operator %>%, you can use the mutate function along with across

hs <- At_data %>%
  dplyr::select(mid, IT_pi, SWE_pi, IT_SWE_fst, IT_SWE_dxy) %>%
  mutate(across(c(IT_SWE_fst, IT_SWE_dxy), ~ ifelse(. < 0, 0, .)))

# use gather to rearrange everything
hs_g <- gather(hs, -mid, key = "stat", value = "value")

#Plot everything together like so:
a <- ggplot(hs_g, aes(mid/10^6, value, colour = stat)) + geom_line()
a <- a + xlab("Position (Mb)")
a + theme_light()

# To take the logarithm of the value variable in your ggplot code, you can use the log10() function
# within the aes() mapping.
hs_g$log_value <- log10(hs_g$value)

a <- ggplot(hs_g, aes(mid/10^6, log_value, colour = stat)) + geom_line()
a <- a + xlab("Position (Mb)")
a + theme_light()

# construct a plot with facets
a <- ggplot(hs_g, aes(mid/10^6, value, colour = stat)) + geom_line()
a <- a + facet_grid(stat~., scales = "free_y")
a <- a + xlab("Position (Mb)")
a + theme_light() + theme(legend.position = "none")

#The facet_grid function allows us to split our data across panels for quick and easy visualization.
#In this case, we split our data by the stat variable - we used stat~. to specify we want this done
#by rows (compare with .~stat for the column equivalent). We also specified that we wanted the scales
#on our y-axes to vary with scales = free_y.

#However, before we examine our plot in detail, it would also be easier if we rearranged everything 
#so FST came at the top, ?? beneath it and then finally, d_XY_. How can we do that? Well we need to
#reorder the stat factor in our hs_g dataset.
# first make a factor
x <- factor(hs_g$stat)
# then reorder the levels
x <- factor(x, levels(x)[c(3, 1, 4, 2)])
# add to data.frame
hs_g$stat <- x
#This looks a little complicated, but in the square brackets above we simply rearranged what order
#our facets are displayed. We can replot our figure to demonstrate this:

# construct a plot with facets
a <- ggplot(hs_g, aes(mid/10^6, value, colour = stat)) + geom_line()
a <- a + facet_grid(stat~., scales = "free_y")
a <- a + xlab("Position (Mb)")
a + theme_light() + theme(legend.position = "none")
