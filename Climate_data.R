library(sp)
library(raster)
library(geodata)
library(terra)

# Don't Foget to change the path to a folder where you stored your climate data and R script
setwd("C:/Users/aless/OneDrive/Desktop/Project_pop&qGen/climate_data_project")

# Load your sample info "coordinates.txt" a tab separated text file with four columns (sample pop lat lon)
# "coordinates.txt" this file is created for sample project, create a tab separated text file consisting of
# four columns (sample pop lat lon) for set of samples you are using for using project
# order of samples should be the same as in your project VCF file.
samples <- read.table("coordinates_project.txt", header = T)
head(samples)
lon<-samples$lon
lat<-samples$lat

# Extract coordinate data
xy <- samples[, c("lon", "lat")]
str(xy)
# Load BioClim data. The following command checks if the data is present at the 
# specified path. If the data is not present, it will be downloaded. Don't Forget to change the path to a folder where you stored your climate data
biodata <- worldclim_global(var = "bio", res = 10, "C:/Users/aless/OneDrive/Desktop/Project_pop&qGen/climate_data_project")
biodata # inspect the data

# Names of all bioclim variables
## BIO1 = Annual Mean Temperature, BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
## BIO3 = Isothermality (BIO2/BIO7) (??100), BIO4 = Temperature Seasonality (standard deviation ??100)
## BIO5 = Max Temperature of Warmest Month, BIO6 = Min Temperature of Coldest Month
## BIO7 = Temperature Annual Range (BIO5-BIO6), BIO8 = Mean Temperature of Wettest Quarter
## BIO9 = Mean Temperature of Driest Quarter, BIO10 = Mean Temperature of Warmest Quarter
## BIO11 = Mean Temperature of Coldest Quarter, BIO12 = Annual Precipitation ####################################################
## BIO13 = Precipitation of Wettest Month, BIO14 = Precipitation of Driest Month
## BIO15 = Precipitation Seasonality (Coefficient of Variation), BIO16 = Precipitation of Wettest Quarter
## BIO17 = Precipitation of Driest Quarter, BIO18 = Precipitation of Warmest Quarter
## BIO19 = Precipitation of Coldest Quarter 

#plot a particular bioclim variable from the list
plot(biodata[[11]], main = "Mean Temperature of Coldest Quarter")
# plot multiple variables
plot(biodata[[1:2]])

# Specify the layer you want to plot (replace 1 with the index of the layer you want to plot)
layer_index <- 11
# Extract the layer data
layer_data <- biodata[[layer_index]]
# Plot the layer
plot(layer_data, col = terrain.colors(255), main = paste("Mean Temperature of Coldest Quarter", layer_index))

# Extract Biolclimatic varaibles using xy coordinates dataframe
biodata_extract <- extract(biodata[[1:19]], xy, df = T)
summary(biodata_extract) 

#Attach it to the original df
samples_bio <- cbind(samples, biodata_extract)

plot(samples_bio$wc2.1_10m_bio_1)
# Extract required Climatic variable, for example we need 
bio1<-samples_bio$wc2.1_10m_bio_1
# Write bio13 data as a text file
write.table(bio1, file = "bio11_v2.txt", row.names = FALSE, col.names = FALSE) 

