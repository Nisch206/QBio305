library(sp)
library(raster)
library(geodata)
library(terra)

# Load your sample info 
samples = read.table("coordinates_project.txt", header = T)
head(samples)
lon<-samples$lon
lat<-samples$lat

# Extract coordinate data
xy <- samples[, c("lon", "lat")]
str(xy)
# Load BioClim data
biodata = worldclim_global(var = "bio", res = 10, "C:/Users/aless/OneDrive/Desktop/Projekt")
biodata # inspect the data

# Bioclim variable
# BIO19 = Precipitation of Coldest Quarter

#plot the bioclim variable 
plot(biodata[[19]], main = "Precipitation of Coldest Quarter")
# plot multiple variables
plot(biodata[[1:2]])

# Plot the index layer 11
layer_index <- 19
# Extract the layer data
layer_data <- biodata[[layer_index]]
# Plot the layer
plot(layer_data, col = terrain.colors(255), main = paste("Mean Temperature of Coldest Quarter", layer_index))

# Extract Biolclimatic varaibles using xy coordinates dataframe
biodata_extract = extract(biodata[[1:19]], xy, df = T)
summary(biodata_extract) 

#Attach it to the original df
samples_bio = cbind(samples, biodata_extract)

plot(samples_bio$wc2.1_10m_bio_1)
# Extract required Climatic variable, for example we need 
bio1<-samples_bio$wc2.1_10m_bio_1
# Write bio13 data as a text file
write.table(bio1, file = "bio19.txt", row.names = FALSE, col.names = FALSE) 

