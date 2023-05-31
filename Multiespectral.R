# This script takes a list of multispectral `.tif` (UAV or satelite imagery) files in a folder and:
# a) computes vegetation indexes
# b) compute texture metrics 
# c) save all outputs in the output folder
# d) load a list of shape files and extract the mean value of the rasters for each polygon
# e) perform an exploratory analysis of the variables  
# Be mindfull of your RAM capacity when loading n raster files at the same time.

# Daniel Augusto da Silva, 01/06/2023

library(raster)
library(moments)
library(vegan)
library(scales)
library(terra)
library(sf)
library(sp)
library(glcm)
library(rgdal)
library(gdata)
library(dplyr)
library(tidyr)
library(ggplot2)
library(GGally)


setwd("Path")
wdname <- "Path"
output_folder <- "Path"

#### ORTOFOTOS ####
  ## Load images ####

raster_files <- list.files(wdname, pattern = "tif$", full.names = TRUE)


# create an empty list to store the loaded stacks
multi_stacks <- list()

# loop through all stack files
for (i in 1:length(raster_files)) {
  # read the raster
  r <- stack(raster_files[i])
  # extract the UA number from the filename
  num <- substr(paste0(basename(raster_files[i])),16,18)
  # add the VI prefix
  nome <- paste0("Multi",num)
  # Rename the object with VI+UA number
  assign(nome, r, envir = .GlobalEnv)
  # Add the loaded stack to the list with the same name as 'nome'
  multi_stacks[[nome]] <- r
}

  ## Compute vegetation indexes ####

# create an empty list to store the output rasters
vi_rasters <- list()

# loop through all raster files
for (i in 1:length(raster_files)) {
  
  # read the raster
  r <- stack(raster_files[i])
  
  # calculate NDVI
  ndvi <- (r[[5]] - r[[3]]) / (r[[5]] + r[[3]])
  ndvi  <- clamp(ndvi, lower = -1, upper = 1)
  # calculate EVI
  # evi <- 2.5 * ((r[[5]] - r[[3]]) / (r[[5]] + 6 * r[[3]] - 7.5 * r[[1]] + 1))
  # 
  # # calculate SAVI
  # savi <- ((r[[5]] - r[[3]]) / (r[[5]] + r[[3]] + 0.5)) * 1.5
  
  # calculate NDRE
  ndre <- (r[[5]] - r[[4]]) / (r[[5]] + r[[4]])
  ndre  <- clamp(ndre, lower = -1, upper = 1)
  
  # stack NDVI, EVI, and SAVI into a new raster
  vi_raster <- stack(ndvi, ndre)
  names(vi_raster) <- c('ndvi', 'ndre')
  names(vi_raster)
  # set output file name
  output_file <- paste0(output_folder, "/VI/", gsub(".tif", "", basename(raster_files[i])), "_VI.tif")
  
  # save the output raster
  raster::writeRaster(vi_raster, filename=output_file, format = "raster", 
                      options=c("INTERLEAVE=BAND","COMPRESS=LZW"), overwrite = TRUE) 
  removeTmpFiles(h=0)
  # add the output raster to the list
  vi_rasters[[i]] <- vi_raster
  
}

# print the list of output rasters
vi_rasters

rm(ndre, ndvi, r, vi_raster, vi_rasters)
gc()

  ## Compute texture metrics ####

# Load VI rasters
vi_folder <- paste0(output_folder,'/VI')
vi_files <- list.files(vi_folder, pattern = "grd$", full.names = TRUE)

# create an empty list to store the loaded stacks
vi_stacks <- list()

# loop through all stack files
for (i in 1:length(vi_files)) {
  # read the raster
  r <- stack(vi_files[i])
  # extract the UA number from the filename
  num <- substr(paste0(basename(vi_files[i])),16,18)
  # add the VI prefix
  nome <- paste0("VI",num)
  # Rename the object with VI+UA number
  assign(nome, r, envir = .GlobalEnv)
    # Add the loaded stack to the list with the same name as 'nome'
  vi_stacks[[nome]] <- r
  }

rm(r)
gc()

# Initialize a list to store mean values
tex_list <- list()

# loop through all layers in raster stacks
for (i in seq_along(vi_stacks)) {
  r <- stack(vi_stacks[i])
  # get the name of the raster
  stack_name <- names(vi_stacks)[i]
  #creates empty list to be used to name the layers of the final texture stack
  l_list <- list()
  
# loop through all layers in raster stack
    for (j in 1:nlayers(r)) {
      # get the name of the layer
      layer_name <- names(r)[j]
      print(stack_name)
      print(layer_name)
      tex_name <- paste0(stack_name,'_', layer_name)
      
      # Calculate texture metrics (GLCM mean)
      tex <- glcm(r[[j]], window = c(3, 3), na_opt = "any", 
                  statistics = c("mean","variance", "homogeneity", "contrast", "dissimilarity", "entropy", 
                                  "second_moment", "correlation")) #
      
      # loop over all layers of the texture stack
      for (k in 1:nlayers(tex)){
        #gets the stack
        s <- stack(tex)
        #gets the name of the metric
        metric_name <- names(s)[k]
        #concatenate the full name including UA and VI
        full_name <- paste0(tex_name,'_',metric_name)
        #append the full name to the list
        l_list <- append(l_list, full_name)
      }
      
      # # Assign the texture metrics to the tex layers
      names(tex) <- l_list
      l_list <- list()
      # # Assign the texture metrics to the concatenated name
      assign(tex_name, tex, envir = .GlobalEnv)

      # # Store the texture metrics in the list
      tex_list[[tex_name]] <- tex
    }
  
    stack_name <-paste0('TEX_',stack_name)
    s <- raster::stack(tex_list)
    # # Assign the texture metrics to the concatenated name
    assign(stack_name, s, envir = .GlobalEnv)
    tex_list <- list()
    rm(tex)
    rm(s)
    rm(list = ls(pattern = "^VI..._*"), envir = globalenv())
  }

rm(r,vi_stacks)


# save texture rasters to output_folder

#re-creates tex_list

tex_list <- list('tex_raster1', 'tex_raster2',...) # I had to re-create and shorten the tex list because of RAM limitations, but you can use the full list from the loop above just comment out the line 176  

# Loop through all rasters in the list
for (i in seq_along(tex_list)) {
  # Get the original name of the raster
  raster_name <- names(tex_list)[i]
  
  # Get the raster object
  raster <- tex_list[[i]]
  
  # Build the output file path with the original name and desired file extension
  output_file <- paste0(output_folder,"/TEX/", raster_name, ".tif")
  
  # Save the raster to the output file
  writeRaster(raster, filename = output_file, format = "raster", 
              options=c("INTERLEAVE=BAND","COMPRESS=LZW"), overwrite = TRUE)
}
  
  


  
  ## Load sample units shapefiles ####
shape_folder <- paste0(output_folder,'/shapes')
shp_files <- list.files(shape_folder, pattern = "shp$", full.names = TRUE)
su <- list()

# loop through all shape files
for (i in 1:length(shp_files)) {
  # read shp 
  s <- readOGR(shp_files[i])
  # extract the UA number and subunit from the filename
  nome <- substr(paste0(basename(shp_files[i])),7,11)
  # Rename the object with UA+sub denomination
  assign(nome, s, envir = .GlobalEnv)
  # add to list
  su [i] <- s
  rm(s)
}
gc()


  ## Extract raster values for the SU ####
# Defines the pattern on the SU
gc()
pattern <- "286" #i.e. SU number

# Select only itens that have the pattern in the name
su_i <- su[grepl(pattern = pattern, su)]


# Load raster with vegetation index
VI <- stack("path_VI.grd")
names(VI)

# Load raster with texture metrics
texture <- stack("path_TEX_VI.grd")
names(texture)

# Create an empty dataframe to receive the VI values from the sample units
VI_values <- data.frame(ID = character(), ndvi = numeric(), ndre = numeric(), stringsAsFactors = FALSE)


# Sample values in the VI raster stack at shape locations
for (i in 1:length(su_i)){
  #subunit ID
  s <- su_i[[i]]
  id <- s$layer
  print(id)
  values <- extract(VI,s,fun = mean, df = T)
  values$ID <- id
  VI_values <- rbind(VI_values, values)
  
  
}


# Create an empty dataframe to receive the texture values from the sample units
TEX_values <- data.frame(ID = character(), stringsAsFactors = FALSE)

tex_cols <- names(texture) #get the names of the variables in the raster stack

tex_cols <- substr(tex_cols,7,30) #remove the id of the UA from the names

tmp <- data.frame(matrix(ncol = 16)) #create a dataframe with the same number of variables
colnames(tmp) <- tex_cols #assign colnames
tmp <- tmp[complete.cases(tmp),] #remove NA line

TEX_values <- cbind(TEX_values, tmp) # bind with the df that contains the ID col
sub_names <- colnames(TEX_values) #Save it in a variable to use after


# Sample values in the VI raster stack at shape locations
for (i in 1:length(su_i)){
  #subunit ID
  s <- su_i[[i]]
  id <- s$layer
  print(id)
  values <- extract(texture,s,fun = mean, df = T)
  values$ID <- id
  TEX_values <- rbind(TEX_values, values)
  
}

# Change the colnames to the standardized one
colnames(TEX_values) <- sub_names

# Bind the VI and TEX dataframes
vitex <- cbind(VI_values, TEX_values)
# Save it to the output_folder
write.csv2(vitex, file = paste0(output_folder,'/vitex_ua1.csv'))
gc()


#### Exploratory analysis of the multispectral data ####
# Read individual dataframes of the sample units
vitex_ua1 <- read.csv2(paste0(output_folder,'/vitex_ua1.csv'))
names(vitex_ua1)
vitex_ua2 <- read.csv2(paste0(output_folder,'/vitex_ua2.csv'))
names(vitex_ua2)
vitex_ua3 <- read.csv2(paste0(output_folder,'/vitex_ua3.csv'))
names(vitex_ua3)

# Unite the sample unit dataframes in a single one 
UAVdf <- rbind(vitex_ua1,vitex_ua2,vitex_ua3) %>% select(-X,-ID.1)
names(UAVdf)

# Read data frame with field work data
fielddf <- read.csv2(paste0(output_folder,'/field_data.csv'))
names(fielddf)

# Create a single df with field and UAV data
df <- cbind(UAVdf, fielddf) %>% select(-UA)
df$AGB <- as.numeric(df$AGB)
df$AB <- as.numeric(df$AB)
df$DAP <- as.numeric(df$DAP)
df$H <- as.numeric(df$H)
df$LAI <- as.numeric(df$LAI)
df$FCOVER <- as.numeric(df$FCOVER)
summary(df)

# Exploratory analysis of the variables

ggcorr(df, label=TRUE) #visualize correlation between variables
dfnum <- df[c(2:10,12:18,20:26)]


# Creates a df in the long format to see correlation between LAI and predictors
lailong = gather(dfnum, key = "var", value = "value", -LAI) %>%
  group_by(var) %>%
  mutate(density = get_density(value, LAI))

# Plot LAI against all predictors
ggplot(lailong, aes(x = value, y = LAI)) + #, color = density
  geom_point(pch = 20, size = 3) +
  geom_smooth(method='lm',color = 'blue',linewidth=0.1,se=T)+
  labs(x = "") +
  facet_wrap(~ var, scales = "free")

# Creates a df in the long format to see correlation between FCOVER and predictors
fclong = gather(dfnum, key = "var", value = "value", -FCOVER) 
 

# Plot FCOVER against all predictors
ggplot(fclong, aes(x = value, y = FCOVER)) + #, color = density
  geom_point(pch = 20, size = 3) +
  geom_smooth(method='lm',color = 'blue',linewidth=0.1,se=T)+
  labs(x = "") +
  facet_wrap(~ var, scales = "free")

# Creates a df in the long format to see correlation between AGB and predictors
agblong = gather(dfnum, key = "var", value = "value", -AGB) 


# Plot AGB against all predictors
ggplot(agblong, aes(x = value, y = AGB)) + #, color = density
  geom_point(pch = 20, size = 3) +
  geom_smooth(method='lm',color = 'blue',linewidth=0.1,se=T)+
  labs(x = "") +
  facet_wrap(~ var, scales = "free")

# Creates a df in the long format to see correlation between H and predictors
hlong = gather(dfnum, key = "var", value = "value", -H) 


# Plot H against all predictors
ggplot(hlong, aes(x = value, y = H)) + #, color = density
  geom_point(pch = 20, size = 3) +
  geom_smooth(method='lm',color = 'blue',linewidth=0.1,se=T)+
  labs(x = "") +
  facet_wrap(~ var, scales = "free")

# Creates a df in the long format to see correlation between DAP and predictors
daplong = gather(dfnum, key = "var", value = "value", -DAP) 


# Plot DAP against all predictors
ggplot(daplong, aes(x = value, y = DAP)) + #, color = density
  geom_point(pch = 20, size = 3) +
  geom_smooth(method='lm',color = 'blue',linewidth=0.1,se=T)+
  labs(x = "") +
  facet_wrap(~ var, scales = "free")

# Creates a df in the long format to see correlation between AB and predictors
ablong = gather(dfnum, key = "var", value = "value", -AB) 


# Plot DAP against all predictors
ggplot(ablong, aes(x = value, y = AB)) + #, color = density
  geom_point(pch = 20, size = 3) +
  geom_smooth(method='lm',color = 'blue',linewidth=0.1,se=T)+
  labs(x = "") +
  facet_wrap(~ var, scales = "free")



