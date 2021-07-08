####Conservatism of invasive leaf beetles####
# Authors of the code: Marianna Simoes, Claudia Nunez-Penichet and Marcus Krull
#Date: 2021
####################################################################################
###############

# packages


# Working directory
setwd('E:\\Expansion')
species <- as.character(read.csv("species_taxmtch.csv")[, 1])

# Functions

source("00Functions.R")

# Downloading
## Directory to save directly
get_gbif_data(species = species, maximum_n = 199800, 
              return = FALSE, save = TRUE)

# Initial cleaning of occurrence data
## Excluding: records with coordinates (0, 0), and records with low precision
D <- list.files(path = "GBIF_data2", pattern = "_georef.csv$", full.names = T)
nam <- list.files(path = "GBIF_data2", pattern = "_georef.csv$", full.names = F)
nam <- gsub("_georef", "", nam)
Finalnam <- paste0("Final2/", nam)
all_species <- list()
i=1
setwd('E:\\Expansion')
i=15
for (i in 1:length(D)){
  occurrences <- read.csv(D[i]) # occurrences 
  colnames(occurrences) <- c("Species", "Longitude", "Latitude")
  
  # Excluding records with (0, 0) coordinates
  occurrences <- occurrences[occurrences$Longitude != 0 & occurrences$Latitude != 0, ]
  
  # Excluding records with low level of precision (<= 2 decimals)
  ## small function to detect precision 
  ## (from https://stackoverflow.com/questions/5173692/how-to-return-number-of-decimal-places-in-r)
  decimalplaces <- function(x) {
    if (abs(x - round(x)) > .Machine$double.eps^0.5) {
      nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
    } else {
      return(0)
    }
  }
  
  occurrences <- occurrences[sapply(occurrences$Longitude, decimalplaces) >= 2 & # keep only the ones with more than 1 decimals
                               sapply(occurrences$Latitude, decimalplaces) >= 2, ]
  
  # saving the new set of occurrences inside continents and area of interest
  write.csv(occurrences, Finalnam[i], row.names = FALSE)
  
  all_species[[i]] <- occurrences
}


#Creating a big table with all species
all_species <- do.call(rbind, all_species)
write.csv(all_species, "all_species2.csv", row.names = FALSE)


#Selecting species
#Number of records
require(spThin)
data <- read.csv('all_species.csv')
head(data)
spvector<-unique(data$Species)
df<-NULL
for (i in 1:length(spvector)) {
  sub <-data[data$Species == spvector[i],]
  row <- nrow(sub) 
  df = rbind(df, data.frame(spvector[i],row))
}  
write.csv(df, 'C:\\Users\\hanie\\Desktop\\summary_data.csv', row.names = F)
colnames(deep) <- c('Latitude','num.records' ,'num.species','chlorophyl','current.velocity','dissolved.oxygen', 'nitrate','salinity','temperature','topo', 'ES15')
head(deep)

### Distribution:Occurence thining and Na omit------- 

#quick fix to delete thin files that were incorrect--------

junk <- list()

for (i in 1:length(spvector)) {
  junk[[i]] <- list.files(path = mydir, pattern="thin.*csv$",, full.names = TRUE, recursive = TRUE)
}

file.remove(junk[[i]])

#Thinning Dataset------
setwd('E:\\Expansion\\Dataset')
spvec <- dir()
spvector <- spvec
final_occ <- list()
final_occ_i <- list()

r <- raster('E:\\Expansion\\WC_2.1\\wc2.1_2.5m_bio_1.tif')

for (i in 1:length(spvector)) {
  # reading species data
  sp <- read.csv(paste(spvector[i], "/", paste(spvector[i], "_n.csv", sep = ""), sep = ""))
  sp_i <- read.csv(paste(spvector[i], "/", paste(spvector[i], "_i.csv", sep = ""), sep = ""))
  
  sp_nodup <- unique(sp)
  sp_nodup_i <- unique(sp_i)
  
  if (dim(sp_nodup)[1] > 1000) {
    sp_nodup <- sp_nodup[sample(1:nrow(sp_nodup), 1000), ]
  }
  if (dim(sp_nodup_i)[1] > 1000) {
    sp_nodup_i <- sp_nodup_i[sample(1:nrow(sp_nodup_i), 1000), ]
  }
  
  # thinning
  thinnedcat_a <- thin( loc.data = sp_nodup, lat.col = "Latitude", 
                        long.col = "Longitude", 
                        spec.col = "Species", thin.par = 10, reps = 5, 
                        locs.thinned.list.return = TRUE, write.files = FALSE, 
                        max.files = 1, out.dir = spvector[i], 
                        out.base = "cat_a_thined", write.log.file = FALSE )
  
  thin <- cbind(as.character(sp_nodup[1, 1]), thinnedcat_a[[5]])
  colnames(thin) <- colnames(sp_nodup)
  correc <- extract(r, thin[2:3]) 
  extract <- cbind(thin[1:3], correc)
  finalOcc <- na.omit(extract)
  occN<- finalOcc[1:3]
  
  write.csv(occN, paste0(spvector[i], "/", paste0(spvector[i], "_n_thin.csv")),
            row.names = FALSE)
  
  thinnedcat_i <- thin( loc.data = sp_nodup_i, lat.col = "Latitude", 
                        long.col = "Longitude", 
                        spec.col = "Species", thin.par = 10, reps = 5, 
                        locs.thinned.list.return = TRUE, write.files = FALSE, 
                        max.files = 1, out.dir = spvector[i], 
                        out.base = "cat_a_thined", write.log.file = FALSE )
  
  thin_i <- cbind(as.character(sp_nodup_i[1, 1]), thinnedcat_i[[5]])
  colnames(thin_i) <- colnames(sp_nodup_i)
  correc_i <- extract(r, thin_i[2:3]) 
  extract_i <- cbind(thin_i[1:3], correc_i)
  finalOcc_i <- na.omit(extract_i)
  occI<- finalOcc_i[1:3]
  write.csv(occI, paste0(spvector[i], "/", paste0(spvector[i], "_i_thin.csv")),
            row.names = FALSE)
  
    # number of final records per species
  final_occ[[i]] <- data.frame(species = spvector[i], n_records = dim(occN)[1])
  final_occ_i[[i]] <- data.frame(species = spvector[i], n_records = dim(occI)[1])
  
  cat("\nSpecies", i, "of", length(spvector), "finished\n")
}

# writing summary of species number of records
final_n <- do.call(rbind, final_occ)
final_i <- do.call(rbind, final_occ_i)

write.csv(final_n, "summary_occ_n.csv", row.names = FALSE)
write.csv(final_i, "summary_occ_i.csv", row.names = FALSE)

#Creating M -----------
install.packages("rangemap")
require(rangemap)
require(maptools)
require(sf)
require(rgdal)
require(raster)
require(rgeos)

setwd('E:\\Expansion\\Dataset2')

spvector <- dir()


var_list <- list(stack(list.files("E:/Expansion/WC_2.1/Vars/", pattern = ".tif",
                                  full.names = TRUE)))

var_names <- list(list.files("E:/Expansion/WC_2.1/Vars/", pattern = ".tif",
                             full.names = FALSE))

var<-c("bio_2","bio_4","bio_5","bio_6","bio_7","bio_13","bio_14", "bio_15")


TEOW <- readOGR(dsn="E:/Expansion/Background/official_teow/official", layer="wwf_terr_ecos")
proj4string(TEOW)  <- sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
#KOPPEN <- readOGR(dsn='E:/Expansion/Background/world_koppen', layer= 'world_koppen')

KOPPEN <- readOGR(dsn='E:/Expansion/Background/WC05_1975H_Koppen_Shapefile/WC05_1975H_Koppen_Shapefile', layer= 'WC05_1975H_Koppen')
proj4string(KOPPEN) <- sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")


for (i in 19:length(spvector)) {
  # reading species data
  #native
  sp <- read.csv(paste(spvector[i], "/", paste(spvector[i], "_n_thin.csv", sep = ""), 
                       sep = ""))
  sp_nodup <- unique(sp)
  
  #invasive
  sp_i <- read.csv(paste(spvector[i], "/", paste(spvector[i], "_i_thin.csv", sep = ""), 
                         sep = ""))
  sp_nodup_i <- unique(sp_i)
  
  #Calibration areas
  WGS84 <- sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  occ_pr <- sp::SpatialPointsDataFrame(coords = sp_nodup[, 2:3], data = sp_nodup,
                                       proj4string = WGS84)
  occ_pri <- sp::SpatialPointsDataFrame(coords = sp_nodup_i[, 2:3], data = sp_nodup_i,
                                        proj4string = WGS84)
  #Buffer
  BUFFAREA_n <-  rgeos::gBuffer(occ_pr, width = 1) 
  projection(BUFFAREA_n) <- CRS('+proj=longlat')
  BUFFAREA_i <- rgeos::gBuffer(occ_pri, width = 1) 
  projection(BUFFAREA_i) <- CRS('+proj=longlat')
  
  #Biomes
  Bio_n <- over(occ_pr, TEOW)
  Bio_n <- unique(Bio_n$BIOME)
  BIOMEPly <- TEOW[TEOW$BIOME %in% Bio_n, ]
  #break down  polygons
  BIOMEPly <- gUnaryUnion(BIOMEPly)
  projection(BIOMEPly) <- CRS('+proj=longlat')
  
  #plot(BIOMEPly)
  Bio_i <- over(occ_pri, TEOW)
  Bio_i <- unique(Bio_i$BIOME)
  BIOMEPlyi <- TEOW[TEOW$BIOME %in% Bio_i, ]
  BIOMEPlyi <- gUnaryUnion(BIOMEPlyi)
  projection(BIOMEPlyi) <- CRS('+proj=longlat')
  #plot(BIOMEPlyi)
  
  #Koppen
  Kop_n <- na.omit(over(occ_pr, KOPPEN))
  Kop_n <- unique(Kop_n$GRIDCODE)
  KOPPENPly <- KOPPEN[KOPPEN$GRIDCODE %in% Kop_n, ]
  KOPPENPly <- gUnaryUnion(KOPPENPly)
  projection(KOPPENPly) <- CRS('+proj=longlat')
  #plot(KOPPENPly)
  #str(Kop_n)
  
  Kop_i <- na.omit(over(occ_pri, KOPPEN))
  Kop_i <- unique(Kop_i$GRIDCODE)
  KOPPENPlyi <- KOPPEN[KOPPEN$GRIDCODE %in% Kop_i, ]
  KOPPENPlyi <- gUnaryUnion(KOPPENPlyi)
  projection(KOPPENPlyi) <- CRS('+proj=longlat')
  
  # folder for Ms
  set <- var_list[[1]]
  infolder <- paste(spvector[i], "M_variables", sep = "/")
  dir.create(infolder)
  
  # mask variables
  masked1 <- mask(crop(set, BUFFAREA_n), BUFFAREA_n)
  mask1 <- paste(infolder,'BUFFAREA_n' , sep = "/")
  dir.create(mask1)
  rnames <- paste0(mask1, '/', names(masked1), ".tif")
  for (k in 1:length(unstack(masked1))){
    writeRaster(masked1[[k]], filename = rnames[k], format = "GTiff",overwrite=TRUE)}
  
  masked2 <- mask(crop(set, BUFFAREA_i), BUFFAREA_i)
  mask2 <- paste(infolder,'BUFFAREA_i' , sep = "/")
  dir.create(mask2)
  rnames2 <- paste0(mask2, '/', names(masked2), ".tif")
  for (k in 1:length(unstack(masked2))) {
    writeRaster(masked2[[k]], filename = rnames2[k], format = "GTiff",overwrite=TRUE)}
  
  
  masked3 <- mask(crop(set, BIOMEPly), BIOMEPly)
  mask3 <- paste(infolder,'BIOMEPly' , sep = "/")
  dir.create(mask3)
  rnames3 <- paste0(mask3, '/', names(masked3), ".tif")
  for (k in 1:length(unstack(masked3))) {
    writeRaster(masked3[[k]], filename = rnames3[k], format = "GTiff",overwrite=TRUE)}
  
  masked4 <- mask(crop(set, BIOMEPlyi), BIOMEPlyi)
  mask4 <- paste(infolder,'BIOMEPlyi' , sep = "/")
  dir.create(mask4)
  rnames4 <- paste0(mask4, '/', names(masked4), ".tif")
  for (k in 1:length(unstack(masked4))) {
    writeRaster(masked4[[k]], filename = rnames4[k], format = "GTiff",overwrite=TRUE)}
  
  masked5 <- mask(crop(set, KOPPENPly), KOPPENPly)
  mask5 <- paste(infolder,'KOPPENPly' , sep = "/")
  dir.create(mask5)
  rnames5 <- paste0(mask5, '/', names(masked5), ".tif")
  for (k in 1:length(unstack(masked5))) {
    writeRaster(masked5[[k]], filename = rnames5[k], format = "GTiff",overwrite=TRUE)}
  
  masked6 <- mask(crop(set, KOPPENPlyi), KOPPENPlyi)
  mask6 <- paste(infolder,'KOPPENPlyi' , sep = "/")
  dir.create(mask6)
  rnames6 <- paste0(mask6, '/', names(masked6), ".tif")
  for (k in 1:length(unstack(masked6))) {
    writeRaster(masked6[[k]], filename = rnames6[k], format = "GTiff",overwrite=TRUE)}
  
  cat(i, "species of", length(spvector), "\n")
}


#Ecospat------------
#I already did the ENM, now Ill check for similarity between Inv and Nativ?
#Then hypervolume: Is the climate niche different between Inv /Nativ?
#Then niche breath to argue for the expansion 

install.packages("devtools")
require(devtools)
require(ecospat)
require(raster)
require(rworldmap)
require(ellipsenm)
require(rgdal)
require(spocc)
require(ade4)

setwd('E:\\Expansion\\Dataset2')

spvector <- dir()

var<-c("bio_2","bio_4","bio_5","bio_6","bio_7","bio_13","bio_14", "bio_15")

clim <- stack(list.files(path = "E:\\Expansion\\WC_2.1\\Vars\\", pattern = ".tif", full.names = TRUE))
crs(clim) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

Results <- NULL

for (i in 1:length(spvector)) {
  # reading species data
  #native
  sp <- read.csv(paste(spvector[i], "/", paste(spvector[i], "_n_thin.csv", sep = ""), 
                       sep = ""))
  sp_nodup <- unique(sp)
  sp_nodup <- sp_nodup[2:3]
  
  #invasive
  sp_i <- read.csv(paste(spvector[i], "/", paste(spvector[i], "_i_thin.csv", sep = ""), 
                         sep = ""))
  sp_nodup_i <- unique(sp_i)
  sp_nodup_i <- sp_nodup_i[2:3]
  
  #Calibration areas
  WGS84 <- sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  occ_pr <- sp::SpatialPointsDataFrame(coords = sp_nodup[, 1:2], data = sp_nodup,
                                       proj4string = WGS84)
  occ_pri <- sp::SpatialPointsDataFrame(coords = sp_nodup_i[, 1:2], data = sp_nodup_i,
                                        proj4string = WGS84)
  
  Buffer_n <- stack(list.files(path = paste(spvector[i], "M_variables", 'BUFFAREA_n/', sep = "/"),
                               pattern = ".tif", full.names = TRUE))
  Buffer_i <- stack(list.files(path = paste(spvector[i], "M_variables", 'BUFFAREA_i/', sep = "/"),
                               pattern = ".tif", full.names = TRUE))
  
  Biome_n <- stack(list.files(path = paste(spvector[i], "M_variables", 'BIOMEPly/', sep = "/"),
                              pattern = ".tif", full.names = TRUE))
  Biome_i <- stack(list.files(path = paste(spvector[i], "M_variables", 'BIOMEPlyi/', sep = "/"),
                              pattern = ".tif", full.names = TRUE))
  
  Koppen_n <- stack(list.files(path = paste(spvector[i], "M_variables", 'KOPPENPly/', sep = "/"),
                               pattern = ".tif", full.names = TRUE))
  Koppen_i <- stack(list.files(path = paste(spvector[i], "M_variables", 'KOPPENPlyi/', sep = "/"),
                               pattern = ".tif", full.names = TRUE))
  
  #Call background calibration------
  #M=Buffer
  env.bkg.Buf_n <- na.exclude(data.frame(rasterToPoints(Buffer_n))) 
  env.bkg.Buf_n <- env.bkg.Buf_n[,3:10]
  env.bkg.Buf_n$Species <- 'Native'
  
  env.bkg.Buf_i <- na.exclude(data.frame(rasterToPoints(Buffer_i))) 
  env.bkg.Buf_i <- env.bkg.Buf_i[,3:10]
  env.bkg.Buf_i$Species <- 'Invasive'
  
  #M=Biome  
  env.bkg.Bio_n <- na.exclude(data.frame(rasterToPoints(Biome_n))) 
  env.bkg.Bio_n <- env.bkg.Bio_n[,3:10]
  env.bkg.Bio_n$Species <- 'Native'
  
  env.bkg.Bio_i <- na.exclude(data.frame(rasterToPoints(Biome_i))) 
  env.bkg.Bio_i <- env.bkg.Bio_i[,3:10]
  env.bkg.Bio_i$Species <- 'Invasive'
  
  #M=Koppen
  env.bkg.Kop_n <- na.exclude(data.frame(rasterToPoints(Koppen_n))) 
  env.bkg.Kop_n <- env.bkg.Kop_n[,3:10]
  env.bkg.Kop_n$Species <- 'Native'
  
  env.bkg.Kop_i <- na.exclude(data.frame(rasterToPoints(Koppen_i))) 
  env.bkg.Kop_i <- env.bkg.Kop_i[,3:10]
  env.bkg.Kop_i$Species <- 'Invasive'
  
  #Occ environemnt: Background----
  env.occ.Bkg_n<- extract(clim, sp_nodup)#Native 
  env.occ.Bkg_n <- na.exclude(cbind(env.occ.Bkg_n,sp_nodup))
  
  env.occ.Bkg_i <-extract(clim,sp_nodup_i)#Invasive
  env.occ.Bkg_i <- na.exclude(cbind(env.occ.Bkg_i,sp_nodup_i))
  
  #calibration of PCA-env 
  #M=Buffer
  pca.env_Buf <-dudi.pca(rbind(env.bkg.Buf_n,env.bkg.Buf_i)[,1:8], center = T, scale = T, scannf = F, nf = 2)
  
  #Buffer Native
  scores.bkg_Bf_n<- pca.env_Buf$li
  scores.bkg.Buf_n<- na.omit(suprow(pca.env_Buf,env.bkg.Buf_n[,1:8])$li)
  #Buffer Invasive
  scores.bkg_Bf_i<- pca.env_Buf$li
  scores.bkg.Buf_i<- na.omit(suprow(pca.env_Buf,env.bkg.Buf_i[,1:8])$li)
  
  #M=Biome  
  pca.env_Bio <-dudi.pca(rbind(env.bkg.Bio_n,env.bkg.Bio_i)[,1:8], center = T, scale = T, scannf = F, nf = 2)
  
  #Biome Native  
  scores.bkg_B_n<- pca.env_Bio$li
  scores.bkg.Bio_n<- na.omit(suprow(pca.env_Bio,env.bkg.Bio_n[,1:8])$li)
  #Biome Invasive
  scores.bkg_B_i<- pca.env_Bio$li
  scores.bkg.Bio_i<- na.omit(suprow(pca.env_Bio,env.bkg.Bio_i[,1:8])$li)  
  
  #M=Koppen 
  pca.env_Kop <-dudi.pca(rbind(env.bkg.Kop_n,env.bkg.Kop_i)[,1:8], center = T, scale = T, scannf = F, nf = 2)
  #Biome Native  
  scores.bkg_K_n<- pca.env_Kop$li
  scores.bkg.Kop_n<- na.omit(suprow(pca.env_Kop,env.bkg.Kop_n[,1:8])$li)
  #Biome Invasive
  scores.bkg_K_i<- pca.env_Kop$li
  scores.bkg.Kop_i<- na.omit(suprow(pca.env_Kop,env.bkg.Kop_i[,1:8])$li)
  
  #Scores for occurences 
  scores.occ.Buff_n<- na.omit(suprow(pca.env_Buf,env.occ.Bkg_n[,1:8])$li)#head(env.occ.nativ)
  scores.occ.Bio_n<- na.omit(suprow(pca.env_Bio,env.occ.Bkg_n[,1:8])$li)
  scores.occ.Kop_n<- na.omit(suprow(pca.env_Kop,env.occ.Bkg_n[,1:8])$li)
  
  scores.occ.Buff_i<- na.omit(suprow(pca.env_Buf,env.occ.Bkg_i[,1:8])$li)
  scores.occ.Bio_i<- na.omit(suprow(pca.env_Bio,env.occ.Bkg_i[,1:8])$li)
  scores.occ.Kop_i<- na.omit(suprow(pca.env_Kop,env.occ.Bkg_i[,1:8])$li)
  
  # calculation of occurence density
  z1_Buf_n<- ecospat.grid.clim.dyn(scores.bkg_Bf_n,scores.bkg.Buf_n,scores.occ.Buff_n, R=100)
  z1_Bio_n<- ecospat.grid.clim.dyn(scores.bkg_B_n,scores.bkg.Bio_n,scores.occ.Bio_n, R=100 )
  z1_Kop_n<- ecospat.grid.clim.dyn(scores.bkg_K_n,scores.bkg.Kop_n,scores.occ.Kop_n,R=100)
  
  z2_Buf_i <- ecospat.grid.clim.dyn(scores.bkg_Bf_i,scores.bkg.Buf_i,scores.occ.Buff_i,R=100)
  z2_Bio_i<- ecospat.grid.clim.dyn(scores.bkg_B_i,scores.bkg.Bio_i,scores.occ.Bio_i,R=100)
  z2_Kop_i<- ecospat.grid.clim.dyn(scores.bkg_K_i,scores.bkg.Kop_i,scores.occ.Kop_i,R=100)
  
  #similarity ---------
  #niche conservatism = greater - p<0.05
  sim1_Buf<-ecospat.niche.similarity.test(z1_Buf_n,z2_Buf_i,rep=1000,alternative = "greater", rand.type = 1) #niches randomly shifted in both area #If rand.type = 1, both z1 and z2 are randomly shifted, if rand.type =2, only z2 is randomly shifted.
  sim1_Bio<-ecospat.niche.similarity.test(z1_Bio_n,z2_Bio_i,rep=1000,alternative = "greater", rand.type = 1) #niches randomly shifted in both area #If rand.type = 1, both z1 and z2 are randomly shifted, if rand.type =2, only z2 is randomly shifted.
  sim1_Kop<-ecospat.niche.similarity.test(z1_Kop_n,z2_Kop_i,rep=1000,alternative = "greater", rand.type = 1) #niches randomly shifted in both area #If rand.type = 1, both z1 and z2 are randomly shifted, if rand.type =2, only z2 is randomly shifted.
  
  
  #Equivalency ----------
  # eq1_Buf<-ecospat.niche.equivalency.test(z1_Buf_n,z2_Buf_i,rep=1000) #niches randomly shifted in both area #If rand.type = 1, both z1 and z2 are randomly shifted, if rand.type =2, only z2 is randomly shifted.
  #  eq1_Bio<-ecospat.niche.equivalency.test(z1_Bio_n,z2_Bio_i,rep=1000) #niches randomly shifted in both area #If rand.type = 1, both z1 and z2 are randomly shifted, if rand.type =2, only z2 is randomly shifted.
  #  eq1_Kop<-ecospat.niche.equivalency.test(z1_Kop_n,z2_Kop_i,rep=1000) #niches randomly shifted in both area #If rand.type = 1, both z1 and z2 are randomly shifted, if rand.type =2, only z2 is randomly shifted.
  
  #Niche overlap corrected by availabilty of background conditions-------------
  NO_Buf <- ecospat.niche.overlap(z1_Buf_n,z2_Buf_i,cor=T) 
  NO_Bio <- ecospat.niche.overlap(z1_Bio_n,z2_Bio_i,cor=T) 
  NO_Kop <- ecospat.niche.overlap(z1_Kop_n,z2_Kop_i,cor=T) 
  
  #Niche dynamics -------------------
  Dy_Buf<-ecospat.niche.dyn.index(z1_Buf_n, z2_Buf_i,intersection=NA)#intersection=NA, the analysis is performed on the whole environmental extent (native and invaded)
  Dy_Bio<-ecospat.niche.dyn.index(z1_Bio_n, z2_Bio_i,intersection=NA)
  Dy_Kop<-ecospat.niche.dyn.index(z1_Kop_n, z2_Kop_i,intersection=NA)
  
  Results <- rbind(Results, data.frame(spvector[i],sim1_Buf$obs$D, sim1_Buf$obs$I, sim1_Buf$p.D, sim1_Buf$p.I,
                                       sim1_Bio$obs$D, sim1_Bio$obs$I, sim1_Bio$p.D, sim1_Bio$p.I, 
                                       sim1_Kop$obs$D, sim1_Kop$obs$I, sim1_Kop$p.D, sim1_Kop$p.I,
                                       NO_Buf$D, NO_Buf$I,NO_Bio$D, NO_Bio$I, NO_Kop$D, 
                                       NO_Kop$I, Dy_Buf$dynamic.index.w, Dy_Bio$dynamic.index.w, Dy_Kop$dynamic.index.w ))  
  cat(i, "of", length(spvector), "\n")
}

colnames(Results) <- c('Species','Similarity: Buffer (D)', 'Similarity: Buffer (I)', 'Similarity: Buffer (p-value: D)', 'Similarity: Buffer (p-value: I)', 
                       'Similarity: Biome (D)', 'Similarity: Biome (I)', 'Similarity: Biome (p-value: D)', 'Similarity: Biome (p-value: I)', 
                       'Similarity: Kopper (D)', 'Similarity: Kopper (I)', 'Similarity: Kopper (p-value: D)', 'Similarity: Kopper (p-value: I)', 
                       'Niche Overlap: Buffer (D)', 'Niche Overlap: Buffer (I)', 'Niche Overlap: Biome (D)', 'Niche Overlap: Biome (I)', 
                       'Niche Overlap: Kopper (D)', 'Niche Overlap: Kopper (I)', 'Niche Dynamic: Buffer', 'Niche Dynamic: Biome', 'Niche Dynamic: Koppe')
write.csv(Results, 'E:/Expansion/Dataset/Dataset1_Ecospat_partial3.csv', row.names = F)

setwd('E:\\Expansion\\Dataset')
spvector <- dir()
i=1

for (i in 1:length(spvector)) {
  
  infolder <- paste('E:/Expansion/Ellipsenm', spvector[i], sep = "/")
  dir.create(infolder)
  
  #native
  sp <- read.csv(paste(spvector[i], "/", paste(spvector[i], "_n_thin.csv", sep = ""), sep = ""))
  write.csv(sp, paste(infolder,'/', spvector[i], "_n_thin.csv", sep = "" ), row.names = F)
  
  #invasive
  sp_i <- read.csv(paste(spvector[i], "/", paste(spvector[i], "_i_thin.csv", sep = ""), sep = ""))
  write.csv(sp_i, paste(infolder,'/', spvector[i], "_i_thin.csv", sep = "" ), row.names = F)
  
}




#Copying occurence for Ellipsenm calculations-------------
setwd('E:\\Expansion\\Dataset2')
spvector <- dir()
i=1

for (i in 1:length(spvector)) {
  
  infolder <- paste('E:/Expansion/Ellipsenm', spvector[i], sep = "/")
  dir.create(infolder)
  
  #native
  sp <- read.csv(paste(spvector[i], "/", paste(spvector[i], "_n_thin.csv", sep = ""), sep = ""))
  write.csv(sp, paste(infolder,'/', spvector[i], "_n_thin.csv", sep = "" ), row.names = F)
  
  #invasive
  sp_i <- read.csv(paste(spvector[i], "/", paste(spvector[i], "_i_thin.csv", sep = ""), sep = ""))
  write.csv(sp_i, paste(infolder,'/', spvector[i], "_i_thin.csv", sep = "" ), row.names = F)
  
}


#Crop layers to G space - work space ----------
require(rgdal)
require(raster)
setwd('E:\\NicheConsv\\LeafBeetles\\Data')
var<-c("bio_2","bio_4","bio_5","bio_6","bio_7","bio_13","bio_14", "bio_15")
var <- list.files(path = "E:\\NicheConsv\\LeafBeetles\\Data\\Layers_8_TIF\\", pattern = ".tif", full.names = TRUE)
var <- stack(var)
shape <- readOGR(dsn = "E:\\NicheConsv\\LeafBeetles\\Data", layer = "ProjectG")
plot(shape)
var <- mask(crop(var, shape), shape)
rnames <- paste0("E:\\NicheConsv\\LeafBeetles\\Data\\Project_G\\", names(var), ".asc") # users select the format
## saving layers in new folder
setwd('E:\\NicheConsv\\LeafBeetles\\Data\\Project_G')
sav <- lapply(1:nlayers(var), function(x) {
  writeRaster(var[[x]], filename = rnames[x], format = "ascii") # change format accordingly
})

#Ellipsenm 
#GAM 
#Bayesian regression


#GLM--------
install.packages('geosphere')
require(geosphere)
require(rgeos)
require(mgcv)
install.packages('dplyr')
require(dplyr)
require(corrplot)
install.packages('splines')
require(splines)
install.packages('ggplot2')
require(ggplot2)
c <- read.csv('E:\\Ecospat\\Occ\\All.csv')
inv <- read.csv('E:\\Ecospat\\Occ\\Invasive.csv')
ntv <- read.csv('E:\\Ecospat\\Occ\\Native.csv')


#Climatic variables-------
AnnTemp <- raster('E:\\Regression\\Vars\\Climate\\bio_1.asc')
MaxTempWarMth <- raster('E:\\Regression\\Vars\\Climate\\bio_4.asc')
TempAnnRange <- raster('E:\\Regression\\Vars\\Climate\\bio_7.asc')

#Humans impact variables -------
#Port distances
prt <- read.csv('E:\\Regression\\Vars\\Ports\\Portsxy.csv')
#Human Influence Index (HII) 1km
hii <- raster('E:\\NicheConsv\\shapes\\hii-global-geo-grid\\hii-global-geo-grid\\hii_v2geo')
#Croplands 1km
cr <- raster('E:\\NicheConsv\\shapes\\gl-croplands-geotif\\cropland.tif')
resample(cr, hii, method="bilinear", filename="E:\\NicheConsv\\shapes\\gl-croplands-geotif\\croplandResampled.tif")
crR<-raster('E:\\Regression\\Vars\\Cropland\\CropsResampled.tif')
#Population density 1km 2010
pop <- raster('E:\\Regression\\Vars\\PopDensity2020\\gpw-v4-population-density-rev11_2020_30_sec_tif\\PopDensity2020.tif')
#Human modification of terrestrial ecosystems
HTE <- raster('E:\\Regression\\Vars\\HUman Modification\\lulc-human-modification-terrestrial-systems-geographic-geotiff\\HMTE.tif')

  
Ports <- SpatialPoints(prt[,5:6])
Species <- SpatialPoints(ntv[,3:4])
ntv$DistPorts <- apply(gDistance(Ports, Species, byid=TRUE), 1, which.min)
hist(ntv$DistPorts)

Atemp <- extract(AnnTemp, ntv[,3:4])
MTmpWrmMth <- extract( MaxTempWarMth, ntv[,3:4])
TmpAnRg <- extract(TempAnnRange, ntv[,3:4])

Hii <- extract( hii, ntv[,3:4])
Crp <- extract( crR,ntv[,3:4])

Pop <- extract( pop, ntv[3:4])
Hte <- extract( HTE,ntv[,3:4])

#Dataset----------
inv.GAM <- na.omit(cbind(inv,Atemp, MTmpWrmMth, TmpAnRg, Hii,Crp, Pop, Hte))
natv.GAM<- na.omit(cbind(ntv, Atemp, MTmpWrmMth, TmpAnRg, Hii,Crp, Pop, Hte))
All.GAM <- rbind(inv.GAM, natv.GAM)
head(All.GAM)

levels(All.GAM$status)[match("invasive",levels(All.GAM$status))] <- "1"
levels(All.GAM$status)[match("native",levels(All.GAM$status))] <- "2"

#Correlation of Variables------
head(All.GAM)
Summary <-select_if(All.GAM, is.numeric)
summary(Summary)
hist(inv.GAM$Expansion2)
hist(inv.GAM$DistPorts)
hist(inv.GAM$Atemp)
hist(inv.GAM$Hii)

analysis.cols <- c('Atemp', 'MTmpWrmMth', 'TmpAnRg', 'Hii','Crp', 'Pop', 'Hte')
cols <- All.GAM[,analysis.cols]
cols <- cols[complete.cases(cols),]
r2 <- corrplot(cor(cols), method='number', number.cex = .7) #high corr ATemp and MTmpWrmMth
res1 <- cor.mtest(cols, conf.level = .95)
res2 <- cor.mtest(cols, conf.level = .99)
corrplot(cor(cols), p.mat = res1$p, insig = "p-value", number.cex = .7)

#GLMs-----

#I wish I could see the difference between Inv x Nat on preferences for the parameters I chose above
#Followed Hill et al 2016. They say they tested GLMs and GAMs to find correlation between expanison ad unfilling of populations

install.packages('glmulti')
require(glmulti)

#head(All.GAM)

#Models with All.GAM--------
global <- glm(Expansion2~DistPorts+Atemp+TmpAnRg+
                Hii+Crp+Pop+Hte,data=All.GAM)
summary(global)
model <- glmulti(global, # use the model with built as a starting point
                      level = 1,  #  just look at main effects
                      crit="aicc")
weightable(model)

global2 <- glm(Unfilling2~DistPorts+Atemp+TmpAnRg+
                 Hii+Crp+Pop+Hte,data=All.GAM)
summary(global2)
model2 <- glmulti(global2, # use the model with built as a starting point
                 level = 1,  #  just look at main effects
                 crit="aicc")
weightable(model2)

#Models with All.GAM$ +status----------

global_ <- glm(Expansion2~DistPorts+Atemp+TmpAnRg+status+
                Hii+Crp+Pop+Hte,data=All.GAM)
summary(global_)
model_ <- glmulti(global_, # use the model with built as a starting point
                 level = 1,  #  just look at main effects
                 crit="aicc")
weightable(model_)

globa2_ <- glm(Unfilling2~DistPorts+Atemp+TmpAnRg+status+
                 Hii+Crp+Pop+Hte,data=All.GAM)
summary(globa2_)
model2_ <- glmulti(globa2_, # use the model with built as a starting point
                  level = 1,  #  just look at main effects
                  crit="aicc")
weightable(model2_)

g1_ <- glm(Expansion2 ~ 1 + status + DistPorts + Atemp + TmpAnRg + Crp + Hte, data = All.GAM)
summary(g1_)
g2_ <- glm(Unfilling2 ~ 1 + status + Atemp + TmpAnRg + Crp, data = All.GAM)
summary(g2_)


#Adding All.GAM$status - for that I had to transform binary catergorical to numeric-------
levels(All.GAM$status)[match("invasive",levels(All.GAM$status))] <- "1"
levels(All.GAM$status)[match("native",levels(All.GAM$status))] <- "2"

g1_ <- glm(Expansion2 ~ 1 + DistPorts + Atemp + TmpAnRg + Hii + Crp + Hte +status, data = All.GAM)
summary(g1)
g2_ <- glm(Unfilling2 ~ 1 + DistPorts + Atemp + Hii + Crp+ status, data = All.GAM)
summary(g2_)


#Models with inv.GAM--------
global3 <- glm(Expansion2~+DistPorts+Atemp+TmpAnRg+
                 Hii+Crp+Pop+Hte,data=inv.GAM)
summary(global3)
model3 <- glmulti(global3, # use the model with built as a starting point
                 level = 1,  #  just look at main effects
                 crit="aicc")
weightable(model3)

global4 <- glm(Unfilling2~+DistPorts+Atemp+TmpAnRg+
                 Hii+Crp+Pop+Hte,data=inv.GAM)
summary(global4)
model4 <- glmulti(global4, # use the model with built as a starting point
                  level = 1,  #  just look at main effects
                  crit="aicc")
weightable(model4)

#Best models with inv.GAM
g3 <- glm(Expansion2 ~ 1 + Atemp + TmpAnRg + Pop + Hte, data = inv.GAM)
summary(g3)
g4 <- glm(Unfilling2 ~ 1 + DistPorts + Atemp + TmpAnRg + Crp + Pop + Hte , data = inv.GAM)
summary(g4)


#Models with ntv.GAM-----------
global5 <- glm(Expansion2~DistPorts+Atemp+TmpAnRg+
                 Hii+Crp+Pop+Hte,data=natv.GAM)
summary(global5)
model4 <- glmulti(global5, # use the model with built as a starting point
                  level = 1,  #  just look at main effects
                  crit="aicc")
weightable(model4)

global6 <- glm(Unfilling2~DistPorts+Atemp+TmpAnRg+
                 Hii+Crp+Pop+Hte,data=natv.GAM)
summary(global6)
model5 <- glmulti(global6, # use the model with built as a starting point
                  level = 1,  #  just look at main effects
                  crit="aicc")
weightable(model5)

#Best models with ntv.GAM
g5 <- glm(Expansion2 ~ 1 + DistPorts + Atemp + TmpAnRg + Crp + Pop + Hte , data = natv.GAM)
summary(g5)
g6 <- glm(  Unfilling2 ~ 1 + DistPorts + Atemp + TmpAnRg + Hii, data = natv.GAM)
summary(g6)

#Visualization of GLMs----------
install.packages('visreg')
library(visreg)
visreg(g1)
visreg(g2)
visreg(g3)
visreg(g4)
visreg(g5)
visreg(g6)
visreg(g1_)
visreg(g2_)
install.packages('bbmle')
library(bbmle)


#GAMs---------
install.packages("devtools")
devtools::install_github("cardiomoon/ggiraphExtra")
install.packages('purr')
install.packages('qpcR')
library(qpcR)
library(mgcv)
library(ggplot2)
install.packages('ggeffects')
library(ggeffects)
install.packages('DHARMa')
library(DHARMa)
install.packages('knitr')
library(knitr)
head(All.GAM)
null <-  gam(Expansion2 ~ 1, data = inv.GAM )
gam_DistPorts <- gam(Expansion2 ~ 1 + s(DistPorts), data = inv.GAM)
gam_Atemp <- gam(Expansion2 ~ 1 + s(Atemp), data = inv.GAM)
gam_TmpAnRg <- gam(Expansion2 ~ 1 + s(TmpAnRg), data = inv.GAM)
gam_Hii <- gam(Expansion2 ~ 1 + s(Hii), data = inv.GAM)
gam_Crp <- gam(Expansion2 ~ 1 + s(Crp), data = inv.GAM)
gam_Pop <- gam(Expansion2 ~ 1 + s(Pop), data = inv.GAM)
gam_Hte <- gam(Expansion2 ~ 1 + s(Hte), data = inv.GAM)

pred.DistP <- ggpredict(gam_DistPorts)
pred.Atemp <- ggpredict(gam_Atemp)
pred.TmpAnRg <- ggpredict(gam_TmpAnRg)
pred.Hii <- ggpredict(gam_Hii)
pred.Crp <- ggpredict(gam_Crp)
pred.Pop <- ggpredict(gam_Pop)
pred.Hte <- ggpredict(gam_Hte)

plot(pred.DistP$DistPorts) + geom_point(aes(x = DistPorts,
                                      y = Expansion2), 
                                  data = inv.GAM)
plot(pred.Atemp$Atemp) + geom_point(aes(x = Atemp,
                                            y = Expansion2), 
                                        data = inv.GAM)
plot(pred.TmpAnRg$TmpAnRg) + geom_point(aes(x = TmpAnRg,
                                        y = Expansion2), 
                                    data = inv.GAM)
plot(pred.Hii$Hii) + geom_point(aes(x = Hii,
                                        y = Expansion2), 
                                    data = inv.GAM)
plot(pred.Crp$Crp) + geom_point(aes(x = Crp,
                                        y = Expansion2), 
                                    data = inv.GAM)

plot(pred.Pop$Pop) + geom_point(aes(x = Pop,
                                           y = Expansion2), 
                                       data = inv.GAM)

simulation.output_DistP <- simulateResiduals(fittedModel = gam_DistPorts, n = 250)
plot(simulation.output_DistP)
simulation.output_Atemp <- simulateResiduals(fittedModel = gam_Atemp, n = 250)
plot(simulation.output_Atemp)
simulation.output_TmpAnRg <- simulateResiduals(fittedModel = gam_TmpAnRg, n = 250)
plot(simulation.output_TmpAnRg)
simulation.output_Hii <- simulateResiduals(fittedModel = gam_Hii, n = 250)
plot(simulation.output_Hii)
simulation.output_crp <- simulateResiduals(fittedModel = gam_Crp, n = 250)
plot(simulation.output_crp)
simulation.output_pop <- simulateResiduals(fittedModel = gam_Pop, n = 250)
plot(simulation.output_pop)
simulation.output_Hte <- simulateResiduals(fittedModel = gam_Hte, n = 250)
plot(simulation.output_Hte)

pdip.models <- list(null = null, DistPort = gam_DistPorts, AnTemp = gam_Atemp,
                    TempAnnRang = gam_TmpAnRg, Hii= gam_Hii, Crops =gam_Crp , Pop=gam_Pop , Hii=gam_Hii , 
                    Hte = gam_Hte)

pdip.aic.df <- data.frame(Model = names(pdip.models), AICc = sapply(pdip.models, function(x) AICc(x)),
                          akaike.weights(sapply(pdip.models, function(x) AICc(x))))

pdip.aic.df <- pdip.aic.df[order(pdip.aic.df$AICc),]
pdip.aic.df$Cumulative.Weight <- cumsum(pdip.aic.df$weights)

kable(pdip.aic.df, row.names = FALSE)


summary(gam_all)
gam.check(gam_all)








head(All.GAM)
fit <- glm(status~Crp,data=All.GAM,family=binomial())
summary(fit)
confint(fit)
predict(fit, type="response") 

fit1 <- glm(Percentage~Crp+DistPorts,data=All.GAM)
summary(fit1)
confint(fit1)

fit2 <- glm(Percentage~Crp+DistPorts+Atemp,data=All.GAM)
summary(fit2)
confint(fit2)

fit3 <- glm(Percentage~Crp+DistPorts+Atemp+Hii,data=All.GAM)
summary(fit3)
confint(fit3)

fit4 <- glm(Percentage~Crp+DistPorts+Atemp+Hii+Crp,data=All.GAM)
summary(fit4)
confint(fit4)

fit5 <- glm(status~Crp+DistPorts+Atemp+Hii+Crp+Pop,data=All.GAM,family=binomial())
summary(fit5)
confint(fit5)

remotes::install_github("MathiasHarrer/dmetar")

data('MVRegressionData')
require(metafor)
mmi = multimodel.inference(TE = 'yi', seTE = 'sei', data = MVRegressionData,
                           predictors = c('pubyear', 'quality',
                                          'reputation', 'continent'))
# Print summary
summary(mmi)



m1 <- glm(Percentage ~ Atemp, family = "poisson", data=All.GAM)
m1pred <- predict(m1, type = "response")
AIC(m1)
par(mfrow = c(1,2))
plot(All.GAM$Atemp, All.GAM$Percentage, main = "DF = 1")
lines(All.GAM$Atemp, exp(mu), lwd = 2, col = "grey")
lines(All.GAM$Atemp, m1pred, col = "orange", lwd = 2)


summary(glm(status~(Percentage),data=All.GAM,family=binomial))
summary(glm(Percentage~bs(Percentage),data=All.GAM,family=binomial))





Expansion<- glm(Percentage ~ DistPorts+ Atemp + Hii+ Crp, family = "poisson", data = All.GAM)
summary(Expansion)

qplot(All.GAM$Atemp, predict(Expansion, All.GAM)) + geom_smooth(method = "glm")

Expansion.Atemp <- glm(Percentage ~ Atemp, family = "poisson", data = All.GAM)
summary(Expansion.Atemp)
Expansion.Atemp <- glm(Percentage ~ DistPorts+Atemp, family = "poisson", data = All.GAM)
summary(Expansion.Atemp)


test <- gam(Percentage ~s(MTmpWrmMth),data=All.GAM, family = "nb", method = "REML", select = TRUE)

summary(test)
gam.check(test)


test1 <- gam(Percentage ~s(MTmpWrmMth),data=All.GAM)
summary(test1)
gam.check(test1)


G1n <- gam(Expansion5 ~s(hii),data=natv.GAM)
G1all <- gam(Expansion5 ~s(hii),data=GAM)



G3 <- gam(Expansion5 ~s(pop),data=inv.GAM)
G4 <- gam(Expansion5 ~s(hii)+s(cr),data=inv.GAM)
G5 <- gam(Expansion5 ~s(hii)+s(cr),data=inv.GAM)
G6 <- gam(Expansion5 ~s(hii)+s(cr)+s(pop),data=inv.GAM, method = "REML", select = TRUE)



test1 <- gam(Percentage ~s(cr),data=inv.GAM, family = gaussian, method = "REML", select = TRUE)
test2 <- gam(Percentage ~s(cr),data=inv.GAM, family = 'nb', method = "REML", select = TRUE)
test3 <- gam(Percentage ~s(cr),data=inv.GAM)
summary(test1)
summary(test2)
summary(test3)
gam.check(test1)
par(mfrow=c(2,2))
plot(test3,scheme=1,unconditional=TRUE) 

plot(G6,pages=1,residuals=TRUE)  ## show partial residuals
plot(G6,pages=1,seWithMean=TRUE)
gam.check(test1)
gam.check(test2)
gam.check(test3, pch=19,cex=.3)


if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
require(ggpubr)
require("dplyr")

ggboxplot(t, x = "species", y = "status", 
          color = "species", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("ctrl", "trt1", "trt2"),
          ylab = "Weight", xlab = "Treatment")

test <- aov(hii ~ status, data = t)

t.test(t$cr~t$status) 



t.test(t$status, t$hii, alternative = c("two.sided"), paired = FALSE, var.equal = FALSE, conf.level = 0.95)
summary(test)
