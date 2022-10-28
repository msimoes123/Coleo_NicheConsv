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

#Extracting occurence csv for claudia ---------
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





#Ellipsenm Models 

#Creating the calibration area with pricipal components 

# PCA------
library(kuenm)

var_folder <- "E:\\Expansion\\WC_2.1\\Vars" # name of folder with variables to be combined in distinct sets
out_folder <- "E:\\Expansion\\WC_2.1\\PCS_Elps" # name of folder that will contain the sets
in_format <- "GTiff" # other options available are "GTiff" and "EHdr" = bil
out_format <- "GTiff" # other options available are "GTiff" and "EHdr" = bil
npcs <-  2 # number of pcs you want as rasters, if not defined all pcs are returned as rasters

# PCA of variables for models
kuenm_rpca(variables  = var_folder, in.format = in_format, var.scale = TRUE, write.result = TRUE, 
           out.format = out_format, out.dir = out_folder, n.pcs = npcs)


#Run Models: Ellipsenm
library(raster)
setwd('E:\\Expansion\\Dataset')

spvector <- dir()

var<-c("bio_2","bio_4","bio_5","bio_6","bio_7","bio_13","bio_14", "bio_15")
out_dir <- paste(spvector, "/", "Calibration_E", sep = "")
i=1
for (i in 1:length(spvector)) {
  # reading species data
  #native
  sp <- read.csv(paste(spvector[i], "/", paste(spvector[i], "_n_thin.csv", sep = ""), 
                       sep = ""))
  sp_nodup <- unique(sp)
  
  colnames(sp_nodup)
  
  # raster layers of environmental data (this ones are masked to the accessible area)
  # users must prepare their layers accordingly if using other data
  #Calibration areas
  WGS84 <- sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  
  Buffer_n <- stack(list.files(path = paste(spvector[i], "M_variables", 'BUFFAREA_n/', sep = "/"),
                               pattern = ".tif", full.names = TRUE))
  
  Biome_n <- stack(list.files(path = paste(spvector[i], "M_variables", 'BIOMEPly/', sep = "/"),
                              pattern = ".tif", full.names = TRUE))
  
  Koppen_n <- stack(list.files(path = paste(spvector[i], "M_variables", 'KOPPENPly/', sep = "/"),
                               pattern = ".tif", full.names = TRUE))
  
 
  # preparing training and testing data
  data_split <- split_data(sp_nodup, method = "random", longitude = "Longitude", 
                           latitude = "Latitude", train_proportion = 0.7, 
                           save = FALSE, name = "Split_Nat")
  
  # sets of variables (example)
  sets <- list(set_1 = c("bio_2","bio_4","bio_5","bio_6","bio_7","bio_13","bio_14", "bio_15")) # change as needed
  
  variable_sets <- prepare_sets(Buffer_n, sets)
 
  dir.create('E_calibration_results')
  
  # methods to create ellipsoids
  methods <- c("mve1")
  # model calibration process
  calib <- ellipsoid_calibration(data_split, species = "Species", longitude = "Longitude", 
                                 latitude = "Latitude", variables = variable_sets,
                                 methods = methods, level = 99, selection_criteria = "S_OR_P",
                                 error = 5, iterations = 100, percentage = 50,
                                 output_directory = out_dir[i], overwrite = TRUE)
 
  cat(i, "species of", length(spvector), "\n")
}
