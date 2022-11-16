####Conservatism of invasive leaf beetles####
# Authors of the code: Marianna Simoes, Claudia Nunez-Penichet and Marcus Krull
#Date: 2021
####################################################################################
###############

# packages


# Working directory
setwd('C:\\Users\\mari1\\Dropbox (Senckenberg)\\Projects\\Beetle_niche\\Data\\Dataset2')

#species <- as.character(read.csv("species_taxmtch.csv")[, 1])

#Creating a big table with all species
#all_species <- do.call(rbind, all_species)
#rite.csv(all_species, "all_species2.csv", row.names = FALSE)

#Selecting species

#Creating M -----------
install.packages("rangemap")
require(rangemap)
require(maptools)
require(sf)
require(rgdal)
require(raster)
require(rgeos)
require(ecospat)
#install.packages("devtools")
require(devtools)
require(ecospat)
require(raster)
require(rworldmap)
require(ellipsenm)
require(rgdal)
require(spocc)
require(ade4)

clim <- stack(list.files(path = "C:\\Users\\mari1\\Dropbox (Senckenberg)\\Projects\\Resampling\\Bootstrap_UE\\Climate\\", pattern = ".tif", full.names = TRUE))
crs(clim) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
spvector <- dir()
Results <- NULL
PCS<-NULL
i=1
#Ecospat------------
#I already did the ENM, now Ill check for similarity between Inv and Nativ?
#Then hypervolume: Is the climate niche different between Inv /Nativ?
#Then niche breath to argue for the expansion 

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
  Biome_n <- stack(list.files(path = paste(spvector[i], "M_variables", 'BIOMEPly/', sep = "/"),
                              pattern = ".tif", full.names = TRUE))
  Biome_i <- stack(list.files(path = paste(spvector[i], "M_variables", 'BIOMEPlyi/', sep = "/"),
                              pattern = ".tif", full.names = TRUE))
 
  #Call background calibration------
  #M=Biome  
  env.bkg.Bio_n <- na.exclude(data.frame(rasterToPoints(Biome_n))) 
  env.bkg.Bio_n <- env.bkg.Bio_n[,3:10]
  env.bkg.Bio_n$Species <- 'Native'
  
  env.bkg.Bio_i <- na.exclude(data.frame(rasterToPoints(Biome_i))) 
  env.bkg.Bio_i <- env.bkg.Bio_i[,3:10]
  env.bkg.Bio_i$Species <- 'Invasive'
  
  #Occ environemnt: Background----
  env.occ.Bkg_n<- extract(clim, sp_nodup)#Native 
  env.occ.Bkg_n <- na.exclude(cbind(env.occ.Bkg_n,sp_nodup))
  
  env.occ.Bkg_i <-extract(clim,sp_nodup_i)#Invasive
  env.occ.Bkg_i <- na.exclude(cbind(env.occ.Bkg_i,sp_nodup_i))
  
  #calibration of PCA-env 
  
  #M=Biome  
  pca.env_Bio <-dudi.pca(rbind(env.bkg.Bio_n,env.bkg.Bio_i)[,1:8], center = T, scale = T, scannf = F, nf = 2)
  ecospat.plot.contrib(contrib=pca.env_Bio$co, eigen=pca.env_Bio$eig)
  jpeg(file=paste("C:\\Users\\mari1\\Desktop\\test\\", spvector[i],".jpeg", sep = ""))
  
  dev.off()
  
  #Biome Native  
  scores.bkg_B_n<- pca.env_Bio$li
  scores.bkg.Bio_n<- na.omit(suprow(pca.env_Bio,env.bkg.Bio_n[,1:8])$li)
  #Biome Invasive
  scores.bkg_B_i<- pca.env_Bio$li
  scores.bkg.Bio_i<- na.omit(suprow(pca.env_Bio,env.bkg.Bio_i[,1:8])$li)  
  
  #Scores for occurences 
 
  scores.occ.Bio_n<- na.omit(suprow(pca.env_Bio,env.occ.Bkg_n[,1:8])$li)
  scores.occ.Bio_i<- na.omit(suprow(pca.env_Bio,env.occ.Bkg_i[,1:8])$li)
 
  # calculation of occurence density

  z1_Bio_n<- ecospat.grid.clim.dyn(scores.bkg_B_n,scores.bkg.Bio_n,scores.occ.Bio_n, R=100 )
  z2_Bio_i<- ecospat.grid.clim.dyn(scores.bkg_B_i,scores.bkg.Bio_i,scores.occ.Bio_i,R=100)
  
  #Niche dynamics -------------------
  Dy_Bio05<-ecospat.niche.dyn.index(z1_Bio_n, z2_Bio_i,intersection=0.05)
  Dy_Bio15<-ecospat.niche.dyn.index(z1_Bio_n, z2_Bio_i,intersection=0.15)
  Dy_Bio25<-ecospat.niche.dyn.index(z1_Bio_n, z2_Bio_i,intersection=0.25)
  
  Results <- rbind(Results, data.frame(spvector[i],Dy_Bio05$dynamic.index.w, Dy_Bio15$dynamic.index.w, Dy_Bio25$dynamic.index.w))  
  PCS <- rbind(PCS, data.frame(spvector[i],pca.env_Bio$eig))
  cat(i, "of", length(spvector), "\n")  
}


colnames(Results) <- c('Species','Niche Dynamic: Biome05', 'Niche Dynamic: Biome15', 'Niche Dynamic: Biome25')

write.csv(Results, 'E:/Expansion/Dataset/Dataset2_Ecospat_partial3.csv', row.names = F)
