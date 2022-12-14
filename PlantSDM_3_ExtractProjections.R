#### Plant Overlay SDMs ####
#renv::deactivate() 

### European Tree SDMS ####
cluster <- FALSE
if(cluster==FALSE){
  #library(neotoma)
  library(dplyr)
  library(data.table)
  library(pbapply)
  library(pbmcapply)
  library(tidyverse)
  library(viridis)
  library(plyr)
  library(raster)
  library(rgeos)
  library(data.table)
  library(magick)
  library(gganimate)
  library(ncdf4)
  library(chron)
  library(foreign)
  library(rgdal)
  library(cowplot)
  library(rcarbon)
  #library(Bchron)
  library(cooccur)
  library(ggalt)
  #library(devtools)
  #install_github("jjvanderwal/climates")
  #library(climates)
  library(biomod2)
  library(sp)
  library(RColorBrewer)
  #library(geometry)
  library(usdm)
  #library(dismo)
  library(ecospat)
  library(sf)
  library(Cairo)
  library(terra)
  library(taxize)
  #devtools::install_github("barnabywalker/kewr")
  library(kewr)
  library(plyr)
  #devtools::install_github("biomodhub/biomod2", dependencies = TRUE)
  
  FP <- "/Users/hannahwauchope/Dropbox/Work/Projects/PaleoSDMs/"
  DataFP <- "/Users/hannahwauchope/Dropbox/Work/Data/"
  BaseFP <- "/Users/hannahwauchope/Dropbox/Work/Projects/"
  Computer <- "MacBookMax"
} else {
  library(sp)
  library(biomod2)
  library(terra)
  library(pbmcapply, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(stringr, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(data.table, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(pbapply, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  FP <- "/gpfs/ts0/projects/Research_Project-T115802/EuropeTrees/"
  DataFP <- "/gpfs/ts0/projects/Research_Project-T115802/Data/"
  
  Computer <- "Cluster"
}

AgeConf <- 1000
PollenThresh <- 5
AgeBin <- "Bin100"

WGSCRS <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
MollCRS <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
RobCRS <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

ncores <- ifelse(Computer=="MacMini", 6, ifelse(Computer=="MacBookMax", 7, ifelse(Computer=="MacBookPro", 6, 16)))
PastFP <- paste0(DataFP, "PalaeoView/Downscale_21000to0BP_20yrInterval_20yrMean/BioClim")

#### Functions ####

source("~/Documents/GitHub/UsefulFunctions/GetLoadGBIF.r")
source("~/Documents/GitHub/UsefulFunctions/SquareRasterCells.R")
source("~/Documents/GitHub/UsefulFunctions/CRUClimateData.R")
source("~/Documents/GitHub/UsefulFunctions/VisualiseYourSpatialData.R")

#### Let's extract points ####

#Load pollen data and gbif data
FosPol <- read.csv(paste0(FP, "PollenData/FosPolCleaned_Harmonised_100Samples_10TimeSteps.csv"))
load(file=paste0(FP, "GBIFDownloads/GBIFReduce.RData"))

#Get the names of all the species for which the SDMs have run through 
FinishedSpec <- gsub("_Finish.csv", "", list.files(paste0(FP, "SDMs/Log"), pattern = "*Finish.csv"))
FinishedSpec <- gsub("[.]", " ", FinishedSpec)

#First, extract the data
PolSDMs <- pblapply(FinishedSpec, function(Spec){
  print(Spec)
  FosPolSpec <- unique(subset(FosPol, Species==Spec)) #Fossil species
  FosPolSpec <- FosPolSpec[,c("Species", "Dataset_Sample", "long", "lat", "age", "x12_5_percent", "x87_5_percent")]
  
  GBIFSpec <- GBIFRed[GBIFRed$species==Spec] #GBIF species
  GBIFSpec <- GBIFSpec[,c("species", "YearGroup", "decimalLongitude", "decimalLatitude")]
  names(GBIFSpec) <- c("Species", "Ages", "x", "y")
  
  #Get all the uncertainty ages for each species in each time step
  FPSAges <- rbindlist(lapply(1:nrow(FosPolSpec), function(x){
    FPS <- FosPolSpec[x,]
    Ages <- seq(round_any(FosPolSpec[x,]$x12_5_percent, 100), round_any(FosPolSpec[x,]$x87_5_percent, 100), 100)
    FPS <-   FPS[rep(seq_len(nrow(FPS)), each = length(Ages)), ]
    FPS$Ages <- Ages
    return(FPS)
  }))
  
  #Extract projection data for all those ages (both the prob of occurrence from the SDM, and the "clamp" returns 0, 1, 2, or 3 for if 0, 1, 2 or 3 PCA axes are outside the training data)
  AgeSDMs <- rbindlist(pblapply(unique(FPSAges$Ages), function(Age){
    AgeSub <- FPSAges[FPSAges$Ages==Age,]
    AgeSub <- vect(AgeSub, geom=c("long", "lat"), crs=WGSCRS)
    AgeSub <- terra::project(AgeSub, MollCRS)
    if(!file.exists(paste0(FP, "/SDMs/", gsub(" ", ".", Spec), "/HindcastRasters/Hindcast_", Age, ".tif"))){
      return(NULL)
    }
    SDM <- rast(paste0(FP, "/SDMs/", gsub(" ", ".", Spec), "/HindcastRasters/Hindcast_", Age, ".tif"))
    Clamp <- rast(paste0(FP, "/SDMs/", gsub(" ", ".", Spec), "/HindcastRasters/Hindcast_", Age, "_Clamp.tif"))
    AgeSub$SDM <- terra::extract(SDM, AgeSub, ID=FALSE)$layer
    AgeSub$Clamp <- terra::extract(Clamp, AgeSub, ID=FALSE)$layer
    return(cbind(as.data.frame(AgeSub), crds(AgeSub)))
  }))
  
  ### Extract current projections from GBIF data
  CurrentSDMs <- rbindlist(pblapply(c("1901_1950", "1951_2000", "2001_2017"), function(Age){
    AgeSub <- GBIFSpec[GBIFSpec$Age==Age]
    if(nrow(AgeSub)==0){
      return(NULL)
    }
    if(!file.exists(paste0(FP, "/SDMs/", gsub(" ", ".", Spec), "/HindcastRasters/Hindcast_", Age, ".tif"))){
      return(NULL)
    }
    AgeSub <- vect(AgeSub, geom=c("x", "y"), crs=WGSCRS)
    AgeSub <- terra::project(AgeSub, MollCRS)
    SDM <- rast(paste0(FP, "/SDMs/", gsub(" ", ".", Spec), "/HindcastRasters/Hindcast_", Age, ".tif"))
    Clamp <- rast(paste0(FP, "/SDMs/", gsub(" ", ".", Spec), "/HindcastRasters/Hindcast_", Age, "_Clamp.tif"))
    AgeSub$SDM <- terra::extract(SDM, AgeSub, ID=FALSE)$layer
    AgeSub$Clamp <- terra::extract(Clamp, AgeSub, ID=FALSE)$layer
    return(cbind(as.data.frame(AgeSub), crds(AgeSub)))
  }))
  CurrentSDMs$Dataset_Sample <- "GBIF"
  
  #Bring together, done
  SDMExtracts <- rbindlist(list(AgeSDMs, CurrentSDMs), fill=T)
  return(SDMExtracts)
}) #, mc.cores=ncores

#Bring together the data for all the species
PolSDMs2 <- rbindlist(PolSDMs)

#Scale prob of occurrence down to proper numebrs (i.e. 0 --> 1)
PolSDMs2$SDMScale <- PolSDMs2$SDM/1000

#Add a value to each row that indicates how for away from the centre the age estimate is (if the age estimate is the central main estimate, it's 1, scales down from there towards 0 for ages above/below)
PolSDMs3 <- rbindlist(pblapply(unique(PolSDMs2$Dataset_Sample), function(DS){
  if(DS=="GBIF"){ #If the sample's from GBIF there's no uncertainty so done  
    PDS <- subset(PolSDMs2, Dataset_Sample==DS)
    PDS$AgeScale <- 1
    return(PDS)
  }
  PDS <- subset(PolSDMs2, Dataset_Sample==DS)
  PDS$Ages <- as.numeric(as.character(PDS$Ages))
  PDSCentre <- round_any(unique(PDS$age), 100)
  
  PDS$AgeScale <- sapply(PDS$Ages, function(y){
    if(y < PDSCentre){
      1 - (abs((y-PDSCentre)) / (abs(min(PDS$Ages)-PDSCentre)+100))
    } else if (y == PDSCentre){
      1
    } else {
      1 - (abs(y-PDSCentre)  / (abs(max(PDS$Ages)-PDSCentre)+100))
    }
  })
  return(PDS)
}))
dir.create(paste0(FP, "SDMs/Output/"), showWarnings=FALSE)
write.csv(PolSDMs3, paste0(FP, "SDMs/Output/SDMSampleOutput.csv"), row.names=FALSE)

#Plot
ggplot(PolSDMs3[PolSDMs3$Dataset_Sample!="GBIF",], aes(x=age, y=SDMScale))+
  geom_point(data=subset(PolSDMs3[PolSDMs3$Dataset_Sample!="GBIF",], AgeScale==1))+
  geom_line(aes(group=Dataset_Sample), alpha=0.5)+
  geom_smooth(method="lm")+
  facet_wrap(~Species)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  theme_classic()

#### Have a go at taking weighting means of probability of occurrence (weigheted based on age)
PolSDMs3 <- read.csv(paste0(FP, "SDMs/Output/SDMSampleOutput.csv"))
PolSDMs3 <- as.data.table(PolSDMs3)

library(dplyr)
library(modi)

FossilPol <- PolSDMs3[PolSDMs3$Dataset_Sample!="GBIF",]
FossilPol <- data.table(FossilPol)

#THIS ISN'T WORKING FOR SOME REASON 
FossilPol <- FossilPol %>%
  group_by(Species, Dataset_Sample, age) %>% 
  mutate(weighted_SDM = weighted.mean(SDMScale, AgeScale)) %>%
  mutate(SDM_var = sqrt(weighted.var(SDMScale, AgeScale)))

#Plot pollen data
ggplot(FossilPol[FossilPol$Dataset_Sample!="GBIF",], aes(x=age, y=weighted_SDM,))+
  geom_point(data=subset(FossilPol[FossilPol$Dataset_Sample!="GBIF",], AgeScale==1))+
  geom_errorbar(aes(ymin=(weighted_SDM - SDM_var), ymax=(weighted_SDM + SDM_var)), alpha=0.5)+
  geom_smooth(method="lm")+
  facet_wrap(~Species)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  theme_classic()

#Plot present day data 
GBIFPol <- PolSDMs3[PolSDMs3$Dataset_Sample=="GBIF",]
GBIFPol$Ages2 <- as.numeric(str_split_fixed(GBIFPol$Ages, "[_]", 2)[,1])
GBIFPol$Ages <- factor(GBIFPol$Ages, levels=c("1901_1950", "1951_2000", "2001_2017"))
ggplot(GBIFPol, aes(x=Ages2, y=SDMScale,))+
  geom_boxplot(data=GBIFPol, aes(group=Ages))+
  #geom_errorbar(aes(ymin=(weighted_SDM - SDM_var), ymax=(weighted_SDM + SDM_var)), alpha=0.5)+
  geom_smooth(method="lm")+
  facet_wrap(~Species)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  theme_classic()

#### Animate points! ####
# Create animations of SDM output, with fossil locations overlayed for the past, and GBIF locations for the present
#Running this for many species
PolSDMs3 <- read.csv(paste0(FP, "SDMs/Output/SDMSampleOutput.csv"))
load(file=paste0(FP, "GBIFDownloads/GBIFReduce.RData"))

SpecToModel <- c("Populus.tremuloides", "Betula.pendula", "Juniperus.communis", "Pinus.contorta", "Alnus.rubra", "Centaurea.scabiosa")
SpecToModel <- gsub("[.]", " ", SpecToModel)
SpecToModel <- data.frame(Species=SpecToModel, Region=c("North America", "Europe Russia Asia Middle East", "Europe Russia North America", "North America", "North America", "Europe Russia Asia Middle East"))

pblapply(1:nrow(SpecToModel), function(STM){
  Spec <- SpecToModel[STM,]$Species
  Region <- SpecToModel[STM,]$Region
  print(Spec)
  #Get the region for the polygon (we're just using one static one)
  if(Region=="Europe Russia Asia Middle East"){
    AnimPolys <- GetThatBasemap(WorldPolyFP, Region, c(-30, 180, -11, 82))
  } else {
    AnimPolys <- GetThatBasemap(WorldPolyFP, Region)
  }
  
  #Get the points for the relevant species, do a bit of cleaning
  AnimPoints <- subset(PolSDMs3, Species==Spec & AgeScale==1)
  AnimPoints <- vect(AnimPoints, geom=c("x", "y"), crs=MollCRS)
  AnimPoints <- project(AnimPoints, WGSCRS)
  AnimPoints <- crop(AnimPoints, AnimPolys)
  AnimPoints <- AnimPoints[!is.na(AnimPoints$SDMScale),]
  names(AnimPoints)[names(AnimPoints)=="Clamp"] <- "NumAxesOutsideRange"
  AnimPoints$NumAxesOutsideRange <- as.factor(AnimPoints$NumAxesOutsideRange)
  AnimPoints$PointVal <- factor(round_any(AnimPoints$SDMScale, 0.25), levels=c(0,0.25, 0.5, 0.75, 1))
  
  #Get the time steps
  TimeSteps <- c(rev(sort(as.numeric(unique(AnimPoints[AnimPoints$Dataset_Sample!="GBIF"]$Ages)))), unique(AnimPoints[AnimPoints$Dataset_Sample=="GBIF"]$Ages))
  TimeSteps <- as.character(TimeSteps)
  
  #Split points vector into a list of vectors based on timesteps. 
  AnimPoints <- pblapply(TimeSteps, function(TS) AnimPoints[AnimPoints$Ages == TS,])
  names(AnimPoints) <- TimeSteps
  
  print("Compile Rasters")
  
  AnimRasts <- pblapply(TimeSteps, function(TS){
    rasty <- rast(paste0(FP, "/SDMs/", gsub(" ", ".", Spec), "/HindcastRasters/Hindcast_", TS, ".tif")) #Read in raster
    rasty <- project(rasty, WGSCRS) #Project
    if(!is.null(AnimPolys)){ #Crop to same region as the polygon
      rasty <- crop(rasty, AnimPolys)
      rasty <- mask(rasty, AnimPolys)
    }
    rasty[!is.na(rasty)] <- round_any(rasty[!is.na(rasty)], 100)/1000 #Convert raster values into 10 breaks (0,1,0.1)
    return(rasty)
  }) #, mc.cores=ncores
  names(AnimRasts) <- TimeSteps
  TimeSteps <- as.character(TimeSteps)
  
  print("Create Images")
  AnimPointsPast <- AnimPoints[!names(AnimPoints) %in% unique(PolSDMs3[PolSDMs3$Dataset_Sample=="GBIF",]$Ages)]
  AnimPointsPresent<- AnimPoints[names(AnimPoints) %in% unique(PolSDMs3[PolSDMs3$Dataset_Sample=="GBIF",]$Ages)]
  AnimRastsPast <- AnimRasts[!names(AnimRasts) %in% unique(PolSDMs3[PolSDMs3$Dataset_Sample=="GBIF",]$Ages)]
  AnimRastsPresent<- AnimRasts[names(AnimRasts) %in% unique(PolSDMs3[PolSDMs3$Dataset_Sample=="GBIF",]$Ages)]
  TimeStepsPast <- TimeSteps[!TimeSteps %in% unique(PolSDMs3[PolSDMs3$Dataset_Sample=="GBIF",]$Ages)]
  TimeStepsPresent <- TimeSteps[TimeSteps %in% unique(PolSDMs3[PolSDMs3$Dataset_Sample=="GBIF",]$Ages)]
  
  names(AnimPointsPast) <- as.character(as.numeric(names(AnimPointsPast))/1000)
  names(AnimRastsPast) <- as.character(as.numeric(names(AnimRastsPast))/1000)
  TimeStepsPast <- as.character(as.numeric(TimeStepsPast)/1000)
  
  CreateThoseImages(AnimPointsPast, AnimPolys, AnimRastsPast, TimeStepsPast, MapOrder=c("Poly", "Rast", "Points"), AnimFP = paste0(FP, "/Animations/", Spec, "/"),
                    PointColVal = "PointVal", 
                    PointColName = "Probability of Occurrence",
                    PointCols = brewer.pal(5, "YlOrRd"))
  CreateThoseImages(AnimPointsPresent, AnimPolys, AnimRastsPresent, TimeStepsPresent, MapOrder=c("Poly", "Rast", "Points"), AnimFP = paste0(FP, "/Animations/", Spec, "/"),
                    PointColVal = "PointVal", 
                    PointColName = "Probability of Occurrence",
                    PointCols = brewer.pal(5, "YlOrRd"),
                    PointSize = 0.5,
                    overwrite=FALSE, 
                    TimeStepTitle=TRUE,
                    SaveReps=5)#
  MakeThatAnimation(AnimFP = paste0(FP, "/Animations/", Spec, "/"), FramesPerSecond=5, SaveName=Spec)
}) #, mc.cores=6





