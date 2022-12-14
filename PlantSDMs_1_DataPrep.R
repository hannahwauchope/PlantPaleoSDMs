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

#### Get taxa info and distribution data ####
##Using Plants of the World Online (Kew) for taxa and distribution information
FosPol <- read.csv("/Users/hannahwauchope/Dropbox/Work/Projects/PaleoSDMs/PollenData/FosPolCleaned.csv")
FosPol$Species <- gsub("_", " ", FosPol$Species)  
FosPol$Species <- sapply(FosPol$Species, function(x) paste(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)), sep=""))
id.powos <- pblapply(unique(FosPol$Species), function(x) search_powo(x))
names(id.powos) <- unique(FosPol$Species)
id.powo <- rbindlist(pblapply(1:length(id.powos), function(x){
  print(x)
  NResults <- length(id.powos[[x]]$results)
  if(NResults == 0){
    Results <- data.frame(PowoSpec = NA, PowoID = NA)
  } else {
    Results <- rbindlist(lapply(1:NResults, function(y){
      Name <- id.powos[[x]]$results[[y]]$name
      ID <- tail(str_split(id.powos[[x]]$results[[y]]$fqId, "[:]")[[1]],1)
      return(data.frame(PowoSpec = Name, PowoID = ID))
    }))
  }
  Results$Species <- names(id.powos)[[x]]
  return(Results)
}))

#For each powo, get the accepted taxon ID
POWOIDsAccepted <- pblapply(1:nrow(id.powo), function(x){
  print(x)
  ID <- id.powo[x,]$PowoID
  if(is.na(ID)){
    Acp <- NA
  } else {
    POWO <- lookup_powo(ID)
    if(POWO$taxonomicStatus == "Accepted"){
      Acp <- ID
    } else if(length(POWO$accepted) == 0 & is.null(POWO$accepted[[1]])){
      Acp <- NA
    } else {
      Acp <- tail(str_split(POWO$accepted$fqId, "[:]")[[1]],1)
    }
  }
  IDDat <- id.powo[x,]
  IDDat$Accepted <- Acp
  return(IDDat)
})
POWOIDsAccepted2 <- rbindlist(POWOIDsAccepted)
POWOIDsAcceptedNAs <- POWOIDsAccepted2[is.na(POWOIDsAccepted2$PowoID),]
POWOIDsAccepted3 <- POWOIDsAccepted2[!POWOIDsAccepted2$Species %in% POWOIDsAcceptedNAs$Species]

#Ok first let's pull out the ones where everything agrees
POWOIDsAcceptedDone <- subset(POWOIDsAccepted3, PowoSpec==Species & PowoID==Accepted)
POWOIDsAccepted4 <- POWOIDsAccepted3[!POWOIDsAccepted3$Species %in% POWOIDsAcceptedDone$Species,]

#Next ones where there's just been a name switch (condition: two rows, ID is consistent for 3/4)
POWOIDsAcceptedSwitch <- POWOIDsAccepted4[POWOIDsAccepted4$Species %in% names(table(POWOIDsAccepted4$Species)[table(POWOIDsAccepted4$Species)==2]),]
POWOIDsAcceptedSwitch <- subset(POWOIDsAcceptedSwitch, PowoID==Accepted)
POWOIDsAcceptedSwitch <- POWOIDsAcceptedSwitch[!POWOIDsAcceptedSwitch$Species %in% POWOIDsAcceptedSwitch[duplicated(POWOIDsAcceptedSwitch$Species),]$Species,]

POWOIDsAccepted5 <- POWOIDsAccepted4[!POWOIDsAccepted4$Species %in% POWOIDsAcceptedSwitch$Species,]

#For the rest, it's gonna be a manual job

WrittenOut <- "DONE"
if(WrittenOut != "DONE"){
  write.csv(POWOIDsAccepted5, paste0(FP, "PollenData/POWOHarmonisation.csv"), row.names=FALSE)
}

Harmonised <- read.csv(paste0(FP, "PollenData/POWOHarmonisation.csv"))

#Ok wtf do we do now. Um. Ok um. We need to........ merge them?
POWOIDsAcceptedManualCheck <- cbind(POWOIDsAccepted5, Harmonised[,c("ManualSpec", "ManualAccepted")])
POWOIDsAcceptedManualCheckNAs <- POWOIDsAcceptedManualCheck[is.na(POWOIDsAcceptedManualCheck$ManualSpec),]
POWOIDsAcceptedManualCheck <- POWOIDsAcceptedManualCheck[!is.na(POWOIDsAcceptedManualCheck$ManualSpec),]
POWOIDsAcceptedManualCheck <- unique(POWOIDsAcceptedManualCheck[,c("Species", "ManualSpec", "ManualAccepted")])

#OKI DOKIE Let's bring it together baby
POWOIDsAcceptedDone <- POWOIDsAcceptedDone[,c("PowoSpec", "Species", "Accepted")]
names(POWOIDsAcceptedDone) <- c("Species", "GivenSpec", "POWOID")

POWOIDsAcceptedSwitch <- POWOIDsAcceptedSwitch[,c("PowoSpec", "Species", "Accepted")]
names(POWOIDsAcceptedSwitch) <- c("Species", "GivenSpec", "POWOID")

POWOIDsAcceptedManualCheck <- POWOIDsAcceptedManualCheck[,c("ManualSpec", "Species", "ManualAccepted")]
names(POWOIDsAcceptedManualCheck) <- c("Species", "GivenSpec", "POWOID")

POWOHarmonised <- rbindlist(list(POWOIDsAcceptedDone, POWOIDsAcceptedSwitch, POWOIDsAcceptedManualCheck))

#Rejected ones
Rej <- length(c(unique(POWOIDsAcceptedNAs$Species), unique(POWOIDsAcceptedManualCheckNAs$Species)))

#Check we haven't lost any species
if(length(unique(FosPol$Species)) != (Rej + length(unique(POWOHarmonised$GivenSpec)))){stop("some species have been lost in harmonisation")}

#### Get harmonised metadata and distributions

###
TaxaDistribution <- pblapply(unique(POWOHarmonised$POWOID), function(x){
  print(x)
  r <- lookup_powo(x, distribution = TRUE)
  if(length(r$distribution) == 0){
    return(NULL)
  }
  native <- sapply(r$distribution$natives, function(x) x$name)
  introduced <- sapply(r$distribution$introduced, function(x) x$name)
  Combine <- list(native, introduced)
  names(Combine) <- c("Native", "Introduced")

  TaxInfo <- data.frame(POWOID=x, Kingdom=r$kingdom, Phylum=r$phylum, Family=r$family)
  
  TaxDist <- list(Combine, TaxInfo)
  names(TaxDist) <- r$name
  return(TaxDist)
})

DistributionInfo <- lapply(TaxaDistribution, function(x) x[[1]])
names(DistributionInfo) <- sapply(TaxaDistribution, function(x) names(x)[[1]])

TaxaInfo <- rbindlist(lapply(TaxaDistribution, function(x) x[[2]]))

POWOHarmonised <- merge(POWOHarmonised, TaxaInfo, by="POWOID")

names(FosPol)[names(FosPol) == "Species"] <- "GivenSpec"
FosPol2 <- merge(FosPol, POWOHarmonised, by="GivenSpec")

write.csv(FosPol2, paste0(FP, "PollenData/FosPolCleaned_Harmonised.csv"))

#Get taxa with at least 100 samples per species
SamplePerSpec <- dcast(unique(FosPol2[,c("Species", "Dataset_Sample")]), Species ~ ., length)
SamplePerSpec <- subset(SamplePerSpec, . >= 100)

FosPol3 <- FosPol2[FosPol2$Species %in% SamplePerSpec$Species,]
write.csv(FosPol3, paste0(FP, "PollenData/FosPolCleaned_Harmonised_100Samples.csv"))

DistributionsSave <- DistributionInfo
tdwg_level3 <- vect(paste0(DataFP, "RandomCrap/tdwg wgsrpd master level3/level3.shp"))
sort(unique(tdwg_level3$LEVEL3_NAM))
#Correct any broken ones
DistributionInfo <- lapply(DistributionInfo, lapply, function(y) gsub("North European Russi", "North European Russia", y))
DistributionInfo <- lapply(DistributionInfo, lapply, function(y) gsub("North European Russiaa", "North European Russia", y))
DistributionInfo <- lapply(DistributionInfo, lapply, function(y) gsub("South European Russi", "South European Russia", y))
DistributionInfo <- lapply(DistributionInfo, lapply, function(y) gsub("Northwest European R", "Northwest European Russia", y))
DistributionInfo <- lapply(DistributionInfo, lapply, function(y) gsub("Northwest Territorie", "Northwest Territories", y))
DistributionInfo <- lapply(DistributionInfo, lapply, function(y) gsub("Central European Rus", "Central European Russia", y))
DistributionInfo <- lapply(DistributionInfo, lapply, function(y) gsub("Kirgizstan", "Kirgizistan", y))
DistributionInfo <- lapply(DistributionInfo, lapply, function(y) gsub("Central African Repu", "Central African Republic", y))
DistributionInfo <- lapply(DistributionInfo, lapply, function(y) gsub("Gambia", "Gambia, The", y, fixed=TRUE))
DistributionInfo <- lapply(DistributionInfo, lapply, function(y) gsub("Mozambique Channel I", "Mozambique Channel Is.", y))
DistributionInfo <- lapply(DistributionInfo, lapply, function(y) gsub("Panamá", "Panama", y))
DistributionInfo <- lapply(DistributionInfo, lapply, function(y) gsub("Suriname", "Surinam", y))
DistributionInfo <- lapply(DistributionInfo, lapply, function(y) gsub("Leeward Is.", "Leeward Is. AB Ant", y, fixed=TRUE))
DistributionInfo <- lapply(DistributionInfo, lapply, function(y) gsub("Amsterdam-St.Paul Is", "Amsterdam-St.Paul Is.", y, fixed=TRUE))
DistributionInfo <- lapply(DistributionInfo, lapply, function(y) gsub("Central American Pac", "C. American Pacific Is.", y, fixed=TRUE))
DistributionInfo <- lapply(DistributionInfo, lapply, function(y) gsub("Marion-Prince Edward", "Marion-Prince Edward Is.", y, fixed=TRUE))
DistributionInfo <- lapply(DistributionInfo, lapply, function(y) gsub("Cocos (Keeling) Is.", "Cocos (Keeling) I.", y, fixed=TRUE))

DistributionSHPs <- pblapply(1:length(DistributionInfo), function(x){
  print(x)
  if(length(DistributionInfo[[x]]$Native)>0){
    NatPol <- tdwg_level3[tdwg_level3$LEVEL3_NAM %in% DistributionInfo[[x]]$Native,]
    if(length(DistributionInfo[[x]]$Native) != length(unique(NatPol$LEVEL3_NAM))){stop("Countries lost")
      DistributionInfo[[x]]$Native[!DistributionInfo[[x]]$Native %in% unique(NatPol$LEVEL3_NAM)]
    }
    NatPol$Status <- "Native"
  } else {NatPol <- NULL}
  if(length(DistributionInfo[[x]]$Introduced)>0){
    IntPol <- tdwg_level3[tdwg_level3$LEVEL3_NAM %in% DistributionInfo[[x]]$Introduced,]
    if(length(DistributionInfo[[x]]$Introduced) != length(unique(IntPol$LEVEL3_NAM))){stop("Countries lost")
      DistributionInfo[[x]]$Introduced[!DistributionInfo[[x]]$Introduced %in% unique(IntPol$LEVEL3_NAM)]
    }
    IntPol$Status <- "Introduced"
  } else {IntPol <- NULL}
  SpecPol <- rbind(NatPol, IntPol)
  SpecPol$Species <- names(DistributionInfo)[[x]]
  
  writeVector(SpecPol, file=paste0(FP, "POWODistributions/", names(DistributionInfo)[[x]], ".shp"), overwrite=TRUE)
  
  return(SpecPol)
})

#### Get and save european tree polgons (for cross checking later) ####
#European tree distributions
#Polygon features (plg) - continuous areas of occupancy
#Point features (pnt) - more fragmented/isolated populations
#Synanthropic occurrences (syn, plg and/or pnt) - Reported synanthropic occurrences outside natural range, i.e. areas where it is introduced and naturalised

#First list the tree polygon types, and make a lil dataframe of them
TreeFiles <- list.files(paste0(DataFP, "EuropeanTreeDistributions/chorological_maps_dataset"), pattern="*.shp", recursive = TRUE, full.names=TRUE)
TreeSpecies <- str_split_fixed(list.files(paste0(DataFP, "EuropeanTreeDistributions/chorological_maps_dataset"), pattern="*.shp", recursive = TRUE), "[/]", 2)[,1]
TreeTypes <- str_split_fixed(list.files(paste0(DataFP, "EuropeanTreeDistributions/chorological_maps_dataset"), pattern="*.shp", recursive = TRUE), "[_]", 3)[,3]
TreeTypes <- gsub(".shp", "", TreeTypes)
TreeMeta <- data.frame(cbind(TreeSpecies, TreeTypes, TreeFiles))
TreeMeta$PointPoly <- ifelse(grepl("plg", TreeMeta$TreeTypes), "Polygon", "Point")
TreeMeta$NatSyn <- ifelse(grepl("syn", TreeMeta$TreeTypes), "Synanthropic", ifelse(grepl("sym", TreeMeta$TreeTypes), "Synanthropic", "Native"))
TreeMeta$SubSpecies <- gsub("sym", "", gsub("_", "", gsub("syn", "", gsub("pnt", "", gsub("plg", "", TreeMeta$TreeTypes)))))
TreeMeta[TreeMeta$SubSpecies == "mugorotundata",]$SubSpecies <- "mugo rotunda"
TreeMeta$TreeTypes <- NULL

TreeMeta <- TreeMeta[TreeMeta$PointPoly == "Polygon",]
TreePolys <- lapply(TreeMeta$TreeFiles, vect)
names(TreePolys) <- TreeMeta$TreeSpecies

TreePolysFinal <- lapply(unique(TreeMeta$TreeSpecies), function(x){
  Polys <- TreePolys[names(TreePolys)==x]
  Polys <- lapply(Polys, function(y){
    MetaDat <- subset(TreeMeta, TreeFiles==sources(y))
    y$Species <- MetaDat$TreeSpecies
    y$Range <- MetaDat$NatSyn
    return(y)
  })
  i <- length(Polys)
  if(i>1){
    PolyUnion <- rbind(Polys[[i]], Polys[[i-1]])
    i <- i - 1
    while(i>1){
      PolyUnion <- rbind(PolyUnion, Polys[[i-1]])
      i <- i - 1
    }
    Polys <- PolyUnion
  } else {
    Polys <- Polys[[1]]
  }
  writeVector(Polys, paste0(DataFP, "EuropeanTreeDistributions/TreePolygons/", x, ".shp"), overwrite=TRUE)
  return("Done")
})

#### Download GBIF data ####
FosPol <- read.csv(paste0(FP, "PollenData/FosPolCleaned_Harmonised_100Samples.csv"))

GBIFuser <- "hannahwauchope"
GBIFpwd <- "abv2xHAaQ4hBXgLePxM6taACFkevVgymAvYK8KuuMepcTDk8aBwj8CkqJgXn3sdDX2"
GBIFemail <- "hannah.wauchope@gmail.com"

#GBIF/R don't respond well to parallelisation, so I have done it the janky way of just opening multiple R sessions
pblapply(tolower(unique(FosPol$Species)), function(x) GetGBIFData(x, paste0(FP, "GBIFDownloads/"), GBIFuser, GBIFpwd, GBIFemail, cluster=TRUE))

LogFiles <- list.files(paste0(FP, "GBIFDownloads/LogBook"))
LogFiles <- gsub(" Begin", "", LogFiles)
LogFiles <- gsub(" Done", "", LogFiles)
LogFiles <- gsub(" NotFoundinGBIF", "", LogFiles)
LogFiles <- names(table(LogFiles)[table(LogFiles)==1])
LogFiles <- gsub(".csv", "", LogFiles)

LogFilesFull <- list.files(paste0(FP, "GBIFDownloads/LogBook"), full.names=TRUE)
LogFilesFull <- sapply(LogFiles, function(x) LogFilesFull[grep(paste0("*", x, "*"), LogFilesFull)], USE.NAMES = FALSE)
unlink(LogFilesFull)

LogDirs <- list.dirs(paste0(FP, "GBIFDownloads"))
LogDirs <- sapply(LogFiles, function(x) LogDirs[grep(paste0("*", x, "*"), LogDirs)])
unlink(LogDirs, recursive = TRUE)

FosPolSpecies <- data.frame(Species = unique(FosPol$Species), Cut = cut(1:length(unique(FosPol$Species)), 20, labels=paste0("Group", 1:20)))

pblapply(unique(FosPolSpecies$Cut)[13:20], function(cutty){
  CutSpec <- unique(subset(FosPolSpecies, Cut==cutty)$Species)
  FosPolSpecGBIF <- LoadGBIFData(CutSpec, paste0(FP, "GBIFDownloads/"), GBIFuser, GBIFpwd, GBIFemail, redownload = TRUE)
  write.csv(FosPolSpecGBIF, paste0(FP, "GBIFDownloads/GBIF", cutty, ".csv"), row.names = FALSE)
  rm(FosPolSpecGBIF)
  gc()
  return("done")
})

#### Read in and collate GBIF data, reduce to 1 occurrence per species/gridcell/year ####
#First, reduce occurrence data to one occurrence per grid cell/time point
#Run it by the 20 groups cos there's SO MUCH DATA

#First, make a unique grid cell raster and convert to mollweide
#First, change rasters to same res as climate data
TemplateRas <- rast(paste0(DataFP, "/Climate_CRU/BioClim/1951_2000BioClim.grd"))[[1]]
crs(TemplateRas) <- "+proj=longlat +datum=WGS84 +no_defs"
TemplateRas <- terra::project(TemplateRas, "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs")
TemplateRas <- squareRastCells(TemplateRas)
TemplateRas[!is.na(TemplateRas)][1:length(TemplateRas[!is.na(TemplateRas)])] <- 1:length(TemplateRas[!is.na(TemplateRas)])
TimeStamps <- list(seq(1901,1950,1), seq(1951,2000,1), seq(2001,2017,1))

#Now, read in occurrence data, and assign each to its gridcell
GBIFRed <- rbindlist(pbmclapply(c(paste0("Group", 1:20)), function(cutty){
  if(!file.exists(paste0(FP, "GBIFDownloads/GBIF", cutty, ".csv"))){
    return(NULL)
  }
  GBIFSpec <- fread(paste0(FP, "GBIFDownloads/GBIF", cutty, ".csv"), header=TRUE)
  GBIFSpecPts <- vect(GBIFSpec, geom=c("decimalLongitude", "decimalLatitude"), crs=WGSCRS)
  GBIFSpecPts <- project(GBIFSpecPts, MollCRS)
  GBIFSpec$Grid <- extract(TemplateRas, GBIFSpecPts, method="simple")$bio_1
  GBIFSpec <- subset(GBIFSpec, year<2018)
  GBIFSpec <- subset(GBIFSpec, year>1900)
  GBIFSpec$YearGroup <- ifelse(GBIFSpec$year < 1951, "1901_1950", ifelse(GBIFSpec$year < 2001, "1951_2000", "2001_2017"))
  
  #Now do some subsetting.
  GBIFSpec <- GBIFSpec[!GBIFSpec$basisOfRecord %in% c("FOSSIL_SPECIMEN", "PRESERVED_SPECIMEN"),]
  GBIFSpec <- subset(GBIFSpec, occurrenceStatus=="PRESENT")
  GBIFSpec <- subset(GBIFSpec, species == GivenSpec)
  
  #Now take one random point per species per grid cell per year group
  set.seed(364)
  GBIFSpec <- GBIFSpec[sample(1:nrow(GBIFSpec), nrow(GBIFSpec), replace = FALSE),]
  GBIFSpec$SpecGridYear <- paste0(GBIFSpec$GivenSpec," ", GBIFSpec$Grid, " ", GBIFSpec$YearGroup)
  GBIFSpec <- GBIFSpec[!duplicated(GBIFSpec$SpecGridYear,),]
  return(GBIFSpec)
}, mc.cores=ncores))
GBIFRedSite <- unique(GBIFRed[,c("decimalLatitude", "decimalLongitude")])
GBIFRedSite$Site <- 1:nrow(GBIFRedSite)
GBIFRed <- merge(GBIFRed, GBIFRedSite)

#Right, think we're good 
save(GBIFRed, file=paste0(FP, "GBIFDownloads/GBIFReduce.RData"))

#### Get current climate data, take 50yr means ####

#Create Bioclim variables from CRU data
AllRasYears <- function(climtype){
  climate_output <- nc_open(list.files(path=paste0(DataFP, "/Climate_CRU/ClimateData"),pattern=paste0("*", climtype, ".dat.nc"), full.names=TRUE))
  
  #Get data
  lon <- ncvar_get(climate_output,"lon")
  lat <- ncvar_get(climate_output,"lat",verbose=F)
  time <- ncvar_get(climate_output,"time")
  tunits <- ncatt_get(climate_output,"time","units")
  fillvalue <- ncatt_get(climate_output, climtype,"_FillValue")
  clim_array <- ncvar_get(climate_output,climtype)
  clim_array[clim_array==fillvalue$value] <- NA #Change NAs to appropriate thing
  dlname <- ncatt_get(climate_output,climtype,"long_name")
  dunits <- ncatt_get(climate_output, climtype,"units")
  
  #Get metadata
  title <- ncatt_get(climate_output,0,"title")
  institution <- ncatt_get(climate_output,0,"institution")
  datasource <- ncatt_get(climate_output,0,"source")
  references <- ncatt_get(climate_output,0,"references")
  history <- ncatt_get(climate_output,0,"history")
  Conventions <- ncatt_get(climate_output,0,"Conventions")
  nc_close(climate_output)
  
  #Get the right times
  tustr <- strsplit(tunits$value, " ")
  tdstr <- strsplit(unlist(tustr)[3], "-")
  tmonth <- as.integer(unlist(tdstr)[2])
  tday <- as.integer(unlist(tdstr)[3])
  tyear <- as.integer(unlist(tdstr)[1])
  dates <- chron(time,origin=c(tmonth, tday, tyear))
  dates <- rbindlist(lapply(1:length(dates), function(x) as.data.frame(t(strsplit(as.character(dates[x]), "/")[[1]][c(1,3)]))))
  names(dates) <- c("Month", "Year")
  dates$RasNum <- rownames(dates)
  
  #Fix up the dates (they only use the second two numbers of the date, meaning currently 1901 and 2001 are the same. Add the 19/20 back in)
  dates[1:(99*12),]$Year <- paste0("19", dates[1:(99*12),]$Year)
  dates[((99*12)+1):nrow(dates),]$Year <- paste0("20", dates[((99*12)+1):nrow(dates),]$Year)
  
  #Extract raster of designated month/year
  AllRas <- pblapply(unique(dates$Year), function(year){
    MonthRass <- lapply(unique(dates$Month), function(month){
      m <- as.numeric(as.character(dates[dates$Month==month & dates$Year==year,]$RasNum))
      clim_slice <- clim_array[,,m]
      #Ok so that gets alll our data out
      # matrix (nlon*nlon rows by 2 cols) of lons and lats
      lonlat <- as.matrix(expand.grid(lon,lat))
      # vector of values
      clim_vec <- as.vector(clim_slice)
      clim_df <- data.frame(cbind(lonlat,clim_vec))
      names(clim_df) <- c("lon","lat",paste(climtype,as.character(m), sep="_"))
      clim_ras <- rasterFromXYZ(clim_df)  #Convert first two columns as lon-lat and third as value
      return(clim_ras)
    })
    names(MonthRass) <- unique(dates$Month)
    return(MonthRass)
  })
  names(AllRas) <- unique(dates$Year)
  return(AllRas)
}

TempMin <- AllRasYears("tmn")
TempMax <- AllRasYears("tmx")
Precip <- AllRasYears("pre")

AllBioClims <- pblapply(names(TempMin), function(Year){
  BioVarss <- biovars(stack(Precip[[Year]]), stack(TempMin[[Year]]), stack(TempMax[[Year]]))
  writeRaster(BioVarss,paste0(DataFP, "/Climate_CRU/BioClim/", Year, "BioClim.grd"), format="raster")
})  

#Now get means for 50 year time periods
TimeStamps <- list(seq(1901,1950,1), seq(1951,2000,1), seq(2001,2017,1))

for(time in TimeStamps){
  BioRasters <- pblapply(time, function(Year){
    stack(paste0(DataFP, "/Climate_CRU/BioClim/", Year, "BioClim.grd"))
  })
  
  MeanBioRasters <- pblapply(c(1:19), function(x){
    BioMean <- stack(lapply(BioRasters, function(y) y[[x]]))
    BioMean <- calc(BioMean, fun = mean, na.rm = T)
    return(BioMean)
  })
  names(MeanBioRasters) <- paste0("bio_", c(1:19))
  MeanBioRasters <- stack(MeanBioRasters)
  writeRaster(MeanBioRasters,paste0(DataFP, "/Climate_CRU/BioClim/", time[1], "_", time[length(time)], "BioClim.grd"), format="raster")
}

#### Check how well distributions are represented in geographic space (sampling of polygons) and environmental space (exdet of climate at occurrences vs climate across polygon) ####
#From this paper: "Here be dragons: a tool for quantifying novelty due to covariate range and correlation change when projecting species distribution models"

DistributionSHPs <- lapply(list.files(paste0(FP, "POWODistributions/"), pattern="*.shp", full.names = TRUE), vect)
names(DistributionSHPs) <- sapply(DistributionSHPs, function(x) unique(x$Species))
ETDs <- list.files(paste0(DataFP, "EuropeanTreeDistributions/TreePolygons/"), pattern="*.shp", full.names=TRUE)

#Get modern temps rasters
ModernTempsAll <- lapply(c("1901_1950", "1951_2000", "2001_2017"), function(yearry){
  ModernTemps <- rast(paste0(DataFP, "/Climate_CRU/BioClim/", yearry, "BioClim.grd"))
  crs(ModernTemps) <- WGSCRS
  ModernTemps <- project(ModernTemps, MollCRS)
  ModernTemps <- squareRastCells(ModernTemps)
  return(ModernTemps)
})
names(ModernTempsAll) <- c("1901_1950", "1951_2000", "2001_2017")

#Make a template
TemplateRas <- ModernTempsAll[[2]][[1]]
TemplateRas[!is.na(TemplateRas)][1:length(TemplateRas[!is.na(TemplateRas)])] <- 1:length(TemplateRas[!is.na(TemplateRas)])

load(file=paste0(FP, "GBIFDownloads/GBIFReduce.RData"))

Spec <- "Abies alba"
CoverageCheck <- function(EnvStack, SpecPolygon, SpecPoints){
  if(class(EnvStack)!="list"){
    EnvStack <- list(EnvStack)
    names(EnvStack) <- "Time"
    SpecPoints$YearGroup <- "Time"
  }
  
  #First create a template raster with unique gridcell IDs
  TemplateRas <- EnvStack[[1]][[1]]
  TemplateRas[!is.na(TemplateRas)][1:length(TemplateRas[!is.na(TemplateRas)])] <- 1:length(TemplateRas[!is.na(TemplateRas)])
  
  #Turn the species polygon into a raster
  SpecPolygon <- project(SpecPolygon, MollCRS)
  SpecRast <- crop(TemplateRas, SpecPolygon)
  SpecRastTemplate <- terra::rasterize(SpecPolygon, SpecRast, getCover=TRUE)
  SpecRastTemplate[SpecRastTemplate==0] <- NA
  SpecRast <- mask(SpecRast, SpecRastTemplate)
  
  #Turn points into a points object, reduce to only points within the occurrence polygon, get the unique grid ID
  SpecPoints <- project(SpecPoints, MollCRS)
  SpecPoints <- SpecPoints[SpecPolygon,]
  SpecPoints$Grid <- extract(TemplateRas, SpecPoints, method="simple")[,2]
  
  #Kill any with less than x points
  if(nrow(SpecPoints) < 20){
    return(data.frame(PropDistribCovered = "TooFewPoints",
                      PropEnvironSampled = "TooFewPoints"))
  }
  
  #Check that there's not too many points outside the template raster
  if(length(SpecPoints$Grid[is.na(SpecPoints$Grid)])/nrow(SpecPoints) > 0.15){stop("there are more than 15% of points falling outside the raster")}
  
  #Now PCA the environmental rasters to check coverage
  EnvStackCropped <- rbindlist(lapply(1:length(EnvStack), function(Rast){
    TempRast <- crop(EnvStack[[Rast]], SpecPolygon)
    TempRastTemplate <- rasterize(SpecPolygon, TempRast, getCover=TRUE)
    TempRastTemplate[TempRastTemplate==0] <- NA
    TempRast <- mask(TempRast, TempRastTemplate)
    TempRastDF <- as.data.frame(TempRast)
    TempRastDF <- cbind(TempRastDF, crds(TempRast))
    TempRastDF <- TempRastDF[complete.cases(TempRastDF[,1:19]),]
    TempRastDF$Time <- names(EnvStack)[[Rast]]
    return(TempRastDF)
  }))
  EnvStackCrdsTime <- EnvStackCropped[,(nlyr(EnvStack[[1]])+1):ncol(EnvStackCropped)]
  EnvStackCropped <- do.call(cbind, lapply(EnvStackCropped[,1:nlyr(EnvStack[[1]])], function(x){
    if(length(unique(x))>1){
      return(x)
    } else {
      return(NULL)
    }
  }))
  pcaEnv <- prcomp(EnvStackCropped, center = T, scale. = T)
  if(summary(pcaEnv)$importance[3,3] < 0.70){stop("PCA is capturing less than 70% of variance")}
  EnvPCA <- as.data.frame(cbind(EnvStackCrdsTime, pcaEnv$x[,c(1:3)]))
  
  EnvPCARas <- pblapply(unique(EnvPCA$Time), function(x){
    PCATime <- EnvPCA[EnvPCA$Time == x,]
    PCATime[,c("Time")] <- NULL
    PCARas <- rast(PCATime, crs=MollCRS, type="xyz")
    return(PCARas)
  })
  names(EnvPCARas) <- unique(EnvPCA$Time)
  
  EnvPCAExtract <- rbindlist(pblapply(names(EnvPCARas), function(x){
    TempRas <- EnvPCARas[[x]]
    SpecTime <- SpecPoints[SpecPoints$YearGroup==x,]
    SpecTimeClim <- as.data.frame(terra::extract(TempRas, SpecTime, method="simple", ID=FALSE))
    SpecTime <- as.data.frame(cbind(SpecTime, SpecTimeClim))
    return(SpecTime)
  }))
  EnvPCAExtract <- EnvPCAExtract[complete.cases(EnvPCAExtract[,c("PC1", "PC2", "PC3")]),]
  
  #Now assess ExDet. So this looks at all the variables in the first entry (in our case, the conditions at the species points)
  #Then it looks at the conditions over the entire range (EnvPCA). It looks to see which cells in EnvPCA are *outside* the climatic space represented by the poitns
  #Values below 0 and above 1 indicate non-analog climates (i.e. we have not covered that climate space with the species points)
  
  EnvPCA$ExDet <- ecospat.climan(EnvPCAExtract[,c("PC1", "PC2", "PC3")], EnvPCA[,c("PC1", "PC2", "PC3")])
  EnvPCA$ExDetCat <- ifelse(EnvPCA$ExDet < 0 | EnvPCA$ExDet > 1, "NonAnalog", "Analog")
  
  #Get proportions, bring it home baby
  
  return(data.frame(PropDistribCovered = length(unique(na.omit(SpecPoints$Grid)))/nrow(na.omit(unique(SpecRast))),
                    PropEnvironSampled = nrow(subset(EnvPCA, ExDetCat == "Analog"))/nrow(EnvPCA)))
}

#Species with <10 points get removed
AssessSpatialClimateCoverage <- pblapply(unique(GBIFRed$species), function(Spec){
  print(Spec)
  if(file.exists(paste0(FP, "Coverages/", Spec, ".csv"))){return(NULL)}
  SpecGBIF <- subset(GBIFRed, species==Spec)
  if(min(SpecGBIF$year)<1901 | max(SpecGBIF$year)>2017){stop(paste0("there's years there should be in ", Spec))}
  if(nrow(SpecGBIF)<10){
    return(NULL)
  }
  SpecPoints <- SpecGBIF
  SpecPoints <- vect(SpecPoints, geom=c("decimalLongitude", "decimalLatitude"), crs=WGSCRS)
  SpecPoints <- project(SpecPoints, MollCRS)
  
  #BEGIN GET HULL
  SpecPointsHull <- SpecPoints
  SpecHull <- convHull(SpecPointsHull)
  SpecHullArea <- expanse(SpecHull)
  SpecPointsHull$HullContrib <- as.numeric(pbmclapply(1:nrow(SpecPointsHull), function(x){
    SpecHullCont <- convHull(SpecPointsHull[-x,])
    return(expanse(SpecHullCont)/SpecHullArea)
  }, mc.cores=ncores))
  while(min(SpecPointsHull$HullContrib) < 0.90){
    SpecPointsHull <- SpecPointsHull[SpecPointsHull$HullContrib > 0.90]
    SpecHull <- convHull(SpecPointsHull)
    SpecHullArea <- expanse(SpecHull)
    SpecPointsHull$HullContrib <- as.numeric(pbmclapply(1:nrow(SpecPointsHull), function(x){
      SpecHullCont <- convHull(SpecPointsHull[-x,])
      return(expanse(SpecHullCont)/SpecHullArea)
    }, mc.cores=ncores))
    if(nrow(SpecPointsHull)<10){break}
    #Let's remove points that contribute more than 10% to the area
  }
  if(nrow(SpecPointsHull)<10){
    return(NULL)
  }
  #END GET HULL
  
  #GET POWO DIST
  POWODist <- DistributionSHPs[[Spec]]
  
  #DOES A EUROPEAN TREE DISTRIBUTION EXIST FOR THE SPECIES?
  if(length(ETDs[grepl(Spec, ETDs)])>0){
    ETD <- vect(ETDs[grepl(Spec, ETDs)])
    DistOptions <- list(SpecHull, POWODist, ETD)
    names(DistOptions) <- c("MCH", "POWO", "EuropeTrees")
  } else {
    DistOptions <- list(SpecHull, POWODist)
    names(DistOptions) <- c("MCH", "POWO")
  }
  
  #POWO Poly
  Covers <- rbindlist(lapply(DistOptions, function(x) CoverageCheck(ModernTempsAll, x, SpecPoints)))
  Covers$Type <- names(DistOptions)
  Covers$Species <- Spec
  write.csv(Covers, file=paste0(FP, "Coverages/", Spec, ".csv"), row.names=FALSE)
  return(Covers)
})

Coverages <- rbindlist(lapply(list.files(path=paste0(FP, "Coverages/"), full.names = TRUE), read.csv))
Coverages$X <- NULL

CoveragesGone <- subset(Coverages, PropEnvironSampled=="TooFewPoints")
Coverages <- subset(Coverages, PropEnvironSampled!="TooFewPoints")
Coverages$PropDistribCovered <- as.numeric(Coverages$PropDistribCovered)
Coverages$PropEnvironSampled <- as.numeric(Coverages$PropEnvironSampled)
Coverages <- subset(Coverages, PropDistribCovered<1)
ggplot(Coverages, aes(x=PropEnvironSampled, y=PropDistribCovered))+
  geom_point()+
  facet_wrap(~Type)+
  xlab("Proportion Climatic Environment Covered")+
  ylab("Proportion Distribution Covered")+
  theme_classic()+
  theme(aspect.rat=1)

#Okay, neverMIND about distribution coverage, it's a bust. That's okay, we can still assess environmental sampling!

ggplot(Coverages, aes(x=Type, y=PropEnvironSampled*100))+
  geom_boxplot()+
  xlab("Range Type")+
  ylab("% Climate Represented")+
  theme_classic()

ggplot(Coverages, aes(x=Type, y=PropDistribCovered*100))+
  geom_boxplot()+
  xlab("Range Type")+
  ylab("% Range Covered")+
  theme_classic()

CoveragesET <- Coverages[Coverages$Species %in% subset(Coverages, Type=="EuropeTrees")$Species,]

library(betareg)

CoveragesET <- dcast(CoveragesET, Species ~ Type, value.var="PropEnvironSampled")
summary(betareg(MCH ~ EuropeTrees, data=CoveragesET))
summary(betareg(POWO ~ EuropeTrees, data=CoveragesET))

CoveragesETPlot <- melt(CoveragesET, id.vars=c("Species", "EuropeTrees"), variable.name="RangeType")

ggplot(CoveragesETPlot, aes(x=EuropeTrees, y=value, fill=RangeType, colour=RangeType))+
  geom_point()+
  geom_smooth(method="lm")+
  xlab("% Climate Represented (European Tree Range)")+
  ylab("% Climate Represented (MCH/POWO)")+
  theme_classic()

#Right, let's go for 70%. Okay. (We can sensitivity test this later? Or just run it for everyone but we have this data. YES that one. )

#### PCA climate vars WOOF THIS IS GONNA BE INTENSE #### 

#OKAY. We have a problem here. The size of the PCA is just too mammoth to make any sense. I'm not sure what this means for what we do.

#PCA on present day occurrences and past rasters
#Right, so we're now doing this for the whole world so, well, fuck. 

#This says use PCAs: Selecting predictors to maximize thetransferability of species distributionmodels: lessons from cross-continentalplant invasions
#This investigates PCA-env vs. PCA-occ (we used PCA env, easier to compare across time steps, and better for multi species): Measuring ecological niche overlap from occurrence and spatial environmental data
#NOPE we've now decided to use PCA-occ of all species for current data, as it better represents sampling bias
#In this one they do the full PCA on both present and projected rasters: Reconstructing the climatic niche breadth of land use for animal production during the African Holocene

BBoxRas <- rast(xmin = -180, xmax= 180, ymin = -90, ymax = 90,nrows=1,ncols=1)
BBoxMoll <- project(BBoxRas, MollCRS)
TimeStamps <- list(seq(1901,1950,1), seq(1951,2000,1), seq(2001,2017,1))
YearGroups <- c("1901_1950", "1951_2000", "2001_2017")

CurrentClimEnv <- lapply(YearGroups, function(YG){
  ModernTemps <- rast(paste0(DataFP, "/Climate_CRU/BioClim/", YG, "BioClim.grd"))
  crs(ModernTemps) <- WGSCRS
  ModernTemps <- project(ModernTemps, MollCRS)
  ModernTemps <- squareRastCells(ModernTemps)
  return(ModernTemps)
})
names(CurrentClimEnv) <- YearGroups

CurrentClimEnvDF <- rbindlist(lapply(1:length(CurrentClimEnv), function(x){
  CurrentClimDF <- as.data.frame(CurrentClimEnv[[x]])
  CurrentClimDF <- cbind(CurrentClimDF, crds(CurrentClimEnv[[x]]))
  CurrentClimDF <- CurrentClimDF[complete.cases(CurrentClimDF[,1:19]),]
  CurrentClimDF$Time <- YearGroups[[x]]
  CurrentClimDF$OccEnv <- "Env"
  return(CurrentClimDF)
}))

#Get past rasters, convert to DF
YearNames <- sapply(list.files(PastFP, pattern = "*.grd"), function(x) gsub("BP", "", str_split_fixed(str_split_fixed(x, "[_]", 3)[,3], "[.]", 2)[,1]))
YearNames <- YearNames[order(as.numeric(YearNames), decreasing = TRUE)]

#Fix downscaling file issue (UGH)
PastFix <- pbmclapply(YearNames, function(YR){
  print(YR)
  if(length(list.files(PastFP, pattern = paste0("Downscaled_BioVars_", YR, "BP.grd"), full.names=TRUE))>1){stop("double raster years!")}
  if(length(list.files(paste0(PastFP,"/", YR, "/"), pattern="*.tif")) >= 19){
    return(NULL)
  }
  BIO_ras <- stack(list.files(PastFP, pattern = paste0("Downscaled_BioVars_", YR, "BP.grd"), full.names=TRUE))
  dir.create(paste0(PastFP,"/", YR), showWarnings = FALSE)
  lapply(c(1:19), function(x) writeRaster(BIO_ras[[x]], paste0(PastFP,"/", YR, "/bio_", x, ".tif")))
  return("Done")
}, mc.cores=ncores-2) #

#Okay, get prep for PCAs (cross products, colmeans, nrows)
#FIRST we're gonna take 100 year means (WOOF)

HundredYears <- gsub("00", "", YearNames[grepl("*00$", YearNames)])
HundredYears <- HundredYears[HundredYears!="210"]
PastFPHundred <- paste0(DataFP, "PalaeoView/Downscale_21000to0BP_20yrInterval_20yrMean/HundredYearAvgsBioclim/")

PastClimHundredYearMeans <- pblapply(HundredYears, function(YR){
  print(YR)
  if(length(list.files(paste0(PastFPHundred, YR, "00")))==19){
    return(NULL)
  }
  mclapply(c(1:19), function(bio){
    PastRast <- rast(paste0(PastFP, "/", YR, c("00", "20", "40", "60", "80"), "/bio_", bio, ".tif"))
    crs(PastRast) <- WGSCRS
    PastRast <- project(PastRast, MollCRS)
    PastRast <- squareRastCells(PastRast)
    PastRast <- mean(PastRast)
    dir.create(paste0(PastFPHundred, YR, "00"), showWarnings = FALSE)
    writeRaster(PastRast, paste0(PastFPHundred, YR, "00/bio_", bio, ".tif"), overwrite=TRUE)
    rm(PastRast)
    return("okay")
  }, mc.cores=ncores)
}) #

#Look and we're just gonna run a separate one for the first 5 chunks. Not elegant but it'll do
pblapply(c(1:19), function(bio){
  PastRast <- rast(paste0(PastFP, "/", c("0", "20", "40", "60", "80"), "/bio_", bio, ".tif"))
  crs(PastRast) <- WGSCRS
  PastRast <- project(PastRast, MollCRS)
  PastRast <- squareRastCells(PastRast)
  PastRast <- mean(PastRast)
  dir.create(paste0(PastFPHundred, "0"), showWarnings = FALSE)
  writeRaster(PastRast, paste0(PastFPHundred, "0/bio_", bio, ".tif"), overwrite=TRUE)
})

HundredYears <- c(paste0(HundredYears, "00"), 0)

#We're gonna crop Antarctica off the rasters
WorldOutline <- vect(paste0(DataFP, "GISLayers/CountryContinentWorldSHPS/World/TM_WORLD_BORDERS-0.3-SSudan.shp"))
WorldOutline <- project(WorldOutline, WGSCRS)
WorldOutline <- project(WorldOutline, MollCRS)
WorldNoAntarctica <- WorldOutline[WorldOutline$REGION!=0,]

#PCA maths is in the notebook and PCA fix code! 

PastClimPCAPrep <- pbmclapply(HundredYears, function(YR){
  PastRast <- rast(list.files(paste0(PastFPHundred, "/", YR), pattern="*.tif", full.names=TRUE))
  names(PastRast) <- gsub(".tif", "", list.files(paste0(PastFPHundred, "/", YR), pattern="*.tif"))
  PastRast <- crop(PastRast, WorldNoAntarctica)
  PastRast <- mask(PastRast, WorldNoAntarctica)
  PastRast <- PastRast[[names(CurrentClimEnv[[1]])]]
  PastRastMat <- as.matrix(as.data.frame(PastRast))
  PastRastMat <- PastRastMat[, !colnames(PastRastMat) %in% c("bio_2", "bio_3")]
  
  CP <- crossprod(sweep(PastRastMat, 2, colMeans(PastRastMat)))
  CM <- colMeans(PastRastMat)
  NR <- nrow(PastRastMat)

  return(list(CP, CM, NR))
}, mc.cores=ncores)

CurrentClimPCAPrep <- pblapply(CurrentClimEnv, function(x){
  RastMat <- as.matrix(as.data.frame(x))
  RastMat <- RastMat[, !colnames(RastMat) %in% c("bio_2", "bio_3")]
  
  CP <- crossprod(sweep(RastMat, 2, colMeans(RastMat)))
  CM <- colMeans(RastMat)
  NR <- nrow(RastMat)
  
  return(list(CP, CM, NR))
})

#Run PCA
TotalCrosProd <- Reduce("+", c(lapply(PastClimPCAPrep, function(x) x[[1]]), lapply(CurrentClimPCAPrep, function(x) x[[1]])))
TotalNRow <- Reduce("+", c(lapply(PastClimPCAPrep, function(x) x[[3]]), lapply(CurrentClimPCAPrep, function(x) x[[3]])))
ClimPCA <- princomp(covmat=TotalCrosProd/TotalNRow)

#Extract average colmeans for reconstructing rasters
AverageColMeans <- Reduce("+", c(lapply(PastClimPCAPrep, function(x) x[[2]]), lapply(CurrentClimPCAPrep, function(x) x[[2]])))/
  (length(PastClimPCAPrep) + length(CurrentClimPCAPrep))

#Reconstruct rasters
CurrentClimPCA <- pblapply(1:length(CurrentClimEnv), function(x){
  RastDF <- as.data.frame(CurrentClimEnv[[x]])
  RastDF <- RastDF[, !colnames(RastDF) %in% c("bio_2", "bio_3")]
  RastDFScale <- sweep(RastDF, 2, AverageColMeans)
  RastDF$PCA1 <- rowSums(sweep(RastDFScale, 2, ClimPCA$loadings[,1], "*"))
  RastDF$PCA2 <- rowSums(sweep(RastDFScale, 2, ClimPCA$loadings[,2], "*"))
  RastDF$PCA3 <- rowSums(sweep(RastDFScale, 2, ClimPCA$loadings[,3], "*"))
  RastDF <- cbind(RastDF, crds(CurrentClimEnv[[x]]))
  Rasty <- rast(RastDF[,c("x", "y", "PCA1", "PCA2", "PCA3")], crs=crs(CurrentClimEnv[[x]]))
  dir.create(paste0(FP, "PCAClims/Current/", names(CurrentClimEnv)[[x]]), showWarnings = FALSE)
  terra::writeRaster(Rasty, paste0(FP, "PCAClims/Current/", names(CurrentClimEnv)[[x]], "/", names(Rasty), ".tif"), overwrite=TRUE)
  return(Rasty)
})

PastClimPCA <- pbmclapply(HundredYears, function(YR){
  if(length(list.files(paste0(FP, "PCAClims/Past/", YR, "/")))==3){
    return(NULL)
  }
  PastRast <- rast(list.files(paste0(PastFPHundred, "/", YR), pattern="*.tif", full.names=TRUE))
  names(PastRast) <- gsub(".tif", "", list.files(paste0(PastFPHundred, "/", YR), pattern="*.tif"))
  PastRast <- crop(PastRast, WorldNoAntarctica)
  PastRast <- mask(PastRast, WorldNoAntarctica)
  PastRast <- PastRast[[names(CurrentClimEnv[[1]])]]
  PastRast <- PastRast[[!names(PastRast) %in% c("bio_2", "bio_3")]]
  
  RastDF <- as.data.frame(PastRast)
  RastDFScale <- sweep(RastDF, 2, AverageColMeans)
  RastDF$PCA1 <- rowSums(sweep(RastDFScale, 2, ClimPCA$loadings[,1], "*"))
  RastDF$PCA2 <- rowSums(sweep(RastDFScale, 2, ClimPCA$loadings[,2], "*"))
  RastDF$PCA3 <- rowSums(sweep(RastDFScale, 2, ClimPCA$loadings[,3], "*"))
  RastDF <- cbind(RastDF, crds(PastRast))
  Rasty <- rast(RastDF[,c("x", "y", "PCA1", "PCA2", "PCA3")], crs=crs(PastRast))
  dir.create(paste0(FP, "PCAClims/Past/", YR), showWarnings = FALSE)
  terra::writeRaster(Rasty, paste0(FP, "PCAClims/Past/", YR, "/", names(Rasty), ".tif"), overwrite=TRUE)
  rm(Rasty)
  return("done")
}, mc.cores=ncores)

#### Attach phylogenetic backbone ####
library(devtools)
install_github("jinyizju/V.PhyloMaker")
library(V.PhyloMaker)

AllSpec <- read.csv(paste0(FP, "PollenData/FosPolCleaned_Harmonised_100Samples.csv"))
dat.phylo <- unique(as.data.frame(cbind(AllSpec$Species, str_split_fixed(AllSpec$Species, "[ ]", 2)[,1], AllSpec$Family)))
colnames(dat.phylo) <- c("species", "genus", "family")

### Make phylo tree
tree1 <- phylo.maker(sp.list=dat.phylo, tree=GBOTB.extended, nodes=nodes.info.1, scenarios="S3")
## Tree information: https://rdrr.io/github/jinyizju/V.PhyloMaker/man/GBOTB.extended.html
unique(tree1$species.list$status) ## all taxa are pruned or bound, none 'failed to bind'.
tree1$scenario.3$node.label ## Most have just "". This causes problems later on.(?)
branching.times(tree1$scenario.3) ## seems fine

write.tree(tree1$scenario.3, "tree1.tre")

### Explore data
plot.phylo(tree1$scenario.3, cex = 0.5, main = "scenario.3")

t2 <- extract.clade(tree1$scenario.3, 504)
plot(t2, edge.width = 2, label.offset = 0.1, type = "cladogram")
nodelabels(t2$node.label)

### Inspect the backbone if needed (e.g. because some taxa don't bind)
bkb <- read.tree(file = "E:/NON_PROJECT/STATISTICS/PHYLOGENY/SEED_PLANT_PHYLOGENIES/v0.1/v0.1/GBOTB.tre")
write.csv(bkb$tip.label, "C:/Users/re259/Downloads/v0.1/v0.1/GBOTB_spnames.csv")

#### Check ExDet (updated MESS) distribution present vs past (This is OLD and hasn't been updated) ####
#From this paper: "Here be dragons: a tool for quantifying novelty due to covariate range and correlation change when projecting species distribution models"

load(file=paste0(FP, "ClimateDat/CurrPastDF.RData"))
load(file=paste0(FP, "ClimateDat/CurrPastStacks.RData"))
load(file=paste0(FP, "ClimateDat/CurrentOcc.RData"))

CurrPastEnv$ExDet <- ecospat.climan(CurrentOcc[,c("PC1", "PC2", "PC3")], CurrPastEnv[,c("PC1", "PC2", "PC3")])
CurrPastEnv$ExDetCat <- ifelse(CurrPastEnv$ExDet > 1 | CurrPastEnv$ExDet < 0, "NoAnalog", 'Analog')

NovelEnvTime <- dcast(as.data.table(CurrPastEnv), Time ~ ExDetCat, length)
NovelEnvTime$PropNoAnalog <- NovelEnvTime$NoAnalog / rowSums(NovelEnvTime[,c(2:3)])
NovelEnvTime <- NovelEnvTime[order(NovelEnvTime$PropNoAnalog, decreasing=TRUE),]

ExDetRasters <- pblapply(unique(CurrPastEnv$Time), function(x){
  ExDetTime <- CurrPastEnv[CurrPastEnv$Time == x, c(1,2, ncol(CurrPastEnv))]
  ExDetTime$ExDetCat <- ifelse(ExDetTime$ExDetCat == "Analog", 0, 1)
  ExDetRes <- if(x %in% c("1901_1950", "1951_2000", "2001_2017")){res(CurrPastStacks[[204]])} else {res(CurrPastStacks[[1]])}
  ExDetRas <- rasterFromXYZ(ExDetTime, crs=MollCRS, res=ExDetRes)
  return(ExDetRas)
})

save(ExDetRasters, file=paste0(FP, "ClimateDat/ExDetRasters.RData"))

#### Hmmm cheeky pollen final update ####
FosPol <- read.csv(paste0(FP, "PollenData/FosPolCleaned_Harmonised_100Samples.csv"))
FosPol$X <- NULL
FosPol <- FosPol[FosPol$Count!=0,]
AlltheAges <- pbmclapply(1:nrow(FosPol), function(x) seq(round_any(FosPol[x,]$x12_5_percent, 100), round_any(FosPol[x,]$x87_5_percent, 100), 100), mc.cores=ncores) #unique(unlist(
FosPol$HundredYrAgeSpan <- sapply(AlltheAges, length)
FosPol <- subset(FosPol, HundredYrAgeSpan <= 10)

SamplePerSpec <- dcast(unique(FosPol[,c("Species", "Dataset_Sample")]), Species ~ ., length)
SamplePerSpec <- subset(SamplePerSpec, . >= 70)

FosPol <- FosPol[FosPol$Species %in% c(SamplePerSpec$Species, "Silene acaulis"),]
write.csv(FosPol, paste0(FP, "PollenData/FosPolCleaned_Harmonised_100Samples_10TimeSteps.csv"), row.names=FALSE)
