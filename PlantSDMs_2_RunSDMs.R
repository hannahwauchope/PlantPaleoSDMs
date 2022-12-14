#renv::deactivate() 

### European Tree SDMS ####
cluster <- TRUE
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
  library(terra)
  library(taxize)
  #devtools::install_github("barnabywalker/kewr")
  library(kewr)
  library(plyr)
  library(stringr)
  #devtools::install_github("biomodhub/biomod2", dependencies = TRUE)
  
  FP <- "/Users/hannahwauchope/Dropbox/Work/Projects/PaleoSDMs/"
  DataFP <- "/Users/hannahwauchope/Dropbox/Work/Data/"
  BaseFP <- "/Users/hannahwauchope/Dropbox/Work/Projects/"
  MaxentPath <- "/Users/hannahwauchope/Dropbox/Work/Software/maxent/"
  Computer <- "MacBookMax"
} else {
  #library(sp)
  library(biomod2)
  library(raster)
  library(terra)
  library(plyr)
  library(pbmcapply, lib.loc="/nobackup/beegfs/home/ISAD/hw656/R/x86_64-pc-linux-gnu-library/3.5/sp/libs")
  library(stringr)
  library(data.table, lib.loc="/nobackup/beegfs/home/ISAD/hw656/R/x86_64-pc-linux-gnu-library/3.5/sp/libs")
  library(pbapply, lib.loc="/nobackup/beegfs/home/ISAD/hw656/R/x86_64-pc-linux-gnu-library/3.5/sp/libs")
  FP <- "/nobackup/beegfs/workspace/hw656/Projects/PaleoSDMs/"
  DataFP <- "/nobackup/beegfs/workspace/hw656/Data/"
  MaxentPath <- "/nobackup/beegfs/workspace/hw656/Data/Software/maxent/"

  Computer <- "Cluster"
}

WGSCRS <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
MollCRS <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
RobCRS <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

ncores <- ifelse(Computer=="MacMini", 6, ifelse(Computer=="MacBookMax", 7, ifelse(Computer=="MacBookPro", 6, 32)))
PastFP <- paste0(DataFP, "PalaeoView/Downscale_21000to0BP_20yrInterval_20yrMean/BioClim")

#### Functions ####
source("~/Documents/GitHub/UsefulFunctions/GetLoadGBIF.r")
source("~/Documents/GitHub/UsefulFunctions/SquareRasterCells.R")
source("~/Documents/GitHub/UsefulFunctions/CRUClimateData.R")

#### Build SDM ####

#On ensemble modelling:"A review of evidence about use and performance of species distribution modelling ensembles like BIOMOD"
#https://onlinelibrary.wiley.com/doi/10.1111/ddi.12892

#On evaluating ensemble models: https://rstudio-pubs-static.s3.amazonaws.com/38564_747d4bbf87704f0394734977bd4905c4.html
#"The best way to address this problem of not fairly evaluated ensemble-models compared with formal-models is to specify an independent data-set to evaluate all models with at 
#BIOMOD_FormatingData step (i.e. filling eval.resp.var , eval.expl.var, eval.resp.xy arguments). 
#Then all models you will produce will be evaluated against this evaluation data.

#Let's set up an independent chunk of data to evaluate everything on (OKAY)

#Read in occurrence data and present day raster stacks
load(file=paste0(FP, "GBIFDownloads/GBIFReduce.RData"))
FosPol <- read.csv(paste0(FP, "PollenData/FosPolCleaned_Harmonised_100Samples_10TimeSteps.csv"))

YearGroups <- c("1901_1950", "1951_2000", "2001_2017")
CurrentOcc <- GBIFRed
rm(GBIFRed)
CurrentOcc <- CurrentOcc[CurrentOcc$GivenSpec %in% FosPol$Species,]

CurrStacks <- lapply(YearGroups, function(YR) rast(list.files(paste0(FP, "PCAClims/Current/", YR), full.names=TRUE)))
names(CurrStacks) <- YearGroups

CurrentOcc$GivenSpec <- gsub(" ", ".", CurrentOcc$GivenSpec)
CurrentOcc$GivenSpec <- gsub("-", "", CurrentOcc$GivenSpec)

PointsPerSpec <- dcast(CurrentOcc, GivenSpec + family~., length)
PointsPerSpec <- PointsPerSpec[order(PointsPerSpec$.),]
PointsPerSpec <- subset(PointsPerSpec, . >= 100)
setwd(paste0(FP, "SDMs/"))

SpecToModel <- sort(unique(PointsPerSpec$GivenSpec))
if(Computer!="Cluster"){
  SpecToModel <- SpecToModel[1:20]
}

#Biomod needs you to do this (then it'll parallelise some bits)
library(doParallel)
registerDoParallel(cores=ncores)

StartTime <- Sys.time()
print(StartTime)

#Right now I've just run this on these 6 species (that were somewhat randomly chosen)
SpecToModel <- c("Populus.tremuloides", "Betula.pendula", "Juniperus.communis", "Pinus.contorta", "Alnus.rubra", "Centaurea.scabiosa")
PointsPerSpec <- PointsPerSpec[sample(1:nrow(PointsPerSpec), nrow(PointsPerSpec))]
Spec <- "Quercus.coccifera"
pbmclapply(PointsPerSpec$GivenSpec, function(Spec){
  dir.create("NoFit", showWarnings = FALSE)
  dir.create(paste0("Log/"), showWarnings=FALSE)
  
  print(Spec)
  
  if(file.exists(paste0("Log/", Spec, "_Begin.csv"))){
    return(NULL)
  }
  write.csv(NULL, paste0("Log/", Spec, "_Begin.csv"))
  
  # If the individual models haven't been run, run them
  if(!file.exists(paste0(FP, "/SDMs/", Spec, "/BioModModels.RData"))){
    #Format presence data and psuedoabsences, split into test and train (psuedoabsences defined as occurrences of other species *not* in the same gridcell)
    Occ <- subset(CurrentOcc, GivenSpec == Spec)
    PseudoA <- subset(CurrentOcc, GivenSpec != Spec)
    PseudoA <- PseudoA[!PseudoA$Grid %in% Occ$Grid,]
    Occ$PA <- 1
    PseudoA$PA <- NA
    
    #We only want to draw pseudoabscences from rouuugghhhlylyyy the same distribution. Let's draw a minimum convex hull, and buffer it a bit
    OccPoints <- vect(Occ, geom=c("decimalLongitude", "decimalLatitude"), crs = WGSCRS)
    OccPoints <- project(OccPoints, MollCRS)
    OccCH <- convHull(OccPoints)
    OccCHBuff <- buffer(OccCH, 250000)
    OccCH <- aggregate(union(OccCH, OccCHBuff))

    #Extract psuedoabsence points from the min convex hull + buffer
    PseudoAPoints <- vect(PseudoA, geom=c("decimalLongitude", "decimalLatitude"), crs = WGSCRS)
    PseudoAPoints <- project(PseudoAPoints, MollCRS)
    PseudoAPoints$CH <- terra::extract(OccCH, PseudoAPoints)[,2]
    PseudoA <- PseudoA[PseudoAPoints$CH==1,]
    
    #Cut down to one pseudoabsence per gridcell per year
    PseudoA$GridYear <- paste0(PseudoA$Grid, " ", PseudoA$YearGroup)
    PseudoA <- PseudoA[!duplicated(PseudoA$GridYear),]
    
    #Bring together occurrences and pseudo absences, extract climate variables at each point in each year
    SpecDat <- rbind(as.data.frame(Occ), as.data.frame(PseudoA[,-c("GridYear")]))
    SpecDat <- rbindlist(pblapply(unique(SpecDat$YearGroup), function(YR){
      CurrClim <- CurrStacks[[YR]]
      SpecYr <- subset(SpecDat, YearGroup==YR)
      SpecYr <- vect(SpecYr, geom=c("decimalLongitude", "decimalLatitude"), crs=WGSCRS)
      SpecYr <- project(SpecYr, MollCRS)
      SpecYr <- cbind(SpecYr, extract(CurrClim, SpecYr, ID=FALSE))
      return(cbind(as.data.frame(SpecYr), crds(SpecYr)))
    }))
    SpecDat <- SpecDat[complete.cases(SpecDat[,c(paste0("PCA", 1:3))]),]
    
    #Format data for BIOMOD
    #Let's pull out a fifth of occurrence records for independent evaluation (biomod doesn't let you run this - weird)
    # set.seed(123)
    # RanSam <- sample(1:nrow(SpecDat), round(nrow(SpecDat)/5), replace = FALSE)
    # SpecDatEval <- SpecDat[RanSam,]
    # SpecDat <- SpecDat[!RanSam,]
    # 
    # 
    # ExplVarEval <- as.matrix(SpecDatEval[,c(paste0("PCA", c(1:3)))])
    
    #Create explanatory var matrix of climate vars for each point  
    ExplVar <- as.matrix(SpecDat[,c(paste0("PCA", c(1:3)))])
    
    #Format input data (Eval. doesn't work with only pseudoabsences)
    InputBioDat <- BIOMOD_FormatingData(resp.var=SpecDat$PA, 
                                        resp.xy=cbind(as.numeric(SpecDat$x), as.numeric(SpecDat$y)),
                                        expl.var=ExplVar,
                                        resp.name=Spec,
                                        # eval.resp.var = SpecDatEval$PA,
                                        # eval.expl.var = ExplVarEval,
                                        # eval.resp.xy = cbind(as.numeric(SpecDatEval$x), as.numeric(SpecDatEval$y)),
                                        PA.nb.rep = 2,
                                        PA.strategy="random",
                                        PA.nb.absences=10*nrow(subset(SpecDat, PA==1)))
    
    #View data (This just histograms the climate distribution of presences vs. pseudoabsences, to see whether/how much they differ)
    #SpecDatShrink <- melt(as.data.table(SpecDat[,c("PA", "PCA1", "PCA2", "PCA3")]), id.vars="PA", value.name = "PCA", variable.name = "Axis")
    # ggplot(SpecDatShrink)+
    #   geom_histogram(aes(x=PCA, y = ..ncount..), bins=100)+
    #   facet_grid(PA ~ Axis)+
    #   theme_classic()
    
    ### Run BioMod Model
    
    #Specify where maxent is stored
    myBiomodOption <- BIOMOD_ModelingOptions(
      MAXENT.Phillips = list( path_to_maxent.jar = MaxentPath,
                              memory_allocated = 2048))
    
    #Delete any old models (if we over writing)
    sapply(list.files(path=Spec, recursive = TRUE, full.names = TRUE), function(x) unlink(x, recursive = TRUE))
    sapply(list.dirs(path=Spec, recursive = TRUE, full.names = TRUE), function(x) unlink(x, recursive = TRUE))
    
    #Run the model (Cut out GBM because it takes FOREVER)
    BioModModel <- BIOMOD_Modeling(InputBioDat, 
                                   bm.options=myBiomodOption, 
                                   models = c("CTA","SRE", "FDA", "RF", "MAXENT.Phillips"),
                                   nb.rep=2, data.split.perc=80, metric.eval=c("TSS"))
    
    #Save model output
    save(BioModModel, file=paste0(FP, "/SDMs/", Spec, "/BioModModels.RData"))
    save(InputBioDat, file=paste0(FP, "/SDMs/", Spec, "/InputBioDat.RData"))
  }
  if(!file.exists(paste0(FP, "/SDMs/", Spec, "/BioModEnsemble.RData"))){
    load(file=paste0(FP, "/SDMs/", Spec, "/BioModModels.RData")) #Load the models
    load(file=paste0(FP, "/SDMs/", Spec, "/InputBioDat.RData")) #Load the input data
    
    ### Evaluate models (We'll remove any with TSS < 0.7) ### 
    #Create a dataframe of all models above 0.7
    BioModEval <- get_evaluations(BioModModel)
    BioModEvalSensible <- as.data.frame(BioModEval["TSS","Testing.data",,,])
    BioModEvalSensible$Model <- rownames(BioModEvalSensible)
    BioModEvalSensible <- melt(as.data.table(BioModEvalSensible), id.vars = c("Model"), variable.name="Run", value.name="TSS")
    BioModEvalSensible <- BioModEvalSensible[complete.cases(BioModEvalSensible),]
    BioModEvalSensible$ModNames <- get_built_models(BioModModel)
    BioModEvalSensible <- subset(BioModEvalSensible, TSS > 0.7)
    
    #If none of them fit, write out "nofit"
    if(nrow(BioModEvalSensible)==0){
      print(c("no fitted models ", Spec))
      write.csv(NULL, paste0("Log/NoFit/", Spec, ".csv"))
      write.csv(NULL, paste0("Log/", Spec, "_Finish.csv"))
      return(NULL)
    }
    
    #How many PAs are there? (It'll determine how we do the ensembl run)
    NPARuns <- length(unique(gsub("PA", "", lapply(str_split(BioModEvalSensible$ModNames, "[_]"), function(x) x[grepl("PA", x)]))))
    em.by.Spec <- ifelse(NPARuns==1, "PA_dataset", "all") 
    
    #Run the ensemble model
    BioModEnsemble <- BIOMOD_EnsembleModeling(bm.mod = BioModModel,
                                              models.chosen = BioModEvalSensible$ModNames,
                                              em.by = em.by.Spec, #
                                              metric.select =  c("TSS"),
                                              metric.eval = c("TSS"),
                                              metric.select.thresh = c(0.7),
                                              #prob.mean = T,
                                              #prob.cv=T,
                                              #prob.ci=T,
                                              #committee.averaging=T,
                                              prob.mean.weight=T)
    
    #Now evaluate the ensemble, and get the final model (it'll be "mean" ['mean'] for when we only had one individual model that fit, and "wmean" ['weighted mean'] for when we have more than one)
    BioModEMEval <- get_evaluations(BioModEnsemble)
    BioModEMEval <- as.data.frame(t(BioModEMEval[1,,]))
    BioModEMEval$Model <- rownames(BioModEMEval)
    if(nrow(BioModEvalSensible)==1){
      BioModEMEval <- BioModEMEval[grepl("EMmean", BioModEMEval$Model),]
    } else {
      BioModEMEval <- BioModEMEval[grepl("EMwmean", BioModEMEval$Model),]
    }
    
    if(BioModEMEval$Testing.data <= 0.7){
      print(c("Ensemble doesn't fit ", Spec))
      write.csv(NULL, paste0("NoFit/", Spec, ".csv"))
      write.csv(NULL, paste0(Spec, "_Finish.csv"))
      return(NULL)
    } 
    modelstokeep <- BioModEMEval$Model
    save(BioModEnsemble, file=paste0(FP, "/SDMs/", Spec, "/BioModEnsemble.RData"))
    save(modelstokeep, file=paste0(FP, "/SDMs/", Spec, "/modelstokeep.RData"))
    save(BioModEvalSensible, file=paste0(FP, "/SDMs/", Spec, "/modelsforensemble.RData"))
  }
  
  ### Hindcast distributions ###
  dir.create(paste0(Spec, "/HindcastRasters/"), showWarnings = FALSE, recursive = TRUE) #Create somewhere for the hindcast rasters to live
  #sapply(list.files(paste0(Spec, "/HindcastRasters/"), full.names=TRUE), unlink)
  
  load(file=paste0(FP, "/SDMs/", Spec, "/BioModEnsemble.RData")) #Load ensemble models
  load(file=paste0(FP, "/SDMs/", Spec, "/modelstokeep.RData")) #Load the name of the ensemble model that fit
  load(file=paste0(FP, "/SDMs/", Spec, "/modelsforensemble.RData")) #Load the names of the individual models
  load(file=paste0(FP, "/SDMs/", Spec, "/BioModModels.RData")) #Load the actual individual models
  load(file=paste0(FP, "/SDMs/", Spec, "/InputBioDat.RData")) #Load initial input data
  
  #Adjust the species name in the fossil data so it matches the gbif data
  FosPolSpec <- subset(FosPol, Species==gsub("[.]", " ", Spec))
  
  #Get years to project to (we will project to every 100 years between the 12.5-87.5 uncertainty brackets for the pollen ages)
  AlltheAges <- unique(unlist(lapply(1:nrow(FosPolSpec), function(x) seq(round_any(FosPolSpec[x,]$x12_5_percent, 100), round_any(FosPolSpec[x,]$x87_5_percent, 100), 100)))) #
  
  #List all the years for which we have past climate data
  PastYears <- list.dirs(paste0(FP, "PCAClims/Past"), full.names=FALSE, recursive=FALSE)
  PastYears <- PastYears[PastYears %in% AlltheAges] #Sub that down to the ages we want for the pollen data
  PastYears <- c(PastYears, names(CurrStacks)) #Add current data to the past year names (Cos we wanna project to current and past)
  
  #Now loop through each time step
  pblapply(PastYears, function (YR){
    print(YR)
    if(file.exists(paste0(FP, "/SDMs/", Spec, "/HindcastRasters/Hindcast_", YR, ".tif"))){
      return("Done")
    }
    
    #Get climate data (past or current)
    if(YR %in% names(CurrStacks)){
      ProjEnv <- CurrStacks[[YR]]
    } else {
      ProjEnv <- rast(paste0(FP, "PCAClims/Past/", YR, "/PCA", c(1:3), ".tif"))
    }
    
    #Project each individual model (of the ones over TSS 0.7)
    myBiomodProj <- BIOMOD_Projection(
      bm.mod = BioModModel, ## The object containing the SDMs. Defined above.
      new.env = stack(ProjEnv), ## Environmental data covering the region where the species' range will be projected
      proj.name = "HindCastTemp",
      #proj.name = paste0(Spec, '_Full'), ## Can project either models made with all presence / pseudo-absence data in a dataset ('_FULL') or with the calibration data only ('_RUN1', '_RUN2' etc...)
      models.chosen = BioModEvalSensible$ModNames, ## Defined above.
      compress="xz",
      build.clamping.mask = T, ## A mask would identify locations where predictions are uncertain because the values of the environmental variables are outside the range used for calibrating the models
      output.format = '.grd' ## Format in which to save the outputted rasters
    )
    
    #Ensemble the projections
    myBiomodEF <- BIOMOD_EnsembleForecasting(
      bm.em = BioModEnsemble,
      bm.proj = myBiomodProj,
      models.chosen = modelstokeep
      #binary.meth = 'KAPPA', ## Whether the continuous projections are thresholded into suitable/unsuitable (atm, not)
    )
    
    #REMEMBER. THE OUTPUT IS GIVEN IN PROBABILITY OF OCCURRENCE*1000 so e.g. 500 = 0.5 probability of occurrence
    ThresholdedMap <- rast(list.files(path = paste0(FP, "/SDMs/", Spec, "/proj_HindcastTemp/"), pattern="*ByTSS_mergedAlgo_mergedRun(.*).grd", recursive=TRUE, full.names = TRUE))
    terra::writeRaster(ThresholdedMap, paste0(FP, "/SDMs/", Spec, "/HindcastRasters/Hindcast_", YR, ".tif"), overwrite=TRUE)
    
    ClampingMask <- rast(list.files(path = paste0(FP, "/SDMs/", Spec), pattern="*_ClampingMask.grd", recursive=TRUE, full.names = TRUE))
    terra::writeRaster(ClampingMask, paste0(FP, "/SDMs/", Spec, "/HindcastRasters/Hindcast_", YR, "_Clamp.tif"), overwrite=TRUE)
    
    return("Done")
  })
  write.csv(NULL, paste0("Log/", Spec, "_Finish.csv"))
}, mc.cores=ncores)
EndTime <- Sys.time()
print(EndTime)
TimeDiff <- EndTime - StartTime
print(TimeDiff)
write.csv(TimeDiff, paste0("TimeDiff", Computer, ".csv"))

