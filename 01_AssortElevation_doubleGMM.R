### General steps for reading in data and running double GMMevents

# Here we have a list of mountain chickadee feeder visitation datasets
# This script runs two Gaussian mixture models at different resolutions 
# output is a group by individual matrix which will later be used to calculate association matrices and run some neat (lengthy and perhaps clunky) analyses

# **** This can be VERY time intensive to run, resulting files are provided in this repository ****

# Libraries used: ----
library(asnipe)

# Double GMM function from VKH:----
dblGMM <- function(meta_events, ids, data) { # Metadata from first GMM (in minutes), global ids, raw observation file
  
  gmmList  <- vector("list", length = nrow(meta_events)) # Store output GMM
  gbiList  <- vector("list", length = nrow(meta_events)) # Extract GBI from GMM and store
  metaList <- vector("list", length = nrow(meta_events)) # Extract Metadata from GMM and store
  
  # We want to loop over all the time slices identified by minute GMM
  for (i in 1:nrow(meta_events)) {
    date_event <- subset(data, data$Date_Loc == meta_events[i,3]) # Date_Loc will be 3rd col in metadata, this will subset obs data by Date_Loc that match this time slice
    date_event <- subset(date_event, date_event$Time_min >= meta_events[i,1] & date_event$Time_min <= meta_events[i,2]) # Further subsetting our observed data to only fit the time window identified by the minute interval GMM (and current slice we are on in our loop! We use start [i,1] and stop [i,2] times to identify this time interval)
    
    if(nrow(date_event) <= 1) { # If the event is too short, you will get an error, so skip the event in this case...
      next()
    }
    
    # Now that everything is set up, we can run our GMM on this slice! (1 s resolution)
    
    gmmE = NULL
    
    try({
      gmmE <- gmmevents(time = date_event$Time_sec, 
                        identity = date_event$ID, 
                        location = date_event$Date_Loc, 
                        global_ids = ids)
    })
    
    gmmList[[i]]  <- gmmE
    gbiList[[i]]  <- gmmE$gbi
    metaList[[i]] <- gmmE$metadata
    
  }
  
  gmmList  <<- gmmList
  gbiList  <<- gbiList
  metaList <<- metaList
}

# It's silly, but I'm going to do this a million times...

# Load in observation data ------

dat.path <- "./icloud/LabThings/MOCHData/AnimalSocialNetworks/Assort/Cleaned/"
dat.file <- "Assort_ObservationData_Pretesting.RData"

load(paste0(dat.path, dat.file))

# let's specify some other paths to store GMM results & GBI files
gbi.path <- "./icloud/LabThings/MOCHData/AnimalSocialNetworks/Assort/Cleaned/GBIs/"
# gmm.path <- "./icloud/LabThings/MOCHData/AnimalSocialNetworks/Assort/Cleaned/GMMs/"

# 2) Use lapply to go over our list of bird feeder visitation data "pretest.obs" ------
# This could be a bit more fuzzy, but I'm going over and editing this script way after the fact and don't want to lose track of objects at this point!

lapply(names(pretest.obs), function(x){
  require(asnipe)
  
  birds          <- pretest.obs[[x]]
  birds$Time_min <- floor(birds$Time_sec/60) + 1 # we need this column to do the first GMM
  
  #Early versions of this script had files save using a convention that I later changed. Here I'd rather just keep the old naming system so code matches output files.
  season    <- strsplit(x, split = '[.]')[[1]][1] # pulling season
  elevation <- strsplit(x, split = '[.]')[[1]][2] # pulling elevation
  
  # singGMMfile <- paste0("minGMM_pre_", elevation, season, ".RData")
  # dubGMMfile  <- paste0("dubGMM_pre_", elevation, season, ".RData")
  dubGBIsfile <- paste0("GBIs_pre_"  , elevation, season, ".RData")
  
  global_ids <- unique(birds$ID)
  
  #) Run single GMM the way you normally would using min col for time.
  gmm <- gmmevents(time     = birds$Time_min, 
                   identity = birds$ID, 
                   location = birds$Date_Loc, 
                   global_ids = global_ids)
  
  # save(gmm, global_ids, file = paste0(gmm.path, singGMMfile)) # definitely don't *need* to save, but I do just in case...
  
  dblGMM(gmm$metadata, global_ids, birds)

  # save(gmmList, global_ids, file = paste0(gmm.path, dubGMMfile)) # definitely don't *need* to save, but I do just in case...
  
  gbi  <- do.call("rbind", gbiList)
  meta <- do.call("rbind", metaList)
  
  row.names(meta) <- 1:nrow(meta)
  row.names(gbi)  <- 1:nrow(gbi)

  tmp <- strsplit(meta$Location, "_")
  tmp <- do.call("rbind", tmp)

  meta$Location  <- tmp[,2]
  meta$Date      <- tmp[,1]
  meta$Season    <- season
  meta$Elevation <- elev
  
  save(gbi, meta, global_ids, file = paste0(gbi.path, dubGBIsfile))
  
})
