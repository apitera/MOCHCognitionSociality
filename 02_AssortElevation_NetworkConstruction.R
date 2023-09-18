# Second script in MOCH spatial cognition & social behavior: subset data given below criteria
# Subsetting will take 2 forms: Daily networks (for double permutation method) & networks calculated from data collected throughout the season (for single, node permutation method)

# ******* All files/datasets created with this script have been provided in the repository:
# "Assort_PreData_Open.RData", "Assort_Dat_1dSP.RData", & "Assort_FullTimeData.RData"

# Daily network subsetting criteria: -----
  # 1-day sampling periods (SPs) but only including days in which # grouping events >= 100 so datastream permutations have something to work with
  # Criteria for inclusion of individuals:
    # Within each sampling period, only include individuals that are observed in at least 5 groups
    # Across SPs, only include individuals present in at least 5 SPs

# Throughout season networks will be built off the same inclusion criteria as the daily networks, but generate a SINGLE network instead of separate networks for each day.

# ====================================================
# Pseudocode: ----
# 1) Compile all the GBI & metadata from the double GMM * Output file available "Assort_PreData_Open.RData"
# 2) Create daily networks, but use all birds present during the period unless they were unbanded during this time... * Output file available "Assort_Dat_1dSP.RData"
# 3) Create networks from all data collected together, but use same inclusion criteria as in step 2. * Output file available "Assort_FullTimeData.RData"

# Libraries used: ------
library(asnipe)
library(igraph)
library(lubridate)
library(sna)

# Functions used: ------

network.cv <- function(x){
  t <- x
  t[lower.tri(t, diag =T)] <- NA # don't want to double up on edges (since we have an undirected network)
  return(sd(t, na.rm = T)/mean(t, na.rm = T))
}

get.ind.ewcv <- function(x, loc.vis){
  t       <- x
  diag(t) <- NA # no self-loops, so make sure to exclude diagonal 
  return(do.call('c', lapply(rownames(t), function(y){
    loc.y <- as.character(loc.vis$locations)[which(loc.vis$RFID == y)]
    if(loc.y == "Both"){
      row <- t[y, ]
    }else{
      locs <- c(loc.y, "Both")
      row <- t[y, which(loc.vis$locations[] %in% locs)]
    }
    b <- (sd(row, na.rm = T)/mean(row, na.rm = T))
    if(is.na(b)){
      b <- 0
    }
    return(b)
  })))
}

# let's make all of our named lists first -----

seas.elev <- list("2015_16.H" = "2015_16.H",
                  "2016_17.H" = "2016_17.H",
                  "2017_18.H" = "2017_18.H",
                  "2018_19.H" = "2018_19.H",
                  "2019_20.H" = "2019_20.H",
                  "2015_16.L" = "2015_16.L",
                  "2016_17.L" = "2016_17.L",
                  "2019_20.L" = "2019_20.L")

sub.dates <- list(
                  '2015_16.H' = c("2015-10-08", "2015-11-10"),
                  '2016_17.H' = c("2016-12-31", "2017-01-10"),
                  '2017_18.H' = c("2017-11-30", "2018-01-16"), 
                  '2018_19.H' = c("2019-03-19", "2019-03-28"),
                  '2019_20.H' = c("2019-12-11", "2019-12-26"),
                  '2015_16.L' = c("2015-10-08", "2015-11-10"),
                  '2016_17.L' = c("2017-02-13", "2017-02-28"),
                  '2019_20.L' = c("2019-12-17", "2019-12-29")
                  )

# 1) Compile all the GBI & metadata from the double GMM * Output file available "Assort_PreData_Open.RData"---------

# will need to unzip the GBIs folder and redefine path if you plan to run this code
gbi.path     <- "./icloud/LabThings/MOCHData/AnimalSocialNetworks/Assort/Cleaned/GBIs/"

gmm_data <- lapply(seas.elev, function(x){
  season    <- strsplit(x, split = '[.]')[[1]][1] # pulling season
  elevation <- strsplit(x, split = '[.]')[[1]][2] # pulling elevation
  
  ele.sea    <- paste(elevation, season, sep = "")
  gbi.filename <- paste("GBIs_pre_", ele.sea, ".RData", sep = "")
  load(file = paste0(gbi.path, gbi.filename))
  
  begin <- sub.dates[[x]][1] # beginning date for timeframe 
  end   <- sub.dates[[x]][2] # ending date for timeframe 
  
  gbi  <- subset(gbi, ymd(meta$Date) >= ymd(begin) & ymd(meta$Date) <= ymd(end))
  meta <- subset(meta, ymd(meta$Date) >= ymd(begin) & ymd(meta$Date) <= ymd(end))
  gbi  <- gbi[ , which(colSums(gbi) > 10)]
  meta <- subset(meta, rowSums(gbi) > 0)
  gbi  <- gbi[which(rowSums(gbi) > 0), ]
  
  outs       <- list()
  outs$gbis  <- gbi
  outs$metas <- meta
})
raw.elev.gbis <- lapply(gmm_data, function(x){
  return(x$gbis)
})

raw.elev.metas <- lapply(gmm_data, function(x){
  return(x$metas)
})

names(raw.elev.gbis)  <- seas.elev
names(raw.elev.metas) <- seas.elev
rm(gmm_data)

# This file has already been created and is available in this repository:
# notice that I have saved this file outside of the GBIs directory
save(seas.elev, raw.elev.gbis, raw.elev.metas, file = "./icloud/LabThings/MOCHData/AnimalSocialNetworks/Assort/Cleaned/Assort_PreData_Open.RData")

rm(gbi.path)

# 2) Create daily networks, but use all birds present during the period unless they were unbanded during this time... * Output file available "Assort_Dat_1dSP.RData"------

# path where the above .RData file lives AND contains the "Assort_IndAtt.RData" file
dat.path <- "./icloud/LabThings/MOCHData/AnimalSocialNetworks/Assort/Cleaned/"

dat      <- "Assort_PreData_Open.RData" # just to read in the file created in step 1
ind.dat  <- "Assort_IndAtt.RData" # available in repository

load(paste0(dat.path, dat))
load(paste0(dat.path, ind.dat))

rm(dat, ind.dat)

# birds$TOIC <- as.Date(birds$TOIC, "%d-%b-%Y") # make this so it can be matched with meta date format

outputs <- lapply(seas_elev, function(x, 
                                      gbi.list          = raw.elev.gbis, 
                                      meta.list         = raw.elev.metas, 
                                      shared.date.list  = shared.dates, 
                                      testing.data      = test.sep, 
                                      locations.visited = movers,
                                      bird.dat          = birds){
  season    <- strsplit(x, split = '[.]')[[1]][1] # pulling season
  elevation <- strsplit(x, split = '[.]')[[1]][2] # pulling elevation
  
  gbis    <- gbi.list[[x]]
  metas   <- meta.list[[x]]
  dates   <- shared.date.list[[x]]
  testing <- testing.data[[x]]
  
  sub.gbi  <- gbis[metas$Date %in% dates, ]
  sub.meta <- metas[metas$Date %in% dates, ]
  sub.gbi  <- sub.gbi[ ,which(colSums(sub.gbi) > 0)]
  sub.meta <- sub.meta[which(rowSums(sub.gbi) > 0), ]
  sub.gbi  <- sub.gbi[which(rowSums(sub.gbi) >0), ]
  
  sub.locs <- locations.visited[[x]]
  
  # just making sure that we only include birds that a) actually learned and b) completed minimum of 20 mean location errors "I.20.LocErr"
  
  testing$I.20.LocErr[which(testing$I.20.LocErr >= 4)]     <- NA
  testing$I.20.LocErr[which(testing$I.TotTrials < 20)]     <- NA
  testing$R.20.LocErr[which(testing$R.TotTrials < 20)]     <- NA
  testing$I.Tot.LocErr[which(is.na(testing$I.20.LocErr))]  <- NA
  testing$R.20.LocErr[which(is.na(testing$I.20.LocErr))]   <- NA
  testing$R.Tot.LocErr[which(is.na(testing$I.Tot.LocErr))] <- NA
  testing$I.TotTrials[which(is.na(testing$I.20.LocErr))]   <- NA
  testing$R.TotTrials[which(is.na(testing$R.20.LocErr))]   <- NA

  # now that we have extracted the GBI & metadata for the entire collection period, subset out each date...
  by.dates <- lapply(dates, function(y){
    sub.g <- sub.gbi[which(sub.meta$Date == y),  ]   # gbi for that day slice
    sub.m <- sub.meta[which(sub.meta$Date == y), ] # metadata for that day slice
    all.sub.g <- apply(sub.g, 2, function(z){
      if(sum(z) < 5) {
        length <- length(z) 
        z <- rep(0, times = length)
      }else{
        z <- z
      }
    })
    sub.g     <- apply(sub.g, 2, function(z){
      if(sum(z) < 5 & sum(z) > 0) { # here we are making a different set of gbi stacks where we single out birds whose daily gbi data have been replaced so we can exclude them from the resulting gbi (but keep birds that just didn't show up)
        length <- length(z) 
        z <- rep(3, times = length)
      }else{
        z <- z
      }
    })
    sub.g     <- sub.g[ , which(colSums(sub.g) <= dim(sub.g)[1])] # throw out all birds whose data have been replaced
    outs <- list()
    outs$tmp.all.g <- all.sub.g
    outs$tmp.g <- sub.g
    outs$tmp.m <- sub.m
    return(outs)
  })
  
  names(by.dates) <- dates

  good.dat.gbi <- do.call('rbind', lapply(by.dates, function(y){
                            return(y[["tmp.all.g"]])
                          }))# making single GBI for subset period data
  
  good.dat.met <- do.call('rbind', lapply(by.dates, function(y){
                            return(y[["tmp.m"]])
                          }))
  
  ids   <- data.frame(ID    = colnames(good.dat.gbi), 
                      n.sps = apply(good.dat.gbi, 2, function(z){length(unique(good.dat.met$Date[which(z == 1)]))}), # number days observed
                      TOIC  = bird.dat$TOIC[match(colnames(good.dat.gbi), bird.dat$RFID)],
                      stringsAsFactors = F)
  
  
  tmp.object.name <- lapply(dates, function(y){
    tmp.gbi  <- by.dates[[y]]$tmp.g
    tmp.meta <- by.dates[[y]]$tmp.m
    tmp.inds <- data.frame(RFID  = colnames(tmp.gbi),
                           n.sps = ids$n.sps[match(colnames(tmp.gbi), ids$ID)], # something funky here?
                           TOIC  = ids$TOIC[match(colnames(tmp.gbi), ids$ID)],
                           stringsAsFactors = F)
    
    tmp.gbi  <- tmp.gbi[ , which(tmp.inds$n.sps >= 5 & ymd(tmp.inds$TOIC) <= ymd(y))] 
    tmp.meta <- tmp.meta[which(rowSums(tmp.gbi) > 0), ]
    tmp.gbi  <- tmp.gbi[which(rowSums(tmp.gbi) > 0), ]
    
    # *** FIX THIS PART...*****
    tmp.net  <- get_network(tmp.gbi)
    
    day.loc <- data.frame(RFID      = rownames(tmp.net),
                          locations = sub.locs$ArraysVisited[match(rownames(tmp.net), sub.locs$RFID)])
    
    b <- data.frame(RFID          = rownames(tmp.net),
                    season        = rep(season, times = dim(tmp.net)[1]),
                    elevation     = rep(elevation, times = dim(tmp.net)[1]), 
                    age           = bird.dat$MaturityTOIC[match(rownames(tmp.net), bird.dat$RFID)],
                    band.yr       = as.numeric(substr(bird.dat$SeasonBanded[match(rownames(tmp.net), bird.dat$RFID)], 1, 4)),
                    diff.yr       = "",
                    ArraysVisited = sub.locs$ArraysVisited[match(rownames(tmp.net), sub.locs$RFID)],
                    in.tot.sp     = ids$n.sps[match(rownames(tmp.net), ids$ID)],
                    tot.sps       = length(dates),
                    in.sps        = apply(tmp.gbi, 2, function(z){ # this way we can differentiate birds that only showed up by themselves from those that weren't present that day
                      if(sum(z) > 0) {
                        return("Y")
                      } else {
                        return("N")
                      }
                    }), 
                    date          = y,
                    n.birds       = length(colnames(tmp.gbi)[which(colSums(tmp.gbi) > 0)]),
                    n.group       = colSums(tmp.gbi),
                    n.events      = rep(dim(tmp.gbi)[1], times = dim(tmp.net)[1]),
                    w.deg         = sna::degree(tmp.net, gmode = "graph", ignore.eval = F),
                    cv.ew         = get.ind.ewcv(tmp.net, loc.vis = day.loc),
                    stringsAsFactors = F)
    
    year                         <- as.numeric(rep(strsplit(season, split = '_')[[1]][1], dim(b)[1]))
    b$diff.yr                    <- year - b$band.yr
    b$age[which(b$diff.yr >= 1)] <- "adult"
    testing$season <- NULL
    testing$elevation <- NULL
    sub.i                        <- dplyr::left_join(b, testing, by = "RFID")
    output <- list()
    output$network <- tmp.net
    output$gbi     <- tmp.gbi
    output$meta    <- tmp.meta
    output$ind.df  <- sub.i
    
    return(output)
  })
  names(tmp.object.name) <- dates
  output <- list()
  output$gbis <- lapply(tmp.object.name, function(y){
    tmp <- y
    return(tmp$gbi)
  })
  
  output$metas <- lapply(tmp.object.name, function(y){
    tmp <- y
    return(tmp$meta)
  })
  
  output$networks <- lapply(tmp.object.name, function(y){
    tmp <- y
    return(tmp$network)
  })
  
  output$daily.dfs <- lapply(tmp.object.name, function(y){
    tmp <- y
    return(tmp$ind.df)
  })
  
  output$inds <- do.call("rbind", lapply(tmp.object.name, function(y){
    tmp <- y
    return(tmp$ind.df)
  }))
  
  return(output)
  
})

# pull everything out...
as.gbi.elev <- lapply(outputs, function(x){
  input <- x
  return(input$gbis)
})

as.met.elev <- lapply(outputs, function(x){
  input <- x
  return(input$metas)
})

as.net.elev <- lapply(outputs, function(x){
  input <- x
  return(input$networks)
})

as.inds.elev <- lapply(outputs, function(x){
  input <- x
  return(input$daily.dfs)
})

as.inds.all.elev <- lapply(outputs, function(x){
  input <- x
  return(input$inds)
})

as.inds.stats.elev <- do.call('rbind', as.inds.all.elev)
rownames(as.inds.stats.elev) <- seq(1, dim(as.inds.stats.elev)[1])

save(as.inds.stats.elev, as.inds.all.elev, as.inds.elev, as.net.elev, as.met.elev, as.gbi.elev, seas.elev,
     file = paste0(dat.path, "Assort_Dat_1dSP.RData"))

# 3) Create networks from all data collected together, but use same inclusion criteria as in step 2. * Output file available "Assort_FullTimeData.RData"------
full.time.elev.dat <- lapply(seas_elev, function(x, 
                                                 gbi.list      = raw.elev.gbis, 
                                                 metadata.list = raw.elev.metas, 
                                                 date.list     = shared.dates, 
                                                 ind.list      = as.inds.all.elev,
                                                 toic.data     = birds){
  season    <- strsplit(x, split = '[.]')[[1]][1] # pulling season
  elevation <- strsplit(x, split = '[.]')[[1]][2] # pulling elevation
  
  gbis  <- gbi.list[[x]]
  metas <- metadata.list[[x]]
  dates <- date.list[[x]]
  inds  <- ind.list[[x]][["RFID"]]
  inds  <- unique(inds)

  sub.gbi  <- gbis [metas$Date %in% dates, ]
  sub.meta <- metas[metas$Date %in% dates, ]
  sub.gbi  <- sub.gbi[ , which(colSums(sub.gbi) > 0)]
  sub.gbi  <- sub.gbi[ , which(colnames(sub.gbi) %in% inds)]
  sub.meta <- sub.meta[which(rowSums(sub.gbi) > 0), ]
  sub.gbi  <- sub.gbi[which(rowSums(sub.gbi) >0),   ]
  
  sub.meta$sub.date <- gsub(pattern = "-", replacement = "", sub.meta$Date)
  sub.meta$sub.date <- as.numeric(sub.meta$sub.date)
  
  toic.dates    <- data.frame(RFID = colnames(sub.gbi),
                              TOIC = as.numeric(toic.data$sub.toic[match(colnames(sub.gbi), toic.data$RFID)]),
                              stringsAsFactors = F)
  
  network <- get_network(association_data = sub.gbi, 
                         times = sub.meta$sub.date, 
                         enter_time = toic.dates$TOIC)
  
  output <- list()
  output$gbi     <- sub.gbi
  output$meta    <- sub.meta
  output$network <- network
  return(output)
  
})

full.time.elev.gbis <- lapply(full.time.elev.dat, function(x){
  return(x$gbi)
})
  
full.time.elev.metas <- lapply(full.time.elev.dat, function(x){
  return(x$meta)
})
  
full.time.elev.networks <- lapply(full.time.elev.dat, function(x){
  return(x$network)
})

# File is readily available in repository
save(full.time.elev.dat, full.time.elev.gbis, full.time.elev.metas, full.time.elev.networks, seas.elev, file = paste0(dat.path, "Assort_FullTimeData.RData"))
