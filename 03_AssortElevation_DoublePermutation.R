# Third script in MOCH spatial cognition & social behavior: 
# Following network construction, run double permutation method of hypothesis testing

# Pseudocode:----
# 1) Run DS permutations 5 times and store 5 reps of DS permutation results 
# 2) Node Permutation Part I: Run node permutations on cognitive testing data
# 3) Node Permutation Part II <- very clunky, essentially appending permuted values to adjusted value DFs...
# 4) Node Permutation Part III: Run LMMs on the permuted data to calculate a gosh darn p-value

# Libraries used: -----
  library(asnipe)

# Functions used: -----

# Datastream permutation accounting for a burn-in period:
# takes in lists of data, corresponding to each of the dates included in our analyses

get.daily.rand.nets <- function(GBIs, NETWORKs, METAs, DATEs, permutations, returns){
  require(asnipe)
  sp_rand <- lapply(DATEs, function(x){
    rand.out <- network_permutation(association_data   = GBIs[[x]],
                                    association_matrix = NETWORKs[[x]],
                                    locations          = METAs[[x]][["Location"]],
                                    within_location    = T,
                                    permutations       = permutations,
                                    returns            = returns)
    rand.out <- rand.out[((500/returns)+1):(permutations/returns),,]
  })
  names(sp_rand) <- DATEs
  return(sp_rand)
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

get.daily.ind.expected.values <- function(GBIs, NETWORKs, NETWORKs_RAND, DATEs, INDs, MOVER_DAT){
  require(sna)
  out <- lapply(as.list(DATEs), function(x){
    locs.date <- data.frame(RFID      = rownames(NETWORKs[[x]]),
                            locations = MOVER_DAT$ArraysVisited[match(rownames(NETWORKs[[x]]), MOVER_DAT$RFID)])
    
    r.wd  <- apply(NETWORKs_RAND[[x]], 1, function(z){sna::degree(z, gmode = "graph", diag = F, ignore.eval = F)})
    r.cv  <- lapply(as.list(1:dim(NETWORKs_RAND[[x]])[1]), function(z){
      network <- NETWORKs_RAND[[x]][z,,]
      rownames(network) <- rownames(NETWORKs[[x]])
      colnames(network) <- colnames(NETWORKs[[x]])
      return(get.ind.ewcv(network, loc.vis = locs.date))
    })
    
    r.cv  <- do.call('cbind', r.cv)
    INDs[[x]]$e.w.deg   <- apply(r.wd,  1, median)
    INDs[[x]]$e.cv.ew   <- apply(r.cv,  1, function(y){median(y, na.rm = T)})
    
    INDs[[x]]$adj.w.deg   <- INDs[[x]]$w.deg   - INDs[[x]]$e.w.deg
    INDs[[x]]$adj.cv.ew   <- INDs[[x]]$cv.ew   - INDs[[x]]$e.cv.ew
    
    INDs[[x]]$cv.ew    [which(INDs[[x]]$in.sps == "N")] <- NA
    INDs[[x]]$e.cv.ew  [which(INDs[[x]]$in.sps == "N")] <- NA
    INDs[[x]]$adj.cv.ew[which(INDs[[x]]$in.sps == "N")] <- NA
    
    return(INDs[[x]])
  })
  names(out)  <- DATEs
  return(out)
}


# Read in files ------
dat.path <- "icloud/LabThings/MOCHData/AnimalSocialNetworks/Assort/Cleaned/"
out.dir <- "FiveRep_LMs/" # directory name to house these files of permuted networks etc... just to keep everything somewhat contained and separated from scripts and other data files
dat      <- "Assort_Dat_1dSP.RData"
ind.dat  <- "Assort_IndAtt.RData"

load(paste0(dat.path, dat))
load(paste0(dat.path, ind.dat))

rm(dat, ind.dat)

# 1) Run DS permutations 5 times and store 5 reps of DS permutation results -----
# *** These files are NOT available because they are far too large to upload to repository.
# *** computing time is rather large, but feel free to run with fewer permutations!
permfile.date <- "2022_0915_BURN-" # this will be added as the prefix of RData files that contain all of the permuted/randomized networks

nums <- list("01", "02", "03", "04", "05")


lapply(nums, function(y, 
                      outfile.date = permfile.date, 
                      out.path = paste0(dat.path, out.dir)){
  elev.rand.daily.nets <- lapply(names(as.gbi.elev), function(x){
                          output <- get.daily.rand.nets(GBIs         = as.gbi.elev[[x]],
                                                        NETWORKs     = as.net.elev[[x]],
                                                        METAs        = as.met.elev[[x]],
                                                        DATEs        = as.list(names(as.gbi.elev[[x]])),
                                                        permutations = 10500,
                                                        returns      = 10)
    return(output)
  })
  names(elev.rand.daily.nets) <- names(as.gbi.elev)
  save(elev.rand.daily.nets, file = paste0(out.path,"AssortElev_RandDailyNets_", outfile.date, y, ".RData"))
  rm(elev.rand.daily.nets)
})

# Yay! the lengthy DS perms are done!

#
# Now we calculate "expected values" based on the values from the permuted networks (previous step)
# Again, this is run 5 times

rand.nets.file <- paste0("AssortElev_RandDailyNets_", permfile.date)
adjfile.name   <- paste0("AssortElev_ADJVals_", permfile.date)

nums <- list("01", "02", "03", "04", "05")

lapply(nums, function(y                                    ,
                      in.path        = paste0(dat.path, out.dir),
                      infile.prefix  = rand.nets.file      ,
                      outfile.prefix = adjfile.name        ,
                      gbis           = as.gbi.elev         ,
                      nets           = as.net.elev         ,
                      rand.nets      = elev.rand.daily.nets,
                      inds           = as.inds.elev        ,
                      movs           = movers              ){
  
  load(paste0(in.path, infile.prefix, y, ".RData")) # opening and closing RData files from the previuos step to calculate expected & ultimately adjusted values
  
  elev.adj.inds <- lapply(names(gbis), function(x){
                     output <- get.daily.ind.expected.values(GBIs          = gbis[[x]],
                                                             NETWORKs      = nets[[x]],
                                                             NETWORKs_RAND = rand.nets[[x]],
                                                             DATEs         = names(gbis[[x]]),
                                                             INDs          = inds[[x]],
                                                             MOVER_DAT     = movs[[x]])
                     return(output)
  })
  
  names(elev.adj.inds) <- names(gbis)
  
  elev.all.adj.inds    <- lapply(elev.adj.inds, function(x){
    return(do.call('rbind', x))
  })
  
  save(elev.adj.inds, elev.all.adj.inds, file = paste0(in.path, adjfile.name, y, ".RData"))
})

# 2) Node Permutation Part I: Run node permutations on cognitive testing data and save results ------

# Come up with a name for a file to save out node permutations *** File provided in Repository
np.file <- "NPDFs_2022_0915"

n.perms <- lapply(as.inds.all.elev, function(dat,
                                             ids              = "RFID"         ,
                                             array.loc        = "ArrayAssigned",
                                             initial.20.var   = "I.20.LocErr"  ,
                                             reversal.20.var  = "R.20.LocErr"  ,
                                             initial.tot.var  = "I.Tot.LocErr" ,
                                             reversal.tot.var = "R.Tot.LocErr" ,
                                             n.perm           = 1000){
  
  require(dplyr)
  
  initial.vars  <- c(initial.20.var , initial.tot.var )
  reversal.vars <- c(reversal.20.var, reversal.tot.var)
  
  tmp.I <- unique(subset(dat, !is.na(dat[[initial.20.var]]), 
                         select = c(ids, array.loc, "ArraysVisited", initial.vars)))
  tmp.R <- unique(subset(dat, !is.na(dat[[reversal.20.var]]) & !is.na(dat[[initial.20.var]]), 
                         select = c(ids, array.loc, "ArraysVisited", reversal.vars)))
  
  tmp.I$all.i <- do.call('paste', c(tmp.I[,c(initial.vars, array.loc)], sep = '_'))
  tmp.R$all.r <- do.call('paste', c(tmp.R[,reversal.vars], sep = '_'))
  
  output <- list()
  
  set.seed(20)
  for (j in 1:n.perm) {
    # swap around the testing performance data 
    #*** find out where birds tested, swap only within testing locations...
    tmp.i <- tmp.I %>% dplyr::group_by(ArrayAssigned) %>% dplyr::mutate(rand.i = sample(all.i))
    tmp.r <- tmp.R %>% dplyr::group_by(ArrayAssigned) %>% dplyr::mutate(rand.r = sample(all.r))
    
    # now need to apply a matching function across all the days sampled so we can do repeated measures

    tmp <- tmp.i
    tmp$rand.r <- tmp.r$rand.r[match(tmp$RFID, tmp.r$RFID)]
    output[[j]] <- tmp
  }
  return(output)
})

save(n.perms, file = paste0(dat.path, np.file, ".RData"))

# 3) Node Permutation Part II: ------
# Output files from this step are far too large to upload to repository! NOT provided!

# let's just get rid of everything that isn't our node perm data
rm(list = ls()[which(ls() != "n.perms")])

# re-establish our path
dat.path <- "icloud/LabThings/MOCHData/AnimalSocialNetworks/Assort/Cleaned/"
out.dir <- "FiveRep_LMs/" # directory name to house these files of permuted networks etc... just to keep everything somewhat contained and separated from scripts and other data files


permfile.date <- "2022_0915_BURN-"

adjfile.name <- paste0("AssortElev_ADJVals_",   permfile.date)
npdf.out     <- paste0("AssortElev_NPadjVals_", permfile.date)

nums <- list("01", "02", "03", "04", "05")

lapply(nums, function(num,
                      np.df          = n.perms  ,
                      in.path        = paste0(dat.path, out.dir),
                      infile.prefix  = adjfile.name,
                      outfile.prefix = npdf.out   ){
  
  load(paste0(in.path, infile.prefix, num, ".RData"))
  
  adj.vals.n.perms <-lapply(names(np.df), function(x, daily.dat = elev.adj.inds){
    rand <- np.df[[x]]
    IND  <- daily.dat[[x]]
    output2 <- lapply(rand, function(y){
      output3 <- do.call("rbind", lapply(IND, function(z){
                  df <- z
                  df$tmp.i <- y$rand.i[match(df$RFID, y$RFID)]
                  df$tmp.r <- y$rand.r[match(df$RFID, y$RFID)]
                  
                  tmp <- strsplit(df$tmp.i, split = "_")
                  tmp <- do.call("rbind", tmp)
                  df$rand.20.i        <- as.numeric(tmp[,1]) # Initial: mean locErr first 20 trials
                  df$rand.tot.i       <- as.numeric(tmp[,2]) # Initial: mean locErr total trials
                  rm(tmp)
                  
                  tmp <- strsplit(df$tmp.r, split = "_")
                  tmp <- do.call("rbind", tmp)
                  df$rand.20.r        <- as.numeric(tmp[,1]) # Reversal: mean locErr first 20 trials
                  df$rand.tot.r       <- as.numeric(tmp[,2]) # Reversal: mean locErr of total trials
                  rm(tmp)
                  
                  return(df)
                })
      )
      return(output3)
    })
    return(output2)
  })
  names(adj.vals.n.perms) <- names(np.df)
  save(adj.vals.n.perms, file = paste0(in.path, outfile.prefix, num, ".RData"))
})


rm(list = ls())

# 4) Node Permutation Part III: Run LMMs on the permuted data to calculate a gosh darn p-value -----

# Output files from this step are available in repository! :)

dat.path <- "icloud/LabThings/MOCHData/AnimalSocialNetworks/Assort/Cleaned/"
out.dir <- "FiveRep_LMs/" # directory name to house these files of permuted networks etc... just to keep everything somewhat contained and separated from scripts and other data files

permfile.date <- "2022_0915_BURN-"

adjfile.name <- paste0("AssortElev_ADJVals_",   permfile.date)
npdf.out     <- paste0("AssortElev_NPadjVals_", permfile.date)
outfiles     <- paste0("AssortElev_DPresults_", permfile.date)

nums <- list("01", "02", "03", "04", "05")

lapply(nums, function(num, 
                      in.path        = paste0(dat.path, out.dir),
                      infile.prefix  = npdf.out,
                      outfile.prefix = outfiles   ){
  require(lme4)
  require(lmerTest)
  
  load(paste0(in.path, infile.prefix, num, ".RData"))
  
  np.lm.results <- lapply(names(adj.vals.n.perms), function(x, rand.dat = adj.vals.n.perms){
    
    # Some functions:
    p_value          <- function(observed, random) {
      sum(observed <= random)/length(random)
    }
    two.tail.p_value <- function(x){
      if(x > 0.5){
        x <- (1-x)*2
      }else{
        x <- x*2
      }
      return(x)
    }
    
    seas <- strsplit(x, split = '[.]')[[1]][1]
    elev <- strsplit(x, split = '[.]')[[1]][2]
    
    rands <- rand.dat[[x]] # could cut this line, but I might get confused down the line....
    
    rand.coefs <- list(do.call("c", lapply(rands, function(y){
      coef(summary(lmer(adj.w.deg ~   rand.20.i + ArraysVisited + (1|RFID) + (1|date), data = subset(y, in.sps == "Y"))))[2,1]
    })), 
                       do.call("c", lapply(rands, function(y){
      coef(summary(lmer(adj.cv.ew ~   rand.20.i + ArraysVisited + (1|RFID) + (1|date), data = subset(y, in.sps == "Y"))))[2,1]
    })),
                       do.call("c", lapply(rands, function(y){
      coef(summary(lmer(adj.w.deg ~   rand.20.r + ArraysVisited + (1|RFID) + (1|date), data = subset(y, in.sps == "Y"))))[2,1]
    })),
                       do.call("c", lapply(rands, function(y){
      coef(summary(lmer(adj.cv.ew ~   rand.20.r + ArraysVisited + (1|RFID) + (1|date), data = subset(y, in.sps == "Y"))))[2,1]
    })),
                       do.call("c", lapply(rands, function(y){
      coef(summary(lmer(adj.w.deg ~   rand.tot.i + ArraysVisited + (1|RFID) + (1|date), data = subset(y, in.sps == "Y"))))[2,1]
    })),
                       do.call("c", lapply(rands, function(y){
      coef(summary(lmer(adj.cv.ew ~   rand.tot.i + ArraysVisited + (1|RFID) + (1|date), data = subset(y, in.sps == "Y"))))[2,1]
    })),
                       do.call("c", lapply(rands, function(y){
      coef(summary(lmer(adj.w.deg ~   rand.tot.r + ArraysVisited + (1|RFID) + (1|date), data = subset(y, in.sps == "Y"))))[2,1]
    })),
                       do.call("c", lapply(rands, function(y){
      coef(summary(lmer(adj.cv.ew ~   rand.tot.r + ArraysVisited + (1|RFID) + (1|date), data = subset(y, in.sps == "Y"))))[2,1]
    })))
    
    dat <- rands[[1]]
    dat <- subset(dat, in.sps == "Y")
    
    # models WITH adjusted values...not using NPerm values here
    adj.lms <- list(lmer(adj.w.deg ~ I.20.LocErr  + ArraysVisited + (1|RFID) + (1|date), data = dat),
                    lmer(adj.cv.ew ~ I.20.LocErr  + ArraysVisited + (1|RFID) + (1|date), data = dat),
                    lmer(adj.w.deg ~ R.20.LocErr  + ArraysVisited + (1|RFID) + (1|date), data = dat),
                    lmer(adj.cv.ew ~ R.20.LocErr  + ArraysVisited + (1|RFID) + (1|date), data = dat),
                    lmer(adj.w.deg ~ I.Tot.LocErr + ArraysVisited + (1|RFID) + (1|date), data = dat),
                    lmer(adj.cv.ew ~ I.Tot.LocErr + ArraysVisited + (1|RFID) + (1|date), data = dat),
                    lmer(adj.w.deg ~ R.Tot.LocErr + ArraysVisited + (1|RFID) + (1|date), data = dat),
                    lmer(adj.cv.ew ~ R.Tot.LocErr + ArraysVisited + (1|RFID) + (1|date), data = dat))
    
    # models WITHOUT adjusted values
    obs.lms <- list(lmer(w.deg ~ I.20.LocErr  + ArraysVisited + (1|RFID) + (1|date), data = dat),
                    lmer(cv.ew ~ I.20.LocErr  + ArraysVisited + (1|RFID) + (1|date), data = dat),
                    lmer(w.deg ~ R.20.LocErr  + ArraysVisited + (1|RFID) + (1|date), data = dat),
                    lmer(cv.ew ~ R.20.LocErr  + ArraysVisited + (1|RFID) + (1|date), data = dat),
                    lmer(w.deg ~ I.Tot.LocErr + ArraysVisited + (1|RFID) + (1|date), data = dat),
                    lmer(cv.ew ~ I.Tot.LocErr + ArraysVisited + (1|RFID) + (1|date), data = dat),
                    lmer(w.deg ~ R.Tot.LocErr + ArraysVisited + (1|RFID) + (1|date), data = dat),
                    lmer(cv.ew ~ R.Tot.LocErr + ArraysVisited + (1|RFID) + (1|date), data = dat))
    
    results <- data.frame(season            = seas,
                          elevation         = elev,
                          dat.type          = rep("1dSP Day Random, Elevation", times = length(obs.lms)),
                          test.type         = rep(c("initial.20" , "reversal.20"  ,
                                                    "initial.tot", "reversal.tot"), each = 2),
                          response          = rep(c("adj.w.deg", "adj.ew.cv"), times = 4),
                          num.perms         = length(rands),
                          first.qu.rand     = do.call('c', lapply(rand.coefs, function(y){summary(y)[[2]]})),
                          med.coef.rand     = do.call('c', lapply(rand.coefs, function(y){median(y)})),
                          third.qu.rand     = do.call('c', lapply(rand.coefs, function(y){summary(y)[[5]]
                          })),
                          adj.coef          = do.call('c', lapply(adj.lms, function(y){
                            coef(summary(y))[2,1]
                          })),
                          tail.p.value      = do.call('c', lapply(as.list(seq(1:length(rand.coefs))), function(y){
                            p_value(coef(summary(adj.lms[[y]]))[2,1], rand.coefs[[y]])
                          })),
                          full.dp.p.val     = "",
                          r.sig             = "",
                          r.sig.level       = "",
                          obs.p.value       = do.call('c', lapply(adj.lms, function(y){coef(summary(y))[2,5]
                          })),
                          obs.p.sig         = "",
                          obs.sig.level     = "",
                          obs.coef.arrvis   = do.call('c', lapply(adj.lms, function(y){
                            coef(summary(y))[3,1]
                          })),
                          obs.p.arrvis      = do.call('c', lapply(adj.lms, function(y){
                            coef(summary(y))[3,5]
                          })),
                          unadj.response    = rep(c("w.deg", "ew.cv"), times = 4),
                          unadj.coef        = do.call('c', lapply(obs.lms, function(y){
                            coef(summary(y))[2,1]
                          })),
                          unadj.t           = do.call('c', lapply(obs.lms, function(y){
                            coef(summary(y))[2,4]
                          })),
                          unadj.se          = do.call('c', lapply(obs.lms, function(y){
                            coef(summary(y))[2,2]
                          })),
                          arrvis.unadj.coef = do.call('c', lapply(obs.lms, function(y){
                            coef(summary(y))[3,1]
                          })),
                          arrvis.unadj.t    = do.call('c', lapply(obs.lms, function(y){
                            coef(summary(y))[3,4]
                          })),
                          arrvis.unadj.se   = do.call('c', lapply(obs.lms, function(y){
                            coef(summary(y))[3,2]
                          })),
                          stringsAsFactors  = F)
    tmp <- as.list(results$tail.p.value)
    tmp <- do.call('c', lapply(tmp, two.tail.p_value))
    results$full.dp.p.val <- tmp
    rm(tmp)
    
    results$r.sig[which(results$tail.p.value <= 0.05)]   <- "."
    results$r.sig[which(results$tail.p.value <= 0.025)]  <- "*"
    results$r.sig[which(results$tail.p.value <= 0.005)]  <- "**"
    results$r.sig[which(results$tail.p.value <= 0.0005)] <- "***"
    
    results$r.sig[which(results$tail.p.value >= 0.95)]   <- "."
    results$r.sig[which(results$tail.p.value >= 0.975)]  <- "*"
    results$r.sig[which(results$tail.p.value >= 0.995)]  <- "**"
    results$r.sig[which(results$tail.p.value >= 0.9995)] <- "***"
    
    
    results$r.sig.level[which(results$tail.p.value <= 0.05)]   <- "0.1"
    results$r.sig.level[which(results$tail.p.value <= 0.025)]  <- "0.05"
    results$r.sig.level[which(results$tail.p.value <= 0.005)]  <- "0.01"
    results$r.sig.level[which(results$tail.p.value <= 0.0005)] <- "0.001"
    
    results$r.sig.level[which(results$tail.p.value >= 0.95)]   <- "0.1"
    results$r.sig.level[which(results$tail.p.value >= 0.975)]  <- "0.05"
    results$r.sig.level[which(results$tail.p.value >= 0.995)]  <- "0.01"
    results$r.sig.level[which(results$tail.p.value >= 0.9995)] <- "0.001"
    
    
    results$obs.p.sig[which(results$obs.p.value <= 0.1)]   <- "."
    results$obs.p.sig[which(results$obs.p.value <= 0.05)]  <- "*"
    results$obs.p.sig[which(results$obs.p.value <= 0.01)]  <- "**"
    results$obs.p.sig[which(results$obs.p.value <= 0.001)] <- "***"
    
    results$obs.sig.level[which(results$obs.p.value <= 0.1)]   <- "0.1"
    results$obs.sig.level[which(results$obs.p.value <= 0.05)]  <- "0.05"
    results$obs.sig.level[which(results$obs.p.value <= 0.01)]  <- "0.01"
    results$obs.sig.level[which(results$obs.p.value <= 0.001)] <- "0.001"
    
    return(results)
  })
  
  
  np.lm.result <- do.call('rbind', np.lm.results)
  save(np.lm.result, np.lm.results, file = paste0(in.path, outfile.prefix, num, ".RData"))
  
})
