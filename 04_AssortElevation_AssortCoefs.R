## Fourth script in MOCH spatial cognition & social behavior: 
# Running basic node permutation to calculate assortment coefficients for assortment by spatial cognitive performance and a separate run by age.

# Pseudocode: -----
# 1) Assortment by spatial cognitive abilities: Calculate assortment coefficients on FULL time network data, evaluate with a node permutation
# 2) Assortment by age class: Calculate assortment coefficients on FULL time network data, evaluate with a node permutation


dat.path    <- "./icloud/LabThings/MOCHData/AnimalSocialNetworks/Assort/Cleaned/"
ft.file     <- "Assort_FullTimeData.RData" # Has full time networks, gbis, metadata files...
ind.file    <- "Assort_IndAtt.RData"
n.perm.file <- "NPDFs_2022_0915"
outfiles    <- "AssortElev_AssortCoefNP_2022_0928"


load(paste0(dat.path, ft.file))
load(paste0(dat.path, ind.file))
load(paste0(dat.path, n.perm.file, ".RData"))

testing.dat <- test.sep # need to rename!
rm(full.time.elev.dat)

# 1) Use previously run node permutations (from script 03) to evaluate any potential patterns of non-random mixing by spatial cognitive performance --------
assort.n.perms <- lapply(names(full.time.elev.networks), function(x){
    require(assortnet)
    
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
    
    np.dfs  <- n.perms[[x]]
    test.df <- testing.dat[[x]]
    
    # Pull IDs of birds that participated in each type of testing and store those -- will need to subset networks based on who actually tested
    testers.i <- np.dfs[[1]]$RFID[which(!is.na(np.dfs[[1]]$rand.i))]
    testers.r <- np.dfs[[1]]$RFID[which(!is.na(np.dfs[[1]]$rand.r))]
    
    # subset our observed network to only include birds that completed each type of testing
    net.i     <- full.time.elev.networks[[x]][rownames(full.time.elev.networks[[x]]) %in% testers.i,
                                              colnames(full.time.elev.networks[[x]]) %in% testers.i]
    
    net.r     <- full.time.elev.networks[[x]][rownames(full.time.elev.networks[[x]]) %in% testers.r,
                                              colnames(full.time.elev.networks[[x]]) %in% testers.r]
    
    # subset our "median" network to only include birds that completed each type of testing

    rm(testers.i, testers.r)

    
    # the below chunk of code will output a list of four lists, each with 1000 elements
    # each element of the four lists will be a deviation score for our assortment coefficients
    # list[[1]] - initial 20 trials location errors; list[[2]] - initial total location errors; list[[3]] - reversal 20 loc.err.; list[[4]] - reversal total loc.err.
    
    rand.r <- list(do.call("c", lapply(np.dfs, function(np){
      
                                   tmp <- strsplit(np$rand.i, split = "_")
                                   tmp <- do.call("rbind", tmp)
                                   np$rand.20.i        <- as.numeric(tmp[,1]) # Initial: mean locErr first 20 trials
                                   rm(tmp)
                                   rands <- data.frame(RFID = rownames(net.i)                               ,
                                                       rand = np$rand.20.i [match(rownames(net.i), np$RFID)])
                                   return(assortment.continuous(graph = net.i, vertex_values = rands$rand, weighted = T, SE = F)$r)
                             })),
                      
                     do.call("c", lapply(np.dfs, function(np){
                                   tmp <- strsplit(np$rand.i, split = "_")
                                   tmp <- do.call("rbind", tmp)
                                   np$rand.tot.i <- as.numeric(tmp[,2]) # Initial: mean locErr total trials
                                   rands <- data.frame(RFID = rownames(net.i)                    ,
                                                       rand = np$rand.tot.i[match(rownames(net.i), np$RFID)])
                                   return(assortment.continuous(graph = net.i, vertex_values = rands$rand, weighted = T, SE = F)$r)
                             })),
                     
                     do.call("c", lapply(np.dfs, function(np){
                                   tmp          <- strsplit(np$rand.r, split = "_")
                                   tmp          <- do.call("rbind", tmp)
                                   np$rand.20.r <- as.numeric(tmp[,1]) # Reversal: mean locErr first 20 trials
                                   rands         <- data.frame(RFID = rownames(net.r)            ,
                                                       rand = np$rand.20.r [match(rownames(net.r), np$RFID)])
                                   return(assortment.continuous(graph = net.r, vertex_values = rands$rand, weighted = T, SE = F)$r)
                             })),
                     
                     do.call("c", lapply(np.dfs, function(np){
                                   tmp <- strsplit(np$rand.r, split = "_")
                                   tmp <- do.call("rbind", tmp)
                                   np$rand.tot.r       <- as.numeric(tmp[,2]) # Reversal: mean locErr of total trials
                                   rm(tmp)
                                   
                                   rands <- data.frame(RFID = rownames(net.r)                               ,
                                                       rand = np$rand.tot.r[match(rownames(net.r), np$RFID)])
                                   return(assortment.continuous(graph = net.r, vertex_values = rands$rand, weighted = T, SE = F)$r)
                    })))
    
    testers.i <- data.frame(RFID         = rownames(net.i)                                           ,
                            I.20.LocErr  = test.df$I.20.LocErr [match(rownames(net.i), test.df$RFID)],
                            I.Tot.LocErr = test.df$I.Tot.LocErr[match(rownames(net.i), test.df$RFID)])
    
    testers.r <- data.frame(RFID         = rownames(net.r)                                           ,
                            R.20.LocErr  = test.df$R.20.LocErr [match(rownames(net.r), test.df$RFID)],
                            R.Tot.LocErr = test.df$R.Tot.LocErr[match(rownames(net.r), test.df$RFID)])
    
    obs.r <- list(assortment.continuous(graph = net.i, vertex_values = testers.i$I.20.LocErr , weighted = T, SE = T),
                       assortment.continuous(graph = net.i, vertex_values = testers.i$I.Tot.LocErr, weighted = T, SE = T),
                       assortment.continuous(graph = net.r, vertex_values = testers.r$R.20.LocErr , weighted = T, SE = T),
                       assortment.continuous(graph = net.r, vertex_values = testers.r$R.Tot.LocErr, weighted = T, SE = T))
    
    
    results <- data.frame(season         = seas,
                          elevation      = elev,
                          dat.type       = "FullTime Elev-Wide",
                          test.stat      = "assortment coef(r)",
                          test.type      = c("initial.20", "initial.tot", "reversal.20", "reversal.tot"),
                          hyp.test      = "NodePermONLY",
                          num.perms      = length(np.dfs),
                          obs.r.val      = do.call("c", lapply(obs.r, function(y){return(y$r)})),
                          obs.se.val     = do.call("c", lapply(obs.r, function(y){return(y$se)})),
                          quant.05 = do.call("c", lapply(rand.r, function(y){quantile(y, probs = 0.05)[[1]]})),
                          first.qu.rand  = do.call("c", lapply(rand.r, function(y){summary(y)[[2]]})),
                          med.coef.rand  = do.call("c", lapply(rand.r, function(y){median(y)})),
                          third.qu.rand  = do.call("c", lapply(rand.r, function(y){summary(y)[[5]]})),
                          quant.95 = do.call("c", lapply(rand.r, function(y){quantile(y, probs = 0.95)[[1]]})),
                          tail.p.value   = do.call('c', lapply(as.list(seq(1:length(rand.r))), function(y){
                                                                return(p_value(obs.r[[y]]$r, rand.r[[y]]))
                          })),
                          full.dp.p.val  = "",
                          r.sig          = "",
                          r.sig.level    = "",
                          stringsAsFactors = F)
    
    tmp <- as.list(results$tail.p.value)
    tmp <- do.call('c', lapply(tmp, two.tail.p_value))
    results$full.dp.p.val <- tmp
    rm(tmp)
    
    results$r.sig[which(results$full.dp.p.val <= 0.1)]   <- "."
    results$r.sig[which(results$full.dp.p.val <= 0.05)]  <- "*"
    results$r.sig[which(results$full.dp.p.val <= 0.01)]  <- "**"
    results$r.sig[which(results$full.dp.p.val <= 0.001)] <- "***"
    
    results$r.sig.level[which(results$full.dp.p.val <= 0.1)]   <- "0.1"
    results$r.sig.level[which(results$full.dp.p.val <= 0.05)]  <- "0.05"
    results$r.sig.level[which(results$full.dp.p.val <= 0.01)]  <- "0.01"
    results$r.sig.level[which(results$full.dp.p.val <= 0.001)] <- "0.001"
    
    return(results)
})
assort.n.perm  <- do.call('rbind', assort.n.perms)
rownames(assort.n.perm) <- seq(1:dim(assort.n.perm)[1])
save(assort.n.perm, file = paste0(dat.path, outfiles, ".RData"))


# 2) Run node permutations on age data to evaluate any potential patterns of non-random mixing by age class ------------
assort.age <- lapply(names(full.time.elev.networks), function(x, 
                                             network.list = full.time.elev.networks, # list of networks
                                             age.data     = as.inds.all.elev, # list of DFs with Individual attributes
                                             age.var      = "age",
                                             id.var       = "RFID",
                                             incl.vars    = c(id.var, "season", "elevation", "sex", age.var),
                                             n.perm       = 1000){
  require(assortnet)

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
  
  # each item in each list is named based on season + elevation, e.g., "2015_16.H"
  seas <- strsplit(x, split = '[.]')[[1]][1]
  elev <- strsplit(x, split = '[.]')[[1]][2]

  network <- network.list[[x]]
  age.dat <- age.data[[x]]
  age.dat <- subset(age.dat, select = incl.vars)
  age.dat <- unique(age.dat)
  
  inds <- data.frame(RFID = rownames(network),
                     age  = age.dat[[age.var]][match(rownames(network), age.dat[[id.var]])],
                     stringsAsFactors = F)
  
  inds$age[which(inds$age == "")] <- NA
  inds <- subset(inds, !is.na(age))
  
  net <- network[rownames(network) %in% inds$RFID, 
                 colnames(network) %in% inds$RFID]
  
  rand.age     <- rep(NA, n.perm)

  # set.seed(20)
  
  for(j in 1:n.perm){
    inds$rand <- sample(inds$age)
    rand.age[j] <- assortment.discrete(graph = net, types = inds$rand, weighted = T, SE = F)$r
  }
  rm(j)
  obs.age <- assortment.discrete(graph = net, types = inds$age, weighted = T, SE = T)

  
  results <- data.frame(season        = seas,
                        elevation     = elev,
                        dat.type      = "FullTime Elev-Wide",
                        test.stat     = "assortment coef(r)",
                        ind.att       = "age",
                        hyp.test      = "NodePermONLY",
                        n.juv         = length(inds$age[which(inds$age == "juv")]),
                        n.adult       = length(inds$age[which(inds$age == "adult")]),
                        num.perms     = n.perm,
                        obs.r.val     = obs.age$r,
                        obs.se.val    = obs.age$se,
                        obs.med.r.val = NA,
                        obs.dev.score = NA,
                        quant.05      = quantile(rand.age, probs = 0.05)[[1]],
                        first.qu.rand = summary(rand.age)[[2]],
                        med.coef.rand = median (rand.age),
                        third.qu.rand = summary(rand.age)[[5]],
                        quant.95      = quantile(rand.age, probs = 0.95)[[1]],
                        tail.p.value  = p_value(obs.age$r, rand.age),
                        full.dp.p.val = two.tail.p_value(p_value(obs.age$r, rand.age)),
                        r.sig         = "",
                        r.sig.level   = "",
                        stringsAsFactors = F)

  
  results$r.sig[which(results$full.dp.p.val <= 0.1)]   <- "."
  results$r.sig[which(results$full.dp.p.val <= 0.05)]  <- "*"
  results$r.sig[which(results$full.dp.p.val <= 0.01)]  <- "**"
  results$r.sig[which(results$full.dp.p.val <= 0.001)] <- "***"
  
  results$r.sig.level[which(results$full.dp.p.val <= 0.1)]   <- "0.1"
  results$r.sig.level[which(results$full.dp.p.val <= 0.05)]  <- "0.05"
  results$r.sig.level[which(results$full.dp.p.val <= 0.01)]  <- "0.01"
  results$r.sig.level[which(results$full.dp.p.val <= 0.001)] <- "0.001"
  
  result <- list()
  
  result$np.results    <- results
  result$mixing.matrix <- obs.age$mixing_matrix
  return(result)
  
})

assort.age.results <- do.call('rbind', lapply(assort.age, function(x){
  return(x$np.results)
}))

age.mix.matrices <- lapply(assort.age, function(x){return(x$mixing.matrix)})
names(age.mix.matrices) <- names(full.time.elev.networks)
# save(age.mix.matrices, assort.age.results, file = paste0(dat.path, "AssortElev_AgeAssortResults_2022_1207.RData"))

print(age.mix.matrices)

