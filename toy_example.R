
rm(list = ls())

library(MoBPS) 
library(optiSel) #version 2.0.9. #nadiv version 2.18.0 on R version >= 4.2


################################################################################
#### you can change settings in this section

# number of individuals
pop.size <- 100

# number of generations to simulate
n.generations <- 11

# initial allele frequencies of the alleles we want to change
af.pick <- c(0.95,0.75,0.5,0.25,0.05)

# desired allele frequency 
af.desired <- c(0.1,0.2,0.3,0.4,0.5)
#af.desired <- c(0.05,0.25,0.5,0.75,0.95)
#af.desired <- c(0.8,0.8,0.8,0.8,0.8)
#af.desired <- c(0,0,0,0,0)

# set 'many.qtl' to TRUE to test scenario with 10,50,100 loci
# If FALSE, the lines in the if-statement have no effect
many.qtl <- FALSE
if(many.qtl){
  n.qtl <- 100 # put number of loci you want to change AF at
  af.pick <- runif(n = n.qtl, min = 0.05, max = 0.95) # drawing random starting AF
  af.desired <- rep(0, n.qtl) # desired frequency for all alleles is 0
}

available.ocs.methods <- c("min.rel", "min.abs.dev.index", "min.sqd.dev.index", "min.ZZ")
# pick a method to test by inserting digit. 
# 1=minimizing average relationship, 2=absolute deviation index, 3=squared deviation index, 4=TAFT
ocs.method <- available.ocs.methods[4]

# inbreeding.rate*100 is percentage point increase in average kinship per generation
# inbreeding.rate=0.01 is 1% point increase in average kinship
# and later in the code, the maximum allowed average kinship is determined as 
# average_kinship_now + 0.01 = average_kinship_next_generation
inbreeding.rate <- 0.01

seed <- 1

################################################################################

set.seed(seed = seed)

# start MoBPS simulation
pop <- creating.diploid(nindi = pop.size
                        , nsnp = 50000
                        , chr.nr = 25
                        , chromosome.length = 1)

# this functions find loci who carry alleles that have the same frequency 
# as specified in parameter 'allele.frequency'
find.loci <- function(population, allele.frequency){
  pop <- population
  af.pick.QTL <- allele.frequency
  
  af <- get.allele.freq(pop, gen = get.ngen(pop))
  repeat{
    
    pos <- c()
    snp.i <- 1
    for(snp.i in 1:length(af.pick.QTL)){
      
      repeat{
        
        # sample a locus whose AF is at most 3% away from desired freq
        if(any(af == af.pick.QTL[snp.i])){
          new.pos <- sample(which(af == af.pick.QTL[snp.i]), size = 1)
        } else{
          new.pos <- sample(which(abs(af - af.pick.QTL[snp.i]) < 0.01), size = 1)
        }
        
        if(!any(new.pos %in% pos)){
          pos <- c(pos, new.pos)
          break
        }
      }
    }
    
    map <- get.map(pop)
    
    Morgan <- as.numeric(map[sort(pos),3])
    
    # stop if the loci are at least 1 cM apart from each other
    if(any(!((Morgan[-1] - Morgan[-length(Morgan)]) < 0.01))){
      break
    }
  }
  
  return(pos)
}

pos <- find.loci(population = pop, allele.frequency = af.pick)

pos.org <- pos
pos <- pos.org

# Object in which later statistics from each generation will be stored
record.all <- c()

gen.i <- 1
for(gen.i in 1:n.generations){
  time.gen.start <- Sys.time()
  set.seed(seed = seed)
  
  kin <- kinship.exp(population = pop, gen = get.ngen(pop), depth.pedigree = get.ngen(pop))
  rownames(kin) <- names(get.id(pop, gen = get.ngen(pop)))
  colnames(kin) <- rownames(kin)
  rel <- kin*2
  av.kin.now <- mean(kin)
  
  av.kin.next <- av.kin.now + inbreeding.rate
  
  av.rel.now <- av.kin.now*2
  # the average relationship we want to achieve in the next generation
  av.rel.next <- av.kin.next*2
  
  af <- get.allele.freq(pop, gen = get.ngen(pop))
  geno <- get.geno(pop, gen = get.ngen(pop))
  
  M <- t(geno[pos,])
  P <- matrix(af.desired*2, nrow = nrow(M), ncol = length(af.desired), byrow = TRUE)
  Z <- M - P
  
  # this stores index values of sum of absolute deviations of allele frequencies 
  # from desired allele frequencies
  abs.dev.index <- rowSums(abs(Z))/2
  # This is identical to: colSums(abs((geno[pos,]/2) - af.desired))
  
  
  # this stores index values of sum of squared deviations of allele frequencies 
  # from desired allele frequencies
  sqd.dev.index <- rowSums((Z/2)^2)
  
  ZZ <- Z %*% t(Z)
  
  # this stores the matrix whose values are the product of the deviations of allele frequencies 
  # from desired allele frequencies
  ZZ.4 <- ZZ/4
  
  # The diagonal of ZZ/4 is identical to the sum of squared deviations
  # identical(round(diag(ZZ/4), 5), round(ZZ.4,5))
  
  # average element of ZZ is identical to 4 times sum of squared deviations 
  # of observed frequency from centering frequency
  #mean(ZZ)
  #4*sum((af[pos] - af.desired)^2)
  
  # division by number of alleles to avoid extremely large values when many
  # alleles are used. This does not change the solutions of OCS. This just 
  # makes it easier to solve in some cases
  ZZ.4 <- ZZ.4/length(af.pick)
  
  ##############################################################################
  # prepare info for optiSel
  
  sex <- get.sex(pop, gen = get.ngen(pop))
  sex[sex==1] <- "male"
  sex[sex==2] <- "female"
  
  isCandidate <- rep(TRUE, times = nrow(rel))
  id <- get.id(pop, gen = get.ngen(pop))
  id <- names(id)
  
  phen <- data.table::data.table(Indiv = id, 
                                 Breed = "BreedA", 
                                 Sex = sex,
                                 
                                 abs.dev.index = abs.dev.index,
                                 sqd.dev.index = sqd.dev.index,
                                 
                                 isCandidate = isCandidate)
  
  snp.k <- 1
  for(snp.k in 1:length(af.desired)){
    phen <- cbind(phen, (geno[pos[snp.k],]/2))
    colnames(phen)[ncol(phen)] <- paste0("AF",snp.k)
  }
  
  cand <- candes(phen = phen
                 , rel = rel
                 , ZZ = ZZ.4
                 , quiet = FALSE)
  
  
  constraints <- list()
  
  # this function adds constraints for the allele frequencies
  # to keep them between 4% and 96% to avoid accidental loss or fixation
  # These constraints are only added if the allele should not be purged or fixed, respectively.
  constraints.af <- function(constraints){
    
    af <- get.allele.freq(pop, gen = get.ngen(pop))
    geno <- get.geno(pop, gen = get.ngen(pop))
    
    con <- list()
    snp.k <- 33
    for(snp.k in 1:length(af.desired)){
      
      if((rowMeans(geno)/2)[pos[snp.k]] %in% c(0,1)){
        # if the allele has been lost or fixed, skip
        next
      }
      
      # all alleles should have frequencies between 4% and 96% to avoid unwanted
      # loss or fixation 
      ub.AF <- 0.96
      lb.AF <- 0.04
      
      # if the desired frequency is higher or lower than that (98% or 2%), the boundary needs 
      # to be change to allow for intentional fixation or purging
      if(af.desired[snp.k] > 0.98){
        ub.AF <- 1.1
      }
      
      if(af.desired[snp.k] < 0.02){
        lb.AF <- -0.1
      }
      
      # optiSel does not allow to specify upper and lower bounds at the same time
      # So we only set that restriction that is more relevant, i.e., closer to the
      # currently observed frequency
      closer.to.1 <- (af[pos[snp.k]]-0.5000001) > 0
      
      if(closer.to.1){
        con[length(con)+1] <- list(ub.AF = ub.AF)
        names(con)[length(con)] <- paste0("ub.AF",snp.k)
      } else{
        con[length(con)+1] <- list(lb.AF = lb.AF)
        names(con)[length(con)] <- paste0("lb.AF",snp.k)
      }
    }
    
    # appending list
    if(length(con) > 0){
      for(i in 1:length(con)){
        constraints[length(constraints)+1] <- con[i]
        names(constraints)[length(constraints)] <- names(con)[i]
      }
    }
    
    return(constraints)
  }
  
  constraints <- constraints.af(constraints = constraints)

  # Add constraint for increase of average kinship.
  constraints[length(constraints)+1] <- list(ub.rel = av.rel.next)
  names(constraints)[length(constraints)] <- "ub.rel"
  
  ##############################################################################
  # TAFT-OCS
  
  # Used solvers for optiSel. cccp2 sometimes gave unreliable results and is thus omitted.
  vec.solvers <- c("slsqp", "cccp", "alabama"
                   #, "cccp2"
  )
  
  # Looping over solvers
  # The solution that is valid and at boundary is taken. So many solvers are 
  # tried until a valid solution is obtained, ideally at the boundary of the solution space.
  # If none of the solvers can come up with such a solution, we use the one from slsqp
  ocs.valid.vec <- c()
  ocs.boundary.vec <- c()
  contr.list <- list()
  ok.solution.found <- FALSE
  solver.i <- 1
  for(solver.i in 1:length(vec.solvers)){
    
    time.solver.start <- Sys.time()
    
    try(
      contr <- optiSel::opticont(method = ocs.method, cand = cand, con = constraints
                                 , quiet = FALSE, solver = vec.solvers[solver.i])
    )
    
    # is the solution valid
    ocs.valid <- contr$info$valid
    
    # is the solution at the boundary of the solution space?
    # This is checked as the absolute difference of the value of the obtained solution
    # from the value used as restriction. If that difference is smaller than 0.0001,
    # it is assumed that the solution is at the boundary of the solution space
    # of the respective dimension.
    # The entire solution is only declared to be a boundary solution if all dimensions
    # are at the boundary.
    ocs.boundary <- !any(!(abs( round(contr$summary$Val - contr$summary$Bound,6) ) < 0.0001), na.rm = TRUE)
    
    contr.list[[solver.i]] <- contr
    ocs.valid.vec <- c(ocs.valid.vec, ocs.valid)
    ocs.boundary.vec <- c(ocs.boundary.vec, ocs.boundary)
    
    time.solver.end <- Sys.time()
    cat(paste0("Solver ", vec.solvers[solver.i]," needed: ",
               as.numeric(round(difftime(time1 = time.solver.end, time2 = time.solver.start, units = "mins"),2)),
               " minutes.\n"
    ))
    
    if(ocs.valid & ocs.boundary){
      ok.solution.found <- TRUE
      break
    }
  }
  
  # If no valid solution at the boundary of the solution space was found, then find the next best one.
  # (i.e. one that is at least valid)
  if(!ok.solution.found){
    if(any(ocs.valid.vec)){
      # select a solution that is at least valid
      sol.pos <- which(ocs.valid.vec)[1]
    } else{
      # if none is valid, use the one from slsqp
      sol.pos <- 1
    }
    
    ocs.valid <- ocs.valid.vec[sol.pos]
    ocs.boundary <- ocs.boundary.vec[sol.pos]
    contr <- contr.list[[sol.pos]]
  }
  
  contr.male <- contr$parent$oc[contr$parent$Sex == "male"]
  contr.female <- contr$parent$oc[contr$parent$Sex == "female"]
  
  ##############################################################################
  
  # Make mate plan for next generation.
  mate.plan <- c()
  for(ind.i in 1:pop.size){
    
    sex <- get.sex(pop, gen = get.ngen(pop))
    id <- get.id(pop, gen = get.ngen(pop))
    id.m <- id[sex == 1]
    id.f <- id[sex == 2]
    
    father <- sample(id.m, size = 1, replace = TRUE, prob = contr.male)
    mother <- sample(id.f, size = 1, replace = TRUE, prob = contr.female)
    
    mate.plan <- rbind(mate.plan,
                       cbind(
                         get.individual.loc(pop, database = get.database(pop, id = father)),
                         get.individual.loc(pop, database = get.database(pop, id = mother))
                       )
    )
  }
  mate.plan <- cbind(mate.plan, sex = c(rep(0, pop.size/2), rep(1, pop.size/2)))
  
  ##############################################################################
  # function to record statistics
  get.stats.QTL.OCS <- function(population){
    pop <- population
    af <- get.allele.freq(pop, gen = get.ngen(pop))
    
    if(exists("contr")){
      valid.ocs <- contr$info$valid
      ocs.boundary <- !any(!(abs( round(contr$summary$Val - contr$summary$Bound,6) ) < 0.0001), na.rm = TRUE)
    } else{
      valid.ocs <- TRUE
      ocs.boundary <- TRUE
    }
    
    n.sel.parents <- sum(round(contr$parent$oc,6) != 0)
    
    af1 <- af[pos[1]]
    af2 <- af[pos[2]]
    af3 <- af[pos[3]]
    af4 <- af[pos[4]]
    af5 <- af[pos[5]]
    
    # Expected average heterozygosity at loci that are to be changed.
    He.qtl <- mean(2*af[pos]*(1-af[pos]))
    frac.qtl.unfavorable.fixed <- mean(af[pos] %in% 1)
    
    af.vec <- c()
    for(snp.k in 1:length(pos)){
      af.vec <- c(af.vec, af[pos[snp.k]])
    }
    
    # I am checking if any allele that show be at 0% is at 100% and thus making 
    # it impossible to achieve the goal
    goal.achievable = !(any((af.vec[af.desired == 0] == 1)) | any((af.vec[af.desired == 1] == 0)))
    
    geno.subset <- get.geno(pop, gen = get.ngen(pop))[pos,]
    # number of individuals that are not carrier for any of the unwanted alleles
    n.free.QTL <- sum(colSums(geno.subset) == 0)
    # number of individuals that are not homozygous for at least one of the unwanted alleles
    n.no.homo.QTL <- sum(colSums(geno.subset == 2) == 0)
    
    ibd.inbreeding <- inbreeding.emp(pop, gen = get.ngen(pop))
    
    pedigree.inbreeding <- diag((kinship.exp(pop, gen = get.ngen(pop), depth.pedigree = get.ngen(pop))*2))-1
    
    ##############################################################################
    
    if(get.ngen(pop) != 1){
      af.prev <- get.allele.freq(pop, gen = get.ngen(pop)-1)
      af.now <- get.allele.freq(pop, gen = get.ngen(pop))
      
      af.change <- af.now[pos] - af.prev[pos]
      
      geno <- get.geno(pop, gen = get.ngen(pop)-1)
      geno.qtl <- geno[pos,]
      
      var.P <- function(vector){
        average <- mean(vector)
        out.var <- sum((vector - average)**2)/length(vector)
        return(out.var)
      }
      
      # var.af.obs is the variance of the allele frequency among individuals.
      # Under Hardy-Weinberg-equilibirum, it is identical to 0.5*p*(1-p)
      var.af.obs <- apply(geno.qtl/2, MARGIN = 1, var.P)
      af.sel.intensity.obs <- af.change / sqrt( var.af.obs )
      
    } else{
      af.change <- rep(NA, times = length(pos))
      af.sel.intensity.obs <- rep(NA, times = length(pos))
    }
    
    ##############################################################################
    
    record.df <- data.frame(
      n.free.QTL = n.free.QTL,
      n.no.homo.QTL = n.no.homo.QTL,
      
      pedigree.inbreeding = mean(round(pedigree.inbreeding,4)),
      ibd.inbreeding = mean(round(ibd.inbreeding,4)),
      
      max.av.rel = ifelse(exists("av.rel.next"),av.rel.next,1/(2*ncol(geno.subset))),
      
      af1 = af1,
      af2 = af2,
      af3 = af3,
      af4 = af4,
      af5 = af5,
      
      ocs.valid = valid.ocs,
      ocs.boundary = ocs.boundary,
      n.sel.parents = n.sel.parents,
      
      av.af = mean(af.vec),
      goal.achievable = goal.achievable,
      
      af1.change = af.change[1],
      af2.change = af.change[2],
      af3.change = af.change[3],
      af4.change = af.change[4],
      af5.change = af.change[5],
      
      af1.sel.intensity.obs = af.sel.intensity.obs[1],
      af2.sel.intensity.obs = af.sel.intensity.obs[2],
      af3.sel.intensity.obs = af.sel.intensity.obs[3],
      af4.sel.intensity.obs = af.sel.intensity.obs[4],
      af5.sel.intensity.obs = af.sel.intensity.obs[5],
      
      He.qtl = He.qtl,
      frac.qtl.unfavorable.fixed = frac.qtl.unfavorable.fixed
      
    )
    rownames(record.df) <- get.ngen(pop)
    
    return(record.df)
  }
  
  record.gen <- get.stats.QTL.OCS(population = pop)
  record.all <- rbind(record.all, record.gen)
  ##############################################################################
  
  # breed next generation according to mate plan
  pop <- breeding.diploid(pop, fixed.breeding = mate.plan
                          , delete.same.origin = TRUE
                          , recombination.minimum.distance = 0.1)
  
  time.gen.end <- Sys.time()
  cat(paste0("Needed: ",
             as.numeric(round(difftime(time1 = time.gen.end, time2 = time.gen.start, units = "mins"),2)),
             " minutes for this generation.\n"
  ))
}

record.all

colosr.all <- 1:5
i <- 1
for(i in 1:5){
  nm.col <- paste0("af",i)
  
  if(i == 1){
    plot(record.all[,nm.col], type = "l", col = colosr.all[i]
         ,xlab = "Generation"
         , xlim = c(1,11)
         , ylab = "Allele frequency"
         , ylim = c(0,1)
         , xaxt = "n"
         , main = paste0("method = ", ocs.method)
    )
    
    u <- par("usr") # The coordinates of the plot area
    rect(u[1], u[3], u[2], u[4], col="gainsboro")
  } 
  
  lines(record.all[,nm.col], col = colosr.all[i], lwd = 2)
  
  if(length(unique(af.desired)) == 1){
    abline(h = af.desired[i], lty = "dashed", col = "black", lwd = 1)
  } else{
    abline(h = af.desired[i], lty = "dashed", col = colosr.all[i], lwd = 1)
  }

}
axis(side = 1, tick = TRUE, at = c(1,3,5,7,9,11), labels = c(1,3,5,7,9,11)-1)

