# nsims: the number of simulation replicates
# here (setting as 10) is tentative, for the code testing purposes only
# set it to be 1000 or more for a complete simulation study

ProcFixW = function(nsims = 10, MdoseH = 9, Ncohorts = 9, dRef = 28, psDoses = c(2, 54),
                 doseH = c(2, 4, 8, 16, 22, 28, 40, 54, 70), CohSize = as.integer(7),
                 MyMixW = 0.5, pTrue = c(0.11, 0.25, 0.35, 0.41, 0.47, 0.52, 0.58, 0.63, 0.70)){

ntrial = as.integer(0)
safestop = as.integer(0)

PjPME = array(0, dim = c(nsims, MdoseH))

selMTD = rep(0, nsims)

TD = 0.25

SIMw = array(0, dim=c(nsims, CohSize))

pMTDhat = array(0, dim = c(nsims, CohSize))
# simTrialInfo = array(0, dim = c(nsims*CohSize, 8))

probEx = array(0, dim = c(nsims, CohSize))

NsubHsum = vector(mode = "numeric", length(doseH))
NtoxHsum = vector(mode = "numeric", length(doseH))

for(i in 1:nsims){
  # initialisation
  dSeq = rep(0, CohSize)
  toxSeq = rep(0, CohSize)
  predSeq = rep(0, CohSize)
  pSeq = rep(0, CohSize)
  p3int = array(0, dim = c(CohSize, 3))
  
  dmax = numeric(length(dSeq))
  dsel = numeric(length(dSeq))
  
  NsubH = c(0, 0, 0, 0, 0, 0, 0, 0, 0)
  NtoxH = c(0, 0, 0, 0, 0, 0, 0, 0, 0)
  
  
  wMix = numeric(2)
  
  nCoh = as.integer(1)
  # force the procedure to begin with the lowest dose
  doseCurr = doseH[2]
  
  while(nCoh <= CohSize){
    dSeq[nCoh] <- doseCurr
    index = match(dSeq[nCoh], doseH)
    # default number of patient as 3 per cohort
    toxSeq[nCoh] <- rbinom(1, size = 3, prob = pTrue[index])
    NtoxH[index] <- NtoxH[index] + toxSeq[nCoh]
    NsubH[index] <- NsubH[index] + 3
    
    wMix[1] = MyMixW
    wMix[2] = 1 - wMix[1]
    
    MCMCSim <- bugs(data, inits, parameters, "MixPrior_PhI.txt", codaPkg = F, bugs.seed = 1, 
                    n.chains = 1, n.burnin = 5000, n.iter = 15000)
    
    if(MCMCSim$mean$pCat[1,3] > 0.25){
      safestop = safestop +1
      break
    }else{
      safestop = safestop +0
      dmax[nCoh] <- max(doseH[MCMCSim$mean$pCat[,3] <= 0.25])
    }
    
    dsel[nCoh] <- min(ifelse(index==9, 70, 
                             ifelse(index==8 | doseH[index+2]/doseH[index] > 2, doseH[index+1],
                                    doseH[index+2])), dmax[nCoh])
    
    probEx[i, nCoh] <- MCMCSim$mean$prob.ex[1]
    
    pMTDhat[i, nCoh] <- MCMCSim$mean$pTox.star[match(TD, pTrue)]
    
    pSeq[nCoh] <- MCMCSim$mean$pTox.star[match(dsel[nCoh], doseH)]
    p3int[nCoh,] <- MCMCSim$mean$pCat[match(dsel[nCoh], doseH),]
    SIMw[i, nCoh] <- wMix[1]
    
    
    doseCurr <- dsel[nCoh]
    nCoh = nCoh + 1
    
  }
  
  for(j in 1:length(doseH)){
    NsubHsum[j] = NsubHsum[j] + NsubH[j]
    NtoxHsum[j] = NtoxHsum[j] + NtoxH[j]
  }
  
  nCoh = nCoh - 1
  ntrial = ntrial + nCoh
  
  # simTrialInfo[((i-1)*CohSize+1):(i*CohSize),] = cbind(i, 1:CohSize, dSeq, toxSeq, pSeq, p3int)
  
  # inconclusive trials will not result in the estimation
  if(nCoh < CohSize)  {
    next
  }
  
  # save the results of the final Bayesian analysis
  PjPME[i,] <- MCMCSim$mean$pTox.star
  selMTD[i] <- min(ifelse(PjPME[i, length(doseH)] == 0, 0, doseH[which.min(abs(PjPME[i,]-TD))]),
                   max(doseH[MCMCSim$mean$pCat[,3] <= 0.25]))
  
}
  return(list(selMTD, data.frame(NsubHsum, NtoxHsum)))

}

