# nsims: the number of simulation replicates
# here (setting as 10) is tentative, for the code testing purposes only
# set it to be 1000 or more for a complete simulation study

ProcA = function(nsims = 10, MdoseH = 9, Ncohorts = 9, dRef = 28, psDoses = c(2, 54),
                 Myu01 = 0.6, doseH = c(2, 4, 8, 16, 22, 28, 40, 54, 70), CohSize = as.integer(7),
                 pTrue = c(0.11, 0.25, 0.35, 0.41, 0.47, 0.52, 0.58, 0.63, 0.70)){

ntrial = as.integer(0)
safestop = as.integer(0)

PjPME = array(0, dim = c(nsims, MdoseH))

selMTD = rep(0, nsims)

TD = 0.25

SIMaveA = array(0, dim=c(nsims, CohSize))
SIMlbd = array(0, dim=c(nsims, CohSize))
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
  
  UtiDose = vector(mode = "numeric", length(doseH))
  maxU = vector(mode = "numeric", length(doseH))
  
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
    
    predSeq[nCoh] = ifelse(doseCurr<psDoses[1], 0, 
                           ifelse(doseCurr>psDoses[2], 1, 
                                  maximeu(priorpredy0=indivyj0(dj = doseCurr, t0 = c(1, 18), n0 = c(30, 30), psdose = psDoses),
                                          priorpredy1=indivyj1(dj = doseCurr, t0 = c(1, 18), n0 = c(30, 30), psdose = psDoses),
                                          u01 = Myu01)
                           )
    )
    
    indUtility = UtiCal(ypred = rep(predSeq[nCoh], 3), 
                        ytrue = c(numeric(3-toxSeq[nCoh]), rep(1, toxSeq[nCoh])),
                        u01 = Myu01)
    
    UtiDose[index] = UtiDose[index] + sum(indUtility)
    maxU[index] = maxU[index] + 3
    pa = UtiDose/maxU
    pa.constraint = pa[(index-1):length(doseH)]
    aveA = mean(pa.constraint[is.finite(pa.constraint)])
    
    if(length(NsubH[NsubH>0]) == 1){
      predij <- indivyj1(dj = doseH, t0 = c(1, 18), n0 = c(30, 30), psdose = psDoses)
    } else{
      predij <- MCMCSim$mean$pTox.star 
    }
    
    # w.p of p1 to observe a DLT from the 1st patient who will receive the lowest dose
    p1 = predij[index]
    # w.p of p0 to observe a no-DLT from
    p0 = 1 - p1
    # the discrete random variable \bar{a}(eta, y) = x1 when y=1; and x0 when y=0
    x1 = UtiCal(ypred = predSeq[nCoh], ytrue = 1, u01 = Myu01)/1
    x0 = UtiCal(ypred = predSeq[nCoh], ytrue = 0, u01 = Myu01)/1
    xbar = p1*x1+p0*x0
    ## 1. the numerator
    sd_yk = sqrt(p1*x1^2 + p0*x0^2 - xbar^2)
    
    minsd = sd_yrmn(N = CohSize, nk = nCoh, doseSet = doseH, ind = index, ppred = predij, 
                    t0=c(1, 18), n0=c(30, 30), psdose = psDoses, u01 = Myu01)
    
    wMix[1] = aveA^(sd_yk/minsd)
    ## An alternative dynamic weight function
    # wMix[1] = sqrt(CohSize/nCoh)
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
    
    SIMaveA[i, nCoh] <- aveA
    SIMlbd[i, nCoh] <- sd_yk/minsd
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

