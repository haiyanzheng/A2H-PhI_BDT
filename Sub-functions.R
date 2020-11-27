#------------------------------ sub functions -------------------------------#

#----------------------------------------------------------------------
## Author: Haiyan Zheng
## Created: Mar 15 2018 (08:47) 
## Version: 
## Last-Updated: Feb 19 2020 (14:10) 
##           By: Haiyan Zheng
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

# m0ptheta_base(dj, t0, u0, psdose, E, pj) is the joint pdf of the human toxicity risk, p_j,
# and the slope parameter, theta_2.
# dj: the human dose, at which we want to predict the toxicity
# t0: vector contains the numbers of animals that experienced DLT at the pseudo doses
# u0: vector contains the numbers of animals that did not experience DLT at the pseudo doses
# psdose: pseudo doses converted from the tested animal doses
# E: exp(theta_2) which will be integrated out 

m0ptheta_base <- function(dj, t0, u0, psdose, E, pj){
  prd.low = (exp(E*log(dj/psdose[1])
                 - log(pj/(1-pj))) + 1)^(-t0[1])*(exp(E*log(psdose[1]/dj)
                                                      + log(pj/(1-pj))) + 1)^(-u0[1])/beta(t0[1], u0[1])
  prd.high = (exp(E*log(dj/psdose[2])
                  - log(pj/(1-pj))) + 1)^(-t0[2])*(exp(E*log(psdose[2]/dj)
                                                       + log(pj/(1-pj))) + 1)^(-u0[2])/beta(t0[2], u0[2])
  return(prd.low*prd.high*abs(log(psdose[1]/psdose[2]))*1/(pj*(1-pj)))
}


## prior predictive probability that DLT=1
indivyj1 <- function(dj, t0, n0, psdose){
  sapply(dj, FUN = function(dj){
    integrate(function(pj){
      sapply(pj, FUN = function(pj){
        integrate(function(E) pj*m0ptheta_base(dj, t0, u0=n0-t0, psdose, E, pj), -Inf, Inf, stop.on.error=F)$value
      })
    }, 0, 1)$value
  })
}

indivyj1(dj=70, t0=c(1, 18), n0=c(30, 30), psdose=c(2, 54))

## prior predictive probability that DLT=0
indivyj0 <- function(dj, t0, n0, psdose){
  sapply(dj, FUN = function(dj){
    integrate(function(pj){
      sapply(pj, FUN = function(pj){
        integrate(function(E) (1-pj)*m0ptheta_base(dj, t0, u0=n0-t0, psdose, E, pj), -Inf, Inf, stop.on.error=F)$value
      })
    }, 0, 1)$value
  })
}

indivyj0(dj=c(28,40), t0=c(1, 18), n0=c(30, 30), psdose=c(2, 54))

## Maximising the expected utility 
# this leads to the optimal prior prediction of human toxicity outcome,
# given the utility values u11, u10, u01 and u00
maximeu <- function(priorpredy0, priorpredy1,
                    u11=1, u10=0, u01=0.6, u00=1){
  # u_{ji} to denote \tilde{y}=j predicted as \eta=i
  exputi = rbind(u00*priorpredy0, u10*priorpredy1, u01*priorpredy0, u11*priorpredy1)
  eta = c(0, 0, 1, 1)
  index = rep(0, length(unique(eta)))
  for(i in 1:length(index))  index[i]=mean(exputi[eta==unique(eta)[i]])
  return(unique(eta)[which.max(index)])
}


maximeu(priorpredy0 = indivyj0(dj=22, t0=c(1, 17), n0=c(30, 30), psdose=c(2, 54)),
        priorpredy1 = indivyj1(dj=22, t0=c(1, 17), n0=c(30, 30), psdose=c(2, 54)))

# To penalise (reward) an incorrect (correct) prediction of the human toxicity outcome
UtiCal <- function(ytrue, ypred, u11=1, u10=0, u01=0.6, u00=1){
  # u_{ji} to denote ytrue=j predicted as \eta=i
  utility = rep(0, length(ypred))
  utility[ytrue==1 & ypred==1] <- u11
  utility[ytrue>ypred] <- u10 
  utility[ytrue<ypred] <- u01
  utility[ypred==0 & ytrue==0] <- u00
  return(utility)
}


## tuning parameter \lambda^{(h)} defined in Model (14) of the manuscript
sd_yrmn <- function(N, nk, doseSet, ind, ppred,
                    t0, n0, psdose, u11=1, u10=0, u01=0.6, u00=1, 
                    sdsims=5000){
  dsel = doseSet[ind]
  nres = N-nk+1
  predprob = rep(ppred[ind], nres)
  optpredy = rep(ifelse(dsel<psdose[1], 0,
                        ifelse(dsel>psdose[2], 1,
                               maximeu(priorpredy0 = indivyj0(dj=dsel, t0, n0, psdose),
                                       priorpredy1 = indivyj1(dj=dsel, t0, n0, psdose),
                                       u11, u10, u01, u00))),
                 nres)
  ## calculate the sd(simaveA)
  simaveA = rep(0, sdsims)
  for(j in 1:sdsims){
    simOut <- sapply(predprob, rbinom, n=1, size=1)
    simpa <- UtiCal(ypred=optpredy, ytrue=simOut)/rep(1, length(simOut))
    simaveA[j] <- mean(simpa)
  }
  return(sd(simaveA))  
}

#------------------------------------------------------------------------#
