#------------------------------ additional sub functions -------------------------------#

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

## To obtain the prior predictive distribution of pj
## given animal data recorded by 
## t0: vector contains the numbers of animals that experienced DLT at the pseudo doses
## n0: vector contains the total number of animals used per pseudo dose
## psdose: the pseudo doses converted from the tested animal doses

probpj <- function(dj, uplim, t0=c(60, 280), n0=300, psdose=c(0.30, 1.00)){
  gd = function(dj, pj){
    sapply(dj, FUN = function(dj){
      sapply(pj, FUN = function(pj){
        integrate(function(E) m0ptheta_base(dj, t0, u0=n0-t0, psdose, E, pj), -Inf, Inf, stop.on.error=F)$value
      })
    })
  }
  
  sapply(dj, FUN=function(dj){
    integrate(function(pj) gd(dj, pj), 0, uplim, stop.on.error=F)$value
  })
}


## To find the key pecentiles, q_ik, defined on page 5 of the manuscript
## taran: the given range, over which the percentiles will be searched
## quant: the target level of percentiles to be found, 
## say, quant = c(0.025, 0.500, 0.975) aims at the 2.5th, 50th and 97.5th percentiles
## dosej: the human dose, at which the risk of toxicity is of interested

quantfind <- function(taran, quant, dosej, t0, n0, psdose){
  p <- NULL
  for(i in taran) p[which(i==taran)] = probpj(dj=dosej, uplim=i, t0, n0, psdose)
  return(list("Fcumu"=p[which.min((p-quant)^2)], "quant"=taran[which.min((p-quant)^2)]))
}

tq <- c(0.025, 0.500, 0.975)


## To obtain a bivariate normal prior for \theta 
## by minimising the defined distance according to the criterion 
## defined on page 5 of the manuscript

ApproxBVN <- function(vparms, qdat, doseSet, dRef){
  mu1 <- vparms[1]
  mu2 <- vparms[2]
  var1 <- vparms[3]^2
  cov12 <- vparms[4]*vparms[3]*vparms[5]
  var2 <- vparms[5]^2
  
  Ez <- mu1 + log(doseSet/dRef)*exp(mu2+1/2*var2)
  Vz <- var1 + 2*log(doseSet/dRef)*exp(mu2+1/2*var2)*cov12 + (log(doseSet))^2*exp(2*mu2+var2)*(exp(var2)-1)
  zL <- Ez - qnorm(1-0.05/2)*sqrt(Vz)
  zU <- Ez + qnorm(1-0.05/2)*sqrt(Vz)
  
  qL = exp(zL)/(1+exp(zL))
  qU = exp(zU)/(1+exp(zU))
  qM = exp(Ez)/(1+exp(Ez))
  
  sum(c((qdat[,1]-qL)^2, (qdat[,2]-qM)^2, (qU-qdat[,3])^2))
}


eps <- .Machine$double.eps  # get a small value for bounding the paramter space to avoid things such as log(0).

# the logistic dose-toxicity model
ddltmod <- function(par, d, dRef){
  pd = sapply(d, FUN=function(d){
    exp(par[1] + exp(par[2])*log(d/dRef))/(1 + exp(par[1] + exp(par[2])*log(d/dRef)))
  })
  return(pd)
}