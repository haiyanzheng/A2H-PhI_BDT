## Load all other sub-functions
source("Sub-functions.R")

## Load the required R package(s)
require(R2OpenBUGS)


#--- Approximating the marginal priors by beta distributions ---#
estBetaParams <- function(mu, var){
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(c(alpha, beta))
}
#---------------------------------------------------------------#

MdoseH = 9
Ncohorts = 9
doseH = c(2, 4, 8, 16, 22, 28, 40, 54, 70)
dRef = 28
psDoses = c(2, 54)

Prior.mw = c(-1.099, 0)
Prior.sw = c(2, 1)
Prior.corr = 0

pTox.cut = c(0.16, 0.33)    

PriorA = c(-0.524, 0.147, 0.151, -0.008, 0.001)

NsubH = c(0, 0, 0, 0, 0, 0, 0, 0, 0)
NtoxH = c(0, 0, 0, 0, 0, 0, 0, 0, 0)

wMix = c(1, 0)

data <- list("MdoseH", "Ncohorts", "doseH", "dRef", "NtoxH", "NsubH", "PriorA", "wMix",
             "Prior.mw", "Prior.sw", "Prior.corr", "pTox.cut")

inits <- function(){
  list(
    mu = c(-0.0142, -0.0919)
  )
}

parameters <- c("pTox.star", "pCat", "prob.ex")

MCMCSim <- bugs(data, inits, parameters, "MixPrior_PhI.txt", codaPkg = F, bugs.seed = 1, 
                n.chains = 1, n.burnin = 5000, n.iter = 15000)


##  matching the first two moments 
estbp = array(0, dim = c(9, 2))
moments <- MCMCSim$summary[1:9,1:2]

for(i in 1:9) estbp[i,] <- estBetaParams(mu=moments[i,1], var=moments[i,2]^2)


##--------------------- Create Table 2 ---------------------##
data.frame(PriorMean = round(moments, 3)[,1], PriorStd = round(moments, 3)[,2])

data.frame(ESS = round(estbp, 1)[,1] + round(estbp, 1)[,2], 
           a = round(estbp, 1)[,1], 
           b = round(estbp, 1)[,2])

