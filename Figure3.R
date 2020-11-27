## Load all other user-defined R functions
source("Sub-functions.R")
source("Approx-functions.R")

## Load the required R package(s)
require(mvtnorm)
require(R2OpenBUGS)
require(ggplot2)
require(coda)
require(lattice)

set.seed(1)
## hypothetical animal data ##
#-----------------------------------------------------------------------------#
#  beta prior: dose0 = 2 - pj ~ Beta(1, 29); dose-1 = 54 - pj ~ Beta(17, 13)  #
#-----------------------------------------------------------------------------#
doseH = c(2, 4, 8, 16, 22, 28, 40, 54, 70)

# probpj(dj=doseH[9], uplim=0.75, t0=c(1, 17), n0=30, psdose=c(2, 54))

qt0025 <- array(0, dim=c(length(doseH), 2))
#taran=seq(0.0008, 0.48, length=100)
qt0500 <- array(0, dim=c(length(doseH), 2))
#taran=seq(0.02, 0.65, length=100)
qt0975 <- array(0, dim=c(length(doseH), 2))
#taran=seq(0.1, 0.85, length=100)

for(j in 1:length(doseH)){
         qt0025[j, ] <- unlist(quantfind(taran=seq(0.0008, 0.48, length=1000), quant=tq[1], 
                                         dosej=doseH[j], t0=c(1, 17), n0=30, psdose=c(2, 54)))
         
         qt0500[j, ] <- unlist(quantfind(taran=seq(0.02, 0.65, length=1000), quant=tq[2], 
                                         dosej=doseH[j], t0=c(1, 17), n0=30, psdose=c(2, 54)))
         
         qt0975[j, ] <- unlist(quantfind(taran=seq(0.1, 0.85, length=1000), quant=tq[3], 
                                         dosej=doseH[j], t0=c(1, 17), n0=30, psdose=c(2, 54)))
         
  print(j)
}


qts1aniPrior <- cbind(qt0025[,2], qt0500[,2], qt0975[,2])

## approximate the marginal priors on pi with a bivariate normal prior of \theta
# Animal Prior 1
bvnopt <- optim(c(1, 0.7, 1, -0.1, 0.1), ApproxBVN, 
                qdat = qts1aniPrior, doseSet = doseH, dRef = 28,
                method = "L-BFGS-B", 
                lower = c(-Inf, -Inf, -Inf, -0.6, -Inf), 
                upper = c(Inf, Inf, Inf, 0.6, Inf))

var.theta1 <- bvnopt$par[3]^2
cov.theta12 <- bvnopt$par[4]*bvnopt$par[3]*bvnopt$par[5]
var.theta2 <- bvnopt$par[5]^2

myparms <- c(bvnopt$par[1], bvnopt$par[2], var.theta1, cov.theta12, var.theta2)
round(myparms, 3)  

# Bivariate normal prior 
mmatrix <- c(-0.524, 0.147)
covmatrix <- matrix(c(0.151, -0.008, -0.008, 0.001), 2, 2)


alphabeta <- rmvnorm(50000, mmatrix, covmatrix)
pd = array(0, dim=c(nrow(alphabeta), length(doseH)))
for(i in 1:nrow(alphabeta)) pd[i,] = ddltmod(par = alphabeta[i,], d = doseH, dRef = 28)

# the implied percentiles
qts = array(rep(0, length(doseH)*3), dim=c(length(doseH), 3))
for(j in 1:length(doseH)) qts[j,] = unname(quantile(pd[,j], probs=c(0.025, 0.500, 0.975), na.rm = T))

##--------------------- Create Figure 3A ---------------------##
predToxS1 <- data.frame(dose= 1:length(doseH), V = qts1aniPrior*100, AM1 = "A")
predToxSF <- data.frame(X = qts*100, AM2 = "F")

predToxSc <- as.data.frame(cbind(predToxS1, predToxSF))

p <- ggplot(predToxSc, aes(x = dose)) + xlab(bquote("Dose ("*mg/m^2*")")) + 
  ylab("Probability of toxicity (%)") + scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = unique(predToxSc$dose), labels = doseH, limits = c(0.75, 9.25)) + theme_bw() 


aPrior <- p + geom_linerange(aes(x = dose, ymin = V.1, ymax = V.3, colour = AM1),
                             size = 1, position = position_dodge(width = 0.65)) +
  geom_point(aes(y =  V.2, colour = AM1), 
             position= position_dodge(width = 0.65), size = 1.5) + 
  theme(legend.position='top', legend.title = element_blank(), panel.grid.minor =  element_blank()) +
  scale_colour_manual(name = element_blank(), labels = "Animal priors",
                      values = "#56B4E9") +
  scale_shape_manual(name = element_blank(), labels = "Animal priors",
                     values = 16) +
  geom_line(aes(y=X.2), colour = "#CC79A7", linetype = "dashed") +
  geom_line(aes(y=X.1), colour = "#CC79A7", linetype = "dashed") +
  geom_line(aes(y=X.3), colour = "#CC79A7", linetype = "dashed") 

aPrior

#---------------------------------------------------------------------------#
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

parameters <- c("pTox.star", "pCat")

MCMCSim <- bugs(data, inits, parameters, "MixPrior_PhI.txt", codaPkg = F, bugs.seed = 1, 
                n.chains = 1, n.burnin = 5000, n.iter = 15000)

##--------------------- Create Figure 3B ---------------------##
# interval probabilities
interim <- data.frame(dose = rep(1:9, 3), pTox = rep(MCMCSim$summary[1:9, 5], 3),
                      pToxL = rep(MCMCSim$summary[1:9, 3], 3),
                      pToxU = rep(MCMCSim$summary[1:9, 7], 3),
                      certainty = c(MCMCSim$mean$pCat[1:9,1], MCMCSim$mean$pCat[1:9,2], MCMCSim$mean$pCat[1:9,3]),
                      type = c(rep("Underdose", 9),
                               rep("Target interval", 9),
                               rep("Overdose", 9)))

pb <- ggplot(interim, aes(x=dose, y=pTox*100)) + scale_x_continuous(breaks=1:9, labels=doseH[1:9]) +
 geom_line(aes(y=pTox*100), colour = "#CC79A7") + 
  geom_point(shape=21, fill="white", colour="#CC79A7") +
  geom_hline(yintercept = 25, linetype = "dashed", colour = "lightgray") +
  xlab(bquote("Dose ("*mg/m^2*")")) + ylab("Prior interval probability (%)")  + ylim(0, 100) +
  theme_bw() + theme(legend.position='top', legend.title = element_blank())

intprob <- pb + geom_bar(aes(y = certainty*100, fill = type), stat = "identity", position = "dodge", width = 0.8) +
  scale_fill_grey(start = 0.6, end = .9) 

intprob 

##--------------------- Create Figure 3C ---------------------##
parameters <- c("pTox.star")

MCMCSim <- bugs(data, inits, parameters, "MixPrior_PhI.txt", codaPkg = T, bugs.seed = 1, 
                n.chains = 1, n.burnin = 5000, n.iter = 15000)

MCMCSim.detail <- read.bugs(MCMCSim)

# densityplot(MCMCSim.detail[, 2:4])

pTox1dat <- data.frame(pTox = round(MCMCSim.detail[[1]][,2], 5)*100, Chain = "1st", 
                        dose = "Dose 2", type = "Robust borrowing")

pTox2dat <- data.frame(pTox = round(MCMCSim.detail[[1]][,3], 5)*100, Chain = "1st", 
                       dose = "Dose 4", type = "Robust borrowing")

pTox3dat <- data.frame(pTox = round(MCMCSim.detail[[1]][,4], 5)*100, Chain = "1st", 
                       dose = "Dose 8", type = "Robust borrowing")

pTox1stdat <- rbind(pTox1dat, pTox2dat, pTox3dat)

m1 <- ggplot(pTox1stdat, aes(x = var1, colour = dose, linetype = dose))

denplot1 <- m1 + geom_line(stat="density", size = 0.8) + labs(x = "Probability of toxicity (%)", 
                                                              y = "Prior density") +
  xlim(0, 50) + ylim(0, 0.3) + theme_bw() + theme(legend.position='top') +
  scale_colour_manual(name = element_blank(),
                      labels = c(expression(paste("2 ", mg/m^2)), 
                                 expression(paste("4 ", mg/m^2)),
                                 expression(paste("8 ", mg/m^2))),
                      values = c("#009E73", "#0072B2", "#D55E00")) + 
  scale_linetype_manual(name = element_blank(),
                        labels = c(expression(paste("2 ", mg/m^2)), 
                                   expression(paste("4 ", mg/m^2)),
                                   expression(paste("8 ", mg/m^2))), values=c(2, 1, 4)) 

denplot1


# save(list = ls(), file = "Animal prior (Figure 3).rds")