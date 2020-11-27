## Load all other user-defined R functions
source("Sub-functions.R")

## Load the required R package(s)
require(R2OpenBUGS)
require(reshape2)
require(gtools)
require(ggplot2)

set.seed(684324)
#-------------------------------------------------------------#
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


data <- list("MdoseH", "Ncohorts", "doseH", "dRef", "NtoxH", "NsubH", "PriorA", "wMix",
             "Prior.mw", "Prior.sw", "Prior.corr", "pTox.cut")


inits <- function(){
  list(
    mu = c(-0.0142, -0.0919)
  )
}

parameters <- c("pTox.star", "pCat", "prob.ex")

# u_01 = 0.6
#---------------------------------- Sequential trial -----------------------------------#

NoRunIn = function(dSeq, toxSeq){

  nCoh = as.integer(1)
  CohSize = 11
  safestop = as.integer(0)
    
  NsubH = c(0, 0, 0, 0, 0, 0, 0, 0, 0)
  NtoxH = c(0, 0, 0, 0, 0, 0, 0, 0, 0)
  
dmax = numeric(length(dSeq))
dsel = numeric(length(dSeq))
predSeq = numeric(length(dSeq))
simpTox = numeric(length(dSeq))
predAc = numeric(length(dSeq))
lambda = numeric(length(dSeq))
stage = numeric(length(dSeq))

wMix = numeric(2)

UtiDose = vector(mode = "numeric", length(doseH))
maxU = vector(mode = "numeric", length(doseH))

probUD = array(0, dim = c(length(dSeq), length(doseH)))
probTI = array(0, dim = c(length(dSeq), length(doseH)))
probOD = array(0, dim = c(length(dSeq), length(doseH)))

priorEx = array(0, dim = c(length(dSeq), 2))
postEx = array(0, dim = c(length(dSeq), 2))

intprob = array(0, dim = c(length(dSeq), 3))

estbp = array(0, c(length(doseH),2))
EHSS = array(0, dim = c(length(dSeq), length(doseH)))

while(nCoh <= length(dSeq)){
  
  doseCurr = dSeq[nCoh]
  index = match(dSeq[nCoh], doseH)
  # default number of patient as 3 each administration
  NsubH[index] <- NsubH[index] + 3
  NtoxH[index] <- NtoxH[index] + toxSeq[nCoh]
  
  predSeq[nCoh] = ifelse(doseCurr<psDoses[1], 0, 
                         ifelse(doseCurr>psDoses[2], 1, 
                                maximeu(priorpredy0=indivyj0(dj = doseCurr, t0 = c(1, 17), n0 = c(30, 30), psdose = psDoses),
                                        priorpredy1=indivyj1(dj = doseCurr, t0 = c(1, 17), n0 = c(30, 30), psdose = psDoses))
                         )
  )
  
  stage[nCoh] = ifelse(all(predSeq[1:nCoh]*3 == toxSeq[1:nCoh]), 1, 2)
  
  indUtility = UtiCal(ypred = rep(predSeq[nCoh], 3), ytrue = c(numeric(3-toxSeq[nCoh]), rep(1, toxSeq[nCoh])))
  
  UtiDose[index] = UtiDose[index] + sum(indUtility)
  maxU[index] = maxU[index] + 3
  pa = UtiDose/maxU
  pa.constraint = pa[(index-1):length(doseH)]
  predAc[nCoh] = mean(pa.constraint[is.finite(pa.constraint)])
  # predAc[nCoh] = sum(pa*NsubH, na.rm=TRUE)/sum(NsubH)
  
  if(length(NsubH[NsubH>0]) == 1){
    predij <- indivyj1(dj = doseH, t0 = c(1, 17), n0 = c(30, 30), psdose = psDoses)
  } else{
    predij <- MCMCSim$mean$pTox.star 
  }
  
  # w.p of p1 to observe a DLT from the 1st patient who will receive the lowest dose
  p1 = predij[index]
  # w.p of p0 to observe a no-DLT from
  p0 = 1 - p1
  # the discrete random variable \bar{a}(eta, y) = x1 when y=1; and x0 when y=0
  x1 = UtiCal(ypred = predSeq[nCoh], ytrue = 1)/1
  x0 = UtiCal(ypred = predSeq[nCoh], ytrue = 0)/1
  xbar = p1*x1+p0*x0
  ## 1. the numerator
  sd_yk = sqrt(p1*x1^2 + p0*x0^2 - xbar^2)
  
  minsd = sd_yrmn(N = CohSize, nk = nCoh, doseSet = doseH, ind = index, ppred = predij, 
                  t0=c(1, 17), n0=c(30, 30), psdose = psDoses)
  
  lambda[nCoh] = sd_yk/minsd
  
  # No run-in period
  wMix[1] = predAc[nCoh]^(sd_yk/minsd)
  
  # # With a run-in period
  # if(stage[nCoh] == 1){
  #   wMix[1] = 0
  # } else{
  #   wMix[1] = predAc[nCoh]^(sd_yk/minsd)
  # }
  
  wMix[2] = 1 - wMix[1]
  
  MCMCSim <- bugs(data, inits, parameters, "MixPrior_PhI.txt", codaPkg = F, bugs.seed = 1, 
                  n.chains = 1, n.burnin = 5000, n.iter = 15000)
  
  dmax[nCoh] <- ifelse(MCMCSim$mean$pCat[1,3]>0.25, 2, max(doseH[MCMCSim$mean$pCat[,3]<= 0.25]))
  dsel[nCoh] <- min(ifelse(index==9, 70, 
                           ifelse(index==8 | doseH[index+2]/doseH[index] > 2, doseH[index+1],
                                  doseH[index+2])), dmax[nCoh])
  simpTox[nCoh] <- MCMCSim$mean$pTox.star[match(dsel[nCoh], doseH)]
  
  interim <- data.frame(dose = doseH, NtoxH, NsubH, MCMCSim$mean$pCat)
  
  probUD[nCoh,] <- interim$X1
  probTI[nCoh,] <- interim$X2
  probOD[nCoh,] <- interim$X3
  priorEx[nCoh,] <- wMix
  postEx[nCoh,] <- MCMCSim$mean$prob.ex 
  
  # double check the numerical results displayed in the table
  intprob[nCoh, ] <- unlist(interim[match(dsel[nCoh], doseH), 4:6])
    
  print(nCoh)
  
  nCoh = nCoh + 1
  
}

tab = as.data.frame(cbind(stage, c(4, dmax[1:(length(dSeq)-1)]), c(4, dsel[1:(length(dSeq)-1)]), simpTox, 
                            round(intprob, 3), round(priorEx, 3), round(postEx, 3), predAc, lambda))
names(tab)[2:3] <- c("dmax", "dsel")
names(tab)[5:7] <- c("probU", "probT", "probO") 
names(tab)[8:9] <- c("Prior.wMixA", "Prior.wMixO")
names(tab)[10:11] <- c("Post.wMixA", "Post.wMixO")

return(tab)

}

#---------------------------------- Data example 1  -----------------------------------#
# Note: the input to toxSeq below was simulated 
# according to the stipulation defined in the manuscript
set.seed(1)
tab1 = NoRunIn(dSeq   = c(4, 8, 16, 22, 28, 22, 16, 16, 16, 22, 22),
               toxSeq = c(0, 0,  0,  0,  2,  2,  1,  0,  0,  1,  0))

MixW1 <- data.frame(prior = c(1, tab1$Prior.wMixA),
                    posterior = c(1, tab1$Post.wMixA))

#---------------------------------- Data example 2  -----------------------------------#
set.seed(1)
tab2 = NoRunIn(dSeq   = c(4, 2, 4, 8, 2, 4, 4, 4, 4, 4, 4),
               toxSeq = c(1, 0, 0, 3, 0, 0, 1, 0, 0, 1, 0))

MixW2 <- data.frame(prior = c(1, tab2$Prior.wMixA),
                    posterior = c(1, tab2$Post.wMixA))

#---------------------------------- Data example 3  -----------------------------------#
set.seed(1)
tab3 = NoRunIn(dSeq   = c(4, 8, 16, 22, 28, 54, 40, 40, 40, 40, 40),
               toxSeq = c(0, 0,  0,  0,  0,  1,  1,  0,  1,  1,  0))

MixW3 <- data.frame(prior = c(1, tab3$Prior.wMixA),
                    posterior = c(1, tab3$Post.wMixA))


##--------------------- Create Figure 4A ---------------------##
esdat1 <- data.frame(cohort = 1:11, 
                     dsel = tab1$dsel,
                     tox = c(0, 0,  0,  0,  2,  2,  1,  0,  0,  1,  0),
                     dsc = "Data example 1")

esdat2 <- data.frame(cohort = 1:11, 
                     dsel = tab2$dsel,
                     tox = c(1, 0, 0, 3, 0, 0, 1, 0, 0, 1, 0),
                     dsc = "Data example 2")

esdat3 <- data.frame(cohort = 1:11, 
                     dsel = tab3$dsel,
                     tox = c(0, 0,  0,  0,  0,  1,  1,  0,  1,  1,  0),
                     dsc = "Data example 3")



# rm(list=mats)
mats <- grep(x= ls(pos=1), pattern="esdat", value=TRUE)
esDFall <- do.call(rbind, mget(mixedsort(mats)))

toxSet = c(0, 1, 2, 3, -1)

p <- ggplot(esDFall, aes(x = cohort, y = match(dsel, doseH))) + theme_bw() + 
  xlab("Cohort number") + ylab(bquote("Dose ("*mg/m^2*")")) + scale_x_continuous(breaks = 0:11, limits = c(0, 11.2)) + 
  scale_y_continuous(breaks = 1:9, label = doseH[1:9], limits = c(0.7, 9.3)) +
  theme(legend.position='top') + 
  theme(axis.text.y=element_text(colour = c(rep("gray", 1), rep("black", 8)))) + 
  theme(axis.text.x=element_text(colour = c(rep("gray", 1), rep("black", 11))))

des <- p + geom_line(colour = "darkgray") + 
  facet_grid(~dsc) +
  geom_point(aes(shape = factor(match(tox, toxSet)), color = factor(match(tox, toxSet))), size = 2.5) + 
  scale_colour_manual(name = "DLT data", labels = c("0/3", "1/3", "2/3", "3/3", "Yet to observe"), 
                      values = c("#66CC99", "#56B4E9", "pink", "red", "darkgray")) +
  scale_shape_manual(name = "DLT data", labels = c("0/3", "1/3", "2/3", "3/3", "Yet to observe"), 
                     values = c(16, 17, 15, 18, 11)) 

# Figure 4A
des + theme(panel.grid.minor =  element_blank())

##--------------------- Create Figure 4B ---------------------##

for(i in 1:3) assign(paste0("MixWeight", i),
                     data.frame(cohort = 0:11,
                                dsc = paste0("Data example ", i), 
                                melt(get(paste0("MixW", i))))
                     )

mats <- c(grep(x= ls(pos=1), pattern="MixWeight", value=TRUE))
probExDFall <- do.call(rbind, mget(mats))

pEx <- ggplot(probExDFall, aes(x = cohort, y = value, color = factor(variable))) + 
  labs(x = "Cohort number", y = "Mixture weight") + 
  scale_x_continuous(breaks = 0:11) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(-0.01, 1.01)) +
  theme_bw() + theme(legend.position='top', legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour = c(rep("gray", 1), rep("black", 11))))

probEx <- pEx + geom_line(aes(y = value, colour=factor(variable), linetype=factor(variable)), size=1) + 
  geom_point(aes(colour=factor(variable), shape = factor(variable)), size=1.6) + guides(colour = guide_legend(nrow = 1)) +
  scale_color_manual(name = element_blank(),
                     labels = c(expression(paste("Prior, ", w^(h), "     ")), expression(paste("Posterior, ", w["*"]^(h)))),
                     values = c("#0072B2", "#999999")) +
  scale_shape_manual(name = element_blank(), 
                     labels = c(expression(paste("Prior, ", w^(h), "     ")), expression(paste("Posterior, ", w["*"]^(h)))),
                     values = c(19, 21)) +
  scale_linetype_manual(name = element_blank(), 
                        labels = c(expression(paste("Prior, ", w^(h), "     ")), expression(paste("Posterior, ", w["*"]^(h)))),
                        values = c(1, 5)) + facet_wrap(~dsc, ncol = 4) 

# Figure 4B
probEx + theme(panel.grid.minor =  element_blank())

## Figure 4 ##
cowplot::plot_grid(des+ theme(panel.grid.minor =  element_blank()), 
                   probEx + theme(panel.grid.minor =  element_blank()), labels = c("A", "B"), ncol = 1, align = "v")

##--------------------- Create Table S1 ---------------------##

commens1 <- data.frame(kappa = tab1$predAc, lambda = tab1$lambda,
                                w = tab1$predAc^tab1$lambda)

commens2 <- data.frame(kappa = tab2$predAc, lambda = tab2$lambda,
                                w = tab2$predAc^tab2$lambda)

commens3 <- data.frame(kappa = tab3$predAc, lambda = tab3$lambda,
                       w = tab3$predAc^tab3$lambda)


## Table S1
rbind(t(round(commens1, 2)), t(round(commens2, 2)), t(round(commens3, 2)))

# save(list = ls(), file = "Simulated data examples (NoRunIn Proc).rds")
