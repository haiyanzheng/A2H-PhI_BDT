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

FixW = function(dSeq, toxSeq, MyMixW = 0){

  nCoh = as.integer(1)
  CohSize = 11
  safestop = as.integer(0)
    
  NsubH = c(0, 0, 0, 0, 0, 0, 0, 0, 0)
  NtoxH = c(0, 0, 0, 0, 0, 0, 0, 0, 0)
  
dmax = numeric(length(dSeq))
dsel = numeric(length(dSeq))
predSeq = numeric(length(dSeq))
simpTox = numeric(length(dSeq))
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


while(nCoh <= length(dSeq)){
  
  doseCurr = dSeq[nCoh]
  index = match(dSeq[nCoh], doseH)
  # default number of patient as 3 each administration
  NsubH[index] <- NsubH[index] + 3
  NtoxH[index] <- NtoxH[index] + toxSeq[nCoh]
  
  wMix[1] = MyMixW
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
                            round(intprob, 3), round(priorEx, 3), round(postEx, 3)))
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
tab1O = FixW(dSeq   = c(4, 8, 16, 28, 22, 22, 16, 16, 16, 16, 16),
            toxSeq = c(0, 0,  0,  2,  0,  2,  1,  0,  0,  0,  1),
            MyMixW = 0)

set.seed(1)
tab1A = FixW(dSeq   = c(4, 8, 16, 22, 28, 22, 16, 16, 16, 22, 22),
             toxSeq = c(0, 0,  0,  0,  2,  2,  1,  0,  0,  1,  0),
             MyMixW = 1)

# check whether column 'dsel' is identical to the argument 'dSeq'
tab1O; tab1A
#---------------------------------- Data example 2  -----------------------------------#
set.seed(1)
tab2O = FixW(dSeq = c(4, 2, 4, 8, 2, 4, 4, 4, 4, 4, 4),
            toxSeq = c(1, 0, 0, 3, 0, 0, 1, 0, 0, 1, 0),
            MyMixW = 0)

set.seed(1)
tab2A = FixW(dSeq = c(4, 8, 4, 8, 4, 4, 4, 4, 8, 4, 4),
             toxSeq = c(1, 3, 0, 2, 0, 1, 0, 0, 2, 1, 0),
             MyMixW = 1)

# check whether column 'dsel' is identical to the argument 'dSeq'
tab2O; tab2A
#---------------------------------- Data example 3  -----------------------------------#
set.seed(1)
tab3O = FixW(dSeq = c(4, 8, 16, 28, 40, 40, 54, 40, 40, 40, 40),
             toxSeq = c(0, 0,  0,  0,  1,  0,  1,  1,  1,  0,  0),
             MyMixW = 0)

set.seed(1)
tab3A = FixW(dSeq = c(4, 8, 16, 22, 28, 28, 28, 28, 28, 28, 28),
             toxSeq = c(0, 0,  0,  0,  0,  1,  0,  0,  1,  0,  0),
             MyMixW = 1)

# check whether column 'dsel' is identical to the argument 'dSeq'
tab3O; tab3A

##--------------------- Create Figure S1 ---------------------##
esdat1O <- data.frame(cohort = 1:11, 
                     dsel = tab1O$dsel,
                     tox = c(0, 0,  0,  2,  0,  2,  1,  0,  0,  0,  1),
                     dsc = "Data example 1",
                     type = "Operational prior")

esdat1A <- data.frame(cohort = 1:11, 
                      dsel = tab1A$dsel,
                      tox = c(0, 0,  0,  0,  2,  2,  1,  0,  0,  1,  0),
                      dsc = "Data example 1",
                      type = "Animal prior")

esdat2O <- data.frame(cohort = 1:11, 
                     dsel = tab2O$dsel,
                     tox = c(1, 0, 0, 3, 0, 0, 1, 0, 0, 1, 0),
                     dsc = "Data example 2",
                     type = "Operational prior")

esdat2A <- data.frame(cohort = 1:11, 
                      dsel = tab2A$dsel,
                      tox = c(1, 3, 0, 2, 0, 1, 0, 0, 2, 1, 0),
                      dsc = "Data example 2",
                      type = "Animal prior")

esdat3O <- data.frame(cohort = 1:11, 
                     dsel = tab3O$dsel,
                     tox = c(0, 0,  0,  0,  1,  0,  1,  1,  1,  0,  0),
                     dsc = "Data example 3",
                     type = "Operational prior")

esdat3A <- data.frame(cohort = 1:11, 
                      dsel = tab3A$dsel,
                      tox = c(0, 0,  0,  0,  0,  1,  0,  0,  1,  0,  0),
                      dsc = "Data example 3",
                      type = "Animal prior")

# rm(list=mats)
mats <- grep(x= ls(pos=1), pattern="esdat", value=TRUE)
esDFall <- do.call(rbind, mget(mixedsort(mats)))

toxSet = c(0, 1, 2, 3, -1)

p <- ggplot(esDFall, aes(x = cohort, y = match(dsel, doseH))) + theme_bw() + 
  xlab("Cohort number") + ylab(bquote("Dose ("*mg/m^2*")")) + scale_x_continuous(breaks = 1:11, limits = c(0.8, 11.2)) + 
  scale_y_continuous(breaks = 1:9, label = doseH[1:9], limits = c(0.7, 9.3)) +
  theme(legend.position='top') + 
  theme(axis.text.y=element_text(colour = c(rep("gray", 1), rep("black", 8)))) 

des <- p + geom_line(colour = "darkgray") + 
  facet_grid(factor(type, levels = c("Operational prior", "Animal prior"))~dsc) +
  geom_point(aes(shape = factor(match(tox, toxSet)), color = factor(match(tox, toxSet))), size = 2.5) + 
  scale_colour_manual(name = "DLT data", labels = c("0/3", "1/3", "2/3", "3/3", "Yet to observe"), 
                      values = c("#66CC99", "#56B4E9", "pink", "red", "darkgray")) +
  scale_shape_manual(name = "DLT data", labels = c("0/3", "1/3", "2/3", "3/3", "Yet to observe"), 
                     values = c(16, 17, 15, 18, 11)) 

# Figure S1
des + theme(panel.grid.minor =  element_blank())

# save(list = ls(), file = "Simulated data examples (Figure S1).rds")