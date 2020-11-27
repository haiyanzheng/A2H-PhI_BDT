## Load all other user-defined R functions
source("Sub-functions.R")
source("Procedure A (No Run In).R")
source("Procedure B (Run In).R")
source("Procedures CDE.R")

## Load the required R package(s)
require(R2OpenBUGS)
require(mvtnorm)
require(gtools)

set.seed(123)
#-------------------------------------------------------------#
# Scenario 1
Sc1 = c(0.11, 0.25, 0.35, 0.41, 0.47, 0.52, 0.58, 0.63, 0.70)
# Scenario 2
Sc2 = c(0.08, 0.16, 0.25, 0.35, 0.42, 0.45, 0.53, 0.60, 0.70)
# Scenario 3
Sc3 = c(0.02, 0.05, 0.14, 0.25, 0.35, 0.42, 0.51, 0.60, 0.68)
# Scenario 4
Sc4 = c(0.03, 0.05, 0.10, 0.16, 0.25, 0.32, 0.40, 0.48, 0.55)
# Scenario 5
Sc5 = c(0.001, 0.005, 0.03, 0.10, 0.16, 0.25, 0.38, 0.50, 0.60)
# Scenario 6
Sc6 = c(0.01, 0.02, 0.05, 0.08, 0.11, 0.14, 0.25, 0.37, 0.47)
# Scenario 7
Sc7 = c(0.35, 0.42, 0.60, 0.75, 0.82, 0.88, 0.91, 0.94, 0.97)
# Scenario 8
Sc8 = c(0.001, 0.005, 0.01, 0.02, 0.04, 0.05, 0.10, 0.16, 0.25)
#-------------------------------------------------------------#
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

#-------------------------- Scenario 1 -----------------------------------#
set.seed(123)
ProcASc1 = ProcA(nsims = 10, pTrue = Sc1, Myu01 = 0.6)
# ProcASc1 = ProcA(nsims = 10, pTrue = Sc1, Myu01 = 0.2)

set.seed(123)
ProcBSc1 = ProcB(nsims = 10, pTrue = Sc1, Myu01 = 0.6)
# ProcBSc1 = ProcB(nsims = 10, pTrue = Sc1, Myu01 = 0.2)

set.seed(123)
ProcCSc1 = ProcFixW(nsims = 10, pTrue = Sc1, MyMixW = 0.5)

set.seed(123)
ProcDSc1 = ProcFixW(nsims = 10, pTrue = Sc1, MyMixW = 1)

set.seed(123)
ProcESc1 = ProcFixW(nsims = 10, pTrue = Sc1, MyMixW = 0)

#-------------------------- Scenario 2 -----------------------------------#
set.seed(123)
ProcASc2 = ProcA(nsims = 10, pTrue = Sc2, Myu01 = 0.6)
# ProcASc2 = ProcA(nsims = 10, pTrue = Sc2, Myu01 = 0.2)

set.seed(123)
ProcBSc2 = ProcB(nsims = 10, pTrue = Sc2, Myu01 = 0.6)
# ProcBSc2 = ProcB(nsims = 10, pTrue = Sc2, Myu01 = 0.2)

set.seed(123)
ProcCSc2 = ProcFixW(nsims = 10, pTrue = Sc2, MyMixW = 0.5)

set.seed(123)
ProcDSc2 = ProcFixW(nsims = 10, pTrue = Sc2, MyMixW = 1)

set.seed(123)
ProcESc2 = ProcFixW(nsims = 10, pTrue = Sc2, MyMixW = 0)

# to generate a compate Figure 6, please run the functions above for scenarios 3 - 8,
# and preserve the results in the object, named as, ProcpSci, 
# where p should be substituted by letters A to E, and i by 3 to 8 (e.g., ProcASc3).

# save(list = ls(), file = "Simulation study (comparison of the procedures).rds")

##--------------------- Create Figure 6/S2/S3 ---------------------##

# Here, we visulise results of two secenarios for illustration
# To reproduce the complete figure 6/S2/S3,
# the user needs to obtain results for additional scenarios using above the R functions,
# namely, ProcA(), ProcB(), ProcFixW(Myu01, MyMixW)

# preserve the operating characteristics, including dose, true toxicity rate (pTox), 
# average number of patients treated per dose (avgPts), and percentage of a dose selected as MTD

for(p in LETTERS[1:5]){
  for(i in 1:2) assign(paste0("Sim", p, "Sc", i),
                     merge(x = data.frame(dose = c(0, 2, 4, 8, 16, 22, 28, 40, 54, 70), 
                                          pTox = c(NA, get(paste0("Sc", i))),
                                          avgPts = c(NA, get(paste0("Proc", p, "Sc", i))[[2]][,1])/10), 
                           y = as.data.frame(table(get(paste0("Proc", p, "Sc", i))[[1]])/10), 
                           by.x = "dose", by.y = "Var1", all = TRUE)
                     )
}


# organise the results with 'dose' ascending in the data frame
# add information of the scenarios and designs used

for(p in LETTERS[1:5]){
  for(i in 1:2) assign(paste0("clSim", p, "Sc", i),
                       data.frame(
                         DoseH = 1:10,
                         get(paste0("Sim", p, "Sc", i))[order(nchar(get(paste0("Sim", p, "Sc", i))$dose)),],
                         Sc = paste0("Scenario ", i),
                         Mod = paste0("Procedure ", p)
                         ))
}

mats <- grep(x= ls(pos=1), pattern="clSim", value=TRUE)
PCSresults <- do.call(rbind, mget(mixedsort(mats)))

PCSresults$hline = ifelse(PCSresults$Sc == "Scenario 1", 3,
                         ifelse(PCSresults$Sc == "Scenario 2", 4, 
                                ifelse(PCSresults$Sc == "Scenario 3", 5,
                                       ifelse(PCSresults$Sc == "Scenario 4", 6, 
                                              ifelse(PCSresults$Sc == "Scenario 5", 7,
                                                     ifelse(PCSresults$Sc == "Scenario 6", 8,
                                                            ifelse(PCSresults$Sc == "Scenario 7", 0, 9)
                                                     )
                                              )
                                       )
                                )
                         )
)

# Subfigure (i) of Figure 6/S2/S3
p1 <- ggplot(PCSresults, aes(x = factor(DoseH), y = Freq*100)) + theme_bw() + 
  scale_y_continuous(limits = c(0, 100)) + 
  ylab("Proportion of times of dose declared as MTD (%)") + 
  xlab(bquote("Doses ("*mg/m^2*")")) + scale_x_discrete(label = c("*", 2, 4, 8, 16, 22, 28, 40, 54, 70)) +
  geom_vline(aes(xintercept = hline)) + theme(legend.position='top') 

PCSplot <- p1  + facet_wrap(~Sc, ncol = 4) +
  geom_point(aes(shape = factor(Mod), color = factor(Mod)), size = 2.5) + 
  scale_colour_manual(name = "Procedure ",
                      labels = c("A", "B", "C", "D", "E"),
                      values = c("pink", "#D55E00", "#56B4E9", "#0072B2", "lightgray")) +
  scale_shape_manual(name = "Procedure ",
                     labels = c("A", "B", "C", "D", "E"),
                     values = c(17, 17, 16, 15, 6)) 

# Subfigure 6(i)
PCSplot


# Subfigure (ii) of Figure 6/S2/S3
p2 <- ggplot(PCSresults, aes(x = factor(DoseH), y = avgPts)) + theme_bw() + 
  scale_y_continuous(limits = c(0, 17)) + 
  ylab("Average number of patients allocated") + 
  xlab(bquote("Doses ("*mg/m^2*")")) + 
  scale_x_discrete(label = c("*", 2, 4, 8, 16, 22, 28, 40, 54, 70)) +
  geom_vline(aes(xintercept = hline)) + theme(legend.position='top') 

Ptsplot <- p2  + facet_wrap(~Sc, ncol = 4) +
  geom_point(aes(shape = factor(Mod), color = factor(Mod)), size = 2.5) + 
  scale_colour_manual(name = "Procedure ",
                      labels = c("A", "B", "C", "D", "E"),
                      values = c("pink", "#D55E00", "#56B4E9", "#0072B2", "lightgray")) +
  scale_shape_manual(name = "Procedure ",
                     labels = c("A", "B", "C", "D", "E"),
                     values = c(17, 17, 16, 15, 6)) 

# Subfigure 6(ii)
Ptsplot

## Figure 6 ##
# cowplot::plot_grid(PCSplot, 
#                    Ptsplot, labels = c("(i)", "(ii)"), ncol = 1, align = "v")

# To re-generate Figures S2 and S3, users need to uncomment line 101 and 108 of 
# [[./Procedure A (No Run In).R]], [[./Procedure B (Run In).R]]
# so as to switch to using wMix[1] = sqrt(CohSize/nCoh)


##--------------------- Create Table S3 ---------------------##
# an example of the table when we investigate Scenarios 1 and 2 only

# Subsection of the table, with results of Scenario 1
rbind(t(SimASc1), t(SimBSc1[, c("avgPts", "Freq")]),
      t(SimCSc1[, c("avgPts", "Freq")]), t(SimDSc1[, c("avgPts", "Freq")]), 
      t(SimESc1[, c("avgPts", "Freq")]))


# Subsection of the table, with results of Scenario 2
rbind(t(SimASc2), t(SimBSc2[, c("avgPts", "Freq")]),
      t(SimCSc2[, c("avgPts", "Freq")]), t(SimDSc2[, c("avgPts", "Freq")]), 
      t(SimESc2[, c("avgPts", "Freq")]))

# Author's note: the NA values especially of the Freq's are 
# because of limited simulation replicates
# the Freq's below dose = 0.0 correspond to 'None' of the Table S2 presented,
# which should be read as proportion of trials stopped early for safety.
